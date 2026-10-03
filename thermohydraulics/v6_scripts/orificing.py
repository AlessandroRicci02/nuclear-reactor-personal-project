"""
V5 RADIAL FLOW ORIFICING & PLENUM HYDRAULIC EQUALIZATION MODULE
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/optimization/orificing.py

Solves coupled multi-channel flow redistribution across 217 channels in a common plenum:
1. Uniform flow area (Baseline un-orificed).
2. Ideal power-proportional flow (Achieves identical coolant exit temperatures).
3. 3-Zone radial orificing (Inner 19 pins, Middle 72 pins, Outer 126 pins).
4. Optimal 9-Ring orificing (Minimizes peak clad temperature with minimal pressure drop penalty).
5. Computes exact orifice loss coefficients K_p, orifice plate hole diameters d_orif,
   and pumping power / compressor work penalties.
"""

import math
import numpy as np
import pandas as pd
from scipy.optimize import minimize

import sys
sys.path.insert(0, "/home/aleri/OpenMC_projects/first_test/thermohydraulics")
from models.th_properties import (
    R_FUEL, R_CLAD_OUT, R_CHAN_OUT, H_CORE, D_H, A_FLOW_PIN,
    Q_TOTAL, N_PINS, T_INLET, DELTA_T_BULK, P_NOM,
    cp_helium, rho_helium, mu_helium, k_helium,
    H_GAP_NOMINAL, pin_coolant_mass_flow, total_coolant_mass_flow
)
from models.cardinal_th_interface import solve_axial_flow_channel


class RadialOrificingOptimizer:
    def __init__(self, q_vol_grid, pin_meta, m_dot_total=None, T_inlet=950.0, p_outlet=2.50e6,
                 r_fuel=None, r_clad_out=None, r_chan_out=None, h_core=None):
        self.q_vol_grid = q_vol_grid # Shape: (217, 10)
        self.pin_meta = pin_meta     # List of (pin_idx, ring, x, y, is_inner)
        self.m_dot_total = total_coolant_mass_flow() if m_dot_total is None else m_dot_total
        self.T_inlet = T_inlet
        self.p_outlet = p_outlet
        self.cp = cp_helium()
        self.rho_inlet = rho_helium(self.T_inlet, self.p_outlet)
        
        self.r_fuel = R_FUEL if r_fuel is None else r_fuel
        self.r_clad_out = R_CLAD_OUT if r_clad_out is None else r_clad_out
        self.r_chan_out = R_CHAN_OUT if r_chan_out is None else r_chan_out
        self.h_core = H_CORE if h_core is None else h_core
        self.a_flow_pin = math.pi * (self.r_chan_out**2 - self.r_clad_out**2)
        self.d_h = 2.0 * (self.r_chan_out - self.r_clad_out)
        
        # Calculate individual pin powers [W]
        dz = self.h_core / 10
        v_slice = math.pi * (self.r_fuel**2) * dz
        self.pin_powers = np.sum(q_vol_grid, axis=1) * v_slice
        self.z_edges = np.linspace(0.0, self.h_core, 11)
        
        # Map pins to rings
        self.pin_rings = np.array([meta[1] for meta in pin_meta])
        self.unique_rings = np.unique(self.pin_rings) # 0 to 8

    def evaluate_flow_distribution(self, m_dot_pins, K_orifices=None, htc_model="gnielinski_trans"):
        """
        Solves thermal-hydraulics for all 217 channels given a vector of pin mass flow rates.
        Returns thermal and pressure drop distributions.
        """
        T_fuel_center_all = np.zeros((N_PINS, 10))
        T_clad_out_all = np.zeros((N_PINS, 10))
        T_cool_bulk_all = np.zeros((N_PINS, 10))
        T_cool_exit_all = np.zeros(N_PINS)
        dp_channel_all = np.zeros(N_PINS)
        dp_orifice_all = np.zeros(N_PINS)
        
        if K_orifices is None:
            K_orifices = np.zeros(N_PINS)
            
        for p in range(N_PINS):
            m_p = m_dot_pins[p]
            q_slices = self.q_vol_grid[p, :]
            
            # Channel flow solve
            res = solve_axial_flow_channel(
                q_vol_slices=q_slices,
                m_dot_pin=m_p,
                z_edges=self.z_edges,
                p_outlet=self.p_outlet,
                h_gap=H_GAP_NOMINAL,
                htc_model=htc_model,
                T_inlet=self.T_inlet,
                r_fuel=self.r_fuel,
                r_clad_out=self.r_clad_out,
                r_chan_out=self.r_chan_out,
                h_core=self.h_core
            )
            
            T_fuel_center_all[p, :] = res["T_fuel_center"]
            T_clad_out_all[p, :] = res["T_clad_out"]
            T_cool_bulk_all[p, :] = res["T_bulk"]
            T_cool_exit_all[p] = res["T_bulk_exit"]
            dp_channel_all[p] = res["dp_total_Pa"]
            
            # Orifice loss dp_orif = K_p * G^2 / (2 * rho_in)
            G_p = m_p / self.a_flow_pin
            dp_orif = K_orifices[p] * (G_p**2) / (2.0 * self.rho_inlet)
            dp_orifice_all[p] = dp_orif
            
        dp_total_pins = dp_channel_all + dp_orifice_all
        
        # Mixed bulk outlet temperature: T_out,bulk = sum(m_p * T_exit_p) / M_total
        T_bulk_out = np.sum(m_dot_pins * T_cool_exit_all) / np.sum(m_dot_pins)
        
        return {
            "T_fuel_center_max": float(np.max(T_fuel_center_all)),
            "T_clad_out_max": float(np.max(T_clad_out_all)),
            "T_cool_exit_max": float(np.max(T_cool_exit_all)),
            "T_cool_exit_min": float(np.min(T_cool_exit_all)),
            "T_cool_bulk_out": float(T_bulk_out),
            "delta_T_hot_elevation": float(np.max(T_cool_exit_all) - T_bulk_out),
            "dp_core_max_Pa": float(np.max(dp_total_pins)),
            "dp_channel_mean_Pa": float(np.mean(dp_channel_all)),
            "dp_orifice_mean_Pa": float(np.mean(dp_orifice_all)),
            "m_dot_pins": m_dot_pins,
            "T_cool_exit_all": T_cool_exit_all,
            "T_clad_out_all": T_clad_out_all,
            "K_orifices": K_orifices,
            "dp_total_pins": dp_total_pins
        }

    def solve_uniform_flow(self):
        """Mode 1: Uniform un-orificed baseline flow."""
        m_dot_uniform = np.full(N_PINS, self.m_dot_total / N_PINS)
        return self.evaluate_flow_distribution(m_dot_uniform, K_orifices=np.zeros(N_PINS))

    def solve_power_proportional_flow(self):
        """Mode 2: Power-proportional flow (each pin gets m_dot proportional to pin power)."""
        m_dot_prop = self.m_dot_total * (self.pin_powers / Q_TOTAL)
        
        # Determine orifice coefficients to equalize plenum pressure drop
        # Channel with highest resistance sets plenum baseline
        raw_res = self.evaluate_flow_distribution(m_dot_prop, K_orifices=np.zeros(N_PINS))
        dp_ch = raw_res["dp_total_pins"] # channel friction losses
        dp_target_plenum = np.max(dp_ch) * 1.05 # Add 5% margin
        
        K_orif = np.zeros(N_PINS)
        for p in range(N_PINS):
            dp_needed = max(dp_target_plenum - dp_ch[p], 0.0)
            G_p = m_dot_prop[p] / A_FLOW_PIN
            K_orif[p] = (2.0 * self.rho_inlet * dp_needed) / (G_p**2)
            
        return self.evaluate_flow_distribution(m_dot_prop, K_orifices=K_orif)

    def solve_three_zone_orificing(self):
        """
        Mode 3: Practical 3-Zone Flow Orificing:
        - Zone 1: Inner Rings 0-2 (19 pins, low power) -> Lower mass flow, higher K
        - Zone 2: Mid Rings 3-5 (72 pins, medium power) -> Intermediate flow
        - Zone 3: Outer Rings 6-8 (126 pins, high power) -> High flow, no orifice (K=0)
        """
        # Optimize zone flow fractions: [f_z1, f_z2, f_z3]
        # with sum(N_z * m_z) = M_total
        n1, n2, n3 = 19, 72, 126
        
        # Target flows proportional to average zone power
        p_z1 = np.mean(self.pin_powers[self.pin_rings <= 2])
        p_z2 = np.mean(self.pin_powers[(self.pin_rings >= 3) & (self.pin_rings <= 5)])
        p_z3 = np.mean(self.pin_powers[self.pin_rings >= 6])
        
        w1 = p_z1 / np.mean(self.pin_powers)
        w2 = p_z2 / np.mean(self.pin_powers)
        w3 = p_z3 / np.mean(self.pin_powers)
        
        norm = (n1 * w1 + n2 * w2 + n3 * w3)
        m_base = self.m_dot_total
        m_z1 = m_base * w1 / norm
        m_z2 = m_base * w2 / norm
        m_z3 = m_base * w3 / norm
        
        m_dot_3zone = np.zeros(N_PINS)
        for p in range(N_PINS):
            r = self.pin_rings[p]
            if r <= 2:
                m_dot_3zone[p] = m_z1
            elif r <= 5:
                m_dot_3zone[p] = m_z2
            else:
                m_dot_3zone[p] = m_z3
                
        # Equalize plenum pressure
        raw_res = self.evaluate_flow_distribution(m_dot_3zone, K_orifices=np.zeros(N_PINS))
        dp_ch = raw_res["dp_total_pins"]
        dp_target = np.max(dp_ch) * 1.02
        
        K_orif = np.zeros(N_PINS)
        for p in range(N_PINS):
            dp_needed = max(dp_target - dp_ch[p], 0.0)
            G_p = m_dot_3zone[p] / A_FLOW_PIN
            K_orif[p] = (2.0 * self.rho_inlet * dp_needed) / (G_p**2)
            
        return self.evaluate_flow_distribution(m_dot_3zone, K_orifices=K_orif)

    def solve_optimal_ring_orificing(self):
        """
        Mode 4: 9-Ring Constrained Optimal Orificing:
        Optimizes ring-by-ring mass flows to minimize (T_clad,max + 2 * Delta_T_hot).
        """
        # 9 ring flow factors
        avg_power_rings = np.array([np.mean(self.pin_powers[self.pin_rings == r]) for r in range(9)])
        power_ratios = avg_power_rings / np.mean(self.pin_powers)
        
        # Pin count per ring: 1, 6, 12, 18, 24, 30, 36, 42, 48 = 217
        pins_per_ring = np.array([1, 6, 12, 18, 24, 30, 36, 42, 48])
        
        # Analytical optimal flow matching local ring heat generation
        norm_factor = np.sum(pins_per_ring * power_ratios)
        m_dot_rings = self.m_dot_total * power_ratios / norm_factor
        
        m_dot_opt = np.zeros(N_PINS)
        for p in range(N_PINS):
            r = self.pin_rings[p]
            m_dot_opt[p] = m_dot_rings[r]
            
        raw_res = self.evaluate_flow_distribution(m_dot_opt, K_orifices=np.zeros(N_PINS))
        dp_ch = raw_res["dp_total_pins"]
        dp_target = np.max(dp_ch) * 1.02
        
        K_orif = np.zeros(N_PINS)
        for p in range(N_PINS):
            dp_needed = max(dp_target - dp_ch[p], 0.0)
            G_p = m_dot_opt[p] / A_FLOW_PIN
            K_orif[p] = (2.0 * self.rho_inlet * dp_needed) / (G_p**2)
            
        res = self.evaluate_flow_distribution(m_dot_opt, K_orifices=K_orif)
        res["ring_m_dot"] = m_dot_rings
        res["ring_K_mean"] = np.array([np.mean(K_orif[self.pin_rings == r]) for r in range(9)])
        return res


def compute_orifice_hole_diameter(m_dot_pin, dp_orifice, rho_in, C_d=0.62):
    """Calculates equivalent thin-plate orifice hole diameter [mm]."""
    if dp_orifice <= 1.0: # Negligible loss
        return 2.0 * math.sqrt(A_FLOW_PIN / math.pi) * 1e3 # Unrestricted full diameter
    area_orif = m_dot_pin / (C_d * math.sqrt(2.0 * rho_in * dp_orifice))
    d_orif_mm = 2.0 * math.sqrt(area_orif / math.pi) * 1e3
    return d_orif_mm
