"""
CLOSED-LOOP PICARD COUPLING DRIVER FOR OPT-17 CAND 17-CLEAN
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/coupling/picard_driver.py

Implements rigorous bidirectional Picard iteration with under-relaxation between
OpenMC 3D neutronics and the annular thermohydraulic solver.
"""

import os
import sys
import time
import json
import math
import numpy as np
import pandas as pd
import openmc

sys.path.insert(0, "/home/aleri/OpenMC_projects/first_test/thermohydraulics")
from models.th_properties import (
    R_GAS, M_HE, P_NOM, N_PINS, H_CORE, R_FUEL, R_CLAD_OUT, R_CHAN_OUT,
    Q_TOTAL, T_INLET, DELTA_T_BULK, cp_helium, pin_coolant_mass_flow,
    total_coolant_mass_flow, rho_helium, h_convection, pressure_drop_channel,
    H_GAP_NOMINAL
)
from models.moose_th_solver import solve_pin_thermal_profile
from inputs.build_opt17_openmc import (
    build_opt17_openmc_model, OPENMC_EXEC, CROSS_SECTIONS_PATH, N_AXIAL_SLICES
)

WORKSPACE_ROOT = "/home/aleri/OpenMC_projects/first_test/thermohydraulics"


class PicardCoupledEngine:
    def __init__(self,
                 work_dir,
                 particles=80000,
                 batches=45,
                 inactive=10,
                 theta_relaxation=0.60,
                 h_gap=H_GAP_NOMINAL,
                 enable_density_feedback=True,
                 enable_moderator_feedback=True,
                 htc_multiplier=1.0,
                 htc_model="gnielinski_trans"):
        self.work_dir = work_dir
        self.particles = particles
        self.batches = batches
        self.inactive = inactive
        self.theta = theta_relaxation
        self.h_gap = h_gap
        self.enable_density_feedback = enable_density_feedback
        self.enable_moderator_feedback = enable_moderator_feedback
        self.htc_multiplier = htc_multiplier
        self.htc_model = htc_model
        
        self.m_dot_pin = pin_coolant_mass_flow()
        self.m_dot_tot = total_coolant_mass_flow()
        self.z_edges = np.linspace(0.0, H_CORE, N_AXIAL_SLICES + 1)
        self.v_fuel_slice = math.pi * (R_FUEL**2) * (H_CORE / N_AXIAL_SLICES)
        self.v_fuel_core = N_PINS * math.pi * (R_FUEL**2) * H_CORE
        
        # Initial thermal state (Uniform 1050 K)
        self.T_fuel_vol = np.full((N_PINS, N_AXIAL_SLICES), 1050.0)
        self.T_fuel_center = np.full((N_PINS, N_AXIAL_SLICES), 1050.0)
        self.T_clad_avg = np.full((N_PINS, N_AXIAL_SLICES), 1050.0)
        self.T_clad_out = np.full((N_PINS, N_AXIAL_SLICES), 1050.0)
        self.T_cool_bulk = np.full((N_PINS, N_AXIAL_SLICES), 1050.0)
        self.T_mod = 800.0

    def run_neutronics_step(self, iter_idx):
        """Runs OpenMC for the current thermal state and extracts 3D power."""
        iter_dir = os.path.join(self.work_dir, f"iter_{iter_idx}")
        os.makedirs(iter_dir, exist_ok=True)
        
        T_cool_grid_to_use = self.T_cool_bulk if self.enable_density_feedback else np.full((N_PINS, N_AXIAL_SLICES), 1050.0)
        T_mod_to_use = self.T_mod if self.enable_moderator_feedback else 800.0
        
        model_info = build_opt17_openmc_model(
            run_dir=iter_dir,
            T_fuel_grid=self.T_fuel_vol,
            T_clad_grid=self.T_clad_avg,
            T_cool_grid=T_cool_grid_to_use,
            T_mod=T_mod_to_use,
            particles=self.particles,
            batches=self.batches,
            inactive=self.inactive
        )
        
        os.environ["OPENMC_CROSS_SECTIONS"] = CROSS_SECTIONS_PATH
        os.environ["OMP_NUM_THREADS"] = "8"
        
        t0 = time.time()
        cwd_orig = os.getcwd()
        os.chdir(iter_dir)
        try:
            openmc.run(openmc_exec=OPENMC_EXEC, output=False)
        finally:
            os.chdir(cwd_orig)
        dt_neut = time.time() - t0
        
        sp_path = os.path.join(iter_dir, f"statepoint.{self.batches}.h5")
        with openmc.StatePoint(sp_path) as sp:
            keff = sp.keff
            k_val = float(keff.nominal_value)
            k_std = float(keff.std_dev)
            tally = sp.get_tally(name="fuel_cell_heating")
            df = tally.get_pandas_dataframe()
            
        df_heat = df[df["score"] == "heating"]
        cell_map = model_info["fuel_cells_map"]
        id_to_coords = {c.id: coords for coords, c in cell_map.items()}
        
        raw_powers = np.zeros((N_PINS, N_AXIAL_SLICES))
        for _, row in df_heat.iterrows():
            cid = int(row["cell"])
            if cid in id_to_coords:
                p_idx, s_idx = id_to_coords[cid]
                raw_powers[p_idx, s_idx] = float(row["mean"])
                
        # Normalization to 400 kWth
        total_raw = np.sum(raw_powers)
        norm_factor = Q_TOTAL / total_raw
        pin_slice_powers = raw_powers * norm_factor
        
        pin_powers = np.sum(pin_slice_powers, axis=1)
        q_avg_pin = Q_TOTAL / N_PINS
        F_r = float(np.max(pin_powers) / q_avg_pin)
        peak_pin_idx = int(np.argmax(pin_powers))
        
        q_vol_3d = pin_slice_powers / self.v_fuel_slice
        q_avg_vol = Q_TOTAL / self.v_fuel_core
        peak_vol_density = float(np.max(q_vol_3d))
        F_q = peak_vol_density / q_avg_vol
        
        core_axial_slice_powers = np.sum(pin_slice_powers, axis=0)
        q_avg_slice = Q_TOTAL / N_AXIAL_SLICES
        axial_profile_norm = core_axial_slice_powers / q_avg_slice
        F_z = float(np.max(axial_profile_norm))
        
        return {
            "k_eff": k_val,
            "k_std": k_std,
            "pin_slice_powers": pin_slice_powers,
            "pin_powers": pin_powers,
            "q_vol_3d": q_vol_3d,
            "peak_pin_idx": peak_pin_idx,
            "F_r": F_r,
            "F_z": F_z,
            "F_q": F_q,
            "peak_vol_density": peak_vol_density,
            "pin_meta": model_info["pin_meta"],
            "runtime_s": dt_neut
        }

    def run_thermal_step(self, neut_res):
        """Solves thermohydraulics across all 217 fuel pins."""
        t0 = time.time()
        q_vol_3d = neut_res["q_vol_3d"]
        
        T_calc_fuel_vol = np.zeros((N_PINS, N_AXIAL_SLICES))
        T_calc_fuel_center = np.zeros((N_PINS, N_AXIAL_SLICES))
        T_calc_clad_out = np.zeros((N_PINS, N_AXIAL_SLICES))
        T_calc_clad_avg = np.zeros((N_PINS, N_AXIAL_SLICES))
        T_calc_cool_bulk = np.zeros((N_PINS, N_AXIAL_SLICES))
        h_matrix = np.zeros((N_PINS, N_AXIAL_SLICES))
        Re_matrix = np.zeros((N_PINS, N_AXIAL_SLICES))
        Nu_matrix = np.zeros((N_PINS, N_AXIAL_SLICES))
        rho_matrix = np.zeros((N_PINS, N_AXIAL_SLICES))
        dp_fric_arr = np.zeros(N_PINS)
        dp_accel_arr = np.zeros(N_PINS)
        
        for p in range(N_PINS):
            q_vol_p = q_vol_3d[p, :]
            sol_p = solve_pin_thermal_profile(
                q_vol_p, self.m_dot_pin, self.z_edges,
                h_gap=self.h_gap,
                htc_multiplier=self.htc_multiplier,
                htc_model=self.htc_model
            )
            T_calc_fuel_vol[p, :] = sol_p["T_fuel_vol"]
            T_calc_fuel_center[p, :] = sol_p["T_fuel_center"]
            T_calc_clad_out[p, :] = sol_p["T_clad_out"]
            T_calc_clad_avg[p, :] = 0.5 * (sol_p["T_clad_out"] + sol_p["T_clad_in"])
            T_calc_cool_bulk[p, :] = sol_p["T_bulk"]
            h_matrix[p, :] = sol_p["h_conv"]
            Re_matrix[p, :] = sol_p["Re"]
            Nu_matrix[p, :] = sol_p["Nu"]
            rho_matrix[p, :] = sol_p["rho_cool"]
            dp_fric_arr[p] = sol_p["dp_fric_total_Pa"]
            dp_accel_arr[p] = sol_p["dp_accel_total_Pa"]
            
        dt_th = time.time() - t0
        
        # Core summary metrics
        max_fuel_temp = float(np.max(T_calc_fuel_center))
        max_clad_temp = float(np.max(T_calc_clad_out))
        coolant_outlet_mixed = float(np.mean(T_calc_cool_bulk[:, -1]))
        peak_coolant_temp = float(np.max(T_calc_cool_bulk[:, -1]))
        
        # Form loss pressure drop
        rho_in = rho_helium(T_INLET)
        rho_out_mean = np.mean(rho_matrix[:, -1])
        u_in = self.m_dot_pin / (rho_in * A_FLOW_PIN)
        u_out_mean = self.m_dot_pin / (rho_out_mean * A_FLOW_PIN)
        dp_form = 0.5 * rho_in * (u_in**2) / 2.0 + 1.0 * rho_out_mean * (u_out_mean**2) / 2.0
        
        dp_total_core = float(np.mean(dp_fric_arr) + np.mean(dp_accel_arr) + dp_form)
        
        cp = cp_helium()
        total_energy_removed = np.sum(self.m_dot_pin * cp * (T_calc_cool_bulk[:, -1] - T_INLET))
        err_abs = abs(total_energy_removed - Q_TOTAL)
        err_rel = err_abs / Q_TOTAL
        
        return {
            "T_fuel_vol": T_calc_fuel_vol,
            "T_fuel_center": T_calc_fuel_center,
            "T_clad_out": T_calc_clad_out,
            "T_clad_avg": T_calc_clad_avg,
            "T_cool_bulk": T_calc_cool_bulk,
            "h_matrix": h_matrix,
            "Re_matrix": Re_matrix,
            "Nu_matrix": Nu_matrix,
            "rho_matrix": rho_matrix,
            "max_fuel_temp": max_fuel_temp,
            "max_clad_temp": max_clad_temp,
            "coolant_outlet_mixed": coolant_outlet_mixed,
            "peak_coolant_temp": peak_coolant_temp,
            "dp_total_core_Pa": dp_total_core,
            "total_energy_removed": total_energy_removed,
            "err_abs": err_abs,
            "err_rel": err_rel,
            "runtime_s": dt_th
        }
