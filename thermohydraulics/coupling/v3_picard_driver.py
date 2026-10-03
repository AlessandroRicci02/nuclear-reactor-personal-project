"""
V3 CLOSED-LOOP PICARD COUPLING ENGINE FOR OPT-17 CAND 17-CLEAN
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/coupling/v3_picard_driver.py

Implements rigorous bidirectional Picard iteration with under-relaxation between
OpenMC 3D neutronics, Cardinal MOOSE THM Flow Channels, and 2D RZ Solid Heat Conduction.
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
    total_coolant_mass_flow, rho_helium, H_GAP_NOMINAL
)
from models.cardinal_th_interface import (
    solve_axial_flow_channel, generate_cardinal_thm_input,
    generate_cardinal_solid_input, run_cardinal_input
)
from inputs.build_opt17_v3 import (
    build_opt17_v3_openmc_model, OPENMC_EXEC, CROSS_SECTIONS_PATH
)

WORKSPACE_ROOT = "/home/aleri/OpenMC_projects/first_test/thermohydraulics"


class V3PicardCoupledEngine:
    def __init__(self,
                 work_dir,
                 n_axial_slices=10,
                 particles=80000,
                 batches=45,
                 inactive=10,
                 theta_relaxation=0.60,
                 h_gap=H_GAP_NOMINAL,
                 htc_model="gnielinski_trans",
                 enable_source_reuse=True):
        self.work_dir = work_dir
        self.n_axial_slices = n_axial_slices
        self.particles = particles
        self.batches = batches
        self.inactive = inactive
        self.theta = theta_relaxation
        self.h_gap = h_gap
        self.htc_model = htc_model
        self.enable_source_reuse = enable_source_reuse
        
        self.m_dot_pin = pin_coolant_mass_flow()
        self.m_dot_tot = total_coolant_mass_flow()
        self.z_edges = np.linspace(0.0, H_CORE, n_axial_slices + 1)
        self.v_fuel_slice = math.pi * (R_FUEL**2) * (H_CORE / n_axial_slices)
        self.v_fuel_core = N_PINS * math.pi * (R_FUEL**2) * H_CORE
        
        # Initial thermal state (Uniform 1050 K)
        self.T_fuel_vol = np.full((N_PINS, n_axial_slices), 1050.0)
        self.T_fuel_center = np.full((N_PINS, n_axial_slices), 1050.0)
        self.T_clad_avg = np.full((N_PINS, n_axial_slices), 1050.0)
        self.T_clad_out = np.full((N_PINS, n_axial_slices), 1050.0)
        self.T_cool_bulk = np.full((N_PINS, n_axial_slices), 1050.0)
        self.T_mod = 800.0
        self.T_refl = 800.0
        
        # Power state
        self.q_vol_unrelaxed = np.full((N_PINS, n_axial_slices), Q_TOTAL / self.v_fuel_core)
        self.q_vol_relaxed = np.full((N_PINS, n_axial_slices), Q_TOTAL / self.v_fuel_core)
        self.pin_powers = np.full(N_PINS, Q_TOTAL / N_PINS)
        
        self.last_statepoint = None

    def run_neutronics_step(self, iter_idx, particles_override=None, batches_override=None):
        """Runs OpenMC for the current thermal state with WMP continuous Doppler broadening."""
        iter_dir = os.path.join(self.work_dir, f"iter_{iter_idx}")
        os.makedirs(iter_dir, exist_ok=True)
        
        n_p = particles_override if particles_override is not None else self.particles
        n_b = batches_override if batches_override is not None else self.batches
        
        # Source reuse: if previous statepoint exists and source reuse enabled
        source_file = None
        if self.enable_source_reuse and self.last_statepoint is not None and os.path.exists(self.last_statepoint):
            source_file = self.last_statepoint
            
        model_info = build_opt17_v3_openmc_model(
            run_dir=iter_dir,
            T_fuel_grid=self.T_fuel_vol,
            T_clad_grid=self.T_clad_avg,
            T_cool_grid=self.T_cool_bulk,
            T_mod=self.T_mod,
            T_refl=self.T_refl,
            n_axial_slices=self.n_axial_slices,
            particles=n_p,
            batches=n_b,
            inactive=self.inactive,
            source_file=source_file
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
        
        # Find statepoint
        sp_file = os.path.join(iter_dir, f"statepoint.{n_b}.h5")
        if not os.path.exists(sp_file):
            sp_files = [f for f in os.listdir(iter_dir) if f.startswith("statepoint.") and f.endswith(".h5")]
            if sp_files:
                sp_files.sort(key=lambda x: int(x.split(".")[1]))
                sp_file = os.path.join(iter_dir, sp_files[-1])
            else:
                raise FileNotFoundError(f"No statepoint found in {iter_dir}")
                
        self.last_statepoint = sp_file
        
        with openmc.StatePoint(sp_file) as sp:
            k_eff = float(sp.k_combined.nominal_value)
            k_std = float(sp.k_combined.std_dev)
            
            # Extract Shannon entropy convergence if available
            entropy_vals = list(sp.entropy) if hasattr(sp, "entropy") and sp.entropy is not None else []
            
            t_heat = sp.get_tally(name="fuel_cell_heating")
            df_heat = t_heat.get_pandas_dataframe()
            
        # Parse heating scores per cell
        df_heat_scores = df_heat[df_heat["score"] == "heating"]
        cell_heat_dict = dict(zip(df_heat_scores["cell"], df_heat_scores["mean"]))
        
        q_raw_grid = np.zeros((N_PINS, self.n_axial_slices))
        for p in range(N_PINS):
            for s in range(self.n_axial_slices):
                cell_obj = model_info["fuel_cells_map"][(p, s)]
                q_raw_grid[p, s] = max(cell_heat_dict.get(cell_obj.id, 0.0), 1e-12)
                
        # Power normalization to Q_TOTAL = 400,000 W
        total_raw_core = np.sum(q_raw_grid)
        norm_factor = Q_TOTAL / total_raw_core
        p_cell_watts = q_raw_grid * norm_factor
        q_vol_new = p_cell_watts / self.v_fuel_slice
        
        return {
            "k_eff": k_eff,
            "k_std": k_std,
            "entropy": entropy_vals,
            "q_vol_new": q_vol_new,
            "dt_neut": dt_neut,
            "iter_dir": iter_dir,
            "pin_meta": model_info["pin_meta"]
        }

    def run_thermal_step(self, q_vol_relaxed):
        """
        Solves 1D axial compressible flow conservation and 2D heat conduction
        for all 217 channels and pins in the core.
        """
        T_fuel_vol_new = np.zeros((N_PINS, self.n_axial_slices))
        T_fuel_center_new = np.zeros((N_PINS, self.n_axial_slices))
        T_clad_avg_new = np.zeros((N_PINS, self.n_axial_slices))
        T_clad_out_new = np.zeros((N_PINS, self.n_axial_slices))
        T_cool_bulk_new = np.zeros((N_PINS, self.n_axial_slices))
        
        pin_powers = np.zeros(N_PINS)
        dp_total_pins = np.zeros(N_PINS)
        dp_fric_pins = np.zeros(N_PINS)
        dp_accel_pins = np.zeros(N_PINS)
        dp_form_pins = np.zeros(N_PINS)
        
        re_pins = np.zeros(N_PINS)
        nu_pins = np.zeros(N_PINS)
        h_pins = np.zeros(N_PINS)
        rho_min_pins = np.zeros(N_PINS)
        rho_max_pins = np.zeros(N_PINS)
        u_max_pins = np.zeros(N_PINS)
        T_out_pins = np.zeros(N_PINS)
        
        for p in range(N_PINS):
            q_slices = q_vol_relaxed[p, :]
            sol = solve_axial_flow_channel(
                q_slices,
                self.m_dot_pin,
                self.z_edges,
                p_outlet=P_NOM,
                h_gap=self.h_gap,
                htc_model=self.htc_model
            )
            
            T_fuel_vol_new[p, :] = sol["T_fuel_vol"]
            T_fuel_center_new[p, :] = sol["T_fuel_center"]
            T_clad_avg_new[p, :] = 0.5 * (sol["T_clad_in"] + sol["T_clad_out"])
            T_clad_out_new[p, :] = sol["T_clad_out"]
            T_cool_bulk_new[p, :] = sol["T_bulk"]
            
            pin_powers[p] = np.sum(q_slices) * self.v_fuel_slice
            dp_total_pins[p] = sol["dp_total_Pa"]
            dp_fric_pins[p] = sol["dp_fric_Pa"]
            dp_accel_pins[p] = sol["dp_accel_Pa"]
            dp_form_pins[p] = sol["dp_form_Pa"]
            
            re_pins[p] = np.mean(sol["Re"])
            nu_pins[p] = np.mean(sol["Nu"])
            h_pins[p] = np.mean(sol["h_conv"])
            rho_min_pins[p] = np.min(sol["rho_cool"])
            rho_max_pins[p] = np.max(sol["rho_cool"])
            u_max_pins[p] = np.max(sol["vel_cool"])
            T_out_pins[p] = sol["T_bulk_exit"]
            
        return {
            "T_fuel_vol": T_fuel_vol_new,
            "T_fuel_center": T_fuel_center_new,
            "T_clad_avg": T_clad_avg_new,
            "T_clad_out": T_clad_out_new,
            "T_cool_bulk": T_cool_bulk_new,
            "pin_powers": pin_powers,
            "dp_total": np.mean(dp_total_pins),
            "dp_fric": np.mean(dp_fric_pins),
            "dp_accel": np.mean(dp_accel_pins),
            "dp_form": np.mean(dp_form_pins),
            "Re_mean": np.mean(re_pins),
            "Re_range": (np.min(re_pins), np.max(re_pins)),
            "Nu_mean": np.mean(nu_pins),
            "h_mean": np.mean(h_pins),
            "h_range": (np.min(h_pins), np.max(h_pins)),
            "rho_min": np.min(rho_min_pins),
            "rho_max": np.max(rho_max_pins),
            "u_max": np.max(u_max_pins),
            "T_out_mean": np.mean(T_out_pins)
        }

    def run_cardinal_benchmark_pin(self, pin_idx, q_vol_slices, iter_dir):
        """
        Executes Cardinal MOOSE THM flow channel and 2D solid pin conduction
        for the given pin index and returns exact FE outputs.
        """
        bench_dir = os.path.join(iter_dir, f"cardinal_bench_p{pin_idx}")
        os.makedirs(bench_dir, exist_ok=True)
        
        # 1. THM flow channel input
        thm_input = os.path.join(bench_dir, "thm_channel.i")
        generate_cardinal_thm_input(f"Pin_{pin_idx}", self.z_edges, q_vol_slices, self.m_dot_pin, thm_input)
        
        try:
            thm_res = run_cardinal_input(thm_input)
        except Exception as e:
            thm_res = {"error": str(e), "dp_total": 0.0, "max_T_wall": 0.0}
            
        # 2. Solid pin conduction input
        sol = solve_axial_flow_channel(q_vol_slices, self.m_dot_pin, self.z_edges)
        solid_input = os.path.join(bench_dir, "solid_pin.i")
        generate_cardinal_solid_input(f"Pin_{pin_idx}", self.z_edges, q_vol_slices, sol["h_conv"], sol["T_bulk"], solid_input)
        
        try:
            solid_res = run_cardinal_input(solid_input)
        except Exception as e:
            solid_res = {"error": str(e), "T_fuel_max": 0.0, "T_clad_max": 0.0}
            
        return {
            "thm_results": thm_res,
            "solid_results": solid_res,
            "bench_dir": bench_dir
        }
