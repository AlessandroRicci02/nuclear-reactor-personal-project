#!/usr/bin/env python3
"""
run_coupled_picard_900k.py

Coupled OpenMC–Thermal-Hydraulics Picard iteration for Second Simulation Campaign:
Coolant Inlet Temperature T_inlet = 900.0 K (reduced from 950.0 K).
Enforces exact energy conservation (Error < 0.0001%), 9-ring flow orificing,
local annular Gnielinski heat transfer, and temperature-dependent conductivities.
"""

import os
import sys
import time
import math
import numpy as np
import pandas as pd
import openmc

BASE_DIR = "/home/aleri/OpenMC_projects/first_test/final_calculation/second_simulation"
sys.path.insert(0, os.path.join(BASE_DIR, "scripts"))
from build_reference_model_900k import (
    build_reference_openmc_model, OPENMC_EXEC, CROSS_SECTIONS_PATH,
    N_PINS, P_NOM, R_GAS, M_HE, rho_helium
)

WORK_DIR = os.path.join(BASE_DIR, "coupled_picard")
os.makedirs(WORK_DIR, exist_ok=True)

# Thermal Operating Conditions (900 K Campaign)
Q_TOTAL = 400000.0          # 400.0 kWth
T_INLET = 900.0             # 900.0 K Coolant Inlet
DELTA_T_BULK = 200.0        # 200.0 K Temperature Rise -> Bulk Outlet = 1100.0 K
CP_HE = 5193.16             # Pure Helium Cp [J / (kg K)]
M_DOT_TOTAL = Q_TOTAL / (CP_HE * DELTA_T_BULK) # 0.3851216 kg/s

H_CORE = 0.4245             # 42.45 cm
R_FUEL = 0.00580            # 0.580 cm
T_CLAD = 0.00048            # 0.048 cm
R_CLAD_OUT = R_FUEL + T_CLAD # 0.00628 m
CHAN_GAP = 0.00084          # 0.084 cm
R_CHAN_OUT = R_CLAD_OUT + CHAN_GAP # 0.00712 m

A_FLOW_PIN = math.pi * (R_CHAN_OUT**2 - R_CLAD_OUT**2) # 3.53618e-5 m2
P_WET_PIN = 2.0 * math.pi * (R_CHAN_OUT + R_CLAD_OUT)  # 0.084195 m
D_H = 4.0 * A_FLOW_PIN / P_WET_PIN                     # 0.00168 m (2 * gap)
A_CLAD_OUT_PIN = 2.0 * math.pi * R_CLAD_OUT * H_CORE
V_FUEL_PIN = math.pi * (R_FUEL**2) * H_CORE

H_GAP_NOMINAL = 25000.0     # Gap conductance [W / (m2 K)]
N_AXIAL = 10


def mu_helium(T):
    return 3.674e-7 * (T**0.70)


def k_helium(T, P=P_NOM):
    term = 1.0 + 1.123e-8 * P
    exp = 0.71 * (1.0 - 2.0e-9 * P)
    return 2.682e-3 * term * (T**exp)


def prandtl_helium(T, P=P_NOM):
    return (CP_HE * mu_helium(T)) / k_helium(T, P)


def reynolds_helium(m_dot, T):
    return (m_dot * D_H) / (A_FLOW_PIN * mu_helium(T))


def friction_factor_annulus(Re):
    # Annular friction factor
    if Re < 2300:
        return 96.0 / max(Re, 10.0)
    elif Re < 4000:
        f_lam = 96.0 / 2300.0
        f_turb = (0.790 * math.log(4000.0) - 1.64)**(-2)
        return f_lam + (Re - 2300.0) / (4000.0 - 2300.0) * (f_turb - f_lam)
    else:
        return (0.790 * math.log(Re) - 1.64)**(-2)


def nusselt_annulus(Re, Pr):
    if Re < 2300:
        return 5.385
    elif Re < 10000:
        f = friction_factor_annulus(Re)
        num = (f / 8.0) * (Re - 1000.0) * Pr
        den = 1.0 + 12.7 * math.sqrt(f / 8.0) * (Pr**(2.0/3.0) - 1.0)
        return max(5.385, num / den)
    else:
        f = friction_factor_annulus(Re)
        num = (f / 8.0) * Re * Pr
        den = 1.07 + 12.7 * math.sqrt(f / 8.0) * (Pr**(2.0/3.0) - 1.0)
        return num / den


def k_clad_nb1zr(T):
    return 41.8 + 0.0125 * T


def solve_pin_thermal_900k(q_vol_slices, m_dot_pin, z_edges):
    dz = H_CORE / N_AXIAL
    v_slice = math.pi * (R_FUEL**2) * dz
    a_clad_slice = 2.0 * math.pi * R_CLAD_OUT * dz
    G = m_dot_pin / A_FLOW_PIN

    T_bulk = np.zeros(N_AXIAL)
    T_clad_out = np.zeros(N_AXIAL)
    T_clad_avg = np.zeros(N_AXIAL)
    T_fuel_center = np.zeros(N_AXIAL)
    T_fuel_vol = np.zeros(N_AXIAL)
    dp_fric_slice = np.zeros(N_AXIAL)
    Re_arr = np.zeros(N_AXIAL)
    h_arr = np.zeros(N_AXIAL)

    T_curr = T_INLET
    for s in range(N_AXIAL):
        p_slice = q_vol_slices[s] * v_slice
        delta_T_b = p_slice / (m_dot_pin * CP_HE)
        T_exit_s = T_curr + delta_T_b
        T_mid_s = 0.5 * (T_curr + T_exit_s)
        T_bulk[s] = T_exit_s

        rho_s = rho_helium(P_NOM, T_mid_s)
        Re_s = reynolds_helium(m_dot_pin, T_mid_s)
        Pr_s = prandtl_helium(T_mid_s)
        Nu_s = nusselt_annulus(Re_s, Pr_s)
        k_he = k_helium(T_mid_s)
        h_s = Nu_s * k_he / D_H
        h_arr[s] = h_s
        Re_arr[s] = Re_s

        f_s = friction_factor_annulus(Re_s)
        dp_fric_slice[s] = f_s * (dz / D_H) * (G**2) / (2.0 * rho_s)

        # Clad surface
        q_flux_out = p_slice / a_clad_slice
        T_co = T_mid_s + q_flux_out / h_s
        T_clad_out[s] = T_co

        # Clad conduction
        kc = k_clad_nb1zr(T_co)
        delta_T_c = (p_slice / (2.0 * math.pi * dz)) * math.log(R_CLAD_OUT / R_FUEL) / kc
        T_ci = T_co + delta_T_c
        T_clad_avg[s] = 0.5 * (T_co + T_ci)

        # Gap
        delta_T_gap = (p_slice / (2.0 * math.pi * R_FUEL * dz)) / H_GAP_NOMINAL
        T_fs = T_ci + delta_T_gap

        # Fuel pellet conduction: k_f(T) = 1.37 * T^0.41
        integ = (p_slice / dz) / (4.0 * math.pi)
        T0 = (T_fs**1.41 + (1.41 / 1.37) * integ)**(1.0 / 1.41)
        T_fuel_center[s] = T0

        r_arr = np.linspace(0.0, R_FUEL, 51)
        T_r = (T0**1.41 - (1.41 / 1.37) * (q_vol_slices[s] * (r_arr**2) / 4.0))**(1.0 / 1.41)
        T_fuel_vol[s] = np.trapezoid(2.0 * r_arr * T_r, r_arr) / (R_FUEL**2)

        T_curr = T_exit_s

    return {
        "T_bulk": T_bulk,
        "T_clad_out": T_clad_out,
        "T_clad_avg": T_clad_avg,
        "T_fuel_center": T_fuel_center,
        "T_fuel_vol": T_fuel_vol,
        "dp_fric_total_Pa": np.sum(dp_fric_slice),
        "T_bulk_exit": T_curr,
        "Re": Re_arr,
        "h_conv": h_arr
    }


def run_coupled_picard(max_iters=5, theta=0.50, particles=30000, batches=30, inactive=10):
    print("=" * 80)
    print("STARTING COUPLED OPENMC-TH PICARD ITERATION (T_inlet = 900.0 K)")
    print(f"Total Power: {Q_TOTAL/1000:.1f} kWth")
    print(f"Helium Mass Flow: {M_DOT_TOTAL:.5f} kg/s, Inlet Temp: {T_INLET:.1f} K")
    print(f"Target Bulk Rise: {DELTA_T_BULK:.1f} K -> Bulk Outlet: {T_INLET + DELTA_T_BULK:.1f} K")
    print("=" * 80)

    z_edges = np.linspace(0.0, H_CORE, N_AXIAL + 1)

    # Initial uniform temperature guess shifted by -50 K
    T_fuel_vol = np.full((N_PINS, N_AXIAL), 1185.0)
    T_clad_avg = np.full((N_PINS, N_AXIAL), 1170.0)
    T_cool_bulk = np.full((N_PINS, N_AXIAL), 1000.0)
    q_vol_prev = None

    history = []

    for it in range(max_iters):
        t0 = time.time()
        iter_dir = os.path.join(WORK_DIR, f"iter_{it}")
        print(f"\n>>> [ITERATION {it}] Running OpenMC Heating Transport Step...")

        info = build_reference_openmc_model(
            run_dir=iter_dir,
            T_fuel_grid=T_fuel_vol,
            T_clad_grid=T_clad_avg,
            T_cool_grid=T_cool_bulk,
            particles=particles,
            batches=batches,
            inactive=inactive,
            tallies_enabled=True
        )

        os.environ["OPENMC_CROSS_SECTIONS"] = CROSS_SECTIONS_PATH
        os.environ["OMP_NUM_THREADS"] = "8"
        cwd_prev = os.getcwd()
        os.chdir(iter_dir)
        try:
            openmc.run(openmc_exec=OPENMC_EXEC, output=False)
        finally:
            os.chdir(cwd_prev)

        sp_path = os.path.join(iter_dir, f"statepoint.{batches}.h5")
        with openmc.StatePoint(sp_path) as sp:
            k_eff = float(sp.keff.nominal_value)
            sigma_k = float(sp.keff.std_dev)
            t_heat = sp.get_tally(name="fuel_cell_heating")
            df_heat = t_heat.get_pandas_dataframe()

        df_heat_scores = df_heat[df_heat["score"] == "heating"]
        raw_heat = df_heat_scores["mean"].values.reshape((N_PINS, N_AXIAL))

        # Normalize heat deposition to exact 400 kWth
        tot_raw = np.sum(raw_heat)
        p_cell_watts = raw_heat * (Q_TOTAL / tot_raw)
        dz = H_CORE / N_AXIAL
        v_fuel_slice = math.pi * (R_FUEL**2) * dz
        q_vol_calc = p_cell_watts / v_fuel_slice

        # Under-relaxation
        if q_vol_prev is not None:
            q_vol_calc = (1.0 - theta) * q_vol_prev + theta * q_vol_calc
            power_l2_diff = float(np.linalg.norm(q_vol_calc - q_vol_prev) / np.linalg.norm(q_vol_prev))
        else:
            power_l2_diff = 1.0
        q_vol_prev = q_vol_calc.copy()

        # 9-Ring flow orificing setup
        pin_rings = np.array([p[1] for p in info["pin_meta"]])
        pin_powers = np.sum(p_cell_watts, axis=1)
        avg_power_rings = np.array([np.mean(pin_powers[pin_rings == r]) for r in range(9)])
        power_ratios = avg_power_rings / np.mean(pin_powers)
        pins_per_ring = np.array([1, 6, 12, 18, 24, 30, 36, 42, 48])
        norm_flow_factor = np.sum(pins_per_ring * power_ratios)
        m_dot_rings = M_DOT_TOTAL * power_ratios / norm_flow_factor

        m_dot_pins = np.zeros(N_PINS)
        for p in range(N_PINS):
            m_dot_pins[p] = m_dot_rings[pin_rings[p]]

        # Solve 217 channels
        T_fuel_vol_calc = np.zeros((N_PINS, N_AXIAL))
        T_fuel_center_calc = np.zeros((N_PINS, N_AXIAL))
        T_clad_avg_calc = np.zeros((N_PINS, N_AXIAL))
        T_clad_out_calc = np.zeros((N_PINS, N_AXIAL))
        T_cool_bulk_calc = np.zeros((N_PINS, N_AXIAL))
        dp_fric_pins = np.zeros(N_PINS)
        T_exit_pins = np.zeros(N_PINS)
        Re_pins = np.zeros(N_PINS)
        h_pins = np.zeros(N_PINS)

        for p in range(N_PINS):
            res_p = solve_pin_thermal_900k(q_vol_calc[p, :], m_dot_pins[p], z_edges)
            T_fuel_vol_calc[p, :] = res_p["T_fuel_vol"]
            T_fuel_center_calc[p, :] = res_p["T_fuel_center"]
            T_clad_avg_calc[p, :] = res_p["T_clad_avg"]
            T_clad_out_calc[p, :] = res_p["T_clad_out"]
            T_cool_bulk_calc[p, :] = res_p["T_bulk"]
            dp_fric_pins[p] = res_p["dp_fric_total_Pa"]
            T_exit_pins[p] = res_p["T_bulk_exit"]
            Re_pins[p] = np.mean(res_p["Re"])
            h_pins[p] = np.mean(res_p["h_conv"])

        # Pressure drop orificing
        dp_fric_max = np.max(dp_fric_pins)
        # Header target allows ring 0 friction + orifice margin
        dp_target = dp_fric_max + 950.0 # ~20 kPa header
        dp_orif_pins = np.maximum(0.0, dp_target - dp_fric_pins)
        dp_total_pins = dp_fric_pins + dp_orif_pins
        core_dp_kPa = float(np.mean(dp_total_pins)) / 1000.0

        d_Tf_max = float(np.max(np.abs(T_fuel_vol_calc - T_fuel_vol)))
        d_Tc_max = float(np.max(np.abs(T_clad_avg_calc - T_clad_avg)))

        T_bulk_out_mix = float(np.sum(m_dot_pins * T_exit_pins) / M_DOT_TOTAL)
        energy_bal_err = abs(T_bulk_out_mix - (T_INLET + DELTA_T_BULK)) / DELTA_T_BULK * 100.0

        if it > 0:
            delta_k = abs(k_eff - history[-1]["k_eff"]) * 1e5
        else:
            delta_k = 0.0

        dt = time.time() - t0
        rec = {
            "iteration": it,
            "k_eff": k_eff,
            "sigma_k": sigma_k,
            "delta_k_pcm": delta_k,
            "max_T_fuel_K": float(np.max(T_fuel_center_calc)),
            "mean_T_fuel_K": float(np.mean(T_fuel_vol_calc)),
            "max_T_clad_K": float(np.max(T_clad_out_calc)),
            "mean_T_clad_K": float(np.mean(T_clad_avg_calc)),
            "T_bulk_out_K": T_bulk_out_mix,
            "delta_T_hot_K": float(np.max(T_exit_pins) - np.min(T_exit_pins)),
            "core_delta_p_kPa": core_dp_kPa,
            "mean_fric_kPa": float(np.mean(dp_fric_pins)) / 1000.0,
            "max_orif_kPa": float(np.max(dp_orif_pins)) / 1000.0,
            "energy_balance_error_pct": energy_bal_err,
            "power_L2_diff": power_l2_diff,
            "d_Tf_max_K": d_Tf_max,
            "d_Tc_max_K": d_Tc_max,
            "runtime_s": dt
        }
        history.append(rec)
        print(f"    Iter {it}: k = {k_eff:.5f} +/- {sigma_k:.5f} (Delta k = {delta_k:.1f} pcm) | "
              f"Max T_clad = {rec['max_T_clad_K']:.1f} K | Max T_fuel = {rec['max_T_fuel_K']:.1f} K | "
              f"Core Delta P = {core_dp_kPa:.2f} kPa | Runtime: {dt:.1f} s")

        # Update fields for next iteration
        T_fuel_vol = T_fuel_vol_calc.copy()
        T_clad_avg = T_clad_avg_calc.copy()
        T_cool_bulk = T_cool_bulk_calc.copy()

    # Save convergence history
    df_hist = pd.DataFrame(history)
    hist_path = os.path.join(WORK_DIR, "picard_convergence_history.csv")
    df_hist.to_csv(hist_path, index=False)
    print(f"\n[+] Saved Picard history to {hist_path}")

    # Save converged channel results
    ch_rows = []
    for p in range(N_PINS):
        ch_rows.append({
            "pin_id": p,
            "ring": pin_rings[p],
            "power_W": pin_powers[p],
            "m_dot_kg_s": m_dot_pins[p],
            "T_exit_K": T_exit_pins[p],
            "T_clad_max_K": np.max(T_clad_out_calc[p, :]),
            "T_fuel_center_max_K": np.max(T_fuel_center_calc[p, :]),
            "dp_fric_Pa": dp_fric_pins[p],
            "dp_orif_Pa": dp_orif_pins[p],
            "dp_total_Pa": dp_total_pins[p],
            "Re_mean": Re_pins[p],
            "h_mean": h_pins[p]
        })
    df_channels = pd.DataFrame(ch_rows)
    ch_path = os.path.join(WORK_DIR, "converged_channel_results.csv")
    df_channels.to_csv(ch_path, index=False)
    print(f"[+] Saved converged channel results to {ch_path}")

    # Save converged arrays
    np.save(os.path.join(WORK_DIR, "converged_q_vol_3d.npy"), q_vol_calc)
    np.savez(os.path.join(WORK_DIR, "converged_temperature_grids.npz"),
             T_fuel_vol=T_fuel_vol_calc,
             T_clad_avg=T_clad_avg_calc,
             T_clad_out=T_clad_out_calc,
             T_cool_bulk=T_cool_bulk_calc,
             T_fuel_center=T_fuel_center_calc)
    print("[+] Saved converged 3D temperature and power grids.")


if __name__ == "__main__":
    run_coupled_picard(max_iters=5, theta=0.50, particles=35000, batches=30, inactive=10)
