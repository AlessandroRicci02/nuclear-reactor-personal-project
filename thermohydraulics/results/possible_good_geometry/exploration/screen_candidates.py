#!/usr/bin/env python3
"""
screen_candidates.py
Parametric screening of candidate NEP reactor geometries to achieve:
1. Lifespan > 15.0 years (keff >= 1.000 at 15.0 yr)
2. T_mod < 950 K
3. Mass < 600 kg
4. Delta P < 35 kPa
"""

import math
import os
import sys
import openmc
import numpy as np

OPENMC_EXEC = "/home/aleri/miniforge3/envs/openmc-env/bin/openmc"
CROSS_SECTIONS_PATH = "/home/aleri/wmp_data/cross_sections_wmp.xml"

# Densities
RHO_UN = 13.50
RHO_CLAD = 8.57
RHO_YH2 = 4.30
RHO_BEO = 2.90
RHO_B4C = 2.52
RHO_AL2O3 = 3.95

P_NOM = 2.50e6
R_GAS = 8314.462618
M_HE = 4.002602
N_PINS = 217
DRUM_COUNT = 10

def calc_geometry_properties(H_core, r_fuel, clad_th, chan_gap, al2o3_th, pitch, refl_ax_th, refl_rad_th, enr_inner, enr_outer):
    r_f = r_fuel
    r_c = r_f + clad_th
    r_ch = r_c + chan_gap
    r_al = r_ch + al2o3_th
    
    # Hydraulic properties
    A_flow_pin = math.pi * (r_ch**2 - r_c**2) * 1e-4 # m2
    P_wet_pin = 2.0 * math.pi * (r_ch + r_c) * 1e-2  # m
    D_h = 4.0 * A_flow_pin / P_wet_pin               # m
    A_flow_core = A_flow_pin * N_PINS               # m2
    
    # Mass calculations
    V_fuel = N_PINS * math.pi * r_f**2 * H_core
    M_un = V_fuel * RHO_UN / 1000.0
    M_hm = M_un * (238.0 / (238.0 + 14.0))
    
    V_clad = N_PINS * math.pi * (r_c**2 - r_f**2) * H_core
    M_clad = V_clad * RHO_CLAD / 1000.0
    
    V_al = N_PINS * math.pi * (r_al**2 - r_ch**2) * H_core
    M_al = V_al * RHO_AL2O3 / 1000.0
    
    A_hex = (math.sqrt(3)/2.0) * pitch**2
    A_chan = math.pi * r_al**2
    A_aux = 6 * math.pi * 0.150**2
    A_mod_cell = max(0.0, A_hex - A_chan - A_aux)
    V_mod = N_PINS * A_mod_cell * H_core
    M_yh2 = V_mod * RHO_YH2 / 1000.0
    
    R_core = 19.45 * (pitch / 2.270)
    R_ref = R_core + refl_rad_th
    
    r_drum = 3.95
    abs_th = 1.50
    abs_ang = 148.0
    A_drum = math.pi * r_drum**2
    A_abs = (abs_ang / 360.0) * math.pi * (r_drum**2 - (r_drum - abs_th)**2)
    A_beo_drum = A_drum - A_abs
    V_abs_tot = 10 * A_abs * H_core
    M_b4c = V_abs_tot * RHO_B4C / 1000.0
    V_beo_drums = 10 * A_beo_drum * H_core
    
    V_rad_annulus = math.pi * (R_ref**2 - R_core**2) * H_core - 10 * A_drum * H_core + V_beo_drums
    V_ax = 2.0 * (math.pi * R_ref**2 * refl_ax_th - N_PINS * math.pi * 0.860**2 * refl_ax_th - N_PINS * 6 * math.pi * 0.150**2 * refl_ax_th)
    V_beo_tot = V_rad_annulus + V_ax
    M_beo = V_beo_tot * RHO_BEO / 1000.0
    
    M_tot = M_un + M_clad + M_al + M_yh2 + M_beo + M_b4c
    
    # Rough Delta P estimate for 920 -> 1130 K, m_dot = 0.36678 kg/s
    m_dot_pin = 0.3667825 / N_PINS
    T_avg = 1025.0
    rho_avg = (P_NOM * M_HE) / (R_GAS * T_avg)
    mu_avg = 3.674e-7 * (T_avg**0.70)
    v_avg = m_dot_pin / (rho_avg * A_flow_pin)
    Re = (rho_avg * v_avg * D_h) / mu_avg
    f = (0.790 * math.log(Re) - 1.64)**(-2) if Re >= 4000 else 96.0 / Re
    dp_fric = f * (H_core * 1e-2 / D_h) * 0.5 * rho_avg * v_avg**2
    dp_form = 1.5 * 0.5 * rho_avg * v_avg**2
    dp_tot_kpa = (dp_fric + dp_form) / 1000.0
    
    # Heavy metal fissile
    # 19 inner pins, 198 outer pins
    u235_kg = M_hm * (19 * enr_inner + 198 * enr_outer) / (217 * 100.0)
    
    return {
        "M_un": M_un,
        "M_hm": M_hm,
        "u235_kg": u235_kg,
        "M_clad": M_clad,
        "M_al": M_al,
        "M_yh2": M_yh2,
        "M_beo": M_beo,
        "M_b4c": M_b4c,
        "M_tot": M_tot,
        "D_h_mm": D_h * 1000.0,
        "A_flow_m2": A_flow_core,
        "dp_kpa": dp_tot_kpa,
        "Re": Re
    }

candidates = {
    "Sim_4B_Reference": {
        "H_core": 42.45, "r_fuel": 0.580, "clad_th": 0.048, "chan_gap": 0.084, "al2o3_th": 0.010,
        "pitch": 2.270, "refl_ax_th": 9.275, "refl_rad_th": 9.80, "enr_inner": 19.30, "enr_outer": 19.75
    },
    "Cand_1_Axial_Extension": {
        # Core lengthened to 45.5 cm (+3.05 cm), slightly trimmed radial reflector to 9.30 cm
        "H_core": 45.50, "r_fuel": 0.580, "clad_th": 0.048, "chan_gap": 0.084, "al2o3_th": 0.010,
        "pitch": 2.270, "refl_ax_th": 9.000, "refl_rad_th": 9.40, "enr_inner": 19.50, "enr_outer": 19.75
    },
    "Cand_2_Dual_Lengthen_Pellet": {
        # Core 45.5 cm, pellet 0.585 cm (+0.05 mm), pitch 2.28 cm, radial refl 9.30 cm
        "H_core": 45.50, "r_fuel": 0.585, "clad_th": 0.048, "chan_gap": 0.084, "al2o3_th": 0.010,
        "pitch": 2.280, "refl_ax_th": 8.800, "refl_rad_th": 9.30, "enr_inner": 19.60, "enr_outer": 19.75
    },
    "Cand_3_Moderator_Enhanced": {
        # Core 46.0 cm, pellet 0.585 cm, pitch 2.29 cm (more YH2 thermalization), radial refl 9.20 cm
        "H_core": 46.00, "r_fuel": 0.585, "clad_th": 0.048, "chan_gap": 0.084, "al2o3_th": 0.010,
        "pitch": 2.290, "refl_ax_th": 8.500, "refl_rad_th": 9.20, "enr_inner": 19.60, "enr_outer": 19.75
    },
    "Cand_4_Maximum_Margin": {
        # Core 46.5 cm, pellet 0.588 cm, pitch 2.29 cm, radial refl 9.20 cm
        "H_core": 46.50, "r_fuel": 0.588, "clad_th": 0.048, "chan_gap": 0.084, "al2o3_th": 0.010,
        "pitch": 2.290, "refl_ax_th": 8.500, "refl_rad_th": 9.20, "enr_inner": 19.65, "enr_outer": 19.75
    }
}

print(f"{'Candidate':<30} | {'M_tot (kg)':<10} | {'M_hm (kg)':<10} | {'U-235 (kg)':<10} | {'DeltaP (kPa)':<12} | {'Margin to 600kg'}")
print("-" * 90)
for name, p in candidates.items():
    res = calc_geometry_properties(**p)
    margin = 600.0 - res["M_tot"]
    print(f"{name:<30} | {res['M_tot']:<10.2f} | {res['M_hm']:<10.2f} | {res['u235_kg']:<10.2f} | {res['dp_kpa']:<12.2f} | {margin:+.2f} kg")

