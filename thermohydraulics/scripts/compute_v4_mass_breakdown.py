"""
SCRIPT: COMPUTE V4 TRANSPARENT COMPONENT MASS BREAKDOWN
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/scripts/compute_v4_mass_breakdown.py
Outputs: /home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v4/v4_mass_breakdown.csv
"""

import math
import numpy as np
import pandas as pd
import os

def compute_mass_breakdown():
    # 1. Geometry parameters [cm]
    pitch = 2.270
    rf = 0.580
    clad_th = 0.050
    rc = rf + clad_th # 0.630 cm
    chan_gap = 0.075
    rch = rc + chan_gap # 0.705 cm
    r_add = 0.150 # aux coolant holes
    
    core_h = 42.45
    refl_ax_th = 9.275
    total_h = core_h + 2.0 * refl_ax_th # 61.00 cm
    core_r = 19.45
    refl_rad_th = 9.80
    refl_rad_outer_r = core_r + refl_rad_th # 29.25 cm
    
    ax_duct_r = 0.860
    drum_count = 10
    drum_r = 3.95
    abs_th = 1.50
    abs_angle = 148.0
    
    # Material densities [g/cm3]
    rho_un = 13.50
    rho_clad = 8.57     # Nb-1Zr
    rho_yh2 = 4.30      # Canned YH2.0
    rho_beo = 2.90      # Sintered BeO
    rho_b4c = 2.52      # B4C absorber
    rho_he = 0.00115    # Helium at 2.5 MPa, ~1050 K (1.15 kg/m3 = 0.00115 g/cm3)
    
    enr_inner = 19.30   # wt%
    enr_outer = 19.75   # wt%
    n_pins = 217
    n_in = 19
    n_out = 198
    
    # 2. Volumes and Masses
    # Fuel UN
    v_f_pin = math.pi * rf**2 * core_h # cm3 per pin
    v_f_tot = n_pins * v_f_pin
    m_f_tot = v_f_tot * rho_un * 1e-3 # kg
    
    # Heavy metal U and U-235
    # UN molecular weight: 238 + 14 = 252 (or 235.04*enr + 238.05*(1-enr) + 14.0)
    w_u_fraction = 238.0 / 252.0
    m_u_tot = m_f_tot * w_u_fraction
    m_u_in = n_in * v_f_pin * rho_un * 1e-3 * w_u_fraction
    m_u_out = n_out * v_f_pin * rho_un * 1e-3 * w_u_fraction
    m_u235_tot = m_u_in * (enr_inner / 100.0) + m_u_out * (enr_outer / 100.0)
    
    # Cladding Nb-1Zr
    v_c_pin = math.pi * (rc**2 - rf**2) * core_h
    v_clad_tot = n_pins * v_c_pin
    m_clad_tot = v_clad_tot * rho_clad * 1e-3
    
    # Helium Coolant in active core channels
    v_ann = math.pi * (rch**2 - rc**2) * core_h
    v_aux_core = 6 * math.pi * r_add**2 * core_h
    v_cool_core = n_pins * (v_ann + v_aux_core)
    
    # Axial reflector coolant channels
    v_duct_ax = math.pi * ax_duct_r**2 * (2.0 * refl_ax_th)
    v_aux_ax = 6 * math.pi * r_add**2 * (2.0 * refl_ax_th)
    v_cool_ax = n_pins * (v_duct_ax + v_aux_ax)
    v_cool_tot = v_cool_core + v_cool_ax
    m_he_tot = v_cool_tot * rho_he * 1e-3 # kg
    
    # Moderator YH2.0
    cell_area = (math.sqrt(3.0) / 2.0) * pitch**2
    v_cell_core = n_pins * cell_area * core_h
    v_mod_tot = v_cell_core - (n_pins * v_f_pin + n_pins * v_c_pin + v_cool_core)
    m_mod_tot = v_mod_tot * rho_yh2 * 1e-3
    
    # BeO Axial Reflector
    v_ax_gross = math.pi * refl_rad_outer_r**2 * (2.0 * refl_ax_th)
    v_beo_ax = v_ax_gross - v_cool_ax
    m_beo_ax = v_beo_ax * rho_beo * 1e-3
    
    # BeO Radial Reflector & Control Drums
    v_cyl_active = math.pi * refl_rad_outer_r**2 * core_h
    v_refl_rad_gross = v_cyl_active - (math.pi * core_r**2 * core_h)
    v_drums = drum_count * math.pi * drum_r**2 * core_h
    v_refl_rad_stat = v_refl_rad_gross - v_drums
    
    # B4C Absorber Sectors (148 deg, 1.5 cm thick)
    v_abs_single = (abs_angle / 360.0) * math.pi * (drum_r**2 - (drum_r - abs_th)**2) * core_h
    v_abs_tot = drum_count * v_abs_single
    m_abs_tot = v_abs_tot * rho_b4c * 1e-3
    
    # BeO in Drums
    v_beo_drums = v_drums - v_abs_tot
    m_beo_drums = v_beo_drums * rho_beo * 1e-3
    m_beo_rad = (v_refl_rad_stat + v_beo_drums) * rho_beo * 1e-3
    m_beo_tot = m_beo_ax + m_beo_rad
    
    # Core Subassembly Total
    m_subassembly = m_f_tot + m_clad_tot + m_mod_tot + m_beo_tot + m_abs_tot + m_he_tot
    
    # Enclosing Vessel & Drum Drive estimate (5-10% allowance)
    # Reactor Vessel Nb-1Zr (e.g. 3 mm thick shell around outer reflector)
    v_vessel = math.pi * ((refl_rad_outer_r + 0.3)**2 - refl_rad_outer_r**2) * total_h + 2.0 * math.pi * (refl_rad_outer_r + 0.3)**2 * 0.3
    m_vessel = v_vessel * rho_clad * 1e-3
    
    # Drum actuation motors & bearings (10 x 2.5 kg = 25.0 kg)
    m_drum_drives = 25.0
    
    # Core support plates & tie rods Nb-1Zr (~15.0 kg)
    m_core_internals = 15.0
    
    m_reactor_module = m_subassembly + m_vessel + m_drum_drives + m_core_internals
    
    # NASA JIMO System Level Allocations for comparison
    # Shadow shield (W/LiH), Brayton turbomachinery, Heat pipe radiators, Power conditioning
    m_shield_est = 750.0      # Scaled for 400 kWth / 100 kWe NEP (vs 2000-3000 kg for 1.5 MWth JIMO)
    m_brayton_est = 350.0     # 2 x 50 kWe TAC units + recuperators
    m_radiator_est = 650.0    # C-C heat pipe radiator ~100 m2 @ ~6.5 kg/m2
    m_power_proc_est = 400.0  # PMAD + cabling
    m_total_nep_est = m_reactor_module + m_shield_est + m_brayton_est + m_radiator_est + m_power_proc_est

    breakdown = [
        {"Component": "UN Fuel (Pellets)", "Material": "UN (19.3 - 19.75% U-235)", "Volume_L": round(v_f_tot * 1e-3, 3), "Density_g_cm3": rho_un, "Mass_kg": round(m_f_tot, 2), "Subassembly_Fraction_%": round(m_f_tot / m_subassembly * 100, 2), "Category": "Active Core"},
        {"Component": "Fissile U-235 Inventory", "Material": "U-235", "Volume_L": "-", "Density_g_cm3": "-", "Mass_kg": round(m_u235_tot, 2), "Subassembly_Fraction_%": round(m_u235_tot / m_subassembly * 100, 2), "Category": "Fissile Subtotal"},
        {"Component": "Total Uranium Metal", "Material": "U (elemental)", "Volume_L": "-", "Density_g_cm3": "-", "Mass_kg": round(m_u_tot, 2), "Subassembly_Fraction_%": round(m_u_tot / m_subassembly * 100, 2), "Category": "Fissile Subtotal"},
        {"Component": "Fuel Cladding", "Material": "Nb-1Zr Alloy", "Volume_L": round(v_clad_tot * 1e-3, 3), "Density_g_cm3": rho_clad, "Mass_kg": round(m_clad_tot, 2), "Subassembly_Fraction_%": round(m_clad_tot / m_subassembly * 100, 2), "Category": "Active Core"},
        {"Component": "Moderator Matrix", "Material": "Canned YH2.0", "Volume_L": round(v_mod_tot * 1e-3, 3), "Density_g_cm3": rho_yh2, "Mass_kg": round(m_mod_tot, 2), "Subassembly_Fraction_%": round(m_mod_tot / m_subassembly * 100, 2), "Category": "Active Core"},
        {"Component": "Helium Coolant (In-Core)", "Material": "4He (2.5 MPa)", "Volume_L": round(v_cool_tot * 1e-3, 3), "Density_g_cm3": round(rho_he, 5), "Mass_kg": round(m_he_tot, 3), "Subassembly_Fraction_%": round(m_he_tot / m_subassembly * 100, 3), "Category": "Active Core"},
        {"Component": "BeO Axial Reflectors", "Material": "Sintered BeO", "Volume_L": round(v_beo_ax * 1e-3, 3), "Density_g_cm3": rho_beo, "Mass_kg": round(m_beo_ax, 2), "Subassembly_Fraction_%": round(m_beo_ax / m_subassembly * 100, 2), "Category": "Reflector System"},
        {"Component": "BeO Radial Reflector (Stationary)", "Material": "Sintered BeO", "Volume_L": round(v_refl_rad_stat * 1e-3, 3), "Density_g_cm3": rho_beo, "Mass_kg": round(v_refl_rad_stat * rho_beo * 1e-3, 2), "Subassembly_Fraction_%": round(v_refl_rad_stat * rho_beo * 1e-3 / m_subassembly * 100, 2), "Category": "Reflector System"},
        {"Component": "BeO Drum Segments", "Material": "Sintered BeO", "Volume_L": round(v_beo_drums * 1e-3, 3), "Density_g_cm3": rho_beo, "Mass_kg": round(m_beo_drums, 2), "Subassembly_Fraction_%": round(m_beo_drums / m_subassembly * 100, 2), "Category": "Control Drums"},
        {"Component": "B4C Absorber Sectors (10 drums)", "Material": "B4C (90% B-10)", "Volume_L": round(v_abs_tot * 1e-3, 3), "Density_g_cm3": rho_b4c, "Mass_kg": round(m_abs_tot, 2), "Subassembly_Fraction_%": round(m_abs_tot / m_subassembly * 100, 2), "Category": "Control Drums"},
        {"Component": "TOTAL SUBASSEMBLY (Opt-17 Verified)", "Material": "Modeled Components", "Volume_L": round(math.pi * refl_rad_outer_r**2 * total_h * 1e-3, 2), "Density_g_cm3": "-", "Mass_kg": round(m_subassembly, 2), "Subassembly_Fraction_%": 100.0, "Category": "Opt-17 Core Subassembly"},
        {"Component": "Reactor Vessel Shell (3 mm Nb-1Zr)", "Material": "Nb-1Zr Alloy", "Volume_L": round(v_vessel * 1e-3, 3), "Density_g_cm3": rho_clad, "Mass_kg": round(m_vessel, 2), "Subassembly_Fraction_%": "-", "Category": "Reactor Structural Allowance"},
        {"Component": "Core Support Plates & Grid", "Material": "Nb-1Zr Alloy", "Volume_L": "-", "Density_g_cm3": "-", "Mass_kg": round(m_core_internals, 2), "Subassembly_Fraction_%": "-", "Category": "Reactor Structural Allowance"},
        {"Component": "Drum Actuation Drives & Bearings", "Material": "Stepper Motors / MoS2", "Volume_L": "-", "Density_g_cm3": "-", "Mass_kg": round(m_drum_drives, 2), "Subassembly_Fraction_%": "-", "Category": "Reactor Control Subsystem"},
        {"Component": "TOTAL REACTOR MODULE ESTIMATE", "Material": "Core + Structure + Drives", "Volume_L": "-", "Density_g_cm3": "-", "Mass_kg": round(m_reactor_module, 2), "Subassembly_Fraction_%": "-", "Category": "Complete Reactor Module"},
        {"Component": "Shadow Radiation Shield (Estimated)", "Material": "W / LiH / B4C", "Volume_L": "-", "Density_g_cm3": "-", "Mass_kg": round(m_shield_est, 2), "Subassembly_Fraction_%": "-", "Category": "NEP Balance-of-Plant"},
        {"Component": "Brayton Turbomachinery & Recuperator", "Material": "Superalloys / TiAl", "Volume_L": "-", "Density_g_cm3": "-", "Mass_kg": round(m_brayton_est, 2), "Subassembly_Fraction_%": "-", "Category": "NEP Balance-of-Plant"},
        {"Component": "Heat Pipe Radiator Subsystem", "Material": "C-C / Titanium-Water/K", "Volume_L": "-", "Density_g_cm3": "-", "Mass_kg": round(m_radiator_est, 2), "Subassembly_Fraction_%": "-", "Category": "NEP Balance-of-Plant"},
        {"Component": "Power Management & Distribution (PMAD)", "Material": "SiC Electronics / Cu", "Volume_L": "-", "Density_g_cm3": "-", "Mass_kg": round(m_power_proc_est, 2), "Subassembly_Fraction_%": "-", "Category": "NEP Balance-of-Plant"},
        {"Component": "TOTAL 100 kWe NEP SYSTEM ESTIMATE", "Material": "Full Spacecraft Power System", "Volume_L": "-", "Density_g_cm3": "-", "Mass_kg": round(m_total_nep_est, 2), "Subassembly_Fraction_%": "-", "Category": "Full Flight System"}
    ]
    
    df = pd.DataFrame(breakdown)
    out_csv = "/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v4/v4_mass_breakdown.csv"
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"Generated {out_csv}")
    print(f"Subassembly Mass: {m_subassembly:.2f} kg (Target < 650 kg, Margin = {650.0 - m_subassembly:+.2f} kg)")
    print(f"Total Reactor Module: {m_reactor_module:.2f} kg")
    print(f"Total NEP System Estimate: {m_total_nep_est:.2f} kg")

if __name__ == "__main__":
    compute_mass_breakdown()
