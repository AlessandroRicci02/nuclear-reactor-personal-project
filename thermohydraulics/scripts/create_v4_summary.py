"""
SCRIPT: ASSEMBLE V4 MACHINE-READABLE SUMMARY JSON
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/scripts/create_v4_summary.py
Outputs: /home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v4/v4_summary.json
"""

import os
import json
import pandas as pd

RESULTS_V4_DIR = "/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v4"

def assemble_summary():
    # Load all generated CSVs
    df_neut = pd.read_csv(os.path.join(RESULTS_V4_DIR, "v4_neutronics.csv"))
    df_fb = pd.read_csv(os.path.join(RESULTS_V4_DIR, "v4_feedback_coefficients.csv"))
    df_th = pd.read_csv(os.path.join(RESULTS_V4_DIR, "v4_thermal_hydraulics.csv"))
    df_mass = pd.read_csv(os.path.join(RESULTS_V4_DIR, "v4_mass_breakdown.csv"))
    df_sens = pd.read_csv(os.path.join(RESULTS_V4_DIR, "v4_sensitivity_studies.csv"))
    df_comp = pd.read_csv(os.path.join(RESULTS_V4_DIR, "v4_reference_comparison.csv"))
    df_iss = pd.read_csv(os.path.join(RESULTS_V4_DIR, "v4_issue_register.csv"))

    summary = {
        "campaign": "Opt-17 Cand 17-Clean V4 Multiphysics Coupling & JIMO Assessment",
        "date": "2026-10-03",
        "architecture": {
            "name": "Opt-17 Cand 17-Clean",
            "lattice_type": "Hexagonal 9-ring",
            "fuel_pin_count": 217,
            "fuel_material": "Sintered UN (19.30% inner / 19.75% outer U-235)",
            "cladding_material": "Nb-1Zr refractory alloy (0.50 mm wall)",
            "moderator_material": "Canned YH2.0 metal hydride",
            "reflector_material": "Sintered BeO (radial 9.80 cm, axial 9.275 cm)",
            "control_system": "10 rotating BeO drums with 148 deg B4C sectors (90% B-10)",
            "thermal_power_kWth": 400.0,
            "electrical_power_kWe": 100.0,
            "brayton_efficiency_nominal": 0.25
        },
        "coupling_status": {
            "picard_iterations": 8,
            "convergence_status": "CONVERGED",
            "delta_k_ratio": 0.592,
            "delta_T_fuel_K": 5.47,
            "power_RMS_change": 0.0048,
            "global_energy_residual_percent": 0.0000,
            "cardinal_fe_peak_fuel_diff_K": 0.27,
            "cardinal_fe_peak_clad_diff_K": 6.17
        },
        "neutronics_results": {
            "k_eff_operating": float(df_neut.loc[df_neut["Parameter"] == "Coupled Operating k_eff", "Value"].values[0]),
            "sigma_k": float(df_neut.loc[df_neut["Parameter"] == "Coupled Operating k_eff", "Uncertainty_1sigma"].values[0]),
            "excess_reactivity_pcm": float(df_neut.loc[df_neut["Parameter"] == "Excess Reactivity", "Value"].values[0]),
            "depletion_lifetime_years": float(df_neut.loc[df_neut["Parameter"] == "Depletion Operational Lifetime", "Value"].values[0]),
            "radial_peaking_Fr": float(df_neut.loc[df_neut["Parameter"] == "Radial Peaking Factor (F_r)", "Value"].values[0]),
            "axial_peaking_Fz": float(df_neut.loc[df_neut["Parameter"] == "Axial Peaking Factor (F_z)", "Value"].values[0]),
            "total_peaking_Fq": float(df_neut.loc[df_neut["Parameter"] == "Total Peaking Factor (F_q)", "Value"].values[0]),
            "peak_volumetric_power_MW_m3": float(df_neut.loc[df_neut["Parameter"] == "Peak Volumetric Power Density", "Value"].values[0]),
            "fissile_u235_kg": float(df_neut.loc[df_neut["Parameter"] == "Core Fissile U-235 Inventory", "Value"].values[0])
        },
        "feedback_coefficients": {
            "doppler_alpha_D_pcm_K": float(df_fb.loc[df_fb["Coefficient"].str.startswith("Fuel Doppler"), "Nominal_Value"].values[0]),
            "doppler_uncertainty_pcm_K": float(df_fb.loc[df_fb["Coefficient"].str.startswith("Fuel Doppler"), "Uncertainty_1sigma"].values[0]),
            "moderator_alpha_mod_pcm_K": float(df_fb.loc[df_fb["Coefficient"].str.startswith("Moderator Temperature"), "Nominal_Value"].values[0]),
            "coolant_void_reactivity_pcm": float(df_fb.loc[df_fb["Coefficient"].str.startswith("Total Coolant Void"), "Nominal_Value"].values[0]),
            "control_drum_worth_pcm": float(df_fb.loc[df_fb["Coefficient"].str.startswith("Total Control Drum"), "Nominal_Value"].values[0]),
            "all_in_sdm_pcm": float(df_fb.loc[df_fb["Coefficient"].str.startswith("All-In Shutdown"), "Nominal_Value"].values[0]),
            "stuck_drum_sdm_pcm": float(df_fb.loc[df_fb["Coefficient"].str.startswith("Single Stuck-Drum"), "Nominal_Value"].values[0])
        },
        "thermal_hydraulics_results": {
            "coolant_fluid": "Pure 4He (2.5 MPa)",
            "total_mass_flow_kg_s": float(df_th.loc[df_th["Parameter"] == "Total Coolant Mass Flow", "Value"].values[0]),
            "T_inlet_K": float(df_th.loc[df_th["Parameter"] == "Coolant Inlet Temperature", "Value"].values[0]),
            "T_outlet_bulk_K": float(df_th.loc[df_th["Parameter"] == "Mixed Coolant Outlet Temperature", "Value"].values[0]),
            "T_cool_exit_max_K": float(df_th.loc[df_th["Parameter"] == "Hot Channel Coolant Exit Temp", "Value"].values[0]),
            "T_fuel_center_max_K": float(df_th.loc[df_th["Parameter"] == "Peak Fuel Centerline Temperature", "Value"].values[0]),
            "T_clad_out_max_K": float(df_th.loc[df_th["Parameter"] == "Peak Cladding Outer Temperature", "Value"].values[0]),
            "core_pressure_drop_kPa": float(df_th.loc[df_th["Parameter"] == "Core Total Pressure Drop", "Value"].values[0]),
            "reynolds_range": str(df_th.loc[df_th["Parameter"] == "Reynolds Number Range", "Value"].values[0]),
            "nusselt_mean": float(df_th.loc[df_th["Parameter"] == "Average Nusselt Number", "Value"].values[0]),
            "htc_mean_W_m2_K": float(df_th.loc[df_th["Parameter"] == "Convective Heat Transfer Coeff", "Value"].values[0]),
            "hydraulic_diameter_mm": float(df_th.loc[df_th["Parameter"] == "Channel Hydraulic Diameter", "Value"].values[0]),
            "gap_width_mm": float(df_th.loc[df_th["Parameter"] == "Annular Coolant Gap Width", "Value"].values[0])
        },
        "mass_accounting": {
            "subassembly_mass_kg": float(df_mass.loc[df_mass["Component"].str.startswith("TOTAL SUBASSEMBLY"), "Mass_kg"].values[0]),
            "subassembly_target_limit_kg": 650.0,
            "subassembly_mass_margin_kg": float(650.0 - df_mass.loc[df_mass["Component"].str.startswith("TOTAL SUBASSEMBLY"), "Mass_kg"].values[0]),
            "total_reactor_module_est_kg": float(df_mass.loc[df_mass["Component"].str.startswith("TOTAL REACTOR MODULE"), "Mass_kg"].values[0]),
            "total_nep_system_est_kg": float(df_mass.loc[df_mass["Component"].str.startswith("TOTAL 100 kWe NEP"), "Mass_kg"].values[0])
        },
        "reference_comparison_status_counts": {
            "within_target": int((df_comp["Status"] == "WITHIN TARGET").sum()),
            "outside_target": int((df_comp["Status"] == "OUTSIDE TARGET").sum()),
            "indeterminate": int((df_comp["Status"] == "INDETERMINATE").sum()),
            "not_assessed": int((df_comp["Status"] == "NOT ASSESSED").sum())
        },
        "ranked_issues_count": len(df_iss),
        "deliverables": {
            "report": "V4_FINAL_REPORT.md",
            "summary_json": "v4_summary.json",
            "iterations_csv": "v4_iterations.csv",
            "neutronics_csv": "v4_neutronics.csv",
            "feedback_coefficients_csv": "v4_feedback_coefficients.csv",
            "thermal_hydraulics_csv": "v4_thermal_hydraulics.csv",
            "mass_breakdown_csv": "v4_mass_breakdown.csv",
            "reference_comparison_csv": "v4_reference_comparison.csv",
            "sensitivity_studies_csv": "v4_sensitivity_studies.csv",
            "issue_register_csv": "v4_issue_register.csv",
            "validation_matrix_csv": "v4_validation_matrix.csv",
            "figures_count": 10,
            "spatial_vtu": "opt17_v4_coupled_field.vtk"
        }
    }

    out_json = os.path.join(RESULTS_V4_DIR, "v4_summary.json")
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Generated {out_json}")

if __name__ == "__main__":
    assemble_summary()
