"""
V5 DESIGN SPACE & BASELINE DESIGN VECTOR DEFINITION
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/optimization/design_space.py

Formulates the constrained multi-objective optimization problem for Opt-17 NEP Reactor:
1. Exact reconstructed Opt-17 baseline design vector (Cand 17-Clean).
2. Parameterized design variables with physical bounds.
3. Hard physical/model constraints.
4. Sourced engineering constraints with acceptance limits.
5. Objective functions definitions.
"""

import os
import json
import math
import numpy as np

RESULTS_V5_DIR = "/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v5"
os.makedirs(RESULTS_V5_DIR, exist_ok=True)

# -----------------------------------------------------------------------------
# 1. EXACT RECONSTRUCTED BASELINE DESIGN VECTOR (OPT-17 CAND 17-CLEAN)
# -----------------------------------------------------------------------------
BASELINE_DESIGN_VECTOR = {
    "identifier": "OPT-17_CAND_17_CLEAN_BASELINE",
    "description": "Validated baseline configuration from V4 campaign",
    
    # Geometric parameters [cm unless specified]
    "fuel_radius_cm": 0.580,
    "cladding_thickness_cm": 0.050,
    "cladding_outer_radius_cm": 0.630,
    "coolant_gap_cm": 0.075,
    "channel_outer_radius_cm": 0.705,
    "pin_pitch_cm": 2.270,
    "number_of_pins": 217,
    "hex_rings": 9,
    "active_core_height_cm": 42.45,
    "active_core_radius_cm": 19.45,
    "axial_reflector_thickness_cm": 9.275,
    "total_assembly_height_cm": 61.00,
    "radial_reflector_thickness_cm": 9.80,
    "assembly_outer_radius_cm": 29.25,
    "number_of_axial_slices": 10,
    
    # Auxiliary flow holes
    "aux_holes_per_pin": 6,
    "aux_hole_radius_cm": 0.150,
    "aux_hole_radial_pos_cm": 0.910,
    "axial_duct_radius_cm": 0.860,
    
    # Control drum system
    "drum_count": 10,
    "drum_radius_cm": 3.95,
    "drum_absorber_thickness_cm": 1.50,
    "drum_absorber_angle_deg": 148.0,
    "drum_clearance_gap_cm": 0.180,
    "drum_center_radial_pos_cm": 23.58,
    "baseline_operating_drum_angle_deg": 0.0,
    
    # Material compositions & densities [g/cm3]
    "fuel_material": "Sintered UN",
    "fuel_density_g_cm3": 13.50,
    "inner_pins_count": 19,
    "inner_pins_enrichment_wt_u235": 19.30,
    "outer_pins_count": 198,
    "outer_pins_enrichment_wt_u235": 19.75,
    "cladding_material": "Nb-1Zr Alloy",
    "cladding_density_g_cm3": 8.57,
    "moderator_material": "Canned YH2.0",
    "moderator_density_g_cm3": 4.30,
    "reflector_material": "Sintered BeO",
    "reflector_density_g_cm3": 2.90,
    "absorber_material": "B4C (90% B-10 enriched)",
    "absorber_density_g_cm3": 2.52,
    
    # Thermal-hydraulic boundary conditions
    "thermal_power_rating_kWth": 400.0,
    "target_electrical_output_kWe": 100.0,
    "nominal_brayton_efficiency": 0.25,
    "coolant_fluid": "Pure 4He",
    "system_operating_pressure_MPa": 2.50,
    "coolant_inlet_temperature_K": 950.0,
    "target_bulk_outlet_temperature_K": 1150.0,
    "target_bulk_temperature_rise_K": 200.0,
    "total_mass_flow_kg_s": 0.385122,
    "pin_mass_flow_kg_s": 0.00177475,
    "gap_conductance_W_m2_K": 25000.0,
    
    # Neutronic operational metrics from V4
    "k_eff_operating": 1.01645,
    "sigma_k_operating": 0.00055,
    "excess_reactivity_pcm": 1618.3,
    "doppler_coefficient_pcm_K": -0.844,
    "moderator_coefficient_pcm_K": -0.811,
    "control_drum_worth_pcm": 4698.1,
    "shutdown_margin_all_in_pcm": 3079.8,
    "shutdown_margin_stuck_pcm": 2763.7,
    
    # Thermal-hydraulic metrics from V4
    "peak_fuel_centerline_temp_K": 1295.56,
    "peak_cladding_outer_temp_K": 1271.40,
    "hot_channel_coolant_exit_temp_K": 1207.80,
    "core_total_pressure_drop_kPa": 23.59,
    "radial_peaking_factor_Fr": 1.327,
    "total_3d_peaking_factor_Fq": 1.707,
    "subassembly_mass_kg": 539.66,
    "total_reactor_module_mass_kg": 622.74
}


# -----------------------------------------------------------------------------
# 2. DESIGN VARIABLES SPECIFICATION & BOUNDS
# -----------------------------------------------------------------------------
DESIGN_VARIABLES_BOUNDS = {
    # 1. Fuel pellet radius [cm] (Baseline: 0.580 cm)
    # Range: 0.500 cm (-13.8%) to 0.620 cm (+6.9%)
    "fuel_radius_cm": {
        "baseline": 0.580,
        "min": 0.520,
        "max": 0.610,
        "step": 0.010,
        "type": "continuous",
        "description": "Fuel pellet radius controlling fuel volume, power density, and clad heat flux."
    },
    
    # 2. Cladding wall thickness [cm] (Baseline: 0.050 cm = 0.50 mm)
    # Range: 0.035 cm (0.35 mm) to 0.065 cm (0.65 mm)
    "cladding_thickness_cm": {
        "baseline": 0.050,
        "min": 0.035,
        "max": 0.065,
        "step": 0.005,
        "type": "continuous",
        "description": "Nb-1Zr clad wall thickness controlling structural creep margin and thermal resistance."
    },
    
    # 3. Annular coolant channel gap [cm] (Baseline: 0.075 cm = 0.75 mm)
    # Range: 0.055 cm (0.55 mm) to 0.105 cm (1.05 mm)
    "coolant_gap_cm": {
        "baseline": 0.075,
        "min": 0.055,
        "max": 0.105,
        "step": 0.005,
        "type": "continuous",
        "description": "Annular coolant gap width determining flow area, velocity, Re, and friction pressure drop."
    },
    
    # 4. Pin lattice pitch [cm] (Baseline: 2.270 cm)
    # Range: 2.150 cm to 2.400 cm
    "pin_pitch_cm": {
        "baseline": 2.270,
        "min": 2.150,
        "max": 2.400,
        "step": 0.020,
        "type": "continuous",
        "description": "Hexagonal lattice pitch controlling moderator volume fraction and core diameter."
    },
    
    # 5. Active core height [cm] (Baseline: 42.45 cm)
    # Range: 38.00 cm to 46.00 cm
    "active_core_height_cm": {
        "baseline": 42.45,
        "min": 38.00,
        "max": 46.00,
        "step": 1.00,
        "type": "continuous",
        "description": "Active core fuel stack height affecting aspect ratio L/D and axial leakage."
    },
    
    # 6. Radial reflector thickness [cm] (Baseline: 9.80 cm)
    # Range: 7.50 cm to 12.00 cm
    "radial_reflector_thickness_cm": {
        "baseline": 9.80,
        "min": 7.50,
        "max": 12.00,
        "step": 0.50,
        "type": "continuous",
        "description": "BeO radial reflector thickness balancing neutron reflection and mass."
    },
    
    # 7. Total helium mass flow rate multiplier (Baseline: 1.00 = 0.38512 kg/s)
    # Range: 0.80 (-20%) to 1.30 (+30%)
    "mass_flow_multiplier": {
        "baseline": 1.00,
        "min": 0.80,
        "max": 1.30,
        "step": 0.05,
        "type": "continuous",
        "description": "Mass flow multiplier adjusting coolant temperature rise and heat transfer coefficient."
    },
    
    # 8. Coolant inlet temperature [K] (Baseline: 950.0 K)
    # Range: 900.0 K to 980.0 K
    "inlet_temperature_K": {
        "baseline": 950.0,
        "min": 900.0,
        "max": 980.0,
        "step": 10.0,
        "type": "continuous",
        "description": "Recuperator return temperature entering the reactor lower plenum."
    },
    
    # 9. Outer zone fuel enrichment [wt% U-235] (Baseline: 19.75 wt%)
    # Range: 18.50 wt% to 19.75 wt% (Strict non-proliferation ceiling: <= 19.75 wt%)
    "outer_enrichment_wt": {
        "baseline": 19.75,
        "min": 18.50,
        "max": 19.75,
        "step": 0.25,
        "type": "continuous",
        "description": "HALEU enrichment in outer 198 pins (Rings 3-8)."
    },
    
    # 10. Inner zone fuel enrichment [wt% U-235] (Baseline: 19.30 wt%)
    # Range: 17.50 wt% to 19.50 wt%
    "inner_enrichment_wt": {
        "baseline": 19.30,
        "min": 17.50,
        "max": 19.50,
        "step": 0.25,
        "type": "continuous",
        "description": "HALEU enrichment in central 19 pins (Rings 0-2) for radial power flattening."
    },
    
    # 11. Radial flow orificing strategy (Baseline: uniform)
    # Options: "uniform", "power_proportional", "two_zone", "optimal_ring"
    "orificing_mode": {
        "baseline": "uniform",
        "options": ["uniform", "power_proportional", "two_zone", "optimal_ring"],
        "type": "discrete",
        "description": "Inlet orifice loss distribution across rings 0 to 8 to equalize exit coolant temperatures."
    }
}


# -----------------------------------------------------------------------------
# 3. ENGINEERING CONSTRAINTS & ACCEPTANCE LIMITS
# -----------------------------------------------------------------------------
ENGINEERING_CONSTRAINTS = {
    "C01_max_cladding_temp": {
        "parameter": "Peak Outer Cladding Temperature (T_c,max)",
        "limit_value": 1350.0,
        "direction": "<=",
        "units": "K",
        "source": "Touloukian (1970), Nb-1Zr thermal creep acceleration limit (~1077 °C)",
        "status": "HARD_ENGINEERING_LIMIT",
        "description": "Operating above 1350 K induces rapid tertiary creep rupture over 10-20 years."
    },
    "C02_max_fuel_temp": {
        "parameter": "Peak Fuel Centerline Temperature (T_f,max)",
        "limit_value": 1800.0,
        "direction": "<=",
        "units": "K",
        "source": "Hayes et al. (1990), UN melting point 3035 K; 1800 K avoids significant fission gas release",
        "status": "HARD_ENGINEERING_LIMIT",
        "description": "Maintains >1200 K margin to UN melting and keeps FGR below 1%."
    },
    "C03_core_pressure_drop": {
        "parameter": "Core Total Coolant Pressure Drop (Delta P)",
        "limit_value": 35.0,
        "direction": "<=",
        "units": "kPa",
        "source": "Galloway et al. (2004) AIAA-2004-5977, CBC pressure ratio allocation (Delta P/P < 1.4%)",
        "status": "HARD_ENGINEERING_LIMIT",
        "description": "Prevents excessive compressor pumping work degradation."
    },
    "C04_subassembly_mass": {
        "parameter": "Core Subassembly Total Mass (m_subassembly)",
        "limit_value": 650.0,
        "direction": "<=",
        "units": "kg",
        "source": "Opt-17 Project Specification / Launch Payload Envelope",
        "status": "HARD_ENGINEERING_LIMIT",
        "description": "Strict launch vehicle mass envelope ceiling for deep space electric propulsion."
    },
    "C05_shutdown_margin_all_in": {
        "parameter": "All-In Shutdown Margin (SDM_all_in)",
        "limit_value": 2000.0,
        "direction": ">=",
        "units": "pcm",
        "source": "NASA / KAPL-4849 Aerospace Reactor Flight Safety Standard",
        "status": "HARD_ENGINEERING_LIMIT",
        "description": "Requires deep subcriticality during launch, ascent, and inadvertent water immersion."
    },
    "C06_shutdown_margin_stuck": {
        "parameter": "Single Stuck-Drum Safety Margin (SDM_stuck)",
        "limit_value": 1000.0,
        "direction": ">=",
        "units": "pcm",
        "source": "Single Failure Safety Criterion (KAPL-4849, p. 42)",
        "status": "HARD_ENGINEERING_LIMIT",
        "description": "Ensures subcriticality even if the highest-worth drum jams fully withdrawn."
    },
    "C07_operating_keff": {
        "parameter": "Coupled Operating Multiplication Factor (k_eff)",
        "limit_value": 1.0120,
        "direction": ">=",
        "units": "-",
        "source": "Burnup Depletion Reserve Criterion (20+ years at 68.2 pcm/year)",
        "status": "HARD_ENGINEERING_LIMIT",
        "description": "Ensures excess reactivity > 1200 pcm for multi-decade NEP flight lifetime."
    },
    "C08_enrichment_ceiling": {
        "parameter": "Maximum Fuel Enrichment (wt% U-235)",
        "limit_value": 19.75,
        "direction": "<=",
        "units": "wt%",
        "source": "Strict Non-Proliferation / HALEU Commercial Licensing Standard",
        "status": "LEGAL_NON_PROLIFERATION_LIMIT",
        "description": "Must not exceed 19.75 wt% under any optimization state."
    },
    "C09_geometric_feasibility": {
        "parameter": "Geometry Constructibility & Clearances",
        "limit_value": 0.0,
        "direction": ">",
        "units": "-",
        "source": "First-Principles Solid Geometry",
        "status": "HARD_PHYSICAL_CONSTRAINT",
        "description": "Requires fuel_radius + clad_th + gap < pitch/sqrt(3) and positive flow areas."
    }
}


# -----------------------------------------------------------------------------
# 4. MULTI-OBJECTIVE FUNCTIONS DEFINITION
# -----------------------------------------------------------------------------
OBJECTIVE_FUNCTIONS = {
    "F1_subassembly_mass": {
        "name": "Minimize Core Subassembly Mass",
        "symbol": "m_sub",
        "units": "kg",
        "direction": "minimize",
        "baseline_val": 539.66,
        "description": "Reduces launch mass and spacecraft inertia."
    },
    "F2_peak_cladding_temp": {
        "name": "Minimize Peak Cladding Temperature",
        "symbol": "T_c,max",
        "units": "K",
        "direction": "minimize",
        "baseline_val": 1271.40,
        "description": "Expands thermal creep margin below the 1350 K limit."
    },
    "F3_core_pressure_drop": {
        "name": "Minimize Core Total Pressure Drop",
        "symbol": "Delta_P",
        "units": "kPa",
        "direction": "minimize",
        "baseline_val": 23.59,
        "description": "Maximizes Brayton cycle efficiency and reduces compressor power."
    },
    "F4_radial_power_peaking": {
        "name": "Minimize Total 3D Power Peaking",
        "symbol": "F_q",
        "units": "-",
        "direction": "minimize",
        "baseline_val": 1.707,
        "description": "Flattens core power profile and eliminates localized thermal hotspots."
    },
    "F5_hot_channel_elevation": {
        "name": "Minimize Hot-Channel Outlet Elevation",
        "symbol": "Delta_T_hot",
        "units": "K",
        "direction": "minimize",
        "baseline_val": 57.80,
        "description": "T_cool,exit,max - T_cool,bulk,out. Eliminates thermal striping into turbine."
    }
}


def export_baseline_design_vector():
    """Exports machine-readable baseline file to results/v5/v5_baseline_design_vector.json."""
    out_file = os.path.join(RESULTS_V5_DIR, "v5_baseline_design_vector.json")
    payload = {
        "baseline_design_vector": BASELINE_DESIGN_VECTOR,
        "design_variables_bounds": DESIGN_VARIABLES_BOUNDS,
        "engineering_constraints": ENGINEERING_CONSTRAINTS,
        "objective_functions": OBJECTIVE_FUNCTIONS
    }
    with open(out_file, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"Generated {out_file}")
    return out_file

if __name__ == "__main__":
    export_baseline_design_vector()
