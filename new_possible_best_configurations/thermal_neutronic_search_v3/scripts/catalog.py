#!/usr/bin/env python3
"""
catalog.py - Catalog of 30 deterministic candidate NEP reactor configurations
(400 kWth, 100 kWe, YHx moderator, He coolant 850-900 K inlet, 1100-1150 K outlet).
"""

from pathlib import Path
import json
import math
import csv

ROOT = Path(__file__).resolve().parents[1]

def get_candidates():
    """
    Generate 30 deterministic, physically well-grounded NEP reactor candidates.
    Organized into 5 distinct architectural families:
    1. Dual-Pass Matrix Block with 6 Aux Channels (217 Pins)
    2. Dual-Annulus Concentric Sleeve (217 Pins, Counter-Current)
    3. Expanded Pitch Matrix with Enhanced Thermal Buffer (217 Pins, Pitch 2.36-2.42 cm)
    4. 271-Pin Core (10 Rings, Lower Linear Rating, 12 Drums)
    5. 169-Pin Core (8 Rings, Deep Moderation, Thick Reflector)
    """
    candidates = []

    # -------------------------------------------------------------------------
    # FAMILY 1: Dual-Pass Hexagonal Matrix Block (217 Pins, 6 Aux Channels/cell)
    # -------------------------------------------------------------------------
    f1_specs = [
        # C01: Baseline evolution of TB35 with slightly wider aux channels (1.65 mm)
        dict(id="C01_matrix_217_yh20_880_1120", family="matrix_dual_pass",
             pins=217, rings=8, pitch_cm=2.28, height_cm=46.0,
             r_fuel_cm=0.585, t_clad_cm=0.048, t_hot_gap_cm=0.084, t_al2o3_cm=0.010,
             moderator_type="YH2.0", y_stoich=2.0, rho_mod=4.30,
             aux_holes=6, r_aux_cm=0.165,
             reflector="BeO", t_refl_rad_cm=9.30, t_refl_ax_cm=8.80,
             drums=10, r_drum_cm=3.95, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=880.0, outlet_K=1120.0, power_kW=400.0, pressure_MPa=2.5),

        # C02: Slightly taller (47.0 cm) with larger aux holes (1.80 mm)
        dict(id="C02_matrix_217_yh20_870_1130", family="matrix_dual_pass",
             pins=217, rings=8, pitch_cm=2.30, height_cm=47.0,
             r_fuel_cm=0.585, t_clad_cm=0.048, t_hot_gap_cm=0.084, t_al2o3_cm=0.010,
             moderator_type="YH2.0", y_stoich=2.0, rho_mod=4.30,
             aux_holes=6, r_aux_cm=0.180,
             reflector="BeO", t_refl_rad_cm=9.30, t_refl_ax_cm=8.80,
             drums=10, r_drum_cm=3.95, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=870.0, outlet_K=1130.0, power_kW=400.0, pressure_MPa=2.5),

        # C03: Sub-stoichiometric YH1.90 with 1.70 mm aux holes
        dict(id="C03_matrix_217_yh19_880_1120", family="matrix_dual_pass",
             pins=217, rings=8, pitch_cm=2.29, height_cm=47.0,
             r_fuel_cm=0.585, t_clad_cm=0.048, t_hot_gap_cm=0.084, t_al2o3_cm=0.010,
             moderator_type="YH1.9", y_stoich=1.9, rho_mod=4.27,
             aux_holes=6, r_aux_cm=0.170,
             reflector="BeO", t_refl_rad_cm=9.30, t_refl_ax_cm=8.80,
             drums=10, r_drum_cm=3.95, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=880.0, outlet_K=1120.0, power_kW=400.0, pressure_MPa=2.5),

        # C04: YH1.84 (often cited in ORNL handbook) with 1.80 mm aux holes
        dict(id="C04_matrix_217_yh184_870_1120", family="matrix_dual_pass",
             pins=217, rings=8, pitch_cm=2.30, height_cm=48.0,
             r_fuel_cm=0.585, t_clad_cm=0.048, t_hot_gap_cm=0.084, t_al2o3_cm=0.010,
             moderator_type="YH1.84", y_stoich=1.84, rho_mod=4.25,
             aux_holes=6, r_aux_cm=0.180,
             reflector="BeO", t_refl_rad_cm=9.30, t_refl_ax_cm=8.80,
             drums=10, r_drum_cm=3.95, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=870.0, outlet_K=1120.0, power_kW=400.0, pressure_MPa=2.5),

        # C05: YH1.80, cold inlet 860 K, outlet 1110 K, 2.0 mm aux holes
        dict(id="C05_matrix_217_yh18_860_1110", family="matrix_dual_pass",
             pins=217, rings=8, pitch_cm=2.32, height_cm=49.0,
             r_fuel_cm=0.585, t_clad_cm=0.048, t_hot_gap_cm=0.084, t_al2o3_cm=0.010,
             moderator_type="YH1.8", y_stoich=1.8, rho_mod=4.24,
             aux_holes=6, r_aux_cm=0.200,
             reflector="BeO", t_refl_rad_cm=9.30, t_refl_ax_cm=8.80,
             drums=10, r_drum_cm=3.95, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=860.0, outlet_K=1110.0, power_kW=400.0, pressure_MPa=2.5),

        # C06: Beryllium metal reflector with YH2.0, inlet 890 K, outlet 1140 K
        dict(id="C06_matrix_217_be_metal_refl", family="matrix_dual_pass",
             pins=217, rings=8, pitch_cm=2.30, height_cm=47.5,
             r_fuel_cm=0.585, t_clad_cm=0.048, t_hot_gap_cm=0.084, t_al2o3_cm=0.010,
             moderator_type="YH2.0", y_stoich=2.0, rho_mod=4.30,
             aux_holes=6, r_aux_cm=0.180,
             reflector="Be", t_refl_rad_cm=9.50, t_refl_ax_cm=9.00,
             drums=10, r_drum_cm=4.10, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=890.0, outlet_K=1140.0, power_kW=400.0, pressure_MPa=2.5),
    ]
    candidates.extend(f1_specs)

    # -------------------------------------------------------------------------
    # FAMILY 2: Dual-Annulus Concentric Sleeve (Counter-Current Cooling)
    # -------------------------------------------------------------------------
    f2_specs = [
        # C07: Dual annulus with YH2.0, inlet 870 K, outlet 1120 K
        dict(id="C07_dual_annulus_yh20_870_1120", family="dual_annulus",
             pins=217, rings=8, pitch_cm=2.35, height_cm=46.0,
             r_fuel_cm=0.585, t_clad_cm=0.048, t_hot_gap_cm=0.080, t_part_cm=0.035,
             t_cold_gap_cm=0.080, t_al2o3_cm=0.010,
             moderator_type="YH2.0", y_stoich=2.0, rho_mod=4.30,
             reflector="BeO", t_refl_rad_cm=9.30, t_refl_ax_cm=8.80,
             drums=10, r_drum_cm=3.95, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=870.0, outlet_K=1120.0, power_kW=400.0, pressure_MPa=2.5),

        # C08: Dual annulus with YH2.0, inlet 880 K, outlet 1130 K
        dict(id="C08_dual_annulus_yh20_880_1130", family="dual_annulus",
             pins=217, rings=8, pitch_cm=2.36, height_cm=47.0,
             r_fuel_cm=0.585, t_clad_cm=0.048, t_hot_gap_cm=0.080, t_part_cm=0.035,
             t_cold_gap_cm=0.085, t_al2o3_cm=0.010,
             moderator_type="YH2.0", y_stoich=2.0, rho_mod=4.30,
             reflector="BeO", t_refl_rad_cm=9.30, t_refl_ax_cm=8.80,
             drums=10, r_drum_cm=3.95, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=880.0, outlet_K=1130.0, power_kW=400.0, pressure_MPa=2.5),

        # C09: Dual annulus with YH1.84, inlet 870 K, outlet 1120 K
        dict(id="C09_dual_annulus_yh184_870_1120", family="dual_annulus",
             pins=217, rings=8, pitch_cm=2.38, height_cm=48.0,
             r_fuel_cm=0.585, t_clad_cm=0.048, t_hot_gap_cm=0.080, t_part_cm=0.035,
             t_cold_gap_cm=0.085, t_al2o3_cm=0.010,
             moderator_type="YH1.84", y_stoich=1.84, rho_mod=4.25,
             reflector="BeO", t_refl_rad_cm=9.30, t_refl_ax_cm=8.80,
             drums=10, r_drum_cm=3.95, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=870.0, outlet_K=1120.0, power_kW=400.0, pressure_MPa=2.5),

        # C10: Dual annulus with YH1.80, inlet 860 K, outlet 1110 K
        dict(id="C10_dual_annulus_yh18_860_1110", family="dual_annulus",
             pins=217, rings=8, pitch_cm=2.40, height_cm=49.0,
             r_fuel_cm=0.585, t_clad_cm=0.048, t_hot_gap_cm=0.085, t_part_cm=0.035,
             t_cold_gap_cm=0.090, t_al2o3_cm=0.010,
             moderator_type="YH1.8", y_stoich=1.8, rho_mod=4.24,
             reflector="BeO", t_refl_rad_cm=9.30, t_refl_ax_cm=8.80,
             drums=10, r_drum_cm=3.95, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=860.0, outlet_K=1110.0, power_kW=400.0, pressure_MPa=2.5),

        # C11: Dual annulus with Be metal reflector
        dict(id="C11_dual_annulus_be_metal", family="dual_annulus",
             pins=217, rings=8, pitch_cm=2.35, height_cm=46.5,
             r_fuel_cm=0.585, t_clad_cm=0.048, t_hot_gap_cm=0.080, t_part_cm=0.035,
             t_cold_gap_cm=0.080, t_al2o3_cm=0.010,
             moderator_type="YH2.0", y_stoich=2.0, rho_mod=4.30,
             reflector="Be", t_refl_rad_cm=9.50, t_refl_ax_cm=9.00,
             drums=10, r_drum_cm=4.10, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=890.0, outlet_K=1140.0, power_kW=400.0, pressure_MPa=2.5),

        # C12: Dual annulus with YH1.90, inlet 880 K, outlet 1120 K
        dict(id="C12_dual_annulus_yh19_880_1120", family="dual_annulus",
             pins=217, rings=8, pitch_cm=2.38, height_cm=48.0,
             r_fuel_cm=0.585, t_clad_cm=0.048, t_hot_gap_cm=0.080, t_part_cm=0.035,
             t_cold_gap_cm=0.085, t_al2o3_cm=0.010,
             moderator_type="YH1.9", y_stoich=1.9, rho_mod=4.27,
             reflector="BeO", t_refl_rad_cm=9.30, t_refl_ax_cm=8.80,
             drums=10, r_drum_cm=3.95, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=880.0, outlet_K=1120.0, power_kW=400.0, pressure_MPa=2.5),
    ]
    candidates.extend(f2_specs)

    # -------------------------------------------------------------------------
    # FAMILY 3: Expanded Pitch Matrix (Enhanced Thermal Buffer & Lifetime)
    # -------------------------------------------------------------------------
    f3_specs = [
        # C13: Pitch 2.36 cm, 8 aux holes (1.40 mm)
        dict(id="C13_exp_pitch_236_yh20", family="matrix_expanded_pitch",
             pins=217, rings=8, pitch_cm=2.36, height_cm=48.0,
             r_fuel_cm=0.585, t_clad_cm=0.048, t_hot_gap_cm=0.084, t_al2o3_cm=0.010,
             moderator_type="YH2.0", y_stoich=2.0, rho_mod=4.30,
             aux_holes=8, r_aux_cm=0.140,
             reflector="BeO", t_refl_rad_cm=9.30, t_refl_ax_cm=8.80,
             drums=10, r_drum_cm=4.00, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=870.0, outlet_K=1120.0, power_kW=400.0, pressure_MPa=2.5),

        # C14: Pitch 2.38 cm, 6 aux holes (1.80 mm)
        dict(id="C14_exp_pitch_238_yh20", family="matrix_expanded_pitch",
             pins=217, rings=8, pitch_cm=2.38, height_cm=49.0,
             r_fuel_cm=0.585, t_clad_cm=0.048, t_hot_gap_cm=0.084, t_al2o3_cm=0.010,
             moderator_type="YH2.0", y_stoich=2.0, rho_mod=4.30,
             aux_holes=6, r_aux_cm=0.180,
             reflector="BeO", t_refl_rad_cm=9.30, t_refl_ax_cm=8.80,
             drums=10, r_drum_cm=4.00, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=880.0, outlet_K=1130.0, power_kW=400.0, pressure_MPa=2.5),

        # C15: Pitch 2.40 cm, YH1.84, 6 aux holes (1.90 mm)
        dict(id="C15_exp_pitch_240_yh184", family="matrix_expanded_pitch",
             pins=217, rings=8, pitch_cm=2.40, height_cm=50.0,
             r_fuel_cm=0.585, t_clad_cm=0.048, t_hot_gap_cm=0.084, t_al2o3_cm=0.010,
             moderator_type="YH1.84", y_stoich=1.84, rho_mod=4.25,
             aux_holes=6, r_aux_cm=0.190,
             reflector="BeO", t_refl_rad_cm=9.30, t_refl_ax_cm=8.80,
             drums=10, r_drum_cm=4.10, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=870.0, outlet_K=1120.0, power_kW=400.0, pressure_MPa=2.5),

        # C16: Pitch 2.42 cm, YH1.80, 6 aux holes (2.00 mm)
        dict(id="C16_exp_pitch_242_yh18", family="matrix_expanded_pitch",
             pins=217, rings=8, pitch_cm=2.42, height_cm=51.0,
             r_fuel_cm=0.585, t_clad_cm=0.048, t_hot_gap_cm=0.084, t_al2o3_cm=0.010,
             moderator_type="YH1.8", y_stoich=1.8, rho_mod=4.24,
             aux_holes=6, r_aux_cm=0.200,
             reflector="BeO", t_refl_rad_cm=9.30, t_refl_ax_cm=8.80,
             drums=10, r_drum_cm=4.10, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=860.0, outlet_K=1110.0, power_kW=400.0, pressure_MPa=2.5),

        # C17: Pitch 2.38 cm, 12 aux holes (1.20 mm)
        dict(id="C17_exp_pitch_238_12aux", family="matrix_expanded_pitch",
             pins=217, rings=8, pitch_cm=2.38, height_cm=48.5,
             r_fuel_cm=0.585, t_clad_cm=0.048, t_hot_gap_cm=0.084, t_al2o3_cm=0.010,
             moderator_type="YH2.0", y_stoich=2.0, rho_mod=4.30,
             aux_holes=12, r_aux_cm=0.120,
             reflector="BeO", t_refl_rad_cm=9.30, t_refl_ax_cm=8.80,
             drums=10, r_drum_cm=4.00, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=890.0, outlet_K=1150.0, power_kW=400.0, pressure_MPa=2.5),

        # C18: Pitch 2.40 cm, YH1.90, 6 aux holes (1.80 mm)
        dict(id="C18_exp_pitch_240_yh19", family="matrix_expanded_pitch",
             pins=217, rings=8, pitch_cm=2.40, height_cm=49.5,
             r_fuel_cm=0.585, t_clad_cm=0.048, t_hot_gap_cm=0.084, t_al2o3_cm=0.010,
             moderator_type="YH1.9", y_stoich=1.9, rho_mod=4.27,
             aux_holes=6, r_aux_cm=0.180,
             reflector="BeO", t_refl_rad_cm=9.30, t_refl_ax_cm=8.80,
             drums=10, r_drum_cm=4.10, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=875.0, outlet_K=1125.0, power_kW=400.0, pressure_MPa=2.5),
    ]
    candidates.extend(f3_specs)

    # -------------------------------------------------------------------------
    # FAMILY 4: High-Burnup 271-Pin Core (10 Rings, 12 Control Drums)
    # -------------------------------------------------------------------------
    f4_specs = [
        # C19: 271 pins, r_fuel = 0.53 cm, pitch = 2.15 cm
        dict(id="C19_271pins_yh20_870_1120", family="pins271_high_burnup",
             pins=271, rings=9, pitch_cm=2.15, height_cm=46.0,
             r_fuel_cm=0.530, t_clad_cm=0.048, t_hot_gap_cm=0.080, t_al2o3_cm=0.010,
             moderator_type="YH2.0", y_stoich=2.0, rho_mod=4.30,
             aux_holes=6, r_aux_cm=0.160,
             reflector="BeO", t_refl_rad_cm=9.00, t_refl_ax_cm=8.50,
             drums=12, r_drum_cm=3.70, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=870.0, outlet_K=1120.0, power_kW=400.0, pressure_MPa=2.5),

        # C20: 271 pins, r_fuel = 0.54 cm, pitch = 2.18 cm
        dict(id="C20_271pins_yh20_880_1130", family="pins271_high_burnup",
             pins=271, rings=9, pitch_cm=2.18, height_cm=47.0,
             r_fuel_cm=0.540, t_clad_cm=0.048, t_hot_gap_cm=0.080, t_al2o3_cm=0.010,
             moderator_type="YH2.0", y_stoich=2.0, rho_mod=4.30,
             aux_holes=6, r_aux_cm=0.165,
             reflector="BeO", t_refl_rad_cm=9.00, t_refl_ax_cm=8.50,
             drums=12, r_drum_cm=3.70, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=880.0, outlet_K=1130.0, power_kW=400.0, pressure_MPa=2.5),

        # C21: 271 pins, YH1.84, pitch = 2.20 cm
        dict(id="C21_271pins_yh184_870_1120", family="pins271_high_burnup",
             pins=271, rings=9, pitch_cm=2.20, height_cm=48.0,
             r_fuel_cm=0.540, t_clad_cm=0.048, t_hot_gap_cm=0.080, t_al2o3_cm=0.010,
             moderator_type="YH1.84", y_stoich=1.84, rho_mod=4.25,
             aux_holes=6, r_aux_cm=0.170,
             reflector="BeO", t_refl_rad_cm=9.00, t_refl_ax_cm=8.50,
             drums=12, r_drum_cm=3.70, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=870.0, outlet_K=1120.0, power_kW=400.0, pressure_MPa=2.5),

        # C22: 271 pins, YH1.80, pitch = 2.22 cm
        dict(id="C22_271pins_yh18_860_1110", family="pins271_high_burnup",
             pins=271, rings=9, pitch_cm=2.22, height_cm=49.0,
             r_fuel_cm=0.550, t_clad_cm=0.048, t_hot_gap_cm=0.080, t_al2o3_cm=0.010,
             moderator_type="YH1.8", y_stoich=1.8, rho_mod=4.24,
             aux_holes=6, r_aux_cm=0.180,
             reflector="BeO", t_refl_rad_cm=9.00, t_refl_ax_cm=8.50,
             drums=12, r_drum_cm=3.70, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=860.0, outlet_K=1110.0, power_kW=400.0, pressure_MPa=2.5),

        # C23: 271 pins, Be metal reflector
        dict(id="C23_271pins_be_metal", family="pins271_high_burnup",
             pins=271, rings=9, pitch_cm=2.16, height_cm=46.5,
             r_fuel_cm=0.530, t_clad_cm=0.048, t_hot_gap_cm=0.080, t_al2o3_cm=0.010,
             moderator_type="YH2.0", y_stoich=2.0, rho_mod=4.30,
             aux_holes=6, r_aux_cm=0.160,
             reflector="Be", t_refl_rad_cm=9.20, t_refl_ax_cm=8.80,
             drums=12, r_drum_cm=3.80, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=890.0, outlet_K=1140.0, power_kW=400.0, pressure_MPa=2.5),

        # C24: 271 pins, YH1.90, pitch = 2.18 cm
        dict(id="C24_271pins_yh19_875_1125", family="pins271_high_burnup",
             pins=271, rings=9, pitch_cm=2.18, height_cm=47.5,
             r_fuel_cm=0.540, t_clad_cm=0.048, t_hot_gap_cm=0.080, t_al2o3_cm=0.010,
             moderator_type="YH1.9", y_stoich=1.9, rho_mod=4.27,
             aux_holes=6, r_aux_cm=0.165,
             reflector="BeO", t_refl_rad_cm=9.00, t_refl_ax_cm=8.50,
             drums=12, r_drum_cm=3.70, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=875.0, outlet_K=1125.0, power_kW=400.0, pressure_MPa=2.5),
    ]
    candidates.extend(f4_specs)

    # -------------------------------------------------------------------------
    # FAMILY 5: High-Moderation 169-Pin Core (8 Rings, Deep Control Margin)
    # -------------------------------------------------------------------------
    f5_specs = [
        # C25: 169 pins, r_fuel = 0.63 cm, pitch = 2.52 cm
        dict(id="C25_169pins_yh20_870_1120", family="pins169_deep_mod",
             pins=169, rings=7, pitch_cm=2.52, height_cm=48.0,
             r_fuel_cm=0.630, t_clad_cm=0.048, t_hot_gap_cm=0.085, t_al2o3_cm=0.010,
             moderator_type="YH2.0", y_stoich=2.0, rho_mod=4.30,
             aux_holes=6, r_aux_cm=0.170,
             reflector="BeO", t_refl_rad_cm=9.80, t_refl_ax_cm=9.00,
             drums=10, r_drum_cm=4.20, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=870.0, outlet_K=1120.0, power_kW=400.0, pressure_MPa=2.5),

        # C26: 169 pins, r_fuel = 0.64 cm, pitch = 2.55 cm
        dict(id="C26_169pins_yh20_880_1130", family="pins169_deep_mod",
             pins=169, rings=7, pitch_cm=2.55, height_cm=49.0,
             r_fuel_cm=0.640, t_clad_cm=0.048, t_hot_gap_cm=0.085, t_al2o3_cm=0.010,
             moderator_type="YH2.0", y_stoich=2.0, rho_mod=4.30,
             aux_holes=6, r_aux_cm=0.180,
             reflector="BeO", t_refl_rad_cm=9.80, t_refl_ax_cm=9.00,
             drums=10, r_drum_cm=4.20, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=880.0, outlet_K=1130.0, power_kW=400.0, pressure_MPa=2.5),

        # C27: 169 pins, YH1.84, pitch = 2.58 cm
        dict(id="C27_169pins_yh184_870_1120", family="pins169_deep_mod",
             pins=169, rings=7, pitch_cm=2.58, height_cm=50.0,
             r_fuel_cm=0.640, t_clad_cm=0.048, t_hot_gap_cm=0.085, t_al2o3_cm=0.010,
             moderator_type="YH1.84", y_stoich=1.84, rho_mod=4.25,
             aux_holes=6, r_aux_cm=0.180,
             reflector="BeO", t_refl_rad_cm=9.80, t_refl_ax_cm=9.00,
             drums=10, r_drum_cm=4.20, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=870.0, outlet_K=1120.0, power_kW=400.0, pressure_MPa=2.5),

        # C28: 169 pins, YH1.80, pitch = 2.60 cm
        dict(id="C28_169pins_yh18_860_1110", family="pins169_deep_mod",
             pins=169, rings=7, pitch_cm=2.60, height_cm=51.0,
             r_fuel_cm=0.650, t_clad_cm=0.048, t_hot_gap_cm=0.085, t_al2o3_cm=0.010,
             moderator_type="YH1.8", y_stoich=1.8, rho_mod=4.24,
             aux_holes=6, r_aux_cm=0.190,
             reflector="BeO", t_refl_rad_cm=9.80, t_refl_ax_cm=9.00,
             drums=10, r_drum_cm=4.20, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=860.0, outlet_K=1110.0, power_kW=400.0, pressure_MPa=2.5),

        # C29: 169 pins, Be metal reflector
        dict(id="C29_169pins_be_metal", family="pins169_deep_mod",
             pins=169, rings=7, pitch_cm=2.54, height_cm=48.5,
             r_fuel_cm=0.630, t_clad_cm=0.048, t_hot_gap_cm=0.085, t_al2o3_cm=0.010,
             moderator_type="YH2.0", y_stoich=2.0, rho_mod=4.30,
             aux_holes=6, r_aux_cm=0.170,
             reflector="Be", t_refl_rad_cm=10.00, t_refl_ax_cm=9.20,
             drums=10, r_drum_cm=4.30, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=890.0, outlet_K=1140.0, power_kW=400.0, pressure_MPa=2.5),

        # C30: 169 pins, YH1.90, pitch = 2.56 cm
        dict(id="C30_169pins_yh19_875_1125", family="pins169_deep_mod",
             pins=169, rings=7, pitch_cm=2.56, height_cm=49.5,
             r_fuel_cm=0.640, t_clad_cm=0.048, t_hot_gap_cm=0.085, t_al2o3_cm=0.010,
             moderator_type="YH1.9", y_stoich=1.9, rho_mod=4.27,
             aux_holes=6, r_aux_cm=0.180,
             reflector="BeO", t_refl_rad_cm=9.80, t_refl_ax_cm=9.00,
             drums=10, r_drum_cm=4.20, t_b4c_cm=1.80, arc_b4c_deg=155.0,
             enrichment_inner=0.1960, enrichment_outer=0.1975,
             inlet_K=875.0, outlet_K=1125.0, power_kW=400.0, pressure_MPa=2.5),
    ]
    candidates.extend(f5_specs)

    # Calculate analytical mass and fuel inventory metrics for all candidates
    for c in candidates:
        n = c['pins']
        h = c['height_cm']
        rf = c['r_fuel_cm']
        v_fuel = n * math.pi * (rf**2) * h
        m_fuel = v_fuel * 13.50 / 1000.0 # kg UN
        m_hm = m_fuel * (238.0 / (238.0 + 14.0)) # approx heavy metal
        # inner 19 pins, outer n-19 pins
        enr_avg = (19 * c['enrichment_inner'] + (n - 19) * c['enrichment_outer']) / n
        m_u235 = m_hm * enr_avg

        # core radius avoiding pin clipping
        r_core = (c['rings'] + 0.57) * c['pitch_cm']
        r_ref = r_core + c['t_refl_rad_cm']
        h_tot = h + 2.0 * c['t_refl_ax_cm']
        rho_ref = 2.90 if c['reflector'] == "BeO" else 1.85
        v_core_cyl = math.pi * (r_core**2) * h
        v_ref_rad = math.pi * (r_ref**2 - r_core**2) * h
        v_ref_ax = 2.0 * math.pi * (r_ref**2) * c['t_refl_ax_cm']
        m_ref = (v_ref_rad + v_ref_ax) * rho_ref / 1000.0

        # Moderator volume approx 50-60% of core
        v_mod = v_core_cyl * 0.45
        m_mod = v_mod * c['rho_mod'] / 1000.0
        m_clad = n * math.pi * ((rf + c['t_clad_cm'])**2 - rf**2) * h * 8.57 / 1000.0
        m_drums_b4c = c['drums'] * (c['arc_b4c_deg']/360.0) * math.pi * (c['r_drum_cm']**2 - (c['r_drum_cm'] - c['t_b4c_cm'])**2) * h * 2.52 / 1000.0

        m_total_est = m_fuel + m_clad + m_mod + m_ref + m_drums_b4c + 30.0 # structure

        c['heavy_metal_kg'] = round(m_hm, 2)
        c['u235_kg'] = round(m_u235, 2)
        c['mass_est_kg'] = round(m_total_est, 1)
        c['core_radius_cm'] = round(r_core, 2)
        c['r_ref_cm'] = round(r_ref, 2)

    return candidates

def export_catalog(dest_dir=None):
    if dest_dir is None:
        dest_dir = ROOT
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    cands = get_candidates()

    json_path = dest_dir / "configurations.json"
    with open(json_path, "w") as f:
        json.dump(cands, f, indent=2)

    csv_path = dest_dir / "configurations.csv"
    keys = ["id", "family", "pins", "pitch_cm", "height_cm", "moderator_type", "y_stoich",
            "inlet_K", "outlet_K", "heavy_metal_kg", "u235_kg", "mass_est_kg", "reflector", "drums"]
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        for c in cands:
            writer.writerow(c)

    print(f"Exported {len(cands)} candidates to {json_path} and {csv_path}")

if __name__ == "__main__":
    export_catalog()
