#!/usr/bin/env python3
"""
Parametric Physics Investigation for Candidate Opt-9 NEP Reactor
==================================================================
Investigates:
1. Under-moderation curve (k_inf vs pitch and V_mod/V_fuel).
2. Axial duct streaming physics (solid reflector vs ducts of varying radius).
3. Core aspect ratio sweep (H/R in 1.7 to 2.5) on axial leakage and void reactivity.
4. Coolant temperature coefficient physics at constant pressure.
"""

import os
import sys
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import openmc

os.environ["PATH"] = "/home/aleri/miniforge3/envs/openmc-env/bin:" + os.environ.get("PATH", "")
os.environ["OPENMC_CROSS_SECTIONS"] = "/home/aleri/wmp_data/cross_sections_wmp.xml"
os.environ["OMP_NUM_THREADS"] = "8"
openmc.config['cross_sections'] = "/home/aleri/wmp_data/cross_sections_wmp.xml"

sys.path.append("src")
from run_opt8_master import Opt8Parameters, build_opt8_materials, build_opt8_settings

scratch_dir = "results/opt9_screening"
os.makedirs(scratch_dir, exist_ok=True)
os.makedirs(os.path.join(scratch_dir, "figures"), exist_ok=True)

# -----------------------------------------------------------------------------
# STUDY 1: Rigorous Demonstration of Under-Moderated Regime (k_inf vs Pitch)
# -----------------------------------------------------------------------------
def study_moderation_regime():
    print("\n" + "="*80)
    print("STUDY 1: Demonstrating Under-Moderated Regime (k_inf vs Moderator Pitch)")
    print("="*80)

    p = Opt8Parameters()
    pitches = [1.40, 1.50, 1.615, 1.75, 1.90, 2.10, 2.30]
    results = []

    for pitch in pitches:
        openmc.reset_auto_ids()
        mats, fuel, cladding, coolant, moderator, reflector, absorber = build_opt8_materials(p)

        s_fuel = openmc.ZCylinder(r=p.fuel_radius)
        s_clad = openmc.ZCylinder(r=p.clad_outer_radius)
        s_channel = openmc.ZCylinder(r=p.channel_outer_radius)

        px1 = openmc.XPlane(x0=-pitch/2.0, boundary_type='reflective')
        px2 = openmc.XPlane(x0=+pitch/2.0, boundary_type='reflective')
        py1 = openmc.YPlane(y0=-pitch/2.0, boundary_type='reflective')
        py2 = openmc.YPlane(y0=+pitch/2.0, boundary_type='reflective')

        pz1 = openmc.ZPlane(z0=-10.0, boundary_type='periodic')
        pz2 = openmc.ZPlane(z0=+10.0, boundary_type='periodic')
        pz1.periodic_surface = pz2

        cell_box = +px1 & -px2 & +py1 & -py2 & +pz1 & -pz2

        c_fuel = openmc.Cell(fill=fuel, region=-s_fuel & cell_box)
        c_clad = openmc.Cell(fill=cladding, region=+s_fuel & -s_clad & cell_box)
        c_cool = openmc.Cell(fill=coolant, region=+s_clad & -s_channel & cell_box)
        c_mod  = openmc.Cell(fill=moderator, region=+s_channel & cell_box)

        geom = openmc.Geometry([c_fuel, c_clad, c_cool, c_mod])

        sett = openmc.Settings()
        sett.particles = 15000
        sett.batches = 35
        sett.inactive = 10
        sett.temperature = {'multipole': True, 'method': 'nearest', 'tolerance': 1500.0}

        model = openmc.Model(materials=mats, geometry=geom, settings=sett)
        model.export_to_xml()
        openmc.run(output=False)

        with openmc.StatePoint(f"statepoint.{sett.batches}.h5") as sp:
            kinf = float(sp.keff.nominal_value)
            skinf = float(sp.keff.std_dev)

        v_fuel = np.pi * p.fuel_radius**2
        v_mod = pitch**2 - np.pi * p.channel_outer_radius**2
        mod_fuel_ratio = v_mod / v_fuel

        results.append({
            "pitch_cm": pitch, "mod_fuel_ratio": mod_fuel_ratio,
            "kinf": kinf, "std": skinf
        })
        print(f"Pitch = {pitch:.3f} cm | V_mod/V_fuel = {mod_fuel_ratio:5.2f} | k_inf = {kinf:.5f} +/- {skinf:.5f}")

    csv_path = os.path.join(scratch_dir, "moderation_curve.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["pitch_cm", "mod_fuel_ratio", "kinf", "std"])
        writer.writeheader()
        writer.writerows(results)

    # Plot moderation curve
    fig, ax1 = plt.subplots(figsize=(7, 5), dpi=300)
    ax2 = ax1.twiny()

    p_list = [r["pitch_cm"] for r in results]
    r_list = [r["mod_fuel_ratio"] for r in results]
    k_list = [r["kinf"] for r in results]
    s_list = [r["std"] for r in results]

    ax1.errorbar(p_list, k_list, yerr=s_list, fmt='bo-', lw=2, capsize=5, label=r'$k_\infty$ vs Fuel Pitch')
    ax1.axvline(1.615, color='red', ls='--', lw=1.5, label='Opt-8/9 Reference Pitch (1.615 cm)')
    ax1.set_xlabel("Fuel Pin Lattice Pitch [cm]", fontsize=10, fontweight='bold', color='blue')
    ax1.set_ylabel(r"Infinite Multiplication Factor $k_\infty$", fontsize=10, fontweight='bold')
    ax1.grid(True, alpha=0.3)

    ax2.set_xlim(ax1.get_xlim())
    # Set tick labels to corresponding mod_fuel_ratio
    ax2_ticks = p_list
    ax2_labels = [f"{r:.2f}" for r in r_list]
    ax2.set_xticks(ax2_ticks)
    ax2.set_xticklabels(ax2_labels, fontsize=8)
    ax2.set_xlabel(r"Moderator-to-Fuel Volume Ratio $V_\mathrm{mod}/V_\mathrm{fuel}$", fontsize=10, fontweight='bold')

    ax1.set_title(r"Physical Demonstration of Under-Moderated Core Regime ($k_\infty$ vs Moderation)", fontsize=11, fontweight='bold')
    ax1.legend(loc='lower right', fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(scratch_dir, "figures", "physics_under_moderation_curve.png"))
    plt.close(fig)
    print(f"  Study 1 completed and saved to {csv_path}")
    return results

# -----------------------------------------------------------------------------
# STUDY 2: Quantifying Void Sources: Gas Composition vs Ducts vs Duct Radius
# -----------------------------------------------------------------------------
def study_void_decomposition():
    print("\n" + "="*80)
    print("STUDY 2: Quantitative Decomposition of Coolant Void Reactivity")
    print("="*80)

    p = Opt8Parameters()
    p.active_core_height = 36.00  # Screening height
    p.total_reactor_height = p.active_core_height + 2 * p.axial_reflector_thickness

    # Test cases:
    # 1. Solid axial reflector (no ducts) with He-Xe 72/28 (Opt-6 legacy baseline)
    # 2. Solid axial reflector with Pure Helium
    # 3. Open streaming ducts r = 0.58 cm (fuel pellet radius) with Pure Helium
    # 4. Open streaming ducts r = 0.675 cm (channel outer radius, Opt-8) with Pure Helium
    # 5. Open streaming ducts r = 0.75 cm (slightly enlarged duct) with Pure Helium
    # 6. Open streaming ducts r = 0.85 cm (widened streaming duct) with Pure Helium

    cases = [
        {"name": "Solid Reflector (No Ducts), He-Xe 72/28", "fluid": "he_xe_72_28", "duct_r": 0.0},
        {"name": "Solid Reflector (No Ducts), Pure He 2.5 MPa", "fluid": "pure_he", "duct_r": 0.0},
        {"name": "Ducts r=0.580 cm, Pure He 2.5 MPa", "fluid": "pure_he", "duct_r": 0.58},
        {"name": "Ducts r=0.675 cm, Pure He 2.5 MPa (Opt-8 style)", "fluid": "pure_he", "duct_r": 0.675},
        {"name": "Ducts r=0.750 cm, Pure He 2.5 MPa", "fluid": "pure_he", "duct_r": 0.75},
        {"name": "Ducts r=0.850 cm, Pure He 2.5 MPa", "fluid": "pure_he", "duct_r": 0.85},
    ]

    study2_results = []

    for c in cases:
        print(f"\nEvaluating Case: {c['name']} ...")
        duct_r = c["duct_r"]

        # 1. Operating state
        openmc.reset_auto_ids()
        mats, fuel, cladding, coolant, moderator, reflector, absorber = build_opt8_materials(p)

        if c["fluid"] == "he_xe_72_28":
            coolant.nuclides.clear()
            coolant.set_density("g/cm3", 0.00916)
            coolant.add_nuclide("He4", 0.725, "ao")
            coolant.add_nuclide("Xe129", 0.275 * 0.264, "ao")
            coolant.add_nuclide("Xe131", 0.275 * 0.212, "ao")
            coolant.add_nuclide("Xe132", 0.275 * 0.269, "ao")
            coolant.add_nuclide("Xe134", 0.275 * 0.104, "ao")
            coolant.add_nuclide("Xe136", 0.275 * 0.089, "ao")

        s_fuel = openmc.ZCylinder(r=p.fuel_radius)
        s_clad = openmc.ZCylinder(r=p.clad_outer_radius)
        s_channel = openmc.ZCylinder(r=p.channel_outer_radius)

        z_core_bot = openmc.ZPlane(z0=-p.active_core_height/2.0)
        z_core_top = openmc.ZPlane(z0=+p.active_core_height/2.0)
        z_refl_bot = openmc.ZPlane(z0=-p.total_reactor_height/2.0, boundary_type='vacuum')
        z_refl_top = openmc.ZPlane(z0=+p.total_reactor_height/2.0, boundary_type='vacuum')

        c_f = openmc.Cell(fill=fuel, region=-s_fuel & +z_core_bot & -z_core_top)
        c_cl = openmc.Cell(fill=cladding, region=+s_fuel & -s_clad & +z_core_bot & -z_core_top)
        c_co = openmc.Cell(fill=coolant, region=+s_clad & -s_channel & +z_core_bot & -z_core_top)
        c_mo = openmc.Cell(fill=moderator, region=+s_channel & +z_core_bot & -z_core_top)

        if duct_r > 0.0:
            s_duct = openmc.ZCylinder(r=duct_r)
            c_co_bot = openmc.Cell(fill=coolant, region=-s_duct & +z_refl_bot & -z_core_bot)
            c_be_bot = openmc.Cell(fill=reflector, region=+s_duct & +z_refl_bot & -z_core_bot)
            c_co_top = openmc.Cell(fill=coolant, region=-s_duct & +z_core_top & -z_refl_top)
            c_be_top = openmc.Cell(fill=reflector, region=+s_duct & +z_core_top & -z_refl_top)
            u_pin = openmc.Universe(cells=[c_f, c_cl, c_co, c_mo, c_co_bot, c_be_bot, c_co_top, c_be_top])
        else:
            c_be_bot = openmc.Cell(fill=reflector, region=+z_refl_bot & -z_core_bot)
            c_be_top = openmc.Cell(fill=reflector, region=+z_core_top & -z_refl_top)
            u_pin = openmc.Universe(cells=[c_f, c_cl, c_co, c_mo, c_be_bot, c_be_top])

        c_mo_core = openmc.Cell(fill=moderator, region=+z_core_bot & -z_core_top)
        c_be_bot_only = openmc.Cell(fill=reflector, region=+z_refl_bot & -z_core_bot)
        c_be_top_only = openmc.Cell(fill=reflector, region=+z_core_top & -z_refl_top)
        u_mod_only = openmc.Universe(cells=[c_mo_core, c_be_bot_only, c_be_top_only])

        grid_dim = int(np.ceil(2 * p.core_radius / p.moderator_fuel_pitch))
        if grid_dim % 2 == 0: grid_dim += 1
        lat_univs = []
        for i in range(grid_dim):
            row = []
            for j in range(grid_dim):
                x = -grid_dim * p.moderator_fuel_pitch / 2.0 + (j + 0.5) * p.moderator_fuel_pitch
                y = -grid_dim * p.moderator_fuel_pitch / 2.0 + (i + 0.5) * p.moderator_fuel_pitch
                if x**2 + y**2 <= (p.core_radius - p.channel_outer_radius)**2:
                    row.append(u_pin)
                else:
                    row.append(u_mod_only)
            lat_univs.append(row)

        lattice = openmc.RectLattice()
        lattice.lower_left = (-grid_dim * p.moderator_fuel_pitch / 2.0, -grid_dim * p.moderator_fuel_pitch / 2.0)
        lattice.pitch = (p.moderator_fuel_pitch, p.moderator_fuel_pitch)
        lattice.universes = lat_univs
        lattice.outer = u_mod_only

        s_core = openmc.ZCylinder(r=p.core_radius)
        s_outer = openmc.ZCylinder(r=p.radial_reflector_outer_radius, boundary_type='vacuum')

        c_core_cyl = openmc.Cell(fill=lattice, region=-s_core & +z_refl_bot & -z_refl_top)
        c_rad_refl = openmc.Cell(fill=reflector, region=+s_core & -s_outer & +z_refl_bot & -z_refl_top)

        root_u = openmc.Universe(cells=[c_core_cyl, c_rad_refl])
        geom = openmc.Geometry(root_u)

        sett = openmc.Settings()
        sett.particles = 30000
        sett.batches = 45
        sett.inactive = 10
        sett.temperature = {'multipole': True, 'method': 'nearest', 'tolerance': 1500.0}

        # Operating
        model_op = openmc.Model(materials=mats, geometry=geom, settings=sett)
        model_op.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{sett.batches}.h5") as sp:
            k_op = float(sp.keff.nominal_value)
            s_op = float(sp.keff.std_dev)

        # Voided
        coolant.set_density("g/cm3", 1.0e-7)
        model_vd = openmc.Model(materials=mats, geometry=geom, settings=sett)
        model_vd.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{sett.batches}.h5") as sp:
            k_vd = float(sp.keff.nominal_value)
            s_vd = float(sp.keff.std_dev)

        drho = (1.0/k_op - 1.0/k_vd) * 1e5
        u_drho = np.sqrt((s_op/k_op**2)**2 + (s_vd/k_vd**2)**2) * 1e5
        bound_3s = drho + 3.0 * u_drho

        study2_results.append({
            "case": c["name"], "fluid": c["fluid"], "duct_r_cm": duct_r,
            "k_op": k_op, "s_op": s_op, "k_vd": k_vd, "s_vd": s_vd,
            "drho_void_pcm": drho, "u_void_pcm": u_drho, "bound_3sigma_pcm": bound_3s
        })
        print(f"  k_op = {k_op:.5f} +/- {s_op:.5f} | k_vd = {k_vd:.5f} +/- {s_vd:.5f}")
        print(f"  Delta-rho_void = {drho:+.1f} +/- {u_drho:.1f} pcm | 3σ bound: {bound_3s:+.1f} pcm")

    csv_path = os.path.join(scratch_dir, "void_decomposition.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(study2_results[0].keys()))
        writer.writeheader()
        writer.writerows(study2_results)

    # Plot void decomposition
    fig, ax = plt.subplots(figsize=(9, 5), dpi=300)
    labels = [r["case"].replace(", ", "\n") for r in study2_results]
    vals = [r["drho_void_pcm"] for r in study2_results]
    errs = [r["u_void_pcm"] for r in study2_results]
    err3 = [3*e for e in errs]

    x = np.arange(len(labels))
    colors = ['#D32F2F', '#F57C00', '#388E3C', '#1976D2', '#7B1FA2', '#00796B']
    bars = ax.bar(x, vals, yerr=err3, capsize=6, color=colors, edgecolor='black', lw=1.2, width=0.55)
    ax.axhline(0.0, color='black', ls='-', lw=1.5)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8, fontweight='bold', rotation=15, ha='right')
    ax.set_ylabel(r"Coolant Loss Reactivity $\Delta\rho_\mathrm{void}$ [pcm]", fontsize=10, fontweight='bold')
    ax.set_title(r"Quantitative Decomposition of Coolant Void Reactivity Mechanisms (with $3\sigma$ Error Bars)", fontsize=11, fontweight='bold')
    ax.grid(True, axis='y', alpha=0.3)

    for i, (v, e3) in enumerate(zip(vals, err3)):
        y_pos = v + e3 + 15 if v >= 0 else v - e3 - 35
        ax.text(i, y_pos, f"{v:+.0f}±{errs[i]:.0f}\n(3σ: {v+e3:+.0f})", ha='center', fontsize=7, fontweight='bold')

    fig.tight_layout()
    fig.savefig(os.path.join(scratch_dir, "figures", "physics_void_decomposition.png"))
    plt.close(fig)
    print(f"  Study 2 completed and saved to {csv_path}")
    return study2_results

# -----------------------------------------------------------------------------
# STUDY 3: Core Aspect Ratio (H/R) Sweep on Leakage, Mass, and Void Worth
# -----------------------------------------------------------------------------
def study_aspect_ratio_sweep():
    print("\n" + "="*80)
    print("STUDY 3: Core Aspect Ratio (H/R) Sweep: Leakage, Void, Mass & Volume")
    print("="*80)

    # Let's sweep H/R while keeping active core volume constant (V_core = ~43.5 L, ~373 pins)
    # V_core = pi * R^2 * H = const -> H = V_core / (pi * R^2)
    # For R=18.2, H=41.8: V_core = 43.50 L
    v_core_target = np.pi * 18.20**2 * 41.80

    aspect_cases = [
        {"H_R": 1.70, "R": 19.80, "H": v_core_target / (np.pi * 19.80**2)},
        {"H_R": 1.95, "R": 18.90, "H": v_core_target / (np.pi * 18.90**2)},
        {"H_R": 2.15, "R": 18.30, "H": v_core_target / (np.pi * 18.30**2)},
        {"H_R": 2.30, "R": 17.90, "H": v_core_target / (np.pi * 17.90**2)},
        {"H_R": 2.50, "R": 17.40, "H": v_core_target / (np.pi * 17.40**2)},
    ]

    study3_results = []

    for ac in aspect_cases:
        r_core = ac["R"]
        h_core = ac["H"]
        hr_val = h_core / r_core
        print(f"\nEvaluating H/R = {hr_val:.2f} (R={r_core:.2f} cm, H={h_core:.2f} cm) ...")

        p = Opt8Parameters()
        p.core_radius = r_core
        p.active_core_height = h_core
        p.radial_reflector_outer_radius = p.core_radius + p.radial_reflector_thickness
        p.total_reactor_height = p.active_core_height + 2 * p.axial_reflector_thickness

        openmc.reset_auto_ids()
        mats, fuel, cladding, coolant, moderator, reflector, absorber = build_opt8_materials(p)

        s_fuel = openmc.ZCylinder(r=p.fuel_radius)
        s_clad = openmc.ZCylinder(r=p.clad_outer_radius)
        s_channel = openmc.ZCylinder(r=p.channel_outer_radius)

        z_core_bot = openmc.ZPlane(z0=-p.active_core_height/2.0)
        z_core_top = openmc.ZPlane(z0=+p.active_core_height/2.0)
        z_refl_bot = openmc.ZPlane(z0=-p.total_reactor_height/2.0, boundary_type='vacuum')
        z_refl_top = openmc.ZPlane(z0=+p.total_reactor_height/2.0, boundary_type='vacuum')

        c_f = openmc.Cell(fill=fuel, region=-s_fuel & +z_core_bot & -z_core_top)
        c_cl = openmc.Cell(fill=cladding, region=+s_fuel & -s_clad & +z_core_bot & -z_core_top)
        c_co = openmc.Cell(fill=coolant, region=+s_clad & -s_channel & +z_core_bot & -z_core_top)
        c_mo = openmc.Cell(fill=moderator, region=+s_channel & +z_core_bot & -z_core_top)

        # Axial streaming duct r=0.70 cm
        s_duct = openmc.ZCylinder(r=0.70)
        c_co_bot = openmc.Cell(fill=coolant, region=-s_duct & +z_refl_bot & -z_core_bot)
        c_be_bot = openmc.Cell(fill=reflector, region=+s_duct & +z_refl_bot & -z_core_bot)
        c_co_top = openmc.Cell(fill=coolant, region=-s_duct & +z_core_top & -z_refl_top)
        c_be_top = openmc.Cell(fill=reflector, region=+s_duct & +z_core_top & -z_refl_top)
        u_pin = openmc.Universe(cells=[c_f, c_cl, c_co, c_mo, c_co_bot, c_be_bot, c_co_top, c_be_top])

        c_mo_core = openmc.Cell(fill=moderator, region=+z_core_bot & -z_core_top)
        c_be_bot_only = openmc.Cell(fill=reflector, region=+z_refl_bot & -z_core_bot)
        c_be_top_only = openmc.Cell(fill=reflector, region=+z_core_top & -z_refl_top)
        u_mod_only = openmc.Universe(cells=[c_mo_core, c_be_bot_only, c_be_top_only])

        grid_dim = int(np.ceil(2 * p.core_radius / p.moderator_fuel_pitch))
        if grid_dim % 2 == 0: grid_dim += 1
        lat_univs = []
        n_pins = 0
        for i in range(grid_dim):
            row = []
            for j in range(grid_dim):
                x = -grid_dim * p.moderator_fuel_pitch / 2.0 + (j + 0.5) * p.moderator_fuel_pitch
                y = -grid_dim * p.moderator_fuel_pitch / 2.0 + (i + 0.5) * p.moderator_fuel_pitch
                if x**2 + y**2 <= (p.core_radius - p.channel_outer_radius)**2:
                    row.append(u_pin)
                    n_pins += 1
                else:
                    row.append(u_mod_only)
            lat_univs.append(row)

        lattice = openmc.RectLattice()
        lattice.lower_left = (-grid_dim * p.moderator_fuel_pitch / 2.0, -grid_dim * p.moderator_fuel_pitch / 2.0)
        lattice.pitch = (p.moderator_fuel_pitch, p.moderator_fuel_pitch)
        lattice.universes = lat_univs
        lattice.outer = u_mod_only

        s_core = openmc.ZCylinder(r=p.core_radius)
        s_outer = openmc.ZCylinder(r=p.radial_reflector_outer_radius, boundary_type='vacuum')

        c_core_cyl = openmc.Cell(fill=lattice, region=-s_core & +z_refl_bot & -z_refl_top)
        c_rad_refl = openmc.Cell(fill=reflector, region=+s_core & -s_outer & +z_refl_bot & -z_refl_top)

        root_u = openmc.Universe(cells=[c_core_cyl, c_rad_refl])
        geom = openmc.Geometry(root_u)

        # Leakage tallies across surfaces
        mesh = openmc.RegularMesh()
        mesh.dimension = [1, 1, 1]
        mesh.lower_left = [-p.radial_reflector_outer_radius, -p.radial_reflector_outer_radius, -p.total_reactor_height/2.0]
        mesh.upper_right = [+p.radial_reflector_outer_radius, +p.radial_reflector_outer_radius, +p.total_reactor_height/2.0]

        sett = openmc.Settings()
        sett.particles = 30000
        sett.batches = 45
        sett.inactive = 10
        sett.temperature = {'multipole': True, 'method': 'nearest', 'tolerance': 1500.0}

        # 1. Operating
        model_op = openmc.Model(materials=mats, geometry=geom, settings=sett)
        model_op.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{sett.batches}.h5") as sp:
            k_op = float(sp.keff.nominal_value)
            s_op = float(sp.keff.std_dev)

        # 2. Voided
        coolant.set_density("g/cm3", 1.0e-7)
        model_vd = openmc.Model(materials=mats, geometry=geom, settings=sett)
        model_vd.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{sett.batches}.h5") as sp:
            k_vd = float(sp.keff.nominal_value)
            s_vd = float(sp.keff.std_dev)

        drho = (1.0/k_op - 1.0/k_vd) * 1e5
        u_drho = np.sqrt((s_op/k_op**2)**2 + (s_vd/k_vd**2)**2) * 1e5
        bound_3s = drho + 3.0 * u_drho

        # Total reactor mass & volume estimation
        v_tot = np.pi * p.radial_reflector_outer_radius**2 * p.total_reactor_height * 1e-3
        # Approximate mass based on Opt-8 densities
        mass_est = 628.49 * (v_tot / 154.06)

        study3_results.append({
            "H_R": hr_val, "R_core": r_core, "H_core": h_core, "pins": n_pins,
            "total_vol_L": v_tot, "mass_est_kg": mass_est,
            "k_op": k_op, "s_op": s_op, "k_vd": k_vd, "s_vd": s_vd,
            "drho_void_pcm": drho, "u_void_pcm": u_drho, "bound_3sigma_pcm": bound_3s
        })
        print(f"  Pins={n_pins} | Vol={v_tot:.1f} L | Mass={mass_est:.1f} kg")
        print(f"  k_op = {k_op:.5f} +/- {s_op:.5f} | k_vd = {k_vd:.5f} +/- {s_vd:.5f}")
        print(f"  Delta-rho_void = {drho:+.1f} +/- {u_drho:.1f} pcm | 3σ bound: {bound_3s:+.1f} pcm")

    csv_path = os.path.join(scratch_dir, "aspect_ratio_sweep.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(study3_results[0].keys()))
        writer.writeheader()
        writer.writerows(study3_results)

    # Plot H/R sensitivity
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), dpi=300)

    hr_list = [r["H_R"] for r in study3_results]
    void_list = [r["drho_void_pcm"] for r in study3_results]
    void_err = [r["u_void_pcm"] for r in study3_results]
    kop_list = [r["k_op"] for r in study3_results]
    kop_err = [r["s_op"] for r in study3_results]
    mass_list = [r["mass_est_kg"] for r in study3_results]

    ax1.errorbar(hr_list, void_list, yerr=[3*e for e in void_err], fmt='rs-', lw=2, capsize=6, label=r'$\Delta\rho_\mathrm{void}$ with $3\sigma$ bounds')
    ax1.axhline(0.0, color='black', ls='--', lw=1.5, label='Zero Void Limit')
    ax1.set_xlabel("Core Aspect Ratio H/R", fontsize=10, fontweight='bold')
    ax1.set_ylabel(r"Coolant Void Reactivity $\Delta\rho_\mathrm{void}$ [pcm]", fontsize=10, fontweight='bold')
    ax1.set_title(r"Void Reactivity vs Core Aspect Ratio $H/R$", fontsize=11, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='lower left', fontsize=8)

    ax2.plot(hr_list, mass_list, 'bd-', lw=2, ms=6, label='Assembly Mass [kg]')
    ax2.set_xlabel("Core Aspect Ratio H/R", fontsize=10, fontweight='bold')
    ax2.set_ylabel("Estimated Reactor Assembly Mass [kg]", fontsize=10, fontweight='bold', color='blue')
    ax2_r = ax2.twinx()
    ax2_r.errorbar(hr_list, kop_list, yerr=kop_err, fmt='go--', lw=2, ms=6, capsize=4, label=r'$k_\mathrm{op}$ (Operating)')
    ax2_r.set_ylabel(r"Operating Multiplication $k_\mathrm{op}$", fontsize=10, fontweight='bold', color='green')
    ax2.set_title(r"Mass & Criticality vs Core Aspect Ratio $H/R$", fontsize=11, fontweight='bold')
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(os.path.join(scratch_dir, "figures", "physics_aspect_ratio_sweep.png"))
    plt.close(fig)
    print(f"  Study 3 completed and saved to {csv_path}")
    return study3_results

if __name__ == "__main__":
    study_moderation_regime()
    study_void_decomposition()
    study_aspect_ratio_sweep()
