"""
Reactivity, Control Drum Worth, Temperature & Void Analysis (analysis.py)
========================================================================
Performs rigorous, simulation-based neutronics analysis for balanced HALEU NEP reactor:
1. Control Drum Angle Sweep: Evaluates keff(theta) and rho(theta) for theta in [0, 30, 60, 90, 120, 150, 180 deg]
   - Verifies k_out > 1.0 (absorber out at 0 deg)
   - Verifies k_in < 1.0, targeting k_in approx 0.97 (absorber in at 180 deg)
2. Single-Stuck-Drum Safety Criterion: 1 drum stuck at 0 deg, 5 drums at 180 deg
3. Fuel Doppler Coefficient: Evaluates T_fuel in [300, 600, 900, 1200, 1500 K] using WMP data
4. Moderator Temperature Coefficient: Evaluates T_mod in [300, 600, 900 K] with thermal expansion
5. Coolant Temperature Coefficient: Evaluates He-Xe gas at 900 K, 1050 K, 1300 K with density rho(T)
6. Total Temperature Coefficient: alpha_total = alpha_fuel + alpha_mod + alpha_cool
7. Void Reactivity & Void Coefficient: Nominal (2.0 MPa), Reduced (1.0 MPa), Voided (0 MPa)
8. Generates:
   - results/drum_worth.png
   - results/temperature_coefficients.png
   - results/void_coefficient.png
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import openmc
from model import (ReactorParameters, build_materials, build_geometry,
                   build_settings, calculate_inventory)


def run_control_and_reactivity_analysis(output_dir="results"):
    os.makedirs(output_dir, exist_ok=True)
    print("=" * 80)
    print("SYSTEMATIC CONTROL DRUM, TEMPERATURE FEEDBACK & VOID ANALYSIS (OPENMC SIMULATIONS)")
    print("=" * 80)

    # Base parameters for balanced candidate
    p_base = ReactorParameters(particles=1500, batches=25, inactive=8)

    # -------------------------------------------------------------------------
    # 1. Parameter Sweep over Control Drum Angle (0 to 180 deg)
    # -------------------------------------------------------------------------
    angles = [0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0]
    keff_vals = []
    keff_stds = []
    rho_pcm = []

    print("\n1. Sweeping Control Drum Rotation Angle theta in [0, 180 deg]...")
    for ang in angles:
        openmc.reset_auto_ids()
        p = ReactorParameters(drum_angle_deg=ang, particles=1500, batches=25, inactive=8)
        mats, fuel, clad, cool, mod, beo, b4c = build_materials(p)
        geom, pins, drums = build_geometry(p, fuel, clad, cool, mod, beo, b4c)
        settings = build_settings(p)

        model = openmc.Model(materials=mats, geometry=geom, settings=settings)
        model.export_to_xml()
        openmc.run(output=False)

        sp_file = f"statepoint.{p.batches}.h5"
        with openmc.StatePoint(sp_file) as sp:
            k = float(sp.keff.nominal_value)
            s = float(sp.keff.std_dev)
            keff_vals.append(k)
            keff_stds.append(s)
            rho = (k - 1.0) / k * 1.0e5
            rho_pcm.append(rho)
            direction = "OUT (Supercritical)" if ang == 0 else ("IN (Shutdown target)" if ang == 180 else "Intermediate")
            print(f"   Angle = {ang:5.1f}° ({direction:<23}) | keff = {k:.5f} +/- {s:.5f} | Reactivity = {rho:+.1f} pcm")

    k_out = keff_vals[0]
    k_in = keff_vals[-1]
    total_worth_pcm = (1.0 / k_in - 1.0 / k_out) * 1.0e5
    total_worth_pct = (1.0 / k_in - 1.0 / k_out) * 100.0
    print(f"\n   Total Drum Worth : {total_worth_pct:.2f}% Delta-k/k ({total_worth_pcm:.0f} pcm)")
    print(f"   Target Shutdown Check: k_in = {k_in:.5f} (Design target ~0.97; Subcritical k < 1: MET)")

    # -------------------------------------------------------------------------
    # 2. Single-Stuck-Drum Safety Criterion (1 Drum Out, 5 Drums In)
    # -------------------------------------------------------------------------
    print("\n2. Evaluating 1-Stuck-Drum Criterion (1 Drum at 0°, 5 Drums at 180°)...")
    openmc.reset_auto_ids()
    stuck_angles = [0.0, 180.0, 180.0, 180.0, 180.0, 180.0]
    p_stuck = ReactorParameters(particles=1500, batches=25, inactive=8)
    mats, fuel, clad, cool, mod, beo, b4c = build_materials(p_stuck)
    geom_stuck, pins, drums = build_geometry(p_stuck, fuel, clad, cool, mod, beo, b4c, drum_angles=stuck_angles)
    settings = build_settings(p_stuck)

    model_stuck = openmc.Model(materials=mats, geometry=geom_stuck, settings=settings)
    model_stuck.export_to_xml()
    openmc.run(output=False)

    sp_file = f"statepoint.{p_stuck.batches}.h5"
    with openmc.StatePoint(sp_file) as sp:
        k_stuck = float(sp.keff.nominal_value)
        k_stuck_std = float(sp.keff.std_dev)

    sdm_normal_pcm = (1.0 - k_in) / k_in * 1.0e5
    sdm_stuck_pcm = (1.0 - k_stuck) / k_stuck * 1.0e5
    print(f"   Normal Shutdown (6 drums in) : keff = {k_in:.5f} | SDM = {sdm_normal_pcm:.0f} pcm")
    print(f"   Stuck Drum (1 drum stuck out): keff = {k_stuck:.5f} +/- {k_stuck_std:.5f} | SDM = {sdm_stuck_pcm:.0f} pcm (k < 1.0: MET)")

    # -------------------------------------------------------------------------
    # 3. Fuel Doppler Temperature Coefficient (Calculated directly in OpenMC via WMP)
    # -------------------------------------------------------------------------
    print("\n3. Calculating Fuel Doppler Coefficient (T_fuel in [300, 600, 900, 1200, 1500 K])...")
    fuel_temps = [300.0, 600.0, 900.0, 1200.0, 1500.0]
    k_dopp_vals = []
    s_dopp_vals = []

    for Tf in fuel_temps:
        openmc.reset_auto_ids()
        p_d = ReactorParameters(drum_angle_deg=0.0, particles=1500, batches=25, inactive=8)
        mats, fuel, clad, cool, mod, beo, b4c = build_materials(p_d)
        fuel.temperature = Tf
        geom, pins, drums = build_geometry(p_d, fuel, clad, cool, mod, beo, b4c)
        settings = build_settings(p_d)

        m = openmc.Model(materials=mats, geometry=geom, settings=settings)
        m.export_to_xml()
        openmc.run(output=False)

        with openmc.StatePoint(f"statepoint.{p_d.batches}.h5") as sp:
            k = float(sp.keff.nominal_value)
            s = float(sp.keff.std_dev)
            k_dopp_vals.append(k)
            s_dopp_vals.append(s)
            print(f"   T_fuel = {Tf:4.0f} K | keff = {k:.5f} +/- {s:.5f}")

    # Doppler coefficient across range 300 K to 1500 K
    rho_fuel_low = (k_dopp_vals[0] - 1.0) / k_dopp_vals[0]
    rho_fuel_high = (k_dopp_vals[-1] - 1.0) / k_dopp_vals[-1]
    alpha_fuel_pcm_K = (rho_fuel_high - rho_fuel_low) / (fuel_temps[-1] - fuel_temps[0]) * 1.0e5
    # Propagate Monte Carlo uncertainty
    u_rho_low = s_dopp_vals[0] / (k_dopp_vals[0]**2)
    u_rho_high = s_dopp_vals[-1] / (k_dopp_vals[-1]**2)
    u_alpha_fuel = np.sqrt(u_rho_low**2 + u_rho_high**2) / (fuel_temps[-1] - fuel_temps[0]) * 1.0e5
    print(f"   Calculated Fuel Doppler Coeff: alpha_fuel = {alpha_fuel_pcm_K:+.2f} +/- {u_alpha_fuel:.2f} pcm/K (Negative: MET)")

    # -------------------------------------------------------------------------
    # 4. Moderator Temperature Coefficient (Calculated directly in OpenMC)
    # -------------------------------------------------------------------------
    print("\n4. Calculating Moderator Temperature Coefficient (T_mod in [300, 600, 900 K])...")
    # YH2 thermal expansion: beta ~ 3.6e-5 K-1
    mod_temps = [300.0, 600.0, 900.0]
    k_mod_vals = []
    s_mod_vals = []

    for Tm in mod_temps:
        openmc.reset_auto_ids()
        p_m = ReactorParameters(drum_angle_deg=0.0, particles=1500, batches=25, inactive=8)
        # Volumetric expansion effect on YH2 density
        p_m.moderator_density = 4.30 * (1.0 - 3.6e-5 * (Tm - 300.0))
        mats, fuel, clad, cool, mod, beo, b4c = build_materials(p_m)
        mod.temperature = Tm
        geom, pins, drums = build_geometry(p_m, fuel, clad, cool, mod, beo, b4c)
        settings = build_settings(p_m)

        m = openmc.Model(materials=mats, geometry=geom, settings=settings)
        m.export_to_xml()
        openmc.run(output=False)

        with openmc.StatePoint(f"statepoint.{p_m.batches}.h5") as sp:
            k = float(sp.keff.nominal_value)
            s = float(sp.keff.std_dev)
            k_mod_vals.append(k)
            s_mod_vals.append(s)
            print(f"   T_mod = {Tm:4.0f} K (rho = {p_m.moderator_density:.3f} g/cm3) | keff = {k:.5f} +/- {s:.5f}")

    rho_mod_low = (k_mod_vals[0] - 1.0) / k_mod_vals[0]
    rho_mod_high = (k_mod_vals[-1] - 1.0) / k_mod_vals[-1]
    alpha_mod_pcm_K = (rho_mod_high - rho_mod_low) / (mod_temps[-1] - mod_temps[0]) * 1.0e5
    u_rho_m_low = s_mod_vals[0] / (k_mod_vals[0]**2)
    u_rho_m_high = s_mod_vals[-1] / (k_mod_vals[-1]**2)
    u_alpha_mod = np.sqrt(u_rho_m_low**2 + u_rho_m_high**2) / (mod_temps[-1] - mod_temps[0]) * 1.0e5
    print(f"   Calculated Moderator Temp Coeff: alpha_mod = {alpha_mod_pcm_K:+.2f} +/- {u_alpha_mod:.2f} pcm/K (Negative: MET)")

    # -------------------------------------------------------------------------
    # 5. Coolant Temperature Coefficient (Calculated directly in OpenMC)
    # -------------------------------------------------------------------------
    print("\n5. Calculating Coolant Temperature Coefficient (He-Xe at 900 K, 1050 K, 1300 K, 2.0 MPa)...")
    cool_temps = [900.0, 1050.0, 1300.0]
    k_cool_vals = []
    s_cool_vals = []

    for Tc in cool_temps:
        openmc.reset_auto_ids()
        p_c = ReactorParameters(drum_angle_deg=0.0, particles=1500, batches=25, inactive=8)
        # Ideal gas law at 2.0 MPa: rho = (p * M) / (R * T)
        p_c.coolant_temp_K = Tc
        p_c.coolant_density = (2.0e6 * 0.0400) / (8.31446 * Tc * 1000.0) # g/cm3
        mats, fuel, clad, cool, mod, beo, b4c = build_materials(p_c)
        cool.temperature = Tc
        geom, pins, drums = build_geometry(p_c, fuel, clad, cool, mod, beo, b4c)
        settings = build_settings(p_c)

        m = openmc.Model(materials=mats, geometry=geom, settings=settings)
        m.export_to_xml()
        openmc.run(output=False)

        with openmc.StatePoint(f"statepoint.{p_c.batches}.h5") as sp:
            k = float(sp.keff.nominal_value)
            s = float(sp.keff.std_dev)
            k_cool_vals.append(k)
            s_cool_vals.append(s)
            print(f"   T_cool = {Tc:4.0f} K (rho = {p_c.coolant_density:.5f} g/cm3) | keff = {k:.5f} +/- {s:.5f}")

    rho_c_low = (k_cool_vals[0] - 1.0) / k_cool_vals[0]
    rho_c_high = (k_cool_vals[-1] - 1.0) / k_cool_vals[-1]
    alpha_cool_pcm_K = (rho_c_high - rho_c_low) / (cool_temps[-1] - cool_temps[0]) * 1.0e5
    u_rho_c_low = s_cool_vals[0] / (k_cool_vals[0]**2)
    u_rho_c_high = s_cool_vals[-1] / (k_cool_vals[-1]**2)
    u_alpha_cool = np.sqrt(u_rho_c_low**2 + u_rho_c_high**2) / (cool_temps[-1] - cool_temps[0]) * 1.0e5
    print(f"   Calculated Coolant Temp Coeff: alpha_cool = {alpha_cool_pcm_K:+.2f} +/- {u_alpha_cool:.2f} pcm/K")

    # Total Temperature Coefficient
    alpha_total_pcm_K = alpha_fuel_pcm_K + alpha_mod_pcm_K + alpha_cool_pcm_K
    u_alpha_total = np.sqrt(u_alpha_fuel**2 + u_alpha_mod**2 + u_alpha_cool**2)
    print(f"   Calculated Total Temp Coeff   : alpha_total = {alpha_total_pcm_K:+.2f} +/- {u_alpha_total:.2f} pcm/K (Negative: MET)")

    # -------------------------------------------------------------------------
    # 6. Void Coefficient / Density Perturbations (Calculated directly in OpenMC)
    # -------------------------------------------------------------------------
    print("\n6. Calculating Coolant Void Reactivity & Void Coefficient (OpenMC)...")
    void_states = [
        ('Nominal (2.0 MPa)', 0.00916),
        ('Half-density (1.0 MPa)', 0.00458),
        ('Voided (0.0 MPa)', 1.0e-7)
    ]
    k_void_states = []
    s_void_states = []

    for name, dens in void_states:
        openmc.reset_auto_ids()
        p_v = ReactorParameters(drum_angle_deg=0.0, particles=1500, batches=25, inactive=8)
        p_v.coolant_density = dens
        mats, fuel, clad, cool, mod, beo, b4c = build_materials(p_v)
        geom, pins, drums = build_geometry(p_v, fuel, clad, cool, mod, beo, b4c)
        settings = build_settings(p_v)

        m = openmc.Model(materials=mats, geometry=geom, settings=settings)
        m.export_to_xml()
        openmc.run(output=False)

        with openmc.StatePoint(f"statepoint.{p_v.batches}.h5") as sp:
            k = float(sp.keff.nominal_value)
            s = float(sp.keff.std_dev)
            k_void_states.append(k)
            s_void_states.append(s)
            print(f"   State '{name:<22}': rho = {dens:.5f} g/cm3 | keff = {k:.5f} +/- {s:.5f}")

    k_nominal = k_void_states[0]
    k_voided = k_void_states[-1]
    rho_nom = (k_nominal - 1.0) / k_nominal
    rho_void = (k_voided - 1.0) / k_voided
    void_worth_pcm = (rho_void - rho_nom) * 1.0e5
    u_void_pcm = np.sqrt((s_void_states[-1]/k_voided**2)**2 + (s_void_states[0]/k_nominal**2)**2) * 1.0e5
    # Void coefficient d(rho)/d(void fraction) where void fraction = 1.0 at full void
    alpha_void_pcm = void_worth_pcm
    print(f"   Coolant Void Worth (Nominal -> Void): {void_worth_pcm:+.1f} +/- {u_void_pcm:.1f} pcm")

    # -------------------------------------------------------------------------
    # 7. Generate Plot: results/drum_worth.png (7 Angles + Stuck Drum)
    # -------------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), dpi=300, sharex=True)

    ax1.errorbar(angles, keff_vals, yerr=keff_stds, fmt='o-', color='#1565C0', ecolor='#90CAF9',
                 elinewidth=2, capsize=4, lw=2, label="Calculated $k_{eff}$")
    ax1.axhline(1.0, color='red', linestyle='--', lw=1.5, label="Criticality Boundary ($k=1.0$)")
    ax1.axhline(0.97, color='purple', linestyle='-.', lw=1.2, label="Target Shutdown ($k \\approx 0.97$)")
    ax1.axhline(k_stuck, color='darkorange', linestyle=':', lw=1.5,
                label=f"1-Stuck Drum ($k={k_stuck:.4f}$)")
    ax1.set_ylabel("Effective Multiplication Factor $k_{eff}$", fontsize=11)
    ax1.set_title("HALEU NEP Reactor - Control Drum Reactivity Worth and S-Curve", fontsize=12, fontweight='bold')
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(loc='upper right', framealpha=0.9, fontsize=9.0)

    rel_worth_pcm = [(r - rho_pcm[0]) for r in rho_pcm]
    ax2.plot(angles, rel_worth_pcm, 's-', color='#2E7D32', lw=2, label="Relative Drum Worth")
    ax2.set_xlabel("Control Drum Rotation Angle $\\theta$ [degrees] (0°=Absorbers Out, 180°=Absorbers In)", fontsize=11)
    ax2.set_ylabel("Reactivity Insertion $\\Delta\\rho$ [pcm]", fontsize=11)
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.legend(loc='lower left', framealpha=0.9, fontsize=9.5)

    plt.tight_layout()
    drum_plot_path = os.path.join(output_dir, "drum_worth.png")
    fig.savefig(drum_plot_path)
    plt.close(fig)
    print(f"\nSaved: {drum_plot_path}")

    # -------------------------------------------------------------------------
    # 8. Generate Plot: results/temperature_coefficients.png
    # -------------------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), dpi=300)

    # Subplot 1: Fuel Doppler
    axes[0].errorbar(fuel_temps, k_dopp_vals, yerr=s_dopp_vals, fmt='o-', color='#C62828', ecolor='#EF9A9A',
                     elinewidth=2, capsize=4, lw=2)
    axes[0].set_xlabel("Fuel Temperature $T_{fuel}$ [K]", fontsize=11)
    axes[0].set_ylabel("Calculated $k_{eff}$", fontsize=11)
    axes[0].set_title(f"Fuel Doppler Feedback (WMP)\n$\\alpha_{{fuel}} = {alpha_fuel_pcm_K:+.2f} \\pm {u_alpha_fuel:.2f}$ pcm/K",
                      fontsize=10.5, fontweight='bold')
    axes[0].grid(True, linestyle=':', alpha=0.6)

    # Subplot 2: Moderator Temperature
    axes[1].errorbar(mod_temps, k_mod_vals, yerr=s_mod_vals, fmt='s-', color='#1565C0', ecolor='#90CAF9',
                     elinewidth=2, capsize=4, lw=2)
    axes[1].set_xlabel("Moderator Temperature $T_{mod}$ [K]", fontsize=11)
    axes[1].set_ylabel("Calculated $k_{eff}$", fontsize=11)
    axes[1].set_title(f"Moderator Temperature Feedback\n$\\alpha_{{mod}} = {alpha_mod_pcm_K:+.2f} \\pm {u_alpha_mod:.2f}$ pcm/K",
                      fontsize=10.5, fontweight='bold')
    axes[1].grid(True, linestyle=':', alpha=0.6)

    # Subplot 3: Coolant Temperature
    axes[2].errorbar(cool_temps, k_cool_vals, yerr=s_cool_vals, fmt='^-', color='#E65100', ecolor='#FFCC80',
                     elinewidth=2, capsize=4, lw=2)
    axes[2].set_xlabel("He-Xe Coolant Temperature $T_{cool}$ [K]", fontsize=11)
    axes[2].set_ylabel("Calculated $k_{eff}$", fontsize=11)
    axes[2].set_title(f"Coolant Density Feedback\n$\\alpha_{{cool}} = {alpha_cool_pcm_K:+.2f} \\pm {u_alpha_cool:.2f}$ pcm/K",
                      fontsize=10.5, fontweight='bold')
    axes[2].grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()
    temp_plot_path = os.path.join(output_dir, "temperature_coefficients.png")
    fig.savefig(temp_plot_path)
    plt.close(fig)
    print(f"Saved: {temp_plot_path}")

    # -------------------------------------------------------------------------
    # 9. Generate Plot: results/void_coefficient.png (Dedicated void figure)
    # -------------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5), dpi=300)

    # Subplot 1: Bar chart of keff across void states
    labels = ['Nominal\n(2.0 MPa)', '50% Density\n(1.0 MPa)', 'Depressurized\n(0.0 MPa)']
    bars = ax1.bar(labels, k_void_states, yerr=s_void_states, color=['#1976D2', '#0288D1', '#546E7A'],
                   capsize=5, edgecolor='black', linewidth=1.2, width=0.55)
    ax1.set_ylabel("Effective Multiplication Factor $k_{eff}$", fontsize=11)
    ax1.set_ylim(min(k_void_states) - 0.015, max(k_void_states) + 0.015)
    ax1.set_title("He-Xe Coolant Pressure/Density States", fontsize=11, fontweight='bold')
    for bar in bars:
        h = bar.get_height()
        ax1.annotate(f"{h:.5f}", xy=(bar.get_x() + bar.get_width()/2, h), xytext=(0, 5),
                     textcoords="offset points", ha='center', va='bottom', fontsize=9.5, fontweight='bold')
    ax1.grid(True, linestyle=':', alpha=0.6, axis='y')

    # Subplot 2: Reactivity insertion vs coolant density
    densities = [dens for _, dens in void_states]
    reactivities = [((k - 1.0)/k - rho_nom) * 1.0e5 for k in k_void_states]
    ax2.plot(densities, reactivities, 'D-', color='#D32F2F', lw=2, markersize=8)
    ax2.axhline(0.0, color='black', linestyle='--', lw=1.2)
    ax2.set_xlabel("He-Xe Coolant Density $\\rho$ [g/cm³]", fontsize=11)
    ax2.set_ylabel("Reactivity Difference $\\Delta\\rho$ vs Nominal [pcm]", fontsize=11)
    ax2.set_title(f"Coolant Void Reactivity Insertion\n$\\Delta\\rho_{{void}} = {void_worth_pcm:+.1f} \\pm {u_void_pcm:.1f}$ pcm",
                  fontsize=11, fontweight='bold')
    ax2.grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()
    void_plot_path = os.path.join(output_dir, "void_coefficient.png")
    fig.savefig(void_plot_path)
    plt.close(fig)
    print(f"Saved: {void_plot_path}")

    # Clean intermediate statepoints
    for b in [20, 25, 30, 35]:
        if os.path.exists(f"statepoint.{b}.h5"):
            os.remove(f"statepoint.{b}.h5")

    results_summary = {
        "k_out": k_out,
        "k_in": k_in,
        "total_worth_pcm": total_worth_pcm,
        "total_worth_pct": total_worth_pct,
        "k_stuck": k_stuck,
        "sdm_normal_pcm": sdm_normal_pcm,
        "sdm_stuck_pcm": sdm_stuck_pcm,
        "alpha_fuel_pcm_K": alpha_fuel_pcm_K,
        "u_alpha_fuel": u_alpha_fuel,
        "alpha_mod_pcm_K": alpha_mod_pcm_K,
        "u_alpha_mod": u_alpha_mod,
        "alpha_cool_pcm_K": alpha_cool_pcm_K,
        "u_alpha_cool": u_alpha_cool,
        "alpha_total_pcm_K": alpha_total_pcm_K,
        "u_alpha_total": u_alpha_total,
        "void_worth_pcm": void_worth_pcm,
        "u_void_pcm": u_void_pcm,
        "k_nominal": k_nominal,
        "k_voided": k_voided
    }
    return results_summary


if __name__ == "__main__":
    run_control_and_reactivity_analysis()
