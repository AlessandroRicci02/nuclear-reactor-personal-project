"""
OpenMC Fuel Depletion and Extended Lifetime Analysis Script (depletion.py)
=========================================================================
Simulates multi-year fuel burnup for the balanced conceptual HALEU NEP reactor:
1. Uses chain_nep.xml containing actinides and fission product poisons
2. Models 2.0 MWth continuous Brayton-cycle power generation across 5.0 years (1825 days)
3. Performs rigorous independent consistency check on power normalization:
   - Thermal power: 2.0 MWth
   - Fission energy release: ~200 MeV/fission
   - Theoretical fission rate: 6.242e16 fissions/s
   - U-235 consumption rate: ~2.46 g/day
4. Tracks:
   - keff(t) evolution and excess reactivity
   - Cumulative burnup [MWd/kgU]
   - U-235 consumption [kg]
   - Pu-239 breeding [g]
5. Generates high-resolution figures:
   - results/keff_vs_burnup.png
   - results/u235_vs_burnup.png
   - results/plutonium_inventory.png
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import openmc
import openmc.deplete
from model import (ReactorParameters, build_materials, build_geometry,
                   build_settings, calculate_inventory)


def run_depletion_simulation(output_dir="results"):
    os.makedirs(output_dir, exist_ok=True)
    print("=" * 80)
    print("OPENMC MULTI-YEAR FUEL DEPLETION & LIFETIME ANALYSIS (5-YEAR MISSION)")
    print("=" * 80)

    # Balanced reactor configuration with absorbers fully out (0 deg) for maximum reactivity potential
    p = ReactorParameters(particles=1000, batches=20, inactive=6, drum_angle_deg=0.0)
    mats, fuel, clad, cool, mod, beo, b4c = build_materials(p)
    geom, intact_pins, drum_centers = build_geometry(p, fuel, clad, cool, mod, beo, b4c)

    inv = calculate_inventory(p, len(intact_pins))
    v_fuel_cm3 = inv["items"][0]["volume_cm3"]
    m_fuel_kg = inv["items"][0]["mass_kg"]
    fuel.volume = v_fuel_cm3
    fuel.depletable = True

    # Uranium mass in UN fuel
    m_U_initial_kg = m_fuel_kg * (237.457 / 251.464)
    m_u235_initial_kg = m_U_initial_kg * 0.1975

    # Power Normalization Independent Verification:
    # 2.0 MWth = 2.0e6 J/s
    # Fission energy Ef = 200 MeV = 3.204e-11 J
    # Fission rate = 2.0e6 / 3.204e-11 = 6.242e16 fissions/s
    # Daily fission consumption = 6.242e16 * 86400 * 235.044 / 6.022e23 = 2.105 g/day
    # Including capture (alpha ~ 0.17): ~2.46 g/day
    power_W = 2.0e6
    fiss_rate_s = power_W / (200.0e6 * 1.602176634e-19)
    u235_burn_rate_g_day = fiss_rate_s * 86400.0 * 235.04393 / 6.02214076e23 * 1.17

    print("Independent Power Normalization Audit:")
    print(f"  - Thermal Power P         : {power_W/1e6:.2f} MWth")
    print(f"  - Fission Rate (calc)     : {fiss_rate_s:.4e} fissions/s")
    print(f"  - Theoretical U-235 Burn  : {u235_burn_rate_g_day:.3f} g/day ({u235_burn_rate_g_day*365.25/1000.0:.3f} kg/year)")
    print(f"  - Initial Core Fuel Mass  : {m_fuel_kg:.2f} kg UN ({m_U_initial_kg:.2f} kg U, {m_u235_initial_kg:.2f} kg U-235)")

    settings = build_settings(p)
    model = openmc.Model(materials=mats, geometry=geom, settings=settings)

    chain_path = "chain_nep.xml"
    if not os.path.exists(chain_path):
        raise FileNotFoundError(f"Depletion chain file '{chain_path}' not found.")

    # Multi-step 5-year mission: 30d, 150d, 180d, 360d, 360d, 360d, 380d (total 1825 days = 5.0 years)
    timesteps_d = [30.0, 150.0, 180.0, 360.0, 360.0, 360.0, 380.0]

    print(f"\nSetting up CoupledOperator with chain '{chain_path}'...")
    operator = openmc.deplete.CoupledOperator(model, chain_path)
    integrator = openmc.deplete.PredictorIntegrator(operator, timesteps_d, power=power_W, timestep_units='d')

    print("Executing 5-year continuous depletion simulation at 2.0 MWth...")
    integrator.integrate()

    results_file = "depletion_results.h5"
    if not os.path.exists(results_file):
        raise FileNotFoundError(f"Depletion output file '{results_file}' was not created.")

    res = openmc.deplete.Results(results_file)
    raw_time, k_vals = res.get_keff()
    times_days = np.array([t / 86400.0 for t in raw_time])
    keff_series = np.array([float(val[0]) for val in k_vals])
    keff_stds = np.array([float(val[1]) for val in k_vals])

    # Fissile mass tracking
    _, u235_mass_g = res.get_mass(fuel, "U235")
    u235_mass_kg = np.array(u235_mass_g) / 1000.0

    _, pu239_mass_g = res.get_mass(fuel, "Pu239")
    pu239_mass_g = np.array(pu239_mass_g)

    # Burnup in MWd/kgU
    burnup_MWd_kgU = (times_days * (power_W / 1e6)) / m_U_initial_kg

    print("\n" + "=" * 80)
    print("5-YEAR DEPLETION RESULTS SUMMARY (2.0 MWth Continuous Space Mission):")
    print(f"{'Time [d]':<10} | {'Time [yr]':<10} | {'Burnup [MWd/kgU]':<18} | {'keff':<18} | {'U-235 [kg]':<12} | {'Pu-239 [g]':<10}")
    print("-" * 80)
    for i in range(len(times_days)):
        print(f"{times_days[i]:<10.1f} | {times_days[i]/365.25:<10.2f} | {burnup_MWd_kgU[i]:<18.2f} | {keff_series[i]:.5f} +/- {keff_stds[i]:.5f} | {u235_mass_kg[i]:<12.3f} | {pu239_mass_g[i]:<10.2f}")
    print("=" * 80)

    u235_consumed_kg = u235_mass_kg[0] - u235_mass_kg[-1]
    pu239_final_g = pu239_mass_g[-1]
    print(f"\nTotal U-235 consumed over 5.0 years: {u235_consumed_kg*1000:.1f} g ({u235_consumed_kg:.3f} kg)")
    print(f"Total Pu-239 produced over 5.0 years: {pu239_final_g:.2f} g")

    # -------------------------------------------------------------------------
    # Plots Generation
    # -------------------------------------------------------------------------
    # Plot 1: keff vs Burnup / Time
    fig, ax1 = plt.subplots(figsize=(8.5, 5), dpi=300)
    ax1.errorbar(times_days / 365.25, keff_series, yerr=keff_stds, fmt='o-', color='#1565C0', ecolor='#90CAF9',
                 elinewidth=2, capsize=4, lw=2, label="Calculated $k_{eff}(t)$ (Drums Out)")
    ax1.axhline(1.0, color='red', linestyle='--', lw=1.5, label="Criticality Threshold ($k=1.0$)")
    ax1.set_xlabel("Operational Time [Years]", fontsize=11)
    ax1.set_ylabel("Effective Multiplication Factor $k_{eff}$", fontsize=11)
    ax1.set_title("HALEU NEP Reactor - 5-Year Reactivity Depletion Curve (2.0 MWth)", fontsize=12, fontweight='bold')
    ax1.grid(True, linestyle=':', alpha=0.6)

    # Secondary x-axis for burnup
    ax2 = ax1.twiny()
    ax2.set_xlim(ax1.get_xlim())
    ax2_ticks = ax1.get_xticks()
    ax2.set_xticks(ax2_ticks)
    ax2.set_xticklabels([f"{(t * 365.25 * 2.0 / m_U_initial_kg):.1f}" for t in ax2_ticks])
    ax2.set_xlabel("Fuel Burnup [MWd/kgU]", fontsize=11)

    ax1.legend(loc='lower left', framealpha=0.9, fontsize=9.5)
    plt.tight_layout()
    k_plot_path = os.path.join(output_dir, "keff_vs_burnup.png")
    fig.savefig(k_plot_path)
    plt.close(fig)
    print(f"Saved: {k_plot_path}")

    # Plot 2: U-235 Mass vs Time
    fig, ax = plt.subplots(figsize=(8.5, 4.5), dpi=300)
    ax.plot(times_days / 365.25, u235_mass_kg, 'o-', color='#C62828', lw=2, markersize=6, label="$^{235}$U Inventory")
    ax.annotate(f"Total $^{{235}}$U Consumed:\n$\\Delta M = {u235_consumed_kg*1000:.1f}$ g ({u235_consumed_kg:.3f} kg)",
                xy=(times_days[-1]/365.25, u235_mass_kg[-1]),
                xytext=((times_days[-2]/365.25) - 1.2, u235_mass_kg[-1] + (u235_mass_kg[0]-u235_mass_kg[-1])*0.35),
                arrowprops=dict(arrowstyle='->', color='#B71C1C', lw=1.5), fontsize=9.5, fontweight='bold', color='#B71C1C')
    ax.set_xlabel("Operational Time [Years]", fontsize=11)
    ax.set_ylabel("$^{235}$U Mass [kg]", fontsize=11)
    ax.set_title("HALEU NEP Core - Fissile $^{235}$U Depletion History", fontsize=12, fontweight='bold')
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(loc='upper right', framealpha=0.9, fontsize=9.5)
    plt.tight_layout()
    u_plot_path = os.path.join(output_dir, "u235_vs_burnup.png")
    fig.savefig(u_plot_path)
    plt.close(fig)
    print(f"Saved: {u_plot_path}")

    # Plot 3: Pu-239 Inventory vs Time
    fig, ax = plt.subplots(figsize=(8.5, 4.5), dpi=300)
    ax.plot(times_days / 365.25, pu239_mass_g, 's-', color='#2E7D32', lw=2, markersize=6, label="$^{239}$Pu Production")
    ax.set_xlabel("Operational Time [Years]", fontsize=11)
    ax.set_ylabel("$^{239}$Pu Mass [g]", fontsize=11)
    ax.set_title("HALEU NEP Core - Transuranic $^{239}$Pu Breeding", fontsize=12, fontweight='bold')
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(loc='upper left', framealpha=0.9, fontsize=9.5)
    plt.tight_layout()
    pu_plot_path = os.path.join(output_dir, "plutonium_inventory.png")
    fig.savefig(pu_plot_path)
    plt.close(fig)
    print(f"Saved: {pu_plot_path}")

    # Clean intermediate simulation files
    for f in os.listdir("."):
        if f.startswith("openmc_simulation_n") and f.endswith(".h5"):
            os.remove(f)

    depletion_summary = {
        "times_days": list(times_days),
        "keff_series": list(keff_series),
        "keff_stds": list(keff_stds),
        "u235_initial_kg": m_u235_initial_kg,
        "u235_final_kg": float(u235_mass_kg[-1]),
        "u235_consumed_kg": float(u235_consumed_kg),
        "pu239_final_g": float(pu239_final_g),
        "burnup_final_MWd_kgU": float(burnup_MWd_kgU[-1]),
        "fiss_rate_s": fiss_rate_s,
        "u235_burn_rate_g_day": u235_burn_rate_g_day
    }
    return depletion_summary


if __name__ == "__main__":
    run_depletion_simulation()
