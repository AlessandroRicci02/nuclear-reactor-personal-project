"""
Criticality-Controlled Depletion & Theoretical Lifetime Runner for Opt-4 at 400 kWth
===================================================================================
Simulates the realistic operating behavior of the Opt-4 NEP Space Reactor:
- System Thermal Power: 400 kWth (100 kWe electrical bus at ~25% Brayton efficiency)
- Criticality condition: At each burnup step, determines the required drum angle theta_crit(t)
  such that k_eff ≈ 1.000 ± 0.002.
- Burnup trajectory: As fuel burns, drums rotate outward (from theta ≈ 100°-120° toward 0°).
- Evaluates:
    1. Operational state (drums at theta_crit, k ≈ 1.000)
    2. Uncontrolled state (drums at 0.0°, k_out)
    3. Normal shutdown state (drums at 180.0°, k_in, SDM)
    4. Single stuck-drum state (1 drum at 0°, 7 at 180°, k_stuck, SDM_stuck)
    5. Fissile U-235 consumption, transuranic Pu-239 breeding, burnup (MWd/kgU)
- Lifetime: Evaluates the 10-year mandatory milestone AND continues beyond 10 years
  until k_out < 1.000 to determine the true Theoretical Neutronic Lifetime.
- Exports publication plots and CSV to results/opt4/depletion/
"""

import os
import sys
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import openmc
import openmc.deplete

os.environ["OPENMC_CROSS_SECTIONS"] = "/home/aleri/wmp_data/cross_sections_wmp.xml"
if "/home/aleri/miniforge3/envs/openmc-env/bin" not in os.environ.get("PATH", ""):
    os.environ["PATH"] = f"/home/aleri/miniforge3/envs/openmc-env/bin:{os.environ.get('PATH', '')}"

from model import ReactorParameters, build_materials, build_geometry, build_settings, calculate_inventory

def update_fuel_in_place(target_fuel, source_fuel):
    """Updates target_fuel with nuclides and density from source_fuel in-place."""
    target_fuel.nuclides.clear()
    for nuc in source_fuel.nuclides:
        target_fuel.add_nuclide(nuc.name, nuc.percent, nuc.percent_type)
    target_fuel.set_density('g/cm3', source_fuel.get_mass_density())

def run_criticality_controlled_depletion_opt4(output_dir="results/opt4/depletion"):
    os.makedirs(output_dir, exist_ok=True)
    print("=" * 80)
    print("CRITICALITY-CONTROLLED DEPLETION AT 400 kWth - OPT-4 SPACE REACTOR")
    print("=" * 80)

    # Opt-4 Reference Reactor Parameters
    p = ReactorParameters(
        fuel_radius=0.58,
        cladding_thickness=0.05,
        coolant_channel_thickness=0.05,
        moderator_fuel_pitch=1.635,
        core_radius=17.35,
        active_core_height=45.60,
        radial_reflector_thickness=10.00,
        axial_reflector_thickness=8.10,
        number_of_control_drums=8,
        control_drum_radius=4.58,
        absorber_thickness=1.45,
        absorber_sector_angular_width=140.0,
        drum_clearance=0.20,
        b10_enrichment=0.90,
        fuel_enrichment=19.75,
        drum_angle_deg=0.0,
        coolant_state='operating',
        particles=4000,
        batches=30,
        inactive=8
    )

    power_W = 400.0e3  # 400 kWth

    # Timestep schedule:
    # 0 -> 30d (0.08y) -> 180d (0.49y) -> 365d (1.0y) -> 730d (2.0y) -> 1460d (4.0y) ->
    # 2190d (6.0y) -> 2920d (8.0y) -> 3650d (10.0y) -> 4380d (12.0y) -> 5475d (15.0y) ->
    # 7300d (20.0y) -> 9125d (25.0y) -> 10950d (30.0y)
    timesteps_d = [
        30.0,    # to 30 d (0.08 y)
        150.0,   # to 180 d (0.49 y)
        185.0,   # to 365 d (1.00 y)
        365.0,   # to 730 d (2.00 y)
        730.0,   # to 1460 d (4.00 y)
        730.0,   # to 2190 d (6.00 y)
        730.0,   # to 2920 d (8.00 y)
        730.0,   # to 3650 d (10.00 y) [MANDATORY 10-YR REQUIREMENT]
        730.0,   # to 4380 d (12.00 y)
        1095.0,  # to 5475 d (15.00 y)
        1825.0,  # to 7300 d (20.00 y)
        1825.0,  # to 9125 d (25.00 y)
        1825.0   # to 10950 d (30.00 y)
    ]

    # Initialize fresh fuel
    openmc.reset_auto_ids()
    mats_init, fuel_init, clad_init, cool_init, mod_init, beo_init, b4c_init = build_materials(p)
    geom_init, pins_init, drums_init = build_geometry(p, fuel_init, clad_init, cool_init, mod_init, beo_init, b4c_init)
    inv = calculate_inventory(p, len(pins_init))
    v_fuel = inv["items"][0]["volume_cm3"]
    m_fuel_kg = inv["items"][0]["mass_kg"]
    m_U_kg = m_fuel_kg * (237.457 / 251.464)
    initial_u235_kg = m_U_kg * 0.1975

    current_fuel = fuel_init
    current_fuel.volume = v_fuel
    current_fuel.depletable = True

    history = []
    cumulative_days = 0.0
    cumulative_burnup = 0.0
    theoretical_lifetime_years = None

    for step_idx in range(len(timesteps_d) + 1):
        years = cumulative_days / 365.25

        # ---------------------------------------------------------------------
        # A. Evaluate Uncontrolled Excess Reactivity (Drums at 0.0° Outward)
        # ---------------------------------------------------------------------
        openmc.reset_auto_ids()
        p.drum_angle_deg = 0.0
        mats_out, f_out, cl_out, co_out, m_out, be_out, b_out = build_materials(p)
        update_fuel_in_place(f_out, current_fuel)
        f_out.volume = v_fuel
        geom_out, _, _ = build_geometry(p, f_out, cl_out, co_out, m_out, be_out, b_out)
        s_out = build_settings(p)
        m_out_model = openmc.Model(materials=mats_out, geometry=geom_out, settings=s_out)
        m_out_model.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
            k_out = float(sp.keff.nominal_value)
            s_out_val = float(sp.keff.std_dev)

        # ---------------------------------------------------------------------
        # B. Evaluate Shutdown Multiplication Factor (Drums at 180.0° Inward)
        # ---------------------------------------------------------------------
        openmc.reset_auto_ids()
        p.drum_angle_deg = 180.0
        mats_in, f_in, cl_in, co_in, m_in, be_in, b_in = build_materials(p)
        update_fuel_in_place(f_in, current_fuel)
        f_in.volume = v_fuel
        geom_in, _, _ = build_geometry(p, f_in, cl_in, co_in, m_in, be_in, b_in)
        s_in = build_settings(p)
        m_in_model = openmc.Model(materials=mats_in, geometry=geom_in, settings=s_in)
        m_in_model.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
            k_in = float(sp.keff.nominal_value)
            s_in_val = float(sp.keff.std_dev)

        drum_worth_pcm = (1.0 / k_in - 1.0 / k_out) * 1e5
        sdm_normal_pcm = (1.0 - k_in) / k_in * 1e5

        # ---------------------------------------------------------------------
        # C. Single-Stuck-Drum Margin (1 Drum at 0°, 7 Drums at 180°)
        # ---------------------------------------------------------------------
        openmc.reset_auto_ids()
        stuck_angles = [0.0] + [180.0] * 7
        mats_stk, f_stk, cl_stk, co_stk, m_stk, be_stk, b_stk = build_materials(p)
        update_fuel_in_place(f_stk, current_fuel)
        f_stk.volume = v_fuel
        geom_stk, _, _ = build_geometry(p, f_stk, cl_stk, co_stk, m_stk, be_stk, b_stk, drum_angles=stuck_angles)
        s_stk = build_settings(p)
        m_stk_model = openmc.Model(materials=mats_stk, geometry=geom_stk, settings=s_stk)
        m_stk_model.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
            k_stuck = float(sp.keff.nominal_value)
            s_stuck_val = float(sp.keff.std_dev)

        sdm_stuck_pcm = (1.0 - k_stuck) / k_stuck * 1e5

        # ---------------------------------------------------------------------
        # D. Determine Critical Drum Angle theta_crit (to hold k ≈ 1.000)
        # ---------------------------------------------------------------------
        if k_out < 1.000:
            # Cannot maintain criticality even with drums fully outward!
            theta_crit = 0.0
            k_crit = k_out
            s_crit = s_out_val
            if theoretical_lifetime_years is None:
                theoretical_lifetime_years = years
                print(f"\n*** THEORETICAL NEUTRONIC LIFETIME LIMIT REACHED AT {years:.2f} YEARS! ***")
        else:
            # Calibrated analytical S-curve:
            fraction = np.clip((k_out - 1.000) / (k_out - k_in), 0.0, 1.0)
            cos_theta = 1.0 - 2.0 * fraction
            cos_theta = np.clip(cos_theta, -1.0, 1.0)
            theta_crit = float(np.degrees(np.arccos(cos_theta)))

            # Verify with simulation at theta_crit
            openmc.reset_auto_ids()
            p.drum_angle_deg = theta_crit
            mats_cr, f_cr, cl_cr, co_cr, m_cr, be_cr, b_cr = build_materials(p)
            update_fuel_in_place(f_cr, current_fuel)
            f_cr.volume = v_fuel
            geom_cr, _, _ = build_geometry(p, f_cr, cl_cr, co_cr, m_cr, be_cr, b_cr)
            s_cr = build_settings(p)
            m_cr_model = openmc.Model(materials=mats_cr, geometry=geom_cr, settings=s_cr)
            m_cr_model.export_to_xml()
            openmc.run(output=False)
            with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
                k_crit = float(sp.keff.nominal_value)
                s_crit = float(sp.keff.std_dev)

        # ---------------------------------------------------------------------
        # E. Extract Fuel Isotopics
        # ---------------------------------------------------------------------
        u235_kg = current_fuel.get_mass(nuclide="U235") / 1000.0
        pu239_g = current_fuel.get_mass(nuclide="Pu239")
        xe135_g = current_fuel.get_mass(nuclide="Xe135")

        print(f"\n--- STEP {step_idx:2d} | t = {cumulative_days:6.1f} d ({years:5.2f} y) | Burnup = {cumulative_burnup:5.2f} MWd/kgU ---")
        print(f"  Drum Critical Angle theta_crit : {theta_crit:5.1f}° (0°=Out, 180°=In)")
        print(f"  Operational Critical k_eff     : {k_crit:.5f} +/- {s_crit:.5f} (Target ~1.000)")
        print(f"  Uncontrolled Excess k_out      : {k_out:.5f} +/- {s_out_val:.5f} (Excess = {(k_out-1)/k_out*1e5:+.0f} pcm)")
        print(f"  Normal Shutdown k_in           : {k_in:.5f} +/- {s_in_val:.5f} (SDM = {sdm_normal_pcm:4.0f} pcm)")
        print(f"  Single Stuck Drum k_stuck      : {k_stuck:.5f} +/- {s_stuck_val:.5f} (SDM_stk = {sdm_stuck_pcm:4.0f} pcm)")
        print(f"  Control Drum Worth             : {drum_worth_pcm:.0f} pcm ({drum_worth_pcm/1000:.2f}% Delta-k/k)")
        print(f"  Fissile U-235 Remaining        : {u235_kg:.3f} kg (Burned: {initial_u235_kg - u235_kg:.3f} kg)")
        print(f"  Transuranic Pu-239 Bred        : {pu239_g:.1f} g")

        record = {
            "step": step_idx, "time_days": cumulative_days, "time_years": years,
            "burnup_MWd_kgU": cumulative_burnup, "theta_crit_deg": theta_crit,
            "k_crit": k_crit, "s_crit": s_crit,
            "k_out": k_out, "s_out": s_out_val,
            "k_in": k_in, "s_in": s_in_val,
            "k_stuck": k_stuck, "s_stuck": s_stuck_val,
            "sdm_normal_pcm": sdm_normal_pcm, "sdm_stuck_pcm": sdm_stuck_pcm,
            "drum_worth_pcm": drum_worth_pcm,
            "u235_kg": u235_kg, "pu239_g": pu239_g, "xe135_g": xe135_g
        }
        history.append(record)

        # Stop depletion if reactor has fallen subcritical
        if k_out < 1.000 and theoretical_lifetime_years is not None:
            print(f"Stopping depletion search: system is subcritical with drums fully outward.")
            break

        # ---------------------------------------------------------------------
        # F. Advance Depletion Step at Power 400 kWth with Drums at theta_crit
        # ---------------------------------------------------------------------
        if step_idx < len(timesteps_d):
            dt = timesteps_d[step_idx]
            print(f"  >> Depleting {dt:.1f} days at 400 kWth with drums at {theta_crit:.1f}°...")

            openmc.reset_auto_ids()
            p.drum_angle_deg = theta_crit
            mats_dep, f_dep, cl_dep, co_dep, m_dep, be_dep, b_dep = build_materials(p)
            update_fuel_in_place(f_dep, current_fuel)
            f_dep.volume = v_fuel
            f_dep.depletable = True
            geom_dep, _, _ = build_geometry(p, f_dep, cl_dep, co_dep, m_dep, be_dep, b_dep)
            s_dep = build_settings(p)
            m_dep_model = openmc.Model(materials=mats_dep, geometry=geom_dep, settings=s_dep)

            op = openmc.deplete.CoupledOperator(m_dep_model, "data/chain_nep.xml")
            integ = openmc.deplete.PredictorIntegrator(op, [dt], power=power_W, timestep_units='d')
            integ.integrate()

            # Load updated burned material from result
            res = openmc.deplete.Results('depletion_results.h5')
            burned_mats = res.export_to_materials(-1)
            current_fuel = burned_mats[0]
            current_fuel.volume = v_fuel
            current_fuel.depletable = True

            cumulative_days += dt
            cumulative_burnup += (power_W * dt / 1e6) / m_U_kg

    # If reactor remained critical across all steps:
    if theoretical_lifetime_years is None:
        theoretical_lifetime_years = cumulative_days / 365.25
        print(f"\nReactor remains critical beyond maximum simulated mission horizon: > {theoretical_lifetime_years:.1f} years!")

    # -------------------------------------------------------------------------
    # G. Export CSV and Generate Publication Figures
    # -------------------------------------------------------------------------
    csv_path = os.path.join(output_dir, "opt4_criticality_depletion_400kw.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Step", "Time [days]", "Time [years]", "Burnup [MWd/kgU]",
            "Critical Drum Angle [deg]", "k_crit", "k_crit_std",
            "k_out (0 deg)", "k_out_std", "k_in (180 deg)", "k_in_std",
            "k_stuck (1-stuck)", "k_stuck_std", "Normal SDM [pcm]",
            "Stuck SDM [pcm]", "Drum Worth [pcm]", "U-235 Mass [kg]",
            "Pu-239 Mass [g]", "Xe-135 Mass [g]"
        ])
        for h in history:
            writer.writerow([
                h["step"], f"{h['time_days']:.1f}", f"{h['time_years']:.3f}", f"{h['burnup_MWd_kgU']:.3f}",
                f"{h['theta_crit_deg']:.2f}", f"{h['k_crit']:.5f}", f"{h['s_crit']:.5f}",
                f"{h['k_out']:.5f}", f"{h['s_out']:.5f}", f"{h['k_in']:.5f}", f"{h['s_in']:.5f}",
                f"{h['k_stuck']:.5f}", f"{h['s_stuck']:.5f}", f"{h['sdm_normal_pcm']:.0f}",
                f"{h['sdm_stuck_pcm']:.0f}", f"{h['drum_worth_pcm']:.0f}", f"{h['u235_kg']:.4f}",
                f"{h['pu239_g']:.2f}", f"{h['xe135_g']:.4e}"
            ])
    print(f"\nSaved CSV: {csv_path}")

    # Plot 1: Criticality Trajectory & Drum Angle Evolution
    t_yr = [h["time_years"] for h in history]
    k_c = [h["k_crit"] for h in history]
    s_c = [h["s_crit"] for h in history]
    k_o = [h["k_out"] for h in history]
    th_c = [h["theta_crit_deg"] for h in history]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7.5), sharex=True, dpi=300,
                                   gridspec_kw={'height_ratios': [1.8, 1.2]})
    ax1.errorbar(t_yr, k_c, yerr=s_c, fmt='o-', color='#1565C0', lw=2.2, capsize=4,
                 label="Operational Criticality $k_{eff}(\\theta_{crit})$")
    ax1.plot(t_yr, k_o, 's--', color='#D32F2F', lw=1.8, label="Uncontrolled Excess $k_{out}$ (Drums Out at 0°)")
    ax1.axhline(1.0, color='black', linestyle=':', lw=1.5, label="Criticality Target ($k=1.000$)")
    ax1.axvline(10.0, color='green', linestyle='-.', lw=1.5, label="Mandatory 10-Yr Requirement")
    ax1.set_ylabel("Eigenvalue $k_{eff}$", fontsize=11)
    ax1.set_title("HALEU NEP Opt-4 - Criticality-Controlled Depletion Trajectory (400 kWth)", fontsize=12, fontweight='bold')
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(loc='upper right', framealpha=0.9, fontsize=9.5)

    ax2.plot(t_yr, th_c, 'D-', color='#E65100', lw=2.2, label="Control Drum Operating Angle $\\theta_{crit}$")
    ax2.axhline(0.0, color='red', linestyle='--', lw=1.2, label="Full Outward Limit (0° - EOL Exhaustion)")
    ax2.axhline(180.0, color='gray', linestyle=':', lw=1.0, label="Full Inward Limit (180° - Shutdown)")
    ax2.axvline(10.0, color='green', linestyle='-.', lw=1.5)
    ax2.set_xlabel("Continuous Operating Time [Years]", fontsize=11)
    ax2.set_ylabel("Drum Angle $\\theta_{crit}$ [degrees]", fontsize=11)
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.legend(loc='center right', framealpha=0.9, fontsize=9.5)
    fig.tight_layout()
    plot1_path = os.path.join(output_dir, "opt4_criticality_controlled_depletion.png")
    fig.savefig(plot1_path)
    plt.close(fig)
    print(f"Saved: {plot1_path}")

    # Plot 2: Safety Margins Along Depletion
    sdm_n = [h["sdm_normal_pcm"] for h in history]
    sdm_s = [h["sdm_stuck_pcm"] for h in history]

    fig, ax = plt.subplots(figsize=(8.5, 4.8), dpi=300)
    ax.plot(t_yr, sdm_n, 's-', color='#2E7D32', lw=2.2, label="Normal Shutdown Margin (8 Drums In)")
    ax.plot(t_yr, sdm_s, 'o-', color='#C2185B', lw=2.2, label="1-Stuck-Drum Margin (1 Drum Stuck Out)")
    ax.axhline(0.0, color='red', linestyle='--', lw=1.5, label="Subcriticality Limit (SDM = 0 pcm)")
    ax.axvline(10.0, color='green', linestyle='-.', lw=1.5, label="10-Yr Requirement")
    ax.set_xlabel("Continuous Operating Time [Years]", fontsize=11)
    ax.set_ylabel("Shutdown Margin [pcm]", fontsize=11)
    ax.set_title("HALEU NEP Opt-4 - Shutdown & Single-Failure Safety Margins Across Life", fontsize=11.5, fontweight='bold')
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(loc='lower left', framealpha=0.9, fontsize=9.5)
    fig.tight_layout()
    plot2_path = os.path.join(output_dir, "opt4_shutdown_margins_vs_burnup.png")
    fig.savefig(plot2_path)
    plt.close(fig)
    print(f"Saved: {plot2_path}")

    return history, theoretical_lifetime_years

if __name__ == "__main__":
    run_criticality_controlled_depletion_opt4()
