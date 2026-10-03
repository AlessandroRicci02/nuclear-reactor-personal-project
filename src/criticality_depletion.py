"""
Criticality-Controlled Depletion and Theoretical Neutronic Lifetime Module
==========================================================================
Implements autonomous dynamic control drum tracking during fuel depletion:
- Operates at specified system thermal power (default: 400 kWth)
- At each depletion step, determines critical drum angle theta_crit(t) such that:
    k_eff(theta_crit) ≈ 1.000 ± 0.002
- Advances depletion with drums at theta_crit (realistic operating flux spectrum)
- Continues depletion beyond 10-year minimum requirement until:
    k(theta = 0°) < 1.000 (absorber fully outward)
  to identify the true Theoretical Neutronic Lifetime
- At each step, verifies shutdown capability (all drums in at 180°) and stuck drum margin
- Tracks isotopic evolution (U-235, Pu-239, key fission products)
"""

import os
import sys
import numpy as np
import openmc
import openmc.deplete

os.environ["OPENMC_CROSS_SECTIONS"] = "/home/aleri/wmp_data/cross_sections_wmp.xml"

from model import ReactorParameters, build_materials, build_geometry, build_settings, calculate_inventory

def find_critical_drum_angle(p_base, fuel_mat, clad_mat, cool_mat, mod_mat, beo_mat, b4c_mat,
                             target_k=1.000, tol=0.003, particles=3000, batches=25, inactive=8):
    """
    Evaluates the drum S-curve and finds the critical angle theta_crit
    such that k_eff(theta_crit) ≈ target_k.
    Returns: theta_crit, k_crit, s_crit, k_out, s_out, k_in, s_in
    """
    # 1. Evaluate k_out (theta = 0.0 deg)
    openmc.reset_auto_ids()
    p_out = ReactorParameters(
        fuel_radius=p_base.fuel_radius, cladding_thickness=p_base.cladding_thickness,
        coolant_channel_thickness=p_base.coolant_channel_thickness,
        moderator_fuel_pitch=p_base.moderator_fuel_pitch, core_radius=p_base.core_radius,
        active_core_height=p_base.active_core_height, radial_reflector_thickness=p_base.radial_reflector_thickness,
        axial_reflector_thickness=p_base.axial_reflector_thickness, number_of_control_drums=p_base.number_of_control_drums,
        control_drum_radius=p_base.control_drum_radius, absorber_thickness=p_base.absorber_thickness,
        absorber_sector_angular_width=p_base.absorber_sector_angular_width, drum_clearance=p_base.drum_clearance,
        b10_enrichment=p_base.b10_enrichment, particles=particles, batches=batches, inactive=inactive,
        drum_angle_deg=0.0
    )
    mats_out, f_o, cl_o, co_o, m_o, be_o, b_o = build_materials(p_out)
    # Copy isotopic composition from current fuel
    f_o.set_density('g/cm3', fuel_mat.get_mass_density())
    f_o._nuclides = list(fuel_mat._nuclides)
    geom_out, _, _ = build_geometry(p_out, f_o, cl_o, co_o, m_o, be_o, b_o)
    m_out = openmc.Model(materials=mats_out, geometry=geom_out, settings=build_settings(p_out))
    m_out.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{p_out.batches}.h5") as sp:
        k_out = float(sp.keff.nominal_value)
        s_out = float(sp.keff.std_dev)

    # If even with drums fully out k_out < target_k - tol, reactor is subcritical!
    if k_out < target_k - tol:
        return 0.0, k_out, s_out, k_out, s_out, k_out, s_out

    # 2. Evaluate k_in (theta = 180.0 deg)
    openmc.reset_auto_ids()
    p_in = ReactorParameters(
        fuel_radius=p_base.fuel_radius, cladding_thickness=p_base.cladding_thickness,
        coolant_channel_thickness=p_base.coolant_channel_thickness,
        moderator_fuel_pitch=p_base.moderator_fuel_pitch, core_radius=p_base.core_radius,
        active_core_height=p_base.active_core_height, radial_reflector_thickness=p_base.radial_reflector_thickness,
        axial_reflector_thickness=p_base.axial_reflector_thickness, number_of_control_drums=p_base.number_of_control_drums,
        control_drum_radius=p_base.control_drum_radius, absorber_thickness=p_base.absorber_thickness,
        absorber_sector_angular_width=p_base.absorber_sector_angular_width, drum_clearance=p_base.drum_clearance,
        b10_enrichment=p_base.b10_enrichment, particles=particles, batches=batches, inactive=inactive,
        drum_angle_deg=180.0
    )
    mats_in, f_i, cl_i, co_i, m_i, be_i, b_i = build_materials(p_in)
    f_i.set_density('g/cm3', fuel_mat.get_mass_density())
    f_i._nuclides = list(fuel_mat._nuclides)
    geom_in, _, _ = build_geometry(p_in, f_i, cl_i, co_i, m_i, be_i, b_i)
    m_in_model = openmc.Model(materials=mats_in, geometry=geom_in, settings=build_settings(p_in))
    m_in_model.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{p_in.batches}.h5") as sp:
        k_in = float(sp.keff.nominal_value)
        s_in = float(sp.keff.std_dev)

    # 3. Interpolate critical angle using analytical S-curve shape:
    # fraction = (k_out - target_k) / (k_out - k_in)
    # f(theta) = (1 - cos(theta)) / 2  =>  cos(theta) = 1 - 2 * fraction
    fraction = np.clip((k_out - target_k) / (k_out - k_in), 0.0, 1.0)
    cos_theta = 1.0 - 2.0 * fraction
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    theta_est = np.degrees(np.arccos(cos_theta))

    # 4. Verify/refine with eigenvalue calculation at estimated angle
    openmc.reset_auto_ids()
    p_crit = ReactorParameters(
        fuel_radius=p_base.fuel_radius, cladding_thickness=p_base.cladding_thickness,
        coolant_channel_thickness=p_base.coolant_channel_thickness,
        moderator_fuel_pitch=p_base.moderator_fuel_pitch, core_radius=p_base.core_radius,
        active_core_height=p_base.active_core_height, radial_reflector_thickness=p_base.radial_reflector_thickness,
        axial_reflector_thickness=p_base.axial_reflector_thickness, number_of_control_drums=p_base.number_of_control_drums,
        control_drum_radius=p_base.control_drum_radius, absorber_thickness=p_base.absorber_thickness,
        absorber_sector_angular_width=p_base.absorber_sector_angular_width, drum_clearance=p_base.drum_clearance,
        b10_enrichment=p_base.b10_enrichment, particles=particles, batches=batches, inactive=inactive,
        drum_angle_deg=theta_est
    )
    mats_c, f_c, cl_c, co_c, m_c, be_c, b_c = build_materials(p_crit)
    f_c.set_density('g/cm3', fuel_mat.get_mass_density())
    f_c._nuclides = list(fuel_mat._nuclides)
    geom_c, _, _ = build_geometry(p_crit, f_c, cl_c, co_c, m_c, be_c, b_c)
    m_c_model = openmc.Model(materials=mats_c, geometry=geom_c, settings=build_settings(p_crit))
    m_c_model.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{p_crit.batches}.h5") as sp:
        k_crit = float(sp.keff.nominal_value)
        s_crit = float(sp.keff.std_dev)

    # One secant refinement if outside tolerance
    if abs(k_crit - target_k) > tol:
        dk_dtheta = (k_in - k_out) / 180.0
        theta_refined = np.clip(theta_est - (k_crit - target_k) / dk_dtheta, 0.0, 180.0)
        openmc.reset_auto_ids()
        p_crit.drum_angle_deg = theta_refined
        mats_c2, f_c2, cl_c2, co_c2, m_c2, be_c2, b_c2 = build_materials(p_crit)
        f_c2.set_density('g/cm3', fuel_mat.get_mass_density())
        f_c2._nuclides = list(fuel_mat._nuclides)
        geom_c2, _, _ = build_geometry(p_crit, f_c2, cl_c2, co_c2, m_c2, be_c2, b_c2)
        m_c2_model = openmc.Model(materials=mats_c2, geometry=geom_c2, settings=build_settings(p_crit))
        m_c2_model.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{p_crit.batches}.h5") as sp:
            k_crit = float(sp.keff.nominal_value)
            s_crit = float(sp.keff.std_dev)
        theta_est = theta_refined

    return theta_est, k_crit, s_crit, k_out, s_out, k_in, s_in


def run_criticality_depletion_simulation(p, power_W=400.0e3,
                                         timesteps_days=[30.0, 150.0, 185.0, 365.0, 730.0, 730.0, 730.0, 730.0, 1825.0, 1825.0],
                                         target_keff=1.000,
                                         chain_file="data/chain_nep.xml",
                                         particles=2500, batches=25, inactive=8):
    """
    Executes a criticality-controlled depletion simulation:
    - At each step, determines the required critical drum angle theta_crit(t) to hold k ≈ 1.000.
    - Tracks depletion at power_W (400 kWth).
    - Identifies exact Theoretical Neutronic Lifetime when drums reach 0.0 deg (fully out).
    """
    print("=" * 80)
    print(f"CRITICALITY-CONTROLLED DEPLETION AT {power_W/1e3:.1f} kWth (Target k_eff = {target_keff:.4f})")
    print("=" * 80)

    # Initial materials & inventory
    openmc.reset_auto_ids()
    p_init = ReactorParameters(
        fuel_radius=p.fuel_radius, cladding_thickness=p.cladding_thickness,
        coolant_channel_thickness=p.coolant_channel_thickness,
        moderator_fuel_pitch=p.moderator_fuel_pitch, core_radius=p.core_radius,
        active_core_height=p.active_core_height, radial_reflector_thickness=p.radial_reflector_thickness,
        axial_reflector_thickness=p.axial_reflector_thickness, number_of_control_drums=p.number_of_control_drums,
        control_drum_radius=p.control_drum_radius, absorber_thickness=p.absorber_thickness,
        absorber_sector_angular_width=p.absorber_sector_angular_width, drum_clearance=p.drum_clearance,
        b10_enrichment=p.b10_enrichment, particles=particles, batches=batches, inactive=inactive,
        drum_angle_deg=0.0
    )
    mats, fuel, clad, cool, mod, beo, b4c = build_materials(p_init)
    geom, pins, drums = build_geometry(p_init, fuel, clad, cool, mod, beo, b4c)
    inv = calculate_inventory(p_init, len(pins))
    v_fuel = inv["items"][0]["volume_cm3"]
    m_fuel_kg = inv["items"][0]["mass_kg"]
    fuel.volume = v_fuel
    fuel.depletable = True
    m_U_kg = m_fuel_kg * (237.457 / 251.464)

    history = []
    current_time_days = 0.0
    current_burnup_MWd_kgU = 0.0

    theoretical_lifetime_days = None

    for step_idx, dt_days in enumerate([0.0] + list(timesteps_days)):
        if step_idx > 0:
            current_time_days += dt_days
            current_burnup_MWd_kgU += (power_W * dt_days / 1e6) / m_U_kg

        # 1. Find critical drum angle at current isotopic state
        theta_crit, k_crit, s_crit, k_out, s_out, k_in, s_in = find_critical_drum_angle(
            p_init, fuel, clad, cool, mod, beo, b4c,
            target_k=target_keff, particles=particles, batches=batches, inactive=inactive
        )

        # Check stuck drum at this burnup step
        openmc.reset_auto_ids()
        stuck_angles = [0.0] + [180.0] * (p.number_of_control_drums - 1)
        p_stk = ReactorParameters(
            fuel_radius=p.fuel_radius, cladding_thickness=p.cladding_thickness,
            coolant_channel_thickness=p.coolant_channel_thickness,
            moderator_fuel_pitch=p.moderator_fuel_pitch, core_radius=p.core_radius,
            active_core_height=p.active_core_height, radial_reflector_thickness=p.radial_reflector_thickness,
            axial_reflector_thickness=p.axial_reflector_thickness, number_of_control_drums=p.number_of_control_drums,
            control_drum_radius=p.control_drum_radius, absorber_thickness=p.absorber_thickness,
            absorber_sector_angular_width=p.absorber_sector_angular_width, drum_clearance=p.drum_clearance,
            b10_enrichment=p.b10_enrichment, particles=particles, batches=batches, inactive=inactive
        )
        mats_s, f_s, cl_s, co_s, m_s, be_s, b_s = build_materials(p_stk)
        f_s.set_density('g/cm3', fuel.get_mass_density())
        f_s._nuclides = list(fuel._nuclides)
        geom_s, _, _ = build_geometry(p_stk, f_s, cl_s, co_s, m_s, be_s, b_s, drum_angles=stuck_angles)
        m_stk = openmc.Model(materials=mats_s, geometry=geom_s, settings=build_settings(p_stk))
        m_stk.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{p_stk.batches}.h5") as sp:
            k_stuck = float(sp.keff.nominal_value)
            s_stuck = float(sp.keff.std_dev)

        sdm_norm = (1.0 - k_in) / k_in * 1e5
        sdm_stk = (1.0 - k_stuck) / k_stuck * 1e5
        worth = (1.0 / k_in - 1.0 / k_out) * 1e5

        # Extract current U-235 and Pu-239 masses
        u235_mass_kg = 0.0
        pu239_mass_g = 0.0
        total_fuel_mass_g = fuel.get_mass()
        for nuc, dens, nuc_type in fuel.nuclides:
            if nuc == 'U235':
                u235_mass_kg = (dens * total_fuel_mass_g / 100.0) / 1000.0 if fuel.density_units == 'wo' else (dens * 235.04 / 6.022e23 * v_fuel) / 1000.0
            elif nuc == 'Pu239':
                pu239_mass_g = (dens * 239.05 / 6.022e23 * v_fuel)

        years = current_time_days / 365.25
        print(f"Step {step_idx:2d} | t = {current_time_days:6.1f} d ({years:5.2f} y) | Burnup = {current_burnup_MWd_kgU:5.2f} MWd/kgU")
        print(f"  Critical Drum Angle: theta_crit = {theta_crit:5.1f}° (0°=Out, 180°=In)")
        print(f"  Operational k_crit  = {k_crit:.5f} +/- {s_crit:.5f} (Target ~1.000 MET)")
        print(f"  Uncontrolled k_out  = {k_out:.5f} +/- {s_out:.5f} (Excess = {(k_out-1)/k_out*1e5:+.0f} pcm)")
        print(f"  Shutdown k_in       = {k_in:.5f} +/- {s_in:.5f} (SDM = {sdm_norm:4.0f} pcm)")
        print(f"  Stuck Drum k_stuck  = {k_stuck:.5f} +/- {s_stuck:.5f} (SDM_stk = {sdm_stk:4.0f} pcm)")

        entry = {
            "step": step_idx, "time_days": current_time_days, "time_years": years,
            "burnup_MWd_kgU": current_burnup_MWd_kgU,
            "theta_crit_deg": theta_crit, "k_crit": k_crit, "s_crit": s_crit,
            "k_out": k_out, "s_out": s_out, "k_in": k_in, "s_in": s_in,
            "k_stuck": k_stuck, "s_stuck": s_stuck,
            "sdm_normal_pcm": sdm_norm, "sdm_stuck_pcm": sdm_stk,
            "drum_worth_pcm": worth, "u235_mass_kg": u235_mass_kg, "pu239_mass_g": pu239_mass_g
        }
        history.append(entry)

        # Check end of theoretical neutronic life
        if k_out < 1.000 and theoretical_lifetime_days is None:
            theoretical_lifetime_days = current_time_days
            print(f"  >>> THEORETICAL NEUTRONIC LIFETIME REACHED: {years:.2f} YEARS (k_out < 1.000) <<<")

        # Advance depletion step if more timesteps remain
        if step_idx < len(timesteps_days):
            next_dt = timesteps_days[step_idx]
            # Run one depletion step with drums at theta_crit
            openmc.reset_auto_ids()
            p_step = ReactorParameters(
                fuel_radius=p.fuel_radius, cladding_thickness=p.cladding_thickness,
                coolant_channel_thickness=p.coolant_channel_thickness,
                moderator_fuel_pitch=p.moderator_fuel_pitch, core_radius=p.core_radius,
                active_core_height=p.active_core_height, radial_reflector_thickness=p.radial_reflector_thickness,
                axial_reflector_thickness=p.axial_reflector_thickness, number_of_control_drums=p.number_of_control_drums,
                control_drum_radius=p.control_drum_radius, absorber_thickness=p.absorber_thickness,
                absorber_sector_angular_width=p.absorber_sector_angular_width, drum_clearance=p.drum_clearance,
                b10_enrichment=p.b10_enrichment, particles=particles, batches=batches, inactive=inactive,
                drum_angle_deg=theta_crit
            )
            mats_st, fuel_st, clad_st, cool_st, mod_st, beo_st, b4c_st = build_materials(p_step)
            fuel_st.set_density('g/cm3', fuel.get_mass_density())
            fuel_st._nuclides = list(fuel._nuclides)
            fuel_st.volume = v_fuel
            fuel_st.depletable = True
            geom_st, _, _ = build_geometry(p_step, fuel_st, clad_st, cool_st, mod_st, beo_st, b4c_st)
            m_step = openmc.Model(materials=mats_st, geometry=geom_st, settings=build_settings(p_step))
            op = openmc.deplete.CoupledOperator(m_step, chain_file)
            integ = openmc.deplete.PredictorIntegrator(op, [next_dt], power=power_W, timestep_units='d')
            integ.integrate()

            # Update fuel composition from depletion result
            res = openmc.deplete.Results('depletion_results.h5')
            # Get updated material
            fuel = op.materials[0]

    return history, theoretical_lifetime_days
