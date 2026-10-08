#!/usr/bin/env python3
"""
thermal.py - Multi-channel Dual-Pass Thermal-Hydraulic Solver for NEP Reactor Candidates
Evaluates helium flow, convective heat transfer, closed Brayton cycle, and peak moderator temperature.
"""

import math
import numpy as np

CP_HE = 5193.16  # J/(kg K)

def mu_helium(T):
    """Dynamic viscosity of helium [Pa s]"""
    return 3.674e-7 * (T**0.70)

def k_helium(T, P=2.5e6):
    """Thermal conductivity of helium [W/(m K)]"""
    term = 1.0 + 1.123e-8 * P
    exp = 0.71 * (1.0 - 2.0e-9 * P)
    return 2.682e-3 * term * (T**exp)

def prandtl_helium(T, P=2.5e6):
    """Prandtl number of helium [-]"""
    return (CP_HE * mu_helium(T)) / k_helium(T, P)

def friction_factor_tube(Re):
    """Darcy friction factor for smooth tube / channel"""
    if Re < 2300:
        return 64.0 / max(Re, 10.0)
    elif Re < 4000:
        f_lam = 64.0 / 2300.0
        f_turb = (0.790 * math.log(4000.0) - 1.64)**(-2)
        return f_lam + (Re - 2300.0) / (4000.0 - 2300.0) * (f_turb - f_lam)
    else:
        return (0.790 * math.log(Re) - 1.64)**(-2)

def nusselt_tube(Re, Pr):
    """Gnielinski correlation for turbulent tube flow"""
    if Re < 2300:
        return 4.364
    f = friction_factor_tube(Re)
    num = (f / 8.0) * (Re - 1000.0) * Pr
    den = 1.0 + 12.7 * math.sqrt(f / 8.0) * (Pr**(2.0/3.0) - 1.0)
    return max(4.364, num / den)

def evaluate_thermal(spec):
    """
    Evaluate the dual-pass thermal-hydraulic performance for candidate specification.
    Returns dictionary with:
      - m_dot_total (kg/s)
      - T_mod_peak (K)
      - T_mod_avg (K)
      - thermal_buffer_K (950.0 - T_mod_peak)
      - passes_thermal (bool: T_mod_peak < 950.0)
      - T_clad_peak (K)
      - T_fuel_peak (K)
    """
    Q_total = spec['power_kW'] * 1000.0 # W
    T_inlet = spec['inlet_K']
    T_outlet = spec['outlet_K']
    delta_T = T_outlet - T_inlet
    P_nom = spec['pressure_MPa'] * 1e6
    m_dot_total = Q_total / (CP_HE * delta_T)

    n_pins = spec['pins']
    H_core = spec['height_cm'] * 1e-2
    r_fuel = spec['r_fuel_cm'] * 1e-2
    t_clad = spec['t_clad_cm'] * 1e-2
    r_clad = r_fuel + t_clad
    t_hot_gap = spec['t_hot_gap_cm'] * 1e-2
    r_hot_chan = r_clad + t_hot_gap
    family = spec.get('family', 'matrix_dual_pass')

    # Moderator thermal conductivity (Hu et al. ORNL 2021: k_YH2 ~ 20.0 - 22.0 W/(m K))
    k_mod = 21.0

    # Fraction of nuclear heat deposited directly in moderator
    f_mod = 0.022
    Q_mod_total = f_mod * Q_total # ~8.8 kW
    Q_fuel_total = (1.0 - f_mod) * Q_total # ~391.2 kW

    # Hot channel peaking factors
    F_z = 1.28
    F_r = 1.24
    F_q = F_z * F_r # ~1.587

    m_dot_pin_avg = m_dot_total / n_pins

    if family in ("matrix_dual_pass", "matrix_expanded_pitch", "pins271_high_burnup", "pins169_deep_mod"):
        # Dual-pass with dedicated aux holes in the YHx matrix
        n_aux = spec.get('aux_holes', 6)
        r_aux = spec.get('r_aux_cm', 0.165) * 1e-2
        
        # Cold gas (Pass 1) flows through the moderator aux holes
        # Orifice distribution: 85% of total flow passes through aux channels
        m_dot_aux_tot = m_dot_total * 0.85
        m_dot_per_aux = m_dot_aux_tot / (n_pins * n_aux)
        
        A_aux = math.pi * (r_aux**2)
        D_h_aux = 2.0 * r_aux
        
        # Bulk temp rise of coolant in Pass 1 across moderator
        delta_T_cool_mod = Q_mod_total / (m_dot_aux_tot * CP_HE)
        T_cool_mod_exit = T_inlet + delta_T_cool_mod
        T_cool_mod_avg = 0.5 * (T_inlet + T_cool_mod_exit)
        
        Re_aux = (m_dot_per_aux * D_h_aux) / (A_aux * mu_helium(T_cool_mod_avg))
        Pr_aux = prandtl_helium(T_cool_mod_avg, P_nom)
        Nu_aux = nusselt_tube(Re_aux, Pr_aux)
        h_aux = Nu_aux * k_helium(T_cool_mod_avg, P_nom) / D_h_aux

        # Convective drop in peak channel
        A_aux_surf = 2.0 * math.pi * r_aux * H_core
        q_mod_peak_pin = (Q_mod_total / n_pins) * F_r
        q_mod_peak_aux = q_mod_peak_pin / n_aux
        q_flux_aux_peak = (q_mod_peak_aux / A_aux_surf) * F_z
        delta_T_film_mod = q_flux_aux_peak / h_aux

        # Conduction drop in moderator matrix cell
        pitch_m = spec['pitch_cm'] * 1e-2
        r_cell_equiv = math.sqrt(math.sqrt(3.0)/(2.0*math.pi)) * pitch_m
        delta_T_cond_mod = (q_mod_peak_pin * F_z / (math.pi * (r_cell_equiv**2 - r_hot_chan**2) * H_core)) * (r_cell_equiv - r_hot_chan)**2 / (4.0 * k_mod)

        # Thermal insulation drop across Al2O3 liner and hot helium boundary
        # Al2O3 layer (100 um) has k ~ 30 W/m K, creates thermal resistance
        t_al2o3 = spec.get('t_al2o3_cm', 0.010) * 1e-2
        delta_T_leak = 14.5 * (0.010 * 1e-2 / max(t_al2o3, 0.005*1e-2))

        T_mod_peak = T_cool_mod_exit + delta_T_film_mod + delta_T_cond_mod + delta_T_leak
        T_mod_avg = T_cool_mod_avg + 0.5 * (delta_T_film_mod + delta_T_cond_mod)

    elif family == "dual_annulus":
        # Dual-annulus concentric counter-current cooling
        t_part = spec.get('t_part_cm', 0.035) * 1e-2
        t_cold_gap = spec.get('t_cold_gap_cm', 0.080) * 1e-2
        r_part_in = r_hot_chan
        r_part_out = r_part_in + t_part
        r_cold_out = r_part_out + t_cold_gap

        A_flow_cold = math.pi * (r_cold_out**2 - r_part_out**2)
        P_wet_cold = 2.0 * math.pi * (r_cold_out + r_part_out)
        D_h_cold = 4.0 * A_flow_cold / P_wet_cold

        m_dot_cold_pin = m_dot_total / n_pins
        delta_T_cool_mod = Q_mod_total / (m_dot_total * CP_HE)
        T_cool_mod_exit = T_inlet + delta_T_cool_mod
        T_cool_mod_avg = 0.5 * (T_inlet + T_cool_mod_exit)

        Re_cold = (m_dot_cold_pin * D_h_cold) / (A_flow_cold * mu_helium(T_cool_mod_avg))
        Pr_cold = prandtl_helium(T_cool_mod_avg, P_nom)
        Nu_cold = nusselt_tube(Re_cold, Pr_cold)
        h_cold = Nu_cold * k_helium(T_cool_mod_avg, P_nom) / D_h_cold

        q_mod_peak_pin = (Q_mod_total / n_pins) * F_r
        A_cold_mod_surf = 2.0 * math.pi * r_cold_out * H_core
        q_flux_cold_peak = (q_mod_peak_pin / A_cold_mod_surf) * F_z
        delta_T_film_mod = q_flux_cold_peak / h_cold

        pitch_m = spec['pitch_cm'] * 1e-2
        r_cell_equiv = math.sqrt(math.sqrt(3.0)/(2.0*math.pi)) * pitch_m
        delta_T_cond_mod = (q_mod_peak_pin * F_z / (math.pi * (r_cell_equiv**2 - r_cold_out**2) * H_core)) * (r_cell_equiv - r_cold_out)**2 / (4.0 * k_mod)

        delta_T_leak = 11.0 # partition tube separates hot and cold streams

        T_mod_peak = T_cool_mod_exit + delta_T_film_mod + delta_T_cond_mod + delta_T_leak
        T_mod_avg = T_cool_mod_avg + 0.5 * (delta_T_film_mod + delta_T_cond_mod)

    # Hot pass (Fuel Pin):
    A_flow_hot = math.pi * (r_hot_chan**2 - r_clad**2)
    P_wet_hot = 2.0 * math.pi * (r_hot_chan + r_clad)
    D_h_hot = 4.0 * A_flow_hot / P_wet_hot
    Re_hot = (m_dot_pin_avg * D_h_hot) / (A_flow_hot * mu_helium(T_outlet))
    Pr_hot = prandtl_helium(T_outlet, P_nom)
    Nu_hot = nusselt_tube(Re_hot, Pr_hot)
    h_hot = Nu_hot * k_helium(T_outlet, P_nom) / D_h_hot

    P_pin_peak = (Q_fuel_total / n_pins) * F_r
    A_clad_surf = 2.0 * math.pi * r_clad * H_core
    q_flux_clad_peak = (P_pin_peak / A_clad_surf) * F_z
    T_clad_peak = T_outlet + q_flux_clad_peak / h_hot
    
    # UN thermal conductivity ~ 18 W/(m K)
    delta_T_fuel = (P_pin_peak * F_z / (math.pi * (r_fuel**2) * H_core)) * (r_fuel**2) / (4.0 * 18.0)
    T_fuel_peak = T_clad_peak + delta_T_fuel + 20.0 # gap

    thermal_buffer = 950.0 - T_mod_peak
    passes_thermal = bool(T_mod_peak < 950.0)

    return {
        "m_dot_total_kg_s": round(m_dot_total, 5),
        "T_mod_peak_K": round(T_mod_peak, 2),
        "T_mod_avg_K": round(T_mod_avg, 2),
        "thermal_buffer_K": round(thermal_buffer, 2),
        "passes_thermal": passes_thermal,
        "T_clad_peak_K": round(T_clad_peak, 2),
        "T_fuel_peak_K": round(T_fuel_peak, 2),
    }

if __name__ == "__main__":
    from catalog import get_candidates
    cands = get_candidates()
    for c in cands[:3]:
        res = evaluate_thermal(c)
        print(c['id'], res)

