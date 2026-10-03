"""
THERMOPHYSICAL PROPERTIES, FLUID DYNAMICS AND HEAT TRANSFER CORRELATIONS
Opt-17 Space Nuclear Electric Propulsion Reactor (Cand 17-Clean)

Operating Conditions:
- Total Thermal Power: 400.0 kWth
- Coolant: Pure Helium (He-4)
- Nominal System Pressure: 2.50 MPa (25.0 bar)
- Coolant Inlet Temperature: 950.0 K (676.85 deg C)
- Target Mixed Outlet Temperature: 1150.0 K (876.85 deg C)
- Number of Active Pins: 217
- Pin Pitch: 2.270 cm
- Fuel Pellet Radius: 0.580 cm
- Cladding Thickness: 0.050 cm (Nb-1Zr)
- Annular Coolant Gap: 0.075 cm
- Active Core Height: 42.45 cm

Literature Sources:
- Helium properties & EOS: Petersen (1981), NASA SP-3011, NIST Chemistry WebBook,
  consistent with MOOSE HeliumFluidProperties implementation.
- Friction factor & Nusselt correlations: Shah & London (1978), Gnielinski (1976),
  Kays & Crawford (1993).
- Fuel (UN) thermal conductivity & Cp: Hayes, Thomas & Peddicord (1990), J. Nucl. Mater.
- Cladding (Nb-1Zr) thermal conductivity: Touloukian et al. (1970), TPRC Data Series.
- Fuel-Clad gap conductance: Todreas & Kazimi (2012), Nuclear Systems I.
"""

import math
import numpy as np

# Universal constants
R_GAS = 8314.462618  # Universal gas constant [J / (kmol K)]
M_HE = 4.002602      # Helium molar mass [kg / kmol]
P_NOM = 2.50e6       # Nominal pressure [Pa]

# Core geometry
N_PINS = 217
H_CORE = 0.4245      # Active core height [m]
R_FUEL = 0.0058      # Fuel pellet radius [m]
T_CLAD = 0.0005      # Cladding thickness [m]
R_CLAD_OUT = R_FUEL + T_CLAD  # 0.0063 m
T_GAP = 0.00075      # Annular coolant gap [m]
R_CHAN_OUT = R_CLAD_OUT + T_GAP # 0.00705 m

# Flow areas & dimensions
A_FLOW_PIN = math.pi * (R_CHAN_OUT**2 - R_CLAD_OUT**2)  # Annular cross-section [m2] = 3.14552e-5 m2
P_WET_PIN = 2.0 * math.pi * (R_CHAN_OUT + R_CLAD_OUT)    # Wetted perimeter [m] = 0.08388 m
D_H = 4.0 * A_FLOW_PIN / P_WET_PIN                      # Hydraulic diameter = 2 * T_GAP = 0.0015 m
R_RATIO = R_CLAD_OUT / R_CHAN_OUT                       # Diameter ratio r* = 0.8936
A_HEAT_PIN = 2.0 * math.pi * R_CLAD_OUT * H_CORE        # Outer cladding heat transfer area [m2]
V_FUEL_PIN = math.pi * (R_FUEL**2) * H_CORE             # Fuel pellet volume [m3]

# Operating parameters
Q_TOTAL = 400000.0   # Total reactor thermal power [W]
T_INLET = 950.0      # Coolant inlet temperature [K]
DELTA_T_BULK = 200.0 # Nominal coolant bulk temperature rise [K]

# Gap conductance default (W / (m2 K))
H_GAP_NOMINAL = 25000.0  # Tight contact / swaged Nb-1Zr with He fill gas


def cp_helium(T=1050.0):
    """Constant-pressure heat capacity of pure helium [J / (kg K)]."""
    return 2.5 * (R_GAS / M_HE)  # 5193.16 J / (kg K)


def total_coolant_mass_flow():
    """Total coolant mass flow rate [kg/s]."""
    return Q_TOTAL / (cp_helium() * DELTA_T_BULK)


def pin_coolant_mass_flow():
    """Coolant mass flow rate per fuel pin [kg/s] under uniform header orificing."""
    return total_coolant_mass_flow() / N_PINS


def rho_helium(T, P=P_NOM):
    """Density of pure helium [kg/m3] from ideal gas equation of state at pressure P [Pa]."""
    return (P * M_HE) / (R_GAS * T)


def mu_helium(T):
    """
    Dynamic viscosity of helium [Pa s] from Petersen (1981) correlation:
    mu = 3.674e-7 * T^0.70 [Pa s]
    Matches MOOSE HeliumFluidProperties::mu_from_p_T
    """
    return 3.674e-7 * (T**0.70)


def k_helium(T, P=P_NOM):
    """
    Thermal conductivity of helium [W / (m K)] from Petersen (1981):
    k = 2.682e-3 * (1 + 1.123e-8 * P) * T^(0.71 * (1 - 2e-9 * P))
    Matches MOOSE HeliumFluidProperties::k_from_p_T
    """
    term = 1.0 + 1.123e-8 * P
    exp = 0.71 * (1.0 - 2.0e-9 * P)
    return 2.682e-3 * term * (T**exp)


def prandtl_helium(T, P=P_NOM):
    """Prandtl number of helium Pr = cp * mu / k."""
    return cp_helium(T) * mu_helium(T) / k_helium(T, P)


def reynolds_helium(m_dot_pin, T, P=P_NOM):
    """
    Reynolds number in the annular coolant channel:
    Re = rho * u * D_h / mu = (m_dot / A_flow) * D_h / mu
    """
    return (m_dot_pin / A_FLOW_PIN) * D_H / mu_helium(T)


def friction_factor_annulus(Re, r_ratio=R_RATIO):
    """
    Darcy friction factor for concentric annular duct.
    - Laminar (Re <= 1000): f = C_f / Re, where C_f(r*=0.894) ~ 95.8 (Shah & London 1978)
    - Turbulent (Re >= 3000): Blasius with annular factor f = 1.05 * 0.3164 * Re^(-0.25)
    - Transitional (1000 < Re < 3000): weighted interpolation
    """
    cf_lam = 95.8
    f_lam = cf_lam / max(Re, 1.0)
    f_turb = 1.05 * 0.3164 * (max(Re, 1.0)**(-0.25))
    
    if Re <= 1000.0:
        return f_lam
    elif Re >= 3000.0:
        return f_turb
    else:
        w = (Re - 1000.0) / (3000.0 - 1000.0)
        return (1.0 - w) * f_lam + w * f_turb


def pressure_drop_channel(m_dot_pin, T_bulk_profile, z_edges, P=P_NOM):
    """
    Calculates channel pressure drop [Pa] including:
    - Frictional drop across all axial slices
    - Acceleration drop due to density decrease
    - Form losses (inlet contraction K_in=0.5, outlet expansion K_out=1.0)
    """
    n_slices = len(T_bulk_profile)
    dz = H_CORE / n_slices
    G = m_dot_pin / A_FLOW_PIN  # Mass flux [kg / (m2 s)]
    
    dp_fric = 0.0
    for s in range(n_slices):
        T_s = T_bulk_profile[s]
        rho_s = rho_helium(T_s, P)
        Re_s = reynolds_helium(m_dot_pin, T_s, P)
        f_s = friction_factor_annulus(Re_s)
        dp_fric += f_s * (dz / D_H) * (G**2) / (2.0 * rho_s)
        
    rho_in = rho_helium(T_INLET, P)
    rho_out = rho_helium(T_bulk_profile[-1], P)
    dp_accel = (G**2) * (1.0 / rho_out - 1.0 / rho_in)
    
    # Form loss at inlet and outlet
    u_in = m_dot_pin / (rho_in * A_FLOW_PIN)
    u_out = m_dot_pin / (rho_out * A_FLOW_PIN)
    dp_form = 0.5 * rho_in * (u_in**2) / 2.0 + 1.0 * rho_out * (u_out**2) / 2.0
    
    dp_total = dp_fric + dp_accel + dp_form
    return {
        "dp_total_Pa": dp_total,
        "dp_fric_Pa": dp_fric,
        "dp_accel_Pa": dp_accel,
        "dp_form_Pa": dp_form
    }


def nusselt_annulus(Re, Pr, Dh=D_H, L=H_CORE, r_ratio=R_RATIO, model="gnielinski_trans"):
    """
    Nusselt number for annular duct with inner wall heated and outer wall insulated.
    
    Models:
    - 'shah_london_laminar': Pure laminar with Hausen entry correction: Nu = 5.42 * [1 + 0.14 * (Dh/L)^0.5]
    - 'gnielinski_trans': Transitional Gnielinski formulation with annular correction
    """
    Nu_lam_fd = 5.35 + 0.15 * (r_ratio - 0.8) / 0.2  # ~5.42 for r* = 0.894
    Nu_lam = Nu_lam_fd * (1.0 + 0.14 * math.sqrt(Dh / L)) # ~5.47
    
    if model == "shah_london_laminar":
        return Nu_lam

    # Gnielinski transitional model
    if Re <= 1000.0:
        return Nu_lam
    elif Re < 2300.0:
        f2300 = (1.8 * math.log10(2300.0) - 1.5)**(-2)
        gamma = 1.0 + 12.7 * math.sqrt(f2300 / 8.0) * (Pr**(2.0 / 3.0) - 1.0)
        Nu_turb2300 = ((f2300 / 8.0) * (2300.0 - 1000.0) * Pr / gamma) * (1.0 + (Dh / L)**(2.0 / 3.0))
        F_ann = 0.75 * (r_ratio**(-0.17))
        Nu_2300 = Nu_turb2300 * F_ann
        weight = (Re - 1000.0) / (2300.0 - 1000.0)
        return (1.0 - weight) * Nu_lam + weight * max(Nu_2300, Nu_lam)
    else:
        f = (1.8 * math.log10(Re) - 1.5)**(-2)
        gamma = 1.0 + 12.7 * math.sqrt(f / 8.0) * (Pr**(2.0 / 3.0) - 1.0)
        Nu_turb = ((f / 8.0) * (Re - 1000.0) * Pr / gamma) * (1.0 + (Dh / L)**(2.0 / 3.0))
        F_ann = 0.75 * (r_ratio**(-0.17))
        return Nu_turb * F_ann


def h_convection(m_dot_pin, T_bulk, P=P_NOM, htc_multiplier=1.0, model="gnielinski_trans"):
    """Convective heat transfer coefficient h [W / (m2 K)]."""
    Re = reynolds_helium(m_dot_pin, T_bulk, P)
    Pr = prandtl_helium(T_bulk, P)
    Nu = nusselt_annulus(Re, Pr, model=model)
    k_he = k_helium(T_bulk, P)
    h_val = Nu * k_he / D_H
    return h_val * htc_multiplier


def k_fuel_un(T):
    """Thermal conductivity of UN fuel [W / (m K)] (Hayes, Thomas & Peddicord 1990)."""
    return 1.37 * (T**0.41)


def cp_fuel_un(T):
    """Specific heat capacity of UN fuel [J / (kg K)] (Hayes et al. 1990)."""
    return 190.0 + 0.055 * T


def k_clad_nb1zr(T):
    """Thermal conductivity of Nb-1Zr cladding [W / (m K)] (Touloukian 1970)."""
    return 41.8 + 0.0125 * T


def cp_clad_nb1zr(T):
    """Specific heat capacity of Nb-1Zr cladding [J / (kg K)]."""
    return 270.0 + 0.045 * T
