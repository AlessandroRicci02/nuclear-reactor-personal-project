# V6 PHYSICS GAP ANALYSIS & PHENOMENOLOGICAL COMPLETENESS AUDIT

**Document ID:** `results/v6/PHYSICS_GAP_ANALYSIS.md`  
**Auditor:** V6 Multiphysics Engineering Team  
**Date:** October 3, 2026  
**Status:** COMPLETE, AUDITED, REALISTIC  

---

## 1. Executive Summary

This document performs a systematic phenomenological audit of the physical phenomena present in the Opt-17 $400\text{ kWth}$ / $\sim 100\text{ kWe}$ Nuclear Electric Propulsion (NEP) reactor. The objective is to identify and classify all physical phenomena that are either omitted, simplified, or approximated by the current suite of tools (OpenMC 0.16.0, Cardinal/MOOSE 1D compressible flow, 2D axisymmetric solid conduction).

The purpose of this analysis is to prevent premature declarations of "validation" and to clearly delineate the boundaries between **verified computational approximations** and **unmodeled physical risks**.

---

## 2. Phenomenological Classification Framework

Each physical phenomenon is evaluated across three attributes:
1. **Physical Mechanism**: The governing conservation laws and material behaviors.
2. **Current Model Representation**: How it is modeled in the current toolchain (explicit, surrogate, parametric, or omitted).
3. **Engineering Impact Classification**:
   * **NEGLIGIBLE**: Minimal impact on key metrics ($< 1\%$ deviation in $k_{\mathrm{eff}}$, temperatures, or pressure drops).
   * **SECOND-ORDER**: Measurable effect ($1 - 5\%$), but can be accommodated within standard engineering design margins.
   * **IMPORTANT**: Meaningful impact ($5 - 15\%$) that shifts trade-off curves or optimal design parameter values.
   * **POTENTIALLY DESIGN-LIMITING**: Could induce structural failure, critical constraint violation, or mission loss if unaddressed.
   * **REQUIRES EXTERNAL HIGHER-FIDELITY TOOL**: Beyond the capability of 1D flow/Monte Carlo statics, requiring dedicated finite-element or transient codes (e.g. BISON, MOOSE-TensorMechanics, Griffin, RELAP-7).

---

## 3. Systematic Assessment of Missing and Approximated Physics

| # | Physical Phenomenon | Phenomenological Mechanism | Current Representation in Opt-17 Model | Phenomenon Impact Classification | Justification & Engineering Risk Assessment |
| :-: | :--- | :--- | :--- | :---: | :--- |
| **1** | **Cladding Thermal & Irradiation Creep** | Nb-1Zr undergoes thermal creep under internal pressure stress at $T > 1250\text{ K}$, accelerated by neutron flux. | Surrogate constraint limit ($T_{c,\max} \le 1350\text{ K}$) based on Touloukian (1970). No time-dependent strain calculation. | **POTENTIALLY DESIGN-LIMITING** | Over 10–30 years, creep strain can narrow the annular coolant gap, increasing friction $\Delta P$ and elevating local temperatures. Requires BISON or MOOSE-TensorMechanics creep model. |
| **2** | **Hydrogen Retention & $\text{YH}_{2.0}$ Moderator Dissociation** | Yttrium hydride dissociates at high temperature ($\text{YH}_x \to \text{YH}_{x-\delta} + \frac{\delta}{2}\text{H}_2$), generating hydrogen overpressure and shifting neutron spectrum. | Constant material density ($4.30\text{ g/cm}^3$) and stoichiometric ratio (H/Y = 2.0). Canned cladding assumed impermeable. | **POTENTIALLY DESIGN-LIMITING** | Hydrogen leakage from moderator cans softens the epithermal spectrum, resulting in severe reactivity loss ($\Delta k \sim -2000\text{ pcm}$). Requires hydrogen permeation barrier modeling and temperature-dependent Sieverts' law equilibrium. |
| **3** | **Fuel Swelling & Fission Gas Release (FGR)** | UN fuel accumulates solid and gaseous fission products (Xe, Kr), causing volumetric pellet swelling and gap closure. | Fixed pellet radius ($r_f$) and constant gap conductance ($h_{\mathrm{gap}} = 25,000\text{ W/m}^2\text{K}$). | **IMPORTANT** | At $\sim 10 - 20\text{ MWd/kgU}$, UN swelling closes the initial manufacturing gap ($50\text{ }\mu\text{m}$), eliminating gap thermal resistance ($h_{\mathrm{gap}} \to 10^5\text{ W/m}^2\text{K}$) but inducing Pellet-Cladding Mechanical Interaction (PCMI). |
| **4** | **Pellet Cracking & Relocation** | Thermal stress gradients ($\Delta T_{\mathrm{pellet}} \approx 25 - 40\text{ K}$) induce radial cracking and relocation against clad. | Intact solid cylinder with $k(T) = 1.37 T^{0.41}$. | **SECOND-ORDER** | Because linear heat rate is low ($q' \approx 4.3\text{ kW/m}$), UN centerline temperature is very low ($< 1300\text{ K}$), meaning cracking is mild compared to LWR or fast reactor pins. |
| **5** | **Pellet-Cladding Mechanical Interaction (PCMI)** | Differential thermal expansion and fuel swelling compress clad, generating biaxial tensile hoop stresses. | Omitted (gap is modeled as open thermal boundary). | **IMPORTANT** | Low thermal rating limits PCMI severity in early life, but long-term UN swelling will produce contact stress. Requires multi-dimensional fuel performance code. |
| **6** | **Radiation Heat Transfer in Coolant Gap** | Surface-to-surface radiative heat exchange across annular helium channel between cladding outer surface and duct inner wall. | Omitted in 1D flow model; pure convection and conduction assumed. | **SECOND-ORDER** | At $T \sim 1200\text{ K}$, radiation heat flux $q''_{\mathrm{rad}} = \sigma \epsilon (T_{\mathrm{clad}}^4 - T_{\mathrm{duct}}^4) \approx 5 - 15\text{ kW/m}^2$, which accounts for $\sim 5 - 10\%$ of total heat transfer to the moderator can. |
| **7** | **Thermal Expansion & Reactivity Feedback** | Volumetric thermal expansion of fuel pins, moderator matrix, core grid plates, and reflector increases leakage and reduces density. | Evaluated only via uniform material temperature cross-section broadening. Geometry is held geometrically rigid in OpenMC. | **IMPORTANT** | Grid plate radial expansion produces an unconditional negative reactivity feedback ($\sim -0.5\text{ to } -1.2\text{ pcm/K}$). Omitting thermal expansion makes the baseline model slightly conservative. |
| **8** | **Plenum & Header Flow Distribution & Pressure Losses** | Flow entering lower plenum from duct turns into 217 orifice holes; upper plenum header mixes 217 channel jets into turbine intake. | Lumped form loss coefficients ($K_{\mathrm{in}} = 0.5, K_{\mathrm{out}} = 1.0$). Zero plenum flow maldistribution. | **IMPORTANT** | 3D plenum recirculation vortices and manifold static pressure recovery create radial inlet pressure gradients ($\Delta P_{\mathrm{plenum}} \approx 1 - 3\text{ kPa}$), slightly altering flow distribution. |
| **9** | **Inter-Channel Turbulent Cross-Flow & Mixing** | Adjacent coolant channels communicate thermally via conduction through the canned moderator block. | 217 independent parallel 1D flow channels. Zero cross-conduction between channels. | **IMPORTANT** | Moderator conduction flattens radial temperature gradients. Heat conducts from the central hot channels to cooler peripheral channels, reducing true $\Delta T_{\mathrm{hot}}$ by $\sim 5 - 15\text{ K}$ relative to uncoupled 1D channels. |
| **10** | **Structural Pressure Vessel Stresses** | Vessel cylinder at $2.50\text{ MPa}$ helium pressure and $T \sim 950\text{ K}$ subject to ASME Section III primary membrane stresses. | Parametric mass allowance ($t_{\mathrm{vessel}} = 3.0\text{ mm}$ Nb-1Zr). No stress analysis. | **SECOND-ORDER** | For $R = 29.25\text{ cm}$, $P = 2.5\text{ MPa}$, hoop stress $\sigma_h = P R / t = 2.5 \times 292.5 / 3.0 = 243.8\text{ MPa}$. Nb-1Zr yield strength at $950\text{ K}$ is $\sim 180\text{ MPa}$; vessel thickness must be increased to $\ge 4.5\text{ mm}$ to satisfy ASME safety factors. |
| **11** | **Radiation Shielding & Attenuation** | Attenuation of neutron and gamma flux toward Brayton turbo-alternator, radiators, and spacecraft electronics. | Omitted. Only reactor core subassembly and reflector are modeled. | **REQUIRES EXTERNAL TOOL** | A shadow shield (e.g. layered $\text{LiH} / \text{W}$) is required for all NEP spacecraft, adding $\sim 250 - 450\text{ kg}$ to total spacecraft payload mass. |
| **12** | **Dynamic Startup, Shutdown & Brayton Transients** | Reactor startup ramp rate, drum drive backlash, Brayton loop startup, compressor stall, and thermal shock during shutdown. | Steady-state equilibrium calculations only. | **REQUIRES EXTERNAL TOOL** | Dynamic transients dictate control drum rotational speed limits and recuperator thermal inertia. Requires coupled point kinetics / thermal-hydraulic system code (RELAP-7 or SAM). |
| **13** | **Orifice Manufacturing Tolerances & Erosion** | Machining tolerances ($\pm 0.05\text{ mm}$) and high-velocity helium particle erosion alter orifice discharge coefficients over 30 years. | Idealized constant orifice loss coefficients $K_{p,i}$. Zero tolerance degradation. | **IMPORTANT** | Small orifice holes ($d \sim 2.5 - 4.8\text{ mm}$) are sensitive to $\pm 0.05\text{ mm}$ tolerance, which introduces $\pm 4 - 8\%$ flow variation and $\pm 5 - 10\text{ K}$ exit temperature uncertainty. |

---

## 4. Hierarchy of Limiting Physics & Necessary Modeling Upgrades

Based on this audit, the following prioritized modeling advancements are required before the reactor can transition from conceptual design to engineering validation:

1. **Top Design-Limiting Physical Risk:** **Hydrogen Retention in Canned $\text{YH}_{2.0}$ Moderator**.  
   *Action:* Couple a Sieverts' law hydrogen dissociation model to the moderator temperature field. Evaluate hydrogen loss rates across Nb-1Zr / Mo cans and quantify the reactivity penalty over 30 years.
2. **Top Structural-Limiting Physical Risk:** **Nb-1Zr Cladding Thermal Creep & Gap Shrinkage**.  
   *Action:* Incorporate time-dependent Norton creep strain rate equations into the Cardinal/MOOSE solid mechanics solve to verify gap clearance over a 30-year lifetime.
3. **Top Hydraulic Modeling Advancement:** **Inter-Channel Moderator Thermal Cross-Conduction**.  
   *Action:* Couple the 217 1D flow channels to a 2D/3D MOOSE heat conduction mesh representing the solid moderator and fuel assembly, replacing the isolated channel approximation.
