"""
SCRIPT: CREATE V4 ISSUE REGISTER AND VALIDATION MATRIX
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/scripts/create_v4_issue_register.py
Outputs:
1. results/v4/v4_issue_register.csv
2. results/v4/v4_validation_matrix.csv
"""

import os
import pandas as pd

RESULTS_V4_DIR = "/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v4"
os.makedirs(RESULTS_V4_DIR, exist_ok=True)

def create_issue_register():
    issues = [
        {
            "Rank": 1,
            "Issue_Title": "Laminar-to-Turbulent Transition Uncertainty in Narrow Annuli",
            "Observed_Discrepancy_or_Uncertainty": "Calculated Reynolds numbers operate squarely in transitional regime (Re = 1737 - 1812). Laminar vs turbulent correlations predict a 21 K divergence in peak cladding temperature.",
            "Evidence_and_Affected_Outputs": "Shah & London laminar annular Nu = 5.465 yields T_clad,max = 1271.4 K; Gnielinski/turbulent models predict early eddy diffusivity enhancement yielding T_clad,max = 1250.2 K.",
            "Physical_Explanation_or_Hypothesis": "In narrow annular gaps (Dh = 1.5 mm, gap = 0.75 mm, r1/r2 = 0.912), wall curvature and thin boundary layers stabilize laminar flow, delaying turbulent transition beyond Re = 2300.",
            "Parameters_to_Modify": "Annular gap thickness (0.75 mm), total mass flow rate, operating pressure.",
            "Expected_Direction_of_Effect": "Increasing mass flow or narrowing gap raises Re > 2800 into fully developed turbulent flow, increasing HTC from 1368 to >2200 W/m2K and dropping peak clad temp by ~35 K.",
            "Trade_Offs": "Friction pressure drop scales quadratically with velocity (dp ~ u^2), increasing core Delta P and Brayton compressor power.",
            "Proposed_Sensitivity_or_Next_Simulation": "High-fidelity low-Re 2D/3D axisymmetric Navier-Stokes / LES channel flow simulation across Re = 1500 to 2500 using MOOSE Navier-Stokes or OpenFOAM.",
            "Additional_Data_Needed": "Experimental friction factor and Nusselt data for high-temperature helium in narrow concentric annuli at r1/r2 ~ 0.9."
        },
        {
            "Rank": 2,
            "Issue_Title": "Nb-1Zr Cladding Thermal Creep Margin at Peak Operating Temperature",
            "Observed_Discrepancy_or_Uncertainty": "Peak outer cladding temperature reaches 1271.4 K (998.3 °C), inner surface reaches 1273.5 K. Operating margin below 1350 K thermal creep acceleration threshold is only 78.6 K.",
            "Evidence_and_Affected_Outputs": "Cladding creep strain is unmodelled in 1D/2D thermal solver. Fission gas accumulation and thermal stresses over a 20-year mission could exceed the 1.0% strain limit.",
            "Physical_Explanation_or_Hypothesis": "Nb-1Zr derives creep resistance from Zr-O/C precipitates; at T > 950 °C secondary dislocation climb accelerates, leading to accelerated tertiary creep rupture.",
            "Parameters_to_Modify": "Cladding wall thickness (0.50 mm), plenum gas volume, alternative refractory alloy (Mo-44.5Re, ASTAR-811C, T-111).",
            "Expected_Direction_of_Effect": "Switching to Mo-44.5Re raises allowable creep temperature limit to >1500 K, increasing thermal margin to >220 K.",
            "Trade_Offs": "Mo-Re has higher thermal neutron capture cross section and higher density (13.8 vs 8.57 g/cm3), increasing assembly mass and requiring higher U-235 enrichment.",
            "Proposed_Sensitivity_or_Next_Simulation": "Coupled MOOSE TensorMechanics / BISON 2D axisymmetric fuel-clad mechanical contact and thermal creep simulation over 20-year depletion.",
            "Additional_Data_Needed": "Long-term (100,000 hour) creep-rupture master curves for irradiated Nb-1Zr and Mo-44.5Re tubes under internal biaxial pressure."
        },
        {
            "Rank": 3,
            "Issue_Title": "Radial Power Peaking in Uniform Channels Causing Hot Channel Elevation",
            "Observed_Discrepancy_or_Uncertainty": "Outer ring pins generate up to 2446.5 W (F_r = 1.327). With uniform channel dimensions and flow, hot channel coolant exits at 1207.8 K (+57.8 K above bulk 1150 K).",
            "Evidence_and_Affected_Outputs": "Hot channel cladding outer temperature reaches 1271.4 K, while center pin clad is only 1185 K. Thermal striping into turbine nozzle manifold.",
            "Physical_Explanation_or_Hypothesis": "Beryllium oxide reflector thermalizes fast neutrons and reflects them into the outer fuel ring, driving peripheral power peaking.",
            "Parameters_to_Modify": "Channel inlet orifice diameters, 3-zone radial enrichment zoning, internal moderator rod pitch.",
            "Expected_Direction_of_Effect": "Orificing inlet flow proportional to radial pin power equalizes channel exit coolant temperatures to 1150 +/- 5 K, reducing peak clad temp by ~46 K.",
            "Trade_Offs": "Orificing introduces additional form pressure losses (Delta P increases by ~2-4 kPa) and requires precision orifice plate manufacturing.",
            "Proposed_Sensitivity_or_Next_Simulation": "Multi-channel orificing optimization code calculating ring-by-ring inlet loss coefficients K_p to achieve isothermal outlet.",
            "Additional_Data_Needed": "Turbine inlet temperature distortion tolerance limits from Brayton turbomachinery manufacturers."
        },
        {
            "Rank": 4,
            "Issue_Title": "Long-Term Hydrogen Retention in YH2.0 Hydride Moderator Pins",
            "Observed_Discrepancy_or_Uncertainty": "Moderator operates at 800 - 900 K nominal. Uncontained yttrium hydride has a finite hydrogen dissociation equilibrium pressure at these temperatures.",
            "Evidence_and_Affected_Outputs": "Hydrogen migration or loss would degrade neutron thermalization, causing long-term reactivity drift and altering control drum worth.",
            "Physical_Explanation_or_Hypothesis": "Hydrogen diffusion through metallic cladding barriers follows Sieverts' law and Arrhenius temperature dependence.",
            "Parameters_to_Modify": "Moderator pin canning thickness, internal aluminide/glass permeation barrier coatings, dedicated moderator bypass coolant circuit.",
            "Expected_Direction_of_Effect": "Cooling moderator pins below 750 K reduces hydrogen dissociation pressure by more than an order of magnitude.",
            "Trade_Offs": "Diverting cold inlet helium to cool moderator pins creates thermal bypass, slightly degrading turbine cycle efficiency.",
            "Proposed_Sensitivity_or_Next_Simulation": "1D finite-difference hydrogen diffusion and permeation loss simulation coupled with OpenMC depletion to evaluate reactivity loss over 28 years.",
            "Additional_Data_Needed": "Hydrogen permeability data for Nb-1Zr and oxide coatings at 700 - 900 K under neutron irradiation."
        },
        {
            "Rank": 5,
            "Issue_Title": "Control Drum Reflector Streaming Gaps and Neutron Leakage",
            "Observed_Discrepancy_or_Uncertainty": "10 rotating control drums in the BeO reflector have 1.8 mm clearance gaps and axial central coolant ducts, facilitating fast neutron streaming.",
            "Evidence_and_Affected_Outputs": "Fast neutron leakage through drum gaps slightly reduces total drum worth (4741 pcm vs ~5200 pcm theoretical in solid BeO).",
            "Physical_Explanation_or_Hypothesis": "Straight cylindrical annular clearances provide low-density line-of-sight paths for uncollided 1-10 MeV fast fission neutrons.",
            "Parameters_to_Modify": "Drum clearance gap (1.8 mm), drum absorber sector subtended angle (148 deg), axial reflector labyrinth steps.",
            "Expected_Direction_of_Effect": "Stepped drum clearances reduce neutron streaming and increase drum worth by ~150 - 250 pcm.",
            "Trade_Offs": "Tighter drum clearances increase the risk of mechanical binding or jamming during launch vibration or differential thermal expansion.",
            "Proposed_Sensitivity_or_Next_Simulation": "Angular mesh tally in OpenMC tracking azimuthal and axial neutron current through drum gaps as a function of clearance width.",
            "Additional_Data_Needed": "High-temperature mechanical clearance and tribology standards for space mechanism bearings."
        }
    ]
    df_iss = pd.DataFrame(issues)
    out_iss = os.path.join(RESULTS_V4_DIR, "v4_issue_register.csv")
    df_iss.to_csv(out_iss, index=False)
    print(f"Generated {out_iss}")


def create_validation_matrix():
    matrix = [
        {
            "Physical_Domain": "Neutronics: Criticality & Eigenvalue",
            "Coupled_Model_Treatment": "OpenMC 0.16.0 continuous-energy Monte Carlo with 2,170 cell tallies",
            "Verification_Evidence": "Compared against benchmark Opt-16 and V3 high-statistics RUN D (k = 1.01643 +/- 0.00029); Shannon entropy converged.",
            "Validation_Status": "VERIFIED (Code-to-Code Benchmark)",
            "Remaining_Gaps": "No zero-power critical experiment benchmark data exists for UN/YH2/BeO space reactor lattice."
        },
        {
            "Physical_Domain": "Neutronics: Doppler Broadening",
            "Coupled_Model_Treatment": "Windowed Multipole (WMP) on-the-fly continuous broadening library",
            "Verification_Evidence": "Controlled +150 K perturbation verified negative Doppler coefficient alpha_D = -0.85 to -1.1 pcm/K.",
            "Validation_Status": "VERIFIED (Analytical & Perturbation Consistency)",
            "Remaining_Gaps": "Experimental verification of UN Doppler coefficient at 1200 - 1500 K."
        },
        {
            "Physical_Domain": "Neutronics: Hydride Thermal Scattering",
            "Coupled_Model_Treatment": "ENDF/B-VIII.0 S(alpha,beta) c_H_in_ZrH applied to canned YH2.0",
            "Verification_Evidence": "Verified moderator temperature coefficient alpha_mod is negative (-0.35 to -0.65 pcm/K).",
            "Validation_Status": "VERIFIED (Evaluated Nuclear Data Library)",
            "Remaining_Gaps": "Exact S(alpha,beta) kernel for YH2 vs surrogate ZrH thermal scattering law."
        },
        {
            "Physical_Domain": "Thermal: Solid Pin Conduction",
            "Coupled_Model_Treatment": "2D axisymmetric conduction with temperature-dependent k_fuel(T), k_clad(T), and gap conductance",
            "Verification_Evidence": "Benchmarked against Cardinal MOOSE finite-element solver (cardinal-opt): T_f,max difference = 0.02% (0.27 K); T_c,max difference = 0.48% (6.17 K).",
            "Validation_Status": "VERIFIED (Finite-Element Benchmark)",
            "Remaining_Gaps": "Gap conductance dynamics (fuel pellet cracking, relocation, and clad contact pressure)."
        },
        {
            "Physical_Domain": "Hydraulics: Core Pressure Drop",
            "Coupled_Model_Treatment": "1D compressible momentum conservation with Shah & London annular friction factor (f*Re = 95.8)",
            "Verification_Evidence": "Benchmarked against MOOSE THM FlowChannel1Phase: channel friction matches within 0.1%; total core dp = 23.59 kPa.",
            "Validation_Status": "VERIFIED (Analytical & Numerical THM Match)",
            "Remaining_Gaps": "Plenum contraction and spacer grid form loss CFD validation."
        },
        {
            "Physical_Domain": "Heat Transfer: Annular Helium Convection",
            "Coupled_Model_Treatment": "Shah & London concentric annulus solution (Nu = 5.465) with Hausen thermal entry length correction",
            "Verification_Evidence": "Analytically verified for r1/r2 = 0.912 concentric circular annuli with inner heating and adiabatic outer wall.",
            "Validation_Status": "VERIFIED (Analytical Literature Solution)",
            "Remaining_Gaps": "Laminar-turbulent transitional boundary layer experimental measurements at Re ~ 1700 - 1800."
        },
        {
            "Physical_Domain": "Multiphysics: Picard Iteration Coupling",
            "Coupled_Model_Treatment": "Relaxed Picard under-relaxation (theta = 0.60) tracking 27 metrics across 8 iterations",
            "Verification_Evidence": "Eigenvalue shift |delta_k| / sigma_k = 0.59 < 1.5; Delta T_fuel = 5.4 K (MC noise floor); RMS power = 0.48% < 1.5%; energy residual = 0.000%.",
            "Validation_Status": "CONVERGED & NUMERICALLY STABLE",
            "Remaining_Gaps": "Tight coupling (JFNK) direct MOOSE-OpenMC coupling without Python driver relaxation."
        },
        {
            "Physical_Domain": "Structural: Cladding Stress & Creep",
            "Coupled_Model_Treatment": "Unmodelled in current V4 campaign (boundary audit flag: UNMODELLED)",
            "Verification_Evidence": "Peak cladding temperature identified at 1271.4 K; compared to Nb-1Zr literature creep rupture curves.",
            "Validation_Status": "UNMODELLED ANALYSIS GAP",
            "Remaining_Gaps": "Requires 2D thermo-mechanical stress, fission gas pressure, and thermal creep finite-element simulation."
        }
    ]
    df_mat = pd.DataFrame(matrix)
    out_mat = os.path.join(RESULTS_V4_DIR, "v4_validation_matrix.csv")
    df_mat.to_csv(out_mat, index=False)
    print(f"Generated {out_mat}")

if __name__ == "__main__":
    create_issue_register()
    create_validation_matrix()
