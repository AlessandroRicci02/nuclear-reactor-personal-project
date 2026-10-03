# FINAL REPORT: PROGRESSIVE NEUTRONICS–THERMOHYDRAULICS COUPLED CAMPAIGN
## Opt-17 / Cand 17-Clean Space Nuclear Electric Propulsion Reactor (400.0 kWth)

**Location**: `/home/aleri/OpenMC_projects/first_test/thermohydraulics/`  
**Date**: October 3, 2026  
**Reference Architecture**: Opt-17 / Cand 17-Clean (Hexagonal Lattice, 217 UN Pins, BeO Reflector, 10 Control Drums)  
**Coupling Suite**: OpenMC 0.16.0 $\longleftrightarrow$ Picard Iteration Engine $\longleftrightarrow$ Cardinal / MOOSE FE Benchmark  
**Cross-Section Data**: Windowed Multipole (WMP) On-the-Fly Doppler Broadening Library (`/home/aleri/wmp_data/cross_sections_wmp.xml`)  

---

## 1. EXECUTIVE SUMMARY

An autonomous, multi-stage coupled campaign was executed for the **Opt-17 / Cand 17-Clean** space NEP reactor across progressive runs (**RUN 0 through RUN 5**). The campaign systematically incorporated:
1. Multi-zone 3D OpenMC power tallies across **2,170 fuel cells** ($217\text{ pins} \times 10\text{ axial slices}$).
2. On-the-fly Windowed Multipole (WMP) Doppler broadening in sintered UN fuel.
3. Annular helium density feedback via the local ideal-gas equation of state $\rho(p, T_{\mathrm{bulk}})$.
4. Metal-hydride $\text{YH}_{2.0}$ moderator feedback with thermal scattering $S(\alpha,\beta) = \text{c\_H\_in\_ZrH}$.
5. Transitional flow thermohydraulics in annular coolant channels ($Re \approx 1600 - 1891$), friction factor, acceleration and form pressure drops ($\Delta P_{\mathrm{tot}} \approx 23.59\text{ kPa}$).
6. Pellet-cladding finite gap conductance ($h_{\mathrm{gap}} = 25,000\text{ W/(m}^2\text{K)}$).
7. Picard under-relaxation ($\theta = 0.60$) tracking 27 quantitative metrics.
8. Finite-element qualification of the converged peak and average pins using the Cardinal / MOOSE executable (`cardinal-opt`).

The converged coupled eigenvalue settled at **$k_{\mathrm{eff}} = 1.01651 \pm 0.00067$**, introducing a net negative feedback Doppler insertion of **$\Delta \rho = -210.2\text{ pcm}$** relative to the cold/uncoupled baseline ($k_{\mathrm{eff}} = 1.01869 \pm 0.00061$). The peak fuel centerline temperature is **$1296.2\text{ K}$**, outer cladding temperature is **$1271.0\text{ K}$**, and mixed coolant outlet temperature is exactly **$1150.0\text{ K}$**.

---

## 2. OPT-17 REFERENCE CONFIGURATION

The model strictly adheres to the qualified **Cand 17-Clean** reference configuration from `/home/aleri/OpenMC_projects/first_test/results_second_architecture/opt17/screening/verify_opt17_clean.py`:
* **Thermal Rating**: $400.0\text{ kWth}$ ($100.0\text{ kWe}$ at $25.0\%$ Brayton conversion efficiency).
* **Fuel Element**: Sintered UN pellets ($\rho = 13.50\text{ g/cm}^3$), pellet radius $r_f = 0.580\text{ cm}$.
* **Two-Zone Enrichment**: Inner 19 pins (Rings 0–2) at $19.30\text{ wt\%}$ $^{235}\text{U}$; Outer 198 pins (Rings 3–8) at $19.75\text{ wt\%}$ $^{235}\text{U}$.
* **Cladding**: Refractory alloy Nb-1Zr ($\rho = 8.57\text{ g/cm}^3$), thickness $t_{\mathrm{clad}} = 0.050\text{ cm}$, outer radius $r_{co} = 0.630\text{ cm}$.
* **Annular Coolant Channel**: Gap $t_{\mathrm{gap}} = 0.075\text{ cm}$, outer radius $r_{\mathrm{cho}} = 0.705\text{ cm}$.
* **Auxiliary Channels**: 6 per pin cell, radius $r_{\mathrm{add}} = 0.150\text{ cm}$ at $r = 0.910\text{ cm}$.
* **Active Core**: 217 fuel pins in a 9-ring hexagonal lattice, pitch $P = 2.270\text{ cm}$, active height $H = 42.45\text{ cm}$, core radius $R = 19.45\text{ cm}$.
* **Moderator**: Canned $\text{YH}_{2.0}$ ($\rho = 4.30\text{ g/cm}^3$) with $S(\alpha, \beta)$ scattering.
* **Reflectors**: Sintered BeO ($\rho = 2.90\text{ g/cm}^3$); radial thickness $9.80\text{ cm}$; axial thickness $9.275\text{ cm}$ top and bottom.
* **Control System**: 10 rotating BeO drums ($r = 3.95\text{ cm}$) with $148.0^\circ$ sectors of $90\%$ $^{10}\text{B}$-enriched $\text{B}_4\text{C}$ absorber ($t = 1.50\text{ cm}$); nominal operating state: $0^\circ$ absorber rotated OUT.

---

## 3. SOFTWARE VERSIONS & CAPABILITIES INSPECTION

* **OpenMC**: Version 0.16.0 (Commit `617d35a506`, Built with OpenMP, 8 threads). Python 3.13 in `/home/aleri/miniforge3/envs/openmc-env/`.
* **Cardinal**: Executable `/home/aleri/projects/cardinal/cardinal-opt`.
  * Object registry inspection confirmed availability of `OpenMCCellAverageProblem`, `HeliumFluidProperties` (from MOOSE `fluid_properties` module), `HeatConduction`, `ConvectiveFluxFunction`, and PETSc Hypre BoomerAMG preconditioning.
  * Verified that MOOSE `HeliumFluidProperties` implements Petersen (1981) transport properties identical to our thermohydraulics module.

---

## 4. NEUTRONICS RESOLUTION & OPENMC MODEL

* **Spatial Discretization**: $217\text{ pins} \times 10\text{ axial slices} = 2,170\text{ fuel cells}$.
* **Tallies**: Cell tally scoring both `heating` and `fission` across all 2,170 cells.
* **Particle Settings**: 80,000 particles/batch, 45 batches (10 inactive, 35 active) $\rightarrow 2.80 \times 10^6$ active particle histories per solve.
* **Doppler Treatment**: Windowed Multipole continuous on-the-fly Doppler evaluation enabled (`<temperature_multipole>true</temperature_multipole>`).

---

## 5. THERMAL-HYDRAULIC MODEL & FORMULATIONS

### 5.1 Fluid Mechanics & Property Evaluation
* **Coolant**: Pure $^{4}\text{He}$, nominal pressure $p = 2.50\text{ MPa}$ ($25.0\text{ bar}$), inlet temperature $T_{\mathrm{in}} = 950.0\text{ K}$, target mixed outlet $T_{\mathrm{out}} = 1150.0\text{ K}$ ($\Delta T = 200.0\text{ K}$).
* **Total Mass Flow**: $\dot{M}_{\mathrm{tot}} = \frac{400\text{ kW}}{c_p \Delta T} = 0.385122\text{ kg/s}$; per pin: $\dot{m}_{\mathrm{pin}} = 1.7748\text{ g/s}$ (`ASSUMPTION`: uniform header orificing).
* **Hydraulic Diameter**: $D_h = 2 (r_{\mathrm{cho}} - r_{co}) = 1.500\text{ mm}$; Flow area: $A_{\mathrm{flow}} = 3.1455 \times 10^{-5}\text{ m}^2$.
* **Dimensionless Regime**:
  * $Re = 1604.8 - 1891.1$ (transitional annular flow regime).
  * $Pr = 0.6455$.
  * $Nu = 5.465$ (Shah & London 1978 concentric annular laminar correlation with Hausen thermal entry length correction).
  * Convective HTC: $h = 1279.1 - 1509.6\text{ W/(m}^2\text{K)}$.
* **Pressure Drop**:
  * Frictional: $\Delta p_{\mathrm{fric}} = 20.84\text{ kPa}$.
  * Acceleration: $\Delta p_{\mathrm{accel}} = 0.67\text{ kPa}$.
  * Form losses: $\Delta p_{\mathrm{form}} = 2.08\text{ kPa}$.
  * Total channel pressure drop: $\Delta P_{\mathrm{tot}} = 23.59\text{ kPa}$ ($0.236\text{ bar}$, $0.94\%$ of system pressure).

### 5.2 Solid Heat Conduction & Gap Conductance
* **UN Fuel**: $k_f(T) = 1.37 T^{0.41}\text{ W/(m}\cdot\text{K)}$ [Hayes et al. 1990].
* **Nb-1Zr Cladding**: $k_c(T) = 41.8 + 0.0125 T\text{ W/(m}\cdot\text{K)}$ [Touloukian 1970].
* **Pellet-Clad Gap**: $h_{\mathrm{gap}} = 25,000\text{ W/(m}^2\text{K)}$ (`COMPUTED`), accounting for thermal expansion closure and residual helium contact, introducing $\Delta T_{\mathrm{gap}} \approx 6.4 - 7.7\text{ K}$.

---

## 6. RUN 0 THROUGH RUN 5 PROGRESSION & COMPARISON

| RUN | Improvement Category | $k_{\mathrm{eff}} \pm \sigma$ | Peak Pin Power | $F_q$ | $T_{f,\max}$ [K] | $T_{c,\max}$ [K] | $T_{\mathrm{out}}$ [K] | $\Delta P$ [kPa] | Status |
| :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **RUN 0** | Baseline (Uncoupled, uniform 1050 K) | $1.01869 \pm 0.00061$ | $2461.3\text{ W}$ | $1.727$ | $1292.6$ | $1274.3$ | $1150.0$ | $23.58$ | `COMPUTED` |
| **RUN 1** | Basic Doppler Feedback | $1.01683 \pm 0.00061$ | $2447.3\text{ W}$ | $1.742$ | $1288.9$ | $1271.5$ | $1150.0$ | $23.59$ | `COMPUTED` |
| **RUN 2** | Doppler + He Density + Mod Feedback | $1.01592 \pm 0.00080$ | $2437.2\text{ W}$ | $1.752$ | $1287.3$ | $1270.8$ | $1150.0$ | $23.58$ | `COMPUTED` |
| **RUN 3** | Transitional Re + $\Delta P$ + Gap Conductance | $1.01673 \pm 0.00056$ | $2445.2\text{ W}$ | $1.740$ | $1295.8$ | $1272.8$ | $1150.0$ | $23.59$ | `COMPUTED` |
| **RUN 4** | Picard Fixed-Point Under-Relaxation Loop | $1.01651 \pm 0.00067$ | $2431.1\text{ W}$ | $1.722$ | $1296.2$ | $1271.0$ | $1150.0$ | $23.59$ | `COMPUTED` |
| **RUN 5** | Final Converged Master Configuration | **$1.01651 \pm 0.00067$** | **$2431.1\text{ W}$** | **$1.722$** | **$1296.2$** | **$1271.0$** | **$1150.0$** | **$23.59$** | `CONVERGED` |

---

## 7. PICARD COUPLING CONVERGENCE HISTORY (27 METRICS)

| Iter | $k_{\mathrm{eff}} \pm \sigma$ | $\Delta k$ | $\Delta \rho$ [pcm] | $P_{\mathrm{peak}}$ [W] | $\Delta P_{\mathrm{peak}}$ | $P_{\mathrm{RMS}}$ | $F_r$ | $F_z$ | $F_q$ | $T_{f,\max}$ [K] | $T_{c,\max}$ [K] | $T_{\mathrm{out}}$ [K] | $\Delta T$ [K] | $\rho_{\mathrm{He},\min}$ | $\rho_{\mathrm{He},\max}$ | $Re_{\min}$ | $Re_{\max}$ | $h_{\min}$ | $h_{\max}$ | $\Delta P$ [Pa] |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **0** | $1.01869 \pm 0.00061$ | $0.0$ | $0.0$ | $2461.3$ | $0.0$ | $0.0$ | $1.335$ | $1.274$ | $1.727$ | $1292.6$ | $1274.3$ | $1150.0$ | $227.7$ | $0.995$ | $1.261$ | $1601.9$ | $1891.0$ | $1279.2$ | $1512.4$ | $23585$ |
| **1** | $1.01698 \pm 0.00071$ | $-0.00171$ | $-165.3$ | $2468.6$ | $+0.30\%$ | $1.02\%$ | $1.339$ | $1.276$ | $1.748$ | $1299.9$ | $1275.2$ | $1150.0$ | $101.0$ | $0.994$ | $1.262$ | $1600.9$ | $1891.2$ | $1279.0$ | $1513.3$ | $23584$ |
| **2** | $1.01682 \pm 0.00063$ | $-0.00016$ | $-180.4$ | $2454.1$ | $-0.59\%$ | $1.07\%$ | $1.331$ | $1.279$ | $1.758$ | $1297.5$ | $1272.3$ | $1150.0$ | $44.2$ | $0.996$ | $1.262$ | $1602.7$ | $1891.2$ | $1279.0$ | $1511.6$ | $23584$ |
| **3** | $1.01684 \pm 0.00054$ | $+0.00002$ | $-178.9$ | $2437.0$ | $-0.70\%$ | $1.24\%$ | $1.322$ | $1.277$ | $1.731$ | $1294.5$ | $1271.0$ | $1150.0$ | $21.0$ | $0.997$ | $1.262$ | $1604.2$ | $1891.1$ | $1279.1$ | $1510.2$ | $23588$ |
| **4** | $1.01699 \pm 0.00054$ | $+0.00015$ | $-164.4$ | $2431.2$ | $-0.24\%$ | $0.92\%$ | $1.319$ | $1.276$ | $1.705$ | $1293.1$ | $1271.0$ | $1150.0$ | $13.3$ | $0.998$ | $1.262$ | $1605.1$ | $1891.1$ | $1279.1$ | $1509.3$ | $23585$ |
| **5** | **$1.01651 \pm 0.00067$** | **$-0.00047$** | **$-210.2$** | **$2431.1$** | **$-0.002\%$** | **$1.21\%$** | **$1.319$** | **$1.278$** | **$1.722$** | **$1296.2$** | **$1271.0$** | **$1150.0$** | **$10.7$** | **$0.998$** | **$1.261$** | **$1604.8$** | **$1891.1$** | **$1279.1$** | **$1509.6$** | **$23589$** |

*Note: In table above, mass flow is constant at $0.385122\text{ kg/s}$, $Nu_{\min}=Nu_{\max}=5.465$, and $T_{\mathrm{cool},\max}=1213.8\text{ K}$.*

---

## 8. CARDINAL / MOOSE FINITE-ELEMENT QUALIFICATION

The converged power distributions were solved directly using Cardinal (`cardinal-opt`) 2D RZ axisymmetric finite elements with PETSc BoomerAMG solver:
* **Peak Pin (#215, $2431.13\text{ W}$)**:
  * Model analytical $T_{f,\max} = 1293.86\text{ K}$; Cardinal FE $T_{f,\max} = 1293.59\text{ K}$ $\rightarrow$ Difference: **$0.27\text{ K}$ ($0.02\%$)**.
  * Model analytical $T_{c,\max} = 1270.97\text{ K}$; Cardinal FE $T_{c,\max} = 1277.14\text{ K}$ $\rightarrow$ Difference: **$6.17\text{ K}$ ($0.48\%$)**.
* **Average Pin (#112, $1850.36\text{ W}$)**:
  * Model analytical $T_{f,\max} = 1213.07\text{ K}$; Cardinal FE $T_{f,\max} = 1212.37\text{ K}$ $\rightarrow$ Difference: **$0.70\text{ K}$ ($0.05\%$)**.
  * Model analytical $T_{c,\max} = 1196.15\text{ K}$; Cardinal FE $T_{c,\max} = 1199.86\text{ K}$ $\rightarrow$ Difference: **$3.71\text{ K}$ ($0.31\%$)**.

---

## 9. SENSITIVITY & ROBUSTNESS ANALYSIS

1. **HTC $+10\%$**: Peak fuel temperature shifts by **$-8.60\text{ K}$** ($1285.26\text{ K}$).
2. **HTC $-10\%$**: Peak fuel temperature shifts by **$+11.63\text{ K}$** ($1305.48\text{ K}$).
3. **Pure Laminar Shah & London**: Peak fuel temperature is identical (**$0.00\text{ K}$** change) due to laminar entry predominance at $Re \approx 1605 - 1891$.
4. **Gap Conductance**: Perfect contact reduces peak fuel temperature by **$-6.43\text{ K}$** ($1287.43\text{ K}$); finite gap ($25,000\text{ W/m}^2\text{K}$) contributes $+6.43\text{ K}$.

---

## 10. EXACT RESULT ARTIFACTS MAP

All generated files reside inside `/home/aleri/OpenMC_projects/first_test/thermohydraulics/`:

| Output Deliverable | Exact Filesystem Path | Format | Description |
| :--- | :--- | :---: | :--- |
| **Run Comparison** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/run_comparison.csv` | CSV | Summary of RUN 0 through RUN 5 progression |
| **Coupling Iterations** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/coupling/coupling_iterations.csv` | CSV | All 27 tracked metrics across Picard iterations |
| **$k_{\mathrm{eff}}$ Convergence** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/coupling/keff_convergence.csv` | CSV | Eigenvalue convergence with statistical uncertainties |
| **Final Pin Powers** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/neutronics/final_pin_power.csv` | CSV | Pin coordinates, ring, power [W], $F_r$, exit temperatures |
| **Final Axial Power** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/neutronics/final_axial_power.csv` | CSV | Sliced axial power [W], $F_z(z)$, linear power, mean $T$ |
| **Final 3D Power Field** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/neutronics/final_3d_power.csv` | CSV | All 2,170 cells: $(x, y, z)$, power, density [$\text{W/m}^3$], $F_q$ |
| **Final 3D Thermal Field** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/thermal/final_thermal_field.csv` | CSV | All 2,170 cells: $T_{\mathrm{bulk}}, T_{\mathrm{clad}}, T_{f,\mathrm{vol}}, T_{f,\mathrm{center}}, h$ |
| **Final Coolant Field** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/thermal/final_coolant_field.csv` | CSV | Local helium density, velocity, $Re, Nu, h$ per slice |
| **Final Neutronics Summary**| `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/final/final_neutronics_summary.json` | JSON | Converged eigenvalue, peaking factors, power |
| **Final Thermal Summary** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/final/final_thermal_summary.json` | JSON | Converged temperatures, pressure drop, Reynolds range |
| **Final Coupled Summary** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/coupling/final_coupled_summary.json` | JSON | Master summary of baseline vs coupled state |
| **Cardinal Benchmark Data** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/validation/cardinal_benchmark_results.json` | JSON | Peak and average pin Cardinal FE comparisons |
| **Sensitivity Data** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/validation/sensitivity_results.json` | JSON | $\pm 10\%$ HTC, laminar, and gap sensitivity results |
| **Plot 1: $k_{\mathrm{eff}}$ vs Iteration** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/figures/keff_vs_iteration.png` | PNG | Eigenvalue convergence with Monte Carlo error bars |
| **Plot 2: Reactivity vs Iteration** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/figures/reactivity_vs_iteration.png` | PNG | Net Doppler and density feedback shift $\Delta \rho$ [pcm] |
| **Plot 3: Peak Power vs Iteration**| `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/figures/peak_power_vs_iteration.png` | PNG | Hot pin power and peaking factors ($F_r, F_q$) |
| **Plot 4: Power Convergence** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/figures/power_convergence.png` | PNG | Power RMS change vs Monte Carlo noise threshold |
| **Plot 5: Fuel $T$ vs Iteration** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/figures/fuel_temperature_vs_iteration.png` | PNG | Peak fuel centerline temperature stabilization curve |
| **Plot 6: Clad $T$ vs Iteration** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/figures/clad_temperature_vs_iteration.png` | PNG | Peak outer cladding temperature stabilization curve |
| **Plot 7: Coolant $T$ vs Iteration** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/figures/coolant_temperature_vs_iteration.png` | PNG | Mixed outlet and hot channel exit temperatures |
| **Plot 8: Density Feedback** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/figures/density_feedback_vs_iteration.png` | PNG | Axial helium density range variation across iterations |
| **Plot 9: Axial Power Profile** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/figures/axial_power_profile_final.png` | PNG | Final core axial cosine power profile $F_z(z)$ |
| **Plot 10: Radial Power Dist** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/figures/radial_power_distribution_final.png` | PNG | 2D scatter map of all 217 fuel pin powers |
| **Plot 11: Fuel Axial Profiles** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/figures/fuel_temperature_axial_profile.png` | PNG | Axial fuel temperature profiles (peak, average, edge) |
| **Plot 12: Coolant Axial Profiles**| `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/figures/coolant_temperature_axial_profile.png` | PNG | Axial bulk coolant and cladding surface curves |
| **Plot 13: Pressure Drop** | `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/figures/pressure_drop.png` | PNG | Core channel pressure drop vs coupling iteration |
