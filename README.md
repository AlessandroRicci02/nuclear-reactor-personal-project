# Conceptual HALEU Space Nuclear Electric Propulsion (NEP) Reactor Model

> [!IMPORTANT]
> ### 🏆 RISULTATO FINALE DEFINITIVO DEL PROGETTO (FINAL DEFINITIVE REACTOR)
> L'architettura nucleare ottimale, consolidata e validata per la missione **NEP da 400 kWth / 136.39 kWe** (configurazione **`C01_matrix_217_yh20_880_1120`**) si trova nella cartella dedicata:
> 
> 📁 **[`definitive_final_geometry/`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry)**
> 
> * 📄 **Relazione Tecnica Ufficiale Consolidata**: [`definitive_final_geometry/report.txt`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/report.txt)
> * 📄 **Guida Tecnica di Dettaglio del Reattore**: [`definitive_final_geometry/README.md`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/README.md)
> * 📊 **Dati Numerici, Tabelle CSV e Benchmark**: [`definitive_final_geometry/data/`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/data/)
> * 🖼️ **Galleria di 13 Figure ad ALTA RISOLUZIONE (300 DPI)**: [`definitive_final_geometry/figures/`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/)
> * ⚛️ **Modello OpenMC CSG Geometrico**: [`definitive_final_geometry/neutronics_openmc/`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/neutronics_openmc/)
> * ⏳ **Simulazione di Depletion 26 Step (0-25 anni, controllo critico)**: [`definitive_final_geometry/depletion_analysis/`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/depletion_analysis/)
> * 🔥 **Multifisica Termoidraulica FEM MOOSE / Cardinal & Ciclo Brayton**: [`definitive_final_geometry/thermohydraulics_moose_cardinal/`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/thermohydraulics_moose_cardinal/)
> 
> *(Vedere la Sezione 8 in fondo a questo documento per la sintesi completa dei parametri del reattore finale).*

An open-source, research and educational computational neutronics model of a **High-Assay Low-Enriched Uranium (HALEU)** space reactor for megawatt- and hundred-kilowatt-class **Nuclear Electric Propulsion (NEP)**, modeled using the **OpenMC Monte Carlo Code** (v0.16.0).

---

## 1. Project Objective & NEP Application

The objective of this project is to develop, verify, and progressively optimize a scientifically defensible, parameterized conceptual space reactor assembly for deep-space NEP cargo and crew transit tugs (e.g. Earth-to-Mars cargo prepositioning).

Key functional characteristics:
* **Reactor Concept**: Epithermal/intermediate space core fueled with Uranium Mononitride (**UN**, $19.75\text{ wt\%}$ $^{235}\text{U}$ HALEU), moderated with high-temperature Yttrium Hydride (**$\text{YH}_2$**), cooled by an inert binary Helium-Xenon gas mixture (**$\text{He-Xe}$**, $M = 40.0\text{ g/mol}$), and reflected by Beryllium Oxide (**BeO**).
* **Reactivity Control**: External rotating BeO control drums with $90\%$ $^{10}\text{B}$-enriched Boron Carbide ($\text{B}_4\text{C}$) sectors embedded in the radial reflector, providing full operational reactivity control and shutdown margin without in-core control rod penetrations.
* **Historical Reference Design**: **Opt-6** (Buckling-Optimal & Under-Moderated Epithermal Core: **$514.21\text{ kg}$**, **$118.79\text{ L}$**, $H/R = 1.955$, $\ge 10\text{ years}$ operational life, **$25.0\text{ years}$** theoretical neutronic lifetime, Total TMC: **$-3.01 \pm 0.57\text{ pcm/K}$** with $3\sigma$ bound **$-1.30\text{ pcm/K} < 0$**).

---

## 2. Power Specifications & System Energy Balance

The design explicitly differentiates between reactor-level neutronic heat generation and vehicle-level power delivery:

* **Reactor Thermal Power ($P_\text{th}$)**: **$400\text{ kWth}$** ($4.00 \times 10^5\text{ W}$).
  * This is the total fission power generated in the active core.
  * In OpenMC depletion, this power governs the fission rate ($\sim 1.25 \times 10^{16}\text{ fissions/s}$) and the rate of fissile $^{235}\text{U}$ consumption.
* **Thermal Power Deposited in Fuel**: $\approx 380\text{--}385\text{ kWth}$ ($\sim 95\text{--}96\%$ of total fission energy deposited directly in the fuel pellets via fission fragment kinetic energy and prompt electrons; remainder deposited in moderator, cladding, and reflectors via gamma rays and fast neutron elastic scattering).
* **Electrical System Power Bus ($P_e$)**: **$100\text{ kWe}$** ($1.00 \times 10^5\text{ W}$).
  * Dedicated electrical power delivered to the Power Processing Units (PPU) and nested Hall effect / gridded ion thrusters.
* **Assumed Thermal-to-Electric Conversion Efficiency ($\eta_e$)**: **$\sim 25\%$**.
  * Supported by closed-loop recuperated Brayton turbomachinery literature for $\text{He-Xe}$ cycles operating between $900\text{ K}$ (compressor inlet) and $1200\text{--}1300\text{ K}$ (turbine inlet).
  * *Important Note*: The $25\%$ conversion efficiency and $400\text{ kWth} \to 100\text{ kWe}$ relationship is a **system-level design input/target**, not a result computed by OpenMC (which solves neutral-particle transport and isotopic depletion).

---

## 3. Computational Software & Nuclear Data

* **Transport Code**: OpenMC Monte Carlo Particle Transport Code (v0.16.0).
* **Multi-threading**: OpenMP shared-memory parallelism ($8\text{ threads}$).
* **Nuclear Data Library**: Evaluated Nuclear Data File ENDF/B-VII.1 (HDF5 format).
* **Doppler Broadening**: On-the-fly continuous-energy **Windowed Multipole (WMP)** methodology, enabling rigorous Doppler evaluation across arbitrary temperatures ($300\text{--}1500\text{ K}$) without multi-temperature library approximations.
* **Thermal Scattering**: $S(\alpha, \beta)$ evaluated with continuous kernels (`c_Be_in_BeO`, `c_O_in_BeO`, and `c_H_in_ZrH` as audited surrogate for $\text{YH}_2$).
* **Depletion Engine**: OpenMC `CoupledOperator` with `PredictorIntegrator` coupled to `data/chain_nep.xml`.

---

## 4. Repository Structure

```
/home/aleri/OpenMC_projects/first_test/
├── definitive_final_geometry/    # 🏆 RISULTATO FINALE DEFINITIVO CONSOLIDATO (C01 NEP 400 kWth, 100+ kWe)
│   ├── README.md                 # Documentazione tecnica completa del reattore finale
│   ├── report.txt                # Relazione tecnica esaustiva (Neutronica, Depletion, MOOSE, Brayton)
│   ├── results.json              # Specifiche e parametri scalari validati
│   ├── figures/                  # 13 Figure ad ALTA RISOLUZIONE (300 DPI)
│   ├── data/                     # Tabelle dati CSV (burnup a step di 1 anno, FEM, ciclo Brayton)
│   ├── neutronics_openmc/        # Script e modelli CSG OpenMC
│   ├── depletion_analysis/       # Script burnup 26 step (0-25 anni) con controllo criticità attivo
│   └── thermohydraulics_moose_cardinal/ # Benchmark FEM e simulatore termoidraulico MOOSE / Cardinal
├── README.md                     # Main project guide, specifications, and execution instructions
├── PROJECT_STATUS.md             # Live tracking of current reference design, metrics, and open items
├── ITERATION_LOG.md              # Chronological engineering log of all design iterations
├── BIBLIOGRAPHY.md -> literature/BIBLIOGRAPHY.md # Complete literature citations & attribution matrix
├── src/                          # Primary Python source code and simulation modules
│   ├── model.py                  # Core parameterized geometry, materials, and settings builder
│   ├── analysis.py               # Physics post-processing, S-curves, Doppler, and void analysis
│   ├── depletion.py              # Depletion setup, chain handling, and isotopic tracking
│   ├── optimize.py               # Geometric sensitivity sweeps and optimization routines
│   ├── run_workflow.py           # End-to-end execution workflow runner
│   ├── run_iterations_master.py  # Historical iteration engine
│   ├── run_opt5_master.py        # Master execution pipeline for reference Opt-5 design
│   └── run_opt6_master.py        # Master execution pipeline for current reference Opt-6 design
├── data/                         # Nuclear data configs, reaction chains, and validated XML inputs
│   ├── chain_nep.xml             # Depletion reaction chain (transuranics & key fission products)
│   ├── materials.xml             # OpenMC materials definition XML
│   ├── geometry.xml              # OpenMC CSG geometry definition XML
│   ├── settings.xml              # OpenMC run settings and source particle distribution XML
│   └── tallies.xml               # OpenMC flux, fission, heating, and mesh tallies XML
├── literature/                   # Scientific literature base and detailed physical justifications
│   ├── literature.md             # Detailed engineering and physics justification document
│   └── BIBLIOGRAPHY.md          # Formal citations, NASA NTRS identifiers, and attribution matrix
├── results/                      # Structured design results organized by iteration
│   ├── baseline/                 # Baseline 0 / Baseline 1 data and inventories
│   ├── opt2/                     # Opt-2 reference data (5.5-year reference baseline)
│   ├── opt3/                     # Opt-3 reference data (initial 10-year baseline, 670.5 kg)
│   ├── opt4/                     # Opt-4 multi-objective core (622.7 kg, 145.2 L)
│   ├── opt5/                     # Opt-5 ultra-compact core (514.1 kg, 120.9 L)
│   ├── opt6/                     # Current Global Multi-Objective Optimum (514.2 kg, 118.8 L)
│   │   ├── geometry/             # Opt-6 CSG geometry plots and pin cell data
│   │   ├── neutronics/           # Material inventory CSV and integral reaction tallies CSV
│   │   ├── depletion/            # Criticality depletion history, isotopic curves, and HDF5 results
│   │   ├── temperature_coefficients/ # Continuous-energy Doppler, MTC, coolant, and total TMC
│   │   ├── void/                 # Coolant void reactivity and composition decomposition
│   │   ├── control_drums/        # 8-drum rotation calibration S-curve and stuck-drum check
│   │   ├── figures/              # High-resolution consolidated publication figures (300 DPI)
│   │   └── reports/              # Full engineering report (report_opt6.txt)
│   ├── design_comparison.csv     # Multi-generational comparison across all candidates
│   └── report.txt                # Symlink/copy of latest reference report
└── archive/                      # Historical scratch scripts, old statepoints, and scratch logs
```

---

## 5. General Methodology & Safety Requirements

The reactor design is governed by strict space nuclear engineering criteria:
1. **Criticality Condition during Mission ($k_\text{eff} \approx 1.0$)**:
   * During normal operation, the reactor operates critical. The 8 control drums are rotated dynamically from their initial inward orientation ($180^\circ$) toward the outward orientation ($0^\circ$) as fuel burnup progresses.
2. **Shutdown Margin (SDM)**:
   * With all 8 drums fully rotated inward ($180^\circ$), the core must achieve a cold shutdown target of $k_\text{in} \le 0.980$ (design target $\approx 0.970$, $\text{SDM} \ge 2,500\text{ pcm}$).
3. **Single Failure Law (1-Stuck-Drum Margin)**:
   * The reactor must remain subcritical ($k_\text{stuck} < 1.000$) even if the single most reactive control drum becomes permanently jammed in the fully outward position ($0^\circ$) while the remaining 7 drums rotate fully inward.
4. **Thermal Stability & Negative Reactivity Feedback**:
   * Prompt fuel Doppler feedback must be strictly negative ($\alpha_\text{fuel} < 0$).
   * Total isothermal temperature coefficient must be negative ($\alpha_\text{total} < -1.5\text{ pcm/K}$).
   * Coolant void reactivity upon depressurization to space vacuum must be non-positive ($\Delta\rho_\text{void} \le 0\text{ pcm}$).
5. **Operational Lifetime**:
   * Minimum required operational lifetime is **10.0 full-power years** ($3,650\text{ days}$).
   * The theoretical neutronic lifetime is determined by depleting until $k(\theta=0^\circ) < 1.000$.

---

## 6. How to Run the Model & Iterations

### Environment Setup
The simulation environment uses a dedicated Miniforge Python environment configured with OpenMC 0.16.0 and WMP nuclear data:
```bash
export PATH="/home/aleri/miniforge3/envs/openmc-env/bin:$PATH"
export OPENMC_CROSS_SECTIONS=/home/aleri/wmp_data/cross_sections_wmp.xml
export OMP_NUM_THREADS=8
```

### Running the Current Reference Model (Opt-6)
To execute the complete simulation and verification pipeline for Opt-6:
```bash
python src/run_opt6_master.py
```
This script runs the nominal eigenvalue simulation, 8-drum S-curve calibration, 1-stuck drum calculation, continuous-energy Doppler evaluation, void reactivity evaluation, depletion, and exports all CSVs and high-resolution plots to `results/opt6/`.

### Running Standalone Model Execution
To generate XML files and run a quick eigenvalue check from `src/model.py`:
```bash
python src/model.py
```

---

## 7. Current Reference Configuration Summary (Opt-6)

* **Assembly Mass**: **$514.21\text{ kg}$** ($-108.47\text{ kg}$ vs Opt-4).
* **Assembly Volume**: **$118.79\text{ L}$** ($-2.07\text{ L}$ vs Opt-5, $-26.44\text{ L}$ vs Opt-4).
* **Core Aspect Ratio ($H/R$)**: **$1.955$** ($H/D = 0.978$ — finite cylinder buckling optimum validated).
* **BOL Multiplication Factor ($0^\circ$ Out)**: $k_\text{out} = 1.02829 \pm 0.00182$ ($+2,751\text{ pcm}$ excess reserve).
* **Shutdown Multiplication Factor ($180^\circ$ In)**: $k_\text{in} = 0.96086 \pm 0.00258$ ($\text{SDM} = 4,073\text{ pcm}$).
* **Single Stuck-Drum Multiplication Factor**: $k_\text{stuck} = 0.97014 \pm 0.00198 < 1.000$ ($\text{SDM}_\text{stk} = 3,077\text{ pcm}$).
* **Total Control Drum Worth**: $6,614\text{ pcm}$ ($6.61\% \Delta k/k$).
* **10-Year Full-Power Operational Lifetime**: Verified across 3,650 days at 400 kWth ($k_\text{out} = 1.02055 \pm 0.00200$, $\theta_\text{crit} = 65.2^\circ$).
* **Theoretical Neutronic Lifetime**: **$25.0\text{ years}$** ($9,125\text{ days}$ of continuous full-power operation).
* **Prompt Fuel Doppler Feedback**: $\alpha_\text{fuel} = -1.96 \pm 0.15\text{ pcm/K}$ (Windowed Multipole).
* **Moderator Feedback**: $\alpha_\text{mod} = -0.54 \pm 0.33\text{ pcm/K}$ (strictly negative via under-moderation).
* **Total Temperature Coefficient**: $\alpha_\text{total} = -3.01 \pm 0.57\text{ pcm/K}$ ($3\sigma$ bound: **$-1.30\text{ pcm/K} < 0$**, inherently stable).
* **Coolant Void Reactivity**: $\Delta\rho_\text{void} = +292.6 \pm 152.8\text{ pcm}$ (He-Xe 72/28 baseline; negative $\approx -125\text{ to } -150\text{ pcm}$ with lighter Brayton mixtures).

See `PROJECT_STATUS.md` for live progress tracking and `ITERATION_LOG.md` for detailed evolutionary history.
 
---

## 8. Posizione e Dettagli del Risultato Finale: `definitive_final_geometry/`

Tutte le simulazioni più recenti e l'architettura definitiva del reattore, che incorpora la tecnologia a **matrice continua di $\text{YH}_{2.0}$ con raffreddamento a doppio passaggio (Dual-Pass)**, barriera in $\text{Al}_2\text{O}_3$ ($100\ \mu\text{m}$), combustibile $\text{UN}$ ($19.75\text{ wt\%}$), ciclo Brayton chiuso da **$136.39\text{ kWe}$** ($\eta = 34.10\%$) e vita operativa verificata a **$25.0\text{ anni}$** (estrapolata a **$44.8\text{ anni}$**), si trovano in:

📁 **[`definitive_final_geometry/`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry)**

### Documenti e File Salienti:
1. **Relazione Tecnica Ufficiale**: [`definitive_final_geometry/report.txt`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/report.txt)
2. **Documentazione Tecnica e Guida Completa**: [`definitive_final_geometry/README.md`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/README.md)
3. **Parametri Macchina Strutturati**: [`definitive_final_geometry/results.json`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/results.json)
4. **Galleria delle 13 Figure a 300 DPI**: [`definitive_final_geometry/figures/`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/)
   * Sezioni geometriche: [`01_reactor_radial_xy.png`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/01_reactor_radial_xy.png), [`02_reactor_axial_xz.png`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/02_reactor_axial_xz.png), [`03_unit_cell_5layer_detailed.png`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/03_unit_cell_5layer_detailed.png), [`04_control_drum_mechanism.png`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/04_control_drum_mechanism.png).
   * Rendering CSG: [`05_openmc_csg_slice_xy.png`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/05_openmc_csg_slice_xy.png), [`06_openmc_csg_slice_xz.png`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/06_openmc_csg_slice_xz.png).
   * Depletion (0-25 anni): [`07_depletion_keff_and_excess_reactivity.png`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/07_depletion_keff_and_excess_reactivity.png), [`08_depletion_control_drum_critical_angle.png`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/08_depletion_control_drum_critical_angle.png), [`09_depletion_isotopic_inventory_evolution.png`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/09_depletion_isotopic_inventory_evolution.png), [`10_depletion_burnup_and_fima.png`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/10_depletion_burnup_and_fima.png).
   * Termica FEM e Brayton: [`11_moose_cardinal_temperature_profiles_pins.png`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/11_moose_cardinal_temperature_profiles_pins.png), [`12_moose_cardinal_moderator_thermal_field.png`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/12_moose_cardinal_moderator_thermal_field.png), [`13_closed_brayton_cycle_ts_diagram.png`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/13_closed_brayton_cycle_ts_diagram.png).
5. **Dati Numerici di Burnup e Termica**: [`definitive_final_geometry/data/`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/data/)
