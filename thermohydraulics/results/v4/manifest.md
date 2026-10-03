# V4 CAMPAIGN EXECUTION MANIFEST & REPRODUCIBILITY AUDIT
## Opt-17 / Cand 17-Clean Multiphysics Coupling & JIMO Assessment Campaign

**Campaign ID**: `OPT17-V4-COUPLING-ASSESSMENT`  
**Execution Date**: October 3, 2026  
**Working Directory**: `/home/aleri/OpenMC_projects/first_test/thermohydraulics/`  
**Deliverables Root**: `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v4/`  

---

## 1. SOFTWARE VERSIONS & EXECUTION ENVIRONMENT

| Component | Version / Identifier | Full Path / Provenance |
| :--- | :--- | :--- |
| **Operating System** | Linux 6.6.137+ x86_64 | GNU/Linux Core Environment |
| **OpenMC** | 0.16.0 (Commit `617d35a5063c57796b43428bc401e627d2011046`) | `/home/aleri/miniforge3/envs/openmc-env/bin/openmc` |
| **Python Environment** | Python 3.13 (Miniforge3 `openmc-env`) | `/home/aleri/miniforge3/envs/openmc-env/bin/python3` |
| **Cardinal / MOOSE** | snapshot-20-10-27-56855-ge6a7269a5c | `/home/aleri/projects/cardinal/cardinal-opt` |
| **MOOSE C++ Wrapper** | OpenMPI `mpicxx` GNU 15.3.0 | `/home/aleri/miniforge3/envs/moose/bin` (`PATH` prepended) |
| **Nuclear Data Library** | Continuous Windowed Multipole (WMP) ENDF/B-VIII.0 | `/home/aleri/wmp_data/cross_sections_wmp.xml` |

---

## 2. INPUT & SCRIPT FILE SHA-256 CHECKSUMS

| File Path | SHA-256 Checksum | Description |
| :--- | :--- | :--- |
| `inputs/build_opt17_v4.py` | `a0ebbd4adee9233dd889b939c0fe9096456684403ff0aaff904ebbee3cb1b156` | OpenMC Cand 17-Clean V4 builder with drum/density/temperature perturbations |
| `coupling/v3_picard_driver.py` | `26ec3ebbff1e3b95ff1e87730d4924e5a668da1a32c6a29f2c6e174f77396933` | 8-iteration Picard under-relaxation multiphysics driver |
| `models/cardinal_th_interface.py` | `427b8dd62caca82f750ba448dde1de63e20b273a65bfa3571745f8a57b6fd92b` | 1D compressible flow channel & 2D conduction solver with Cardinal interface |
| `models/th_properties.py` | `79750e540b8432f8869ebb8aa890ea6f2a54073fee0695ea6d0b4ded9b932237` | Helium (Petersen 1981), UN, Nb-1Zr, BeO thermophysical property evaluations |
| `scripts/compute_v4_mass_breakdown.py` | `776de3543ac2564c987c8ca0a823efb235f72690c9a0f820cc808a329627a8b8` | Transparent component mass accounting script |
| `scripts/run_v4_perturbations.py` | `d8acce97032e919a49e39101ceb9ab378dbd84b29b44ca42e98889bbb7746774` | 8-case neutronic perturbation & feedback coefficient solver |
| `scripts/run_v4_sensitivities.py` | `ce12525b1e399f13e3796b0036cddd7a5f2ae0e0a5749b635d153e17c250426c` | Thermal-hydraulic sensitivity sweeps and local derivative calculations |
| `scripts/create_v4_reference_comparison.py` | `754800381f8b26200cc51ae0e6cf773ff493d19ce7d35521e732239165ac29d9` | NASA Prometheus / JIMO requirements matrix generator |
| `scripts/create_v4_issue_register.py` | `779eafb4d0d53b74778d04e998401c7f3fdd7441c4d3b19ea84a22ee898feb08` | Issue register and validation matrix generator |
| `scripts/create_v4_summary.py` | `22b661f4aa0ab2488a6929d399180c52c892424caf7dd9af3403afaffa581929` | Summary JSON assembly routine |
| `scripts/generate_v4_plots.py` | `baf5bd5e8ef197202e328b82953bf3086c58a45d789f92f4a6ac6af9796bff41` | Publication figure generator (10 plots) & ParaView VTK export |

---

## 3. EXECUTION COMMANDS RECORD

1. **Mass Breakdown Accounting**:
   ```bash
   /home/aleri/miniforge3/envs/openmc-env/bin/python3 scripts/compute_v4_mass_breakdown.py
   ```
2. **Neutronic Perturbations & Coefficients**:
   ```bash
   /home/aleri/miniforge3/envs/openmc-env/bin/python3 scripts/run_v4_perturbations.py
   ```
3. **Thermal-Hydraulic Sensitivities & Derivatives**:
   ```bash
   /home/aleri/miniforge3/envs/openmc-env/bin/python3 scripts/run_v4_sensitivities.py
   ```
4. **Reference Comparison Dataset**:
   ```bash
   /home/aleri/miniforge3/envs/openmc-env/bin/python3 scripts/create_v4_reference_comparison.py
   ```
5. **Issue Register & Validation Matrix**:
   ```bash
   /home/aleri/miniforge3/envs/openmc-env/bin/python3 scripts/create_v4_issue_register.py
   ```
6. **Machine-Readable Summary**:
   ```bash
   /home/aleri/miniforge3/envs/openmc-env/bin/python3 scripts/create_v4_summary.py
   ```
7. **Visualization Suite & VTK Export**:
   ```bash
   /home/aleri/miniforge3/envs/openmc-env/bin/python3 scripts/generate_v4_plots.py
   ```

---

## 4. VERIFIED DELIVERABLES DIRECTORY MAP

All deliverables are verified to exist in `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v4/`:

| Output Deliverable | Format | Status | Description |
| :--- | :---: | :---: | :--- |
| `V4_FINAL_REPORT.md` | Markdown | Verified | Comprehensive auditable engineering report |
| `v4_summary.json` | JSON | Verified | Machine-readable metrics and status flags |
| `v4_iterations.csv` | CSV | Verified | 27 tracked metrics across 8 Picard iterations |
| `v4_neutronics.csv` | CSV | Verified | Converged k_eff, excess reactivity, lifetime, peaking factors |
| `v4_feedback_coefficients.csv` | CSV | Verified | Doppler, moderator, void, drum worth, SDM with 1-sigma errors |
| `v4_thermal_hydraulics.csv` | CSV | Verified | Complete TH field, margins, Reynolds range, boundary audits |
| `v4_mass_breakdown.csv` | CSV | Verified | Transparent component mass breakdown and system estimates |
| `v4_reference_comparison.csv` | CSV | Verified | 24-parameter comparison matrix vs NASA Prometheus / JIMO |
| `v4_sensitivity_studies.csv` | CSV | Verified | Parametric sweeps across mass flow, gap, inlet T, and HTC |
| `v4_issue_register.csv` | CSV | Verified | Ranked-by-priority technical issue register |
| `v4_validation_matrix.csv` | CSV | Verified | V&V domain evidence and remaining modeling gaps |
| `opt17_v4_coupled_field.vtk` | ASCII VTK | Verified | ParaView 3D dataset of all 2,170 cells (power, T_f, T_c, T_cool) |
| `figures/` (10 files) | PNG (300 DPI)| Verified | Complete suite of high-resolution engineering figures |
