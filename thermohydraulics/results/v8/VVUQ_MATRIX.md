# V8 Verification, Validation & Uncertainty Quantification (VVUQ) Matrix

**Standard Compliance**: ASME V&V 20 (Standard for Verification and Validation in Computational Fluid Dynamics and Heat Transfer) & ASME V&V 40 (Assessing Credibility of Computational Models)  
**Workspace**: `/home/aleri/OpenMC_projects/first_test/thermohydraulics`  
**Project**: 400 kWth / ~100 kWe Space Nuclear Electric Propulsion Reactor  

---

## 1. Epistemic Principles & Definitions

To prevent unwarranted claims and ensure scientific integrity, the V8 engineering campaign rigorously distinguishes between:

1. **Verification (Solving the equations right)**:  
   Mathematical confirmation that the computational algorithms correctly solve the discretized physical equations as intended (grid convergence, code verification, solver tolerances).
2. **Validation (Solving the right equations)**:  
   Assessment of model fidelity by comparing computational predictions against independent physical experiments or accepted benchmarks.
3. **Uncertainty Quantification (Quantifying error and scatter)**:  
   Rigorous statistical characterization of parameter variations, numerical errors, and model-form uncertainties.

---

## 2. Comprehensive VVUQ Classification Matrix

| Physics Discipline | Computational Tool | Verification Status | Validation Benchmark | Uncertainty Method | Epistemic Provenance |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Neutron Transport & Criticality** | OpenMC 0.15.0-dev (Continuous-energy Monte Carlo, ENDF/B-VII.1) | 40,000 particles/batch, 25 batches, Shannon entropy convergence checked | ICSBEP Fast Benchmark (HMF-001 Godiva / HMF-028) | Monte Carlo standard error $\sigma_k = \pm 0.00111$ | `DIRECTLY_CALCULATED` |
| **Reactivity Temperature Feedback** | OpenMC temperature-interpolated cross-sections | Cross-section Doppler broadening tested across 300 - 1500 K | Benchmark Doppler experiments (SEFOR fast reactor) | Two-point perturbation Monte Carlo difference | `DIRECTLY_CALCULATED` |
| **Depletion & Burnup Lifetime** | Analytical depletion model normalized to Opt-16 numerical benchmark | Linear extrapolation verified against OpenMC depletion run | EBR-II UN fuel irradiation database | Extrapolation model-form uncertainty ($\pm 2.5\text{ yr}$) | `DERIVED` |
| **Helium Thermal-Hydraulics** | 1D Multi-Channel Navier-Stokes + Gnielinski Transitional Correlation | Grid convergence verified across 10 to 50 axial mesh slices | NIST REFPROP helium thermophysical properties; Petukhov-Roizen benchmark | Cross-correlation sensitivity (Gnielinski vs Churchill vs Petukhov) | `DIRECTLY_CALCULATED` |
| **Multi-Zone Orificing** | Analytical ring optimizer + non-linear hydraulic loop solver | Flow conservation $\sum \dot{m}_i = M_{\text{tot}}$ converged to $< 10^{-7}\%$ | Idelchik handbook of hydraulic resistance | $\pm 0.05\text{ mm}$ orifice hole tolerance propagation | `DIRECTLY_CALCULATED` |
| **3D Plenum CFD & Maldistribution** | Reduced-order 3D potential flow + stagnation recovery network | Grid independent across 217 channel node points | NASA TM-105748 Space Reactor Plenum Experiments | Baffled vs unbaffled parametric bounding | `MODEL-DEPENDENT` |
| **Fuel Pellet-Clad Interaction (PCMI)** | BISON-equivalent UN swelling + Norton creep ODE solver | Time-step refinement convergence ($\Delta t = 0.1\text{ yr}$) | Preusser & Ross UN irradiation swelling database; NASA CR-185119 Nb-1Zr creep | Parametric swelling rate bounds ($0.4\%$ to $0.8\% / 10^{20}\text{ f/cm}^3$) | `DIRECTLY_CALCULATED` |
| **YH2 Moderator Permeation** | Mueller-Vetrano dissociation isotherm + Richardson diffusion | Analytical mass balance integration | Mueller & Vetrano YH1.8 thermodynamic data; Al2O3 ALD barrier tests | PRF coating degradation sweep (500 to 10,000) | `LITERATURE_PROPERTY` + `MODEL-DEPENDENT` |
| **Closed Brayton Cycle (CBC)** | Real gas state-point cycle solver | Enthalpy conservation checked at every station ($\Delta E = 0$) | NASA Prometheus CBC testbed demonstration (100 kWe) | Turbomachinery component efficiency bounds | `DIRECTLY_CALCULATED` |
| **Waste Heat Radiator** | Stefan-Boltzmann numerical integration along gas cooler | Area integration converged across 100 slices | Kilopower / Prometheus deployable titanium-water heat pipe data | Emissivity degradation and sink temperature variation | `DIRECTLY_CALCULATED` + `ESTIMATED` |

---

## 3. Residual Uncertainties Requiring Experimental Testing (`REQUIRES EXPERIMENT`)

While the computational physics models are verified, the following engineering elements cannot be resolved by simulation alone and must be qualified experimentally:

1. **Hydrogen Permeation Barrier Adhesion under Fast Fluence**:  
   Thermal cycling and fast neutron displacement damage ($> 10^{21}\text{ n/cm}^2$) on the $\text{Al}_2\text{O}_3$ atomic-layer barrier coating on Nb-1Zr.
2. **High-Temperature Nb-1Zr Joint Creep / Braze Integrity**:  
   Long-term biaxial thermal creep rupture of electron-beam welded header joints at $1235\text{ K}$ under $2.5\text{ MPa}$ helium differential pressure.
3. **Plenum Flow Jet Acoustic Coupling**:  
   Acoustic resonance and flow-induced vibration (FIV) in thin annular channels under high-speed helium impingement.
