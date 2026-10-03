# Autonomous Research & Qualification Plan: Opt-14 Hexagonal Architecture Exploration

## 1. Project Objectives & Mandatory Design Thresholds

The Opt-14 research cycle builds upon the audited benchmark of Candidate Opt-13, expanding into radially symmetric hexagonal lattices with dual-pass / counter-current gas coolant topologies.

### 1.1 Mandatory Performance Constraints
- **Thermal Power**: 400 kWth ($P_\text{elec} = 100$ kWe at 25% net Brayton efficiency).
- **Core Assembly Mass**: $< 650.0$ kg (evaluated at core & reflector boundary).
- **Full-Power Lifetime**: $\ge 10.0$ years verified by multi-step depletion solving.
- **Control Authority**:
  - All-in Shutdown Margin (SDM): $\ge 2000$ pcm.
  - Single Stuck-Drum Margin (SDM$_\text{stk}$): $\ge 1000$ pcm.
  - Total Control Drum Worth: $> 4000$ pcm.
- **Reactivity Feedback Coefficients**:
  - Fuel Doppler ($\alpha_\text{fuel}$): $< 0$ at $3\sigma$ confidence.
  - Moderator Feedback ($\alpha_\text{mod}$): $< 0$ at $3\sigma$ confidence.
  - Coolant Feedback ($\alpha_\text{cool}$): $\le 0$ nominal and $1\sigma$.
  - Total Core Isothermal TMC ($\alpha_\text{total}$): $< 0$ at $3\sigma$ confidence.
- **Coolant-Void Reactivity ($\Delta\rho_\text{void}$)**:
  - Strictly negative nominally AND targeting negative upper confidence bounds ($1\sigma$, $2\sigma$, $3\sigma$).
  - Never claim $2\sigma$ or $3\sigma$ negativity if the interval crosses zero.
- **Metallurgical Integrity & Fabricability**:
  - Zero unverified or physically unstable composites (no co-sintered BeO-hydride mixtures).
  - Hermetically canned delta-phase YH$_{2.0}$ moderator with Nb-1Zr cladding and sintered BeO reflectors.

---

## 2. Experimental Architecture Exploration Matrix

The exploration encompasses three distinct geometric channel paradigms across pitch and material dimensions:

### 2.1 The Three Hexagonal Topologies
1. **Architecture A (Shared Vertex Channels)**:
   - Auxiliary channels placed at cell corners (shared by 3 cells, effective $N_\text{aux} = 2$/cell).
   - Minimizes hydride moderator displacement.
2. **Architecture B (Internal Vertex-Aligned Channels)**:
   - Auxiliary channels placed inside the cell along radii towards the 6 corners ($N_\text{aux} = 6$/cell).
   - Maximizes moderator ligament thickness along the flat faces.
3. **Architecture C (Internal Face-Aligned Channels, Shifted 30°)**:
   - Auxiliary channels placed inside the cell along rays towards the 6 flat faces ($N_\text{aux} = 6$/cell).
   - Places coolant directly in the shortest inter-pin moderation corridors, maximizing axial streaming when voided.

### 2.2 Material & Coolant Variants
- **Fuels**: UN (reference, 13.50 g/cm³), UC (13.63 g/cm³), UO$_2$ (10.40 g/cm³).
- **Coolants**: Pure Helium (2.5 MPa, zero neutron capture) vs. He-Xe binary mixtures (72/28 mole%, Brayton-optimized but parasitically absorbing).

---

## 3. Qualification Protocol for the Downselected Candidate

Once the fine sweep confirms the superior architecture and geometry, the candidate will undergo comprehensive qualification:
1. **Criticality & Operating State**: $k_\text{op}$ (drums OUT, 0°), 200k particles $\times$ 80 batches.
2. **Control Drum S-Curve Calibration**: 7 drum angles ($0^\circ, 30^\circ, 60^\circ, 90^\circ, 120^\circ, 150^\circ, 180^\circ$) to map differential and integral worth.
3. **Stuck-Drum Safety Margin**: Criticality with 1 most-reactive drum stuck OUT at 0° while 9 drums rotate to 180°.
4. **Precision Coolant-Void Analysis**: High-statistics Monte Carlo comparison (200k $\times$ 80 batches) of nominal 2.5 MPa state vs. 0.0 MPa vacuum state.
5. **Reactivity Feedback Matrix**:
   - Prompt Fuel Doppler: Perturbation to 1500 K.
   - Moderator Temperature Coefficient: Perturbation to 1100 K.
   - Coolant Temperature Feedback: Perturbation from 800 K to 1200 K.
   - Total Isothermal Core TMC: Uniform core temperature rise.
6. **25-Year Full-Power Depletion Simulation**: Depletion at 400 kWth across 11 discrete timesteps (0, 0.5, 1, 2, 5, 8, 10, 12, 15, 20, 25 years), tracking $^{235}\text{U}$, $^{238}\text{U}$, $^{239}\text{Pu}$, $^{135}\text{Xe}$, $^{149}\text{Sm}$, and critical drum angle $\theta_\text{crit}(t)$.
7. **Spatial Flux & Power Density Tallies**: 3-group neutron flux mesh tallies (thermal, epithermal, fast) and 2D fission heat deposition mesh.
8. **Documentation & Deliverables**:
   - 12 high-resolution publication-quality PNG figures (300 DPI).
   - Plain-text comprehensive report (`opt14/reports/opt14_comprehensive_report.txt`).
   - Markdown summary (`opt14/OPT14_FINAL_REPORT.md`).
   - Comparison database update (`results/design_comparison.csv`).
