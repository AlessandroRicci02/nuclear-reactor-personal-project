# Scientific Literature Review & Thermochemical Feasibility Analysis: Opt-14 Space NEP Reactor

## 1. Introduction and Architectural Motivation

Space Nuclear Electric Propulsion (NEP) missions operating in the 100 kWe / 400 kWth regime require ultra-compact reactor cores capable of continuous, autonomous operation for 10 to 20+ years within stringent launch mass constraints (<650 kg total assembly mass). In earlier iterations (Opt-1 through Opt-13), square-pitch lattices with single-pass annular gas channels demonstrated feasibility, culminating in Candidate Opt-13.

However, three fundamental engineering and physical challenges motivate the exploration of **radially symmetric hexagonal fuel-cell architectures with dual-pass / counter-current gas cooling**:
1. **Thermal-Hydraulic Azimuthal Symmetry**: In square lattices, the moderator volume surrounding each cylindrical pin possesses 4-fold symmetry, creating non-uniform azimuthal temperature distributions across the fuel cladding and moderator blocks. A hexagonal lattice exhibits 6-fold symmetry, yielding significantly more uniform heat transfer and lower thermal stress.
2. **Dual-Pass / Counter-Current Flow Mechanics**: Integrating auxiliary coolant channels into each hexagonal moderator cell enables a two-pass coolant flow. Cold gas (e.g. He at ~800 K) first enters the core via auxiliary channels, keeping the solid hydride moderator cool, before reversing in the bottom plenum to ascend through the high-temperature fuel annular channels (reaching ~1200 K). This isolates the moderator from peak fuel temperatures.
3. **Coolant Void Reactivity Management**: Gas-cooled fast and epithermal space reactors often suffer from small or positive void coefficients if parasitic absorption or moderation effects dominate. A dual-pass hexagonal design introduces dedicated axial streaming paths that enhance neutron leakage upon depressurization, driving the coolant-void reactivity strongly negative.

---

## 2. Thermochemical Stability & Fabricability of Moderator Materials

A core directive of the user requires rigorous evaluation of material fabricability, specifically examining whether multi-phase composite mixes (such as BeO + ZrH$_2$) are metallurgically and chemically viable.

### 2.1 The Myth of Co-Sintered BeO-ZrH$_2$ or BeO-YH$_2$ Composites
In powder metallurgy and nuclear materials engineering, mixing metal hydrides with oxide ceramics to form a single "composite" block faces severe thermodynamic and chemical barriers:
1. **Hydrogen Dissociation vs. Ceramic Sintering Temperatures**:
   - Sintering high-density monolithic BeO (>98% theoretical density) requires thermal processing between **$1750^\circ\text{C}$ and $1950^\circ\text{C}$** in reducing or inert atmospheres.
   - Zirconium hydride ($\delta$-phase ZrH$_{1.6}$) begins substantial hydrogen outgassing above **$500^\circ\text{C}$** and dissociates completely before reaching $900^\circ\text{C}$ at atmospheric pressure.
   - Yttrium hydride ($\delta$-phase YH$_{1.8-2.0}$) exhibits superior thermal stability compared to ZrH$_2$, but decomposes rapidly above **$950^\circ\text{C} - 1050^\circ\text{C}$** at equilibrium hydrogen overpressures below 1 atm.
   - **Conclusion**: Co-sintering a BeO powder compact containing hydride grains at $\sim 1800^\circ\text{C}$ is physically impossible; violent hydrogen evolution would disintegrate the ceramic matrix.
2. **High-Temperature Solid-State Redox Reactions**:
   - Thermodynamically, zirconium and yttrium metals possess high oxygen affinities ($\Delta G^\circ_f$ of ZrO$_2$ and Y$_2$O$_3$ are strongly negative).
   - At elevated temperatures, contact between BeO and metal hydrides promotes reduction of BeO by the metallic constituent:
     $$\text{Zr} + 2\text{BeO} \longrightarrow \text{ZrO}_2 + 2\text{Be}$$
     $$\text{Y} + \frac{3}{2}\text{BeO} \longrightarrow \frac{1}{2}\text{Y}_2\text{O}_3 + \frac{3}{2}\text{Be}$$
   - This reaction degrades both the moderator hydrogen content and the integrity of the oxide matrix, producing low-melting Be metal ($T_m \approx 1287^\circ\text{C}$).

### 2.2 Qualified Engineering Solutions from Historical & Modern Space Programs
Aerospace nuclear history (SNAP-10A, SNAP-8, SP-100, NASA Prometheus, SAFE-400) provides two distinct, fully qualified architectural paradigms:

#### Paradigm A: Hermetically Canned Metal Hydride (The YH$_{2.0}$ Reference)
- **Material**: Delta-phase Yttrium Hydride (YH$_{1.8} - \text{YH}_{2.0}$).
- **Why YH$_2$ over ZrH$_2$**: Yttrium hydride maintains an order-of-magnitude lower hydrogen dissociation pressure than zirconium hydride at equivalent temperatures (Mueller et al., *Metal Hydrides*, Academic Press; General Atomics reports GA-8600). At 1000 K, the equilibrium hydrogen dissociation pressure of YH$_{1.9}$ is $\sim 10^{-2}$ bar, compared to $>10$ bar for ZrH$_{1.8}$.
- **Canning & Containment**: YH$_2$ is fabricated as pre-hydrided solid blocks, hermetically encapsulated in thin Nb-1Zr or Mo-47.5Re canning sleeves with an internal hydrogen barrier (e.g. Al$_2$O$_3$ or glass-ceramic coating).
- **Core Role**: Canned YH$_2$ serves as the internal core moderator block, while pure sintered BeO is utilized exclusively in the axial and radial reflectors where temperatures are moderate ($<1000$ K) and zero hydrogen is present.

#### Paradigm B: Monolithic Sintered BeO Matrix (Pure Ceramic Core)
- **Material**: High-purity sintered Beryllium Oxide (BeO, density 2.86–2.95 g/cm³).
- **Advantages**: 100% chemically inert, stable up to $2200^\circ\text{C}$, zero hydrogen dissociation risk, excellent thermal conductivity ($k \approx 40-70$ W/m·K at 1000 K).
- **Neutronic Limitation**: Beryllium has a neutron slowing-down power $\xi \Sigma_s$ an order of magnitude lower than hydrogen ($\xi \Sigma_s \approx 0.125\text{ cm}^{-1}$ for BeO vs. $1.35\text{ cm}^{-1}$ for YH$_2$). Achieving criticality in a BeO-moderated thermal core requires a large moderator-to-fuel volume ratio ($V_\text{mod}/V_\text{fuel} > 5$), expanding the active core volume to $>200$ L and assembly mass to $>800$ kg.
- **Synthesis**: In compact space reactors (<650 kg limit), canned YH$_{2.0}$ remains the only moderation technology capable of satisfying both the mass envelope and multi-decade lifetime requirements simultaneously. Pure BeO is the premier reflector material.

---

## 3. High-Density Fuel Candidates: UN vs. UC vs. UO$_2$

| Property | Uranium Mononitride (UN) | Uranium Monocarbide (UC) | Uranium Dioxide (UO$_2$) |
| :--- | :--- | :--- | :--- |
| **Heavy Metal Density [g/cm³]** | **13.50** (95% TD) | **12.97** (95% TD) | **9.65** (95% TD) |
| **Uranium Volume Fraction** | Highest (0.944) | High (0.952) | Moderate (0.881) |
| **Thermal Conductivity at 1200 K** | **23.5 W/m·K** | **21.0 W/m·K** | **2.8 W/m·K** |
| **Max Centerline Temp. Margin** | High ($T_m \approx 2850^\circ\text{C}$) | High ($T_m \approx 2400^\circ\text{C}$) | Lower ($T_m \approx 2865^\circ\text{C}$) |
| **Cladding Compatibility** | Excellent with Nb-1Zr | Risk of carburization | Risk of oxidation |
| **Parasitic Neutron Absorption** | $^{14}\text{N}(n,p)^{14}\text{C}$ requires $^{15}\text{N}$ enrichment or thin spectrum | Negligible carbon capture | Moderate oxygen slowing |

- **Selection**: UN enriched in $^{15}\text{N}$ (or naturally bound in hard epithermal spectra) provides the highest fissile packing density, superior thermal conductivity, and optimal chemical compatibility with refractory Nb-1Zr cladding.

---

## 4. Coolant Fluid Selection: Pure Helium vs. Helium-Xenon Binary Mixtures

The choice of coolant gas plays a profound role in both thermal-hydraulic design and coolant-void reactivity:
1. **Helium-Xenon Mixtures (e.g. He-Xe 72/28 mole%, MW ~ 40 g/mol)**:
   - *Thermal-Hydraulics*: Excellent for closed-loop Brayton turbomachinery because the higher molecular weight reduces the number of compressor/turbine stages.
   - *Neutronics*: **Severe disadvantage**. Natural Xenon contains isotopes ($^{129}\text{Xe}$, $^{131}\text{Xe}$) with non-trivial thermal/epithermal absorption cross sections. At 2.5 MPa, Xenon acts as a parasitic absorber in the core. Upon complete coolant depressurization (voiding), the removal of this absorber introduces a **large positive reactivity insertion** ($\Delta\rho_\text{void} \approx +256$ pcm in screening)! This violates the mandatory negative void reactivity criterion.
2. **Pure Helium ($^4\text{He}$, MW = 4.00 g/mol)**:
   - *Thermal-Hydraulics*: Highest specific heat capacity ($c_p = 5193$ J/kg·K) and highest thermal conductivity ($k \approx 0.35$ W/m·K at 1000 K).
   - *Neutronics*: Zero thermal absorption ($\sigma_\gamma < 0.001$ mb). Helium acts solely as an elastic scatterer. When evacuated, the removal of He eliminates small-angle scattering, allowing neutrons to stream axially out of the annular and auxiliary channels into the vacuum plenums.
   - *Result*: **Consistently negative void reactivity** ($\Delta\rho_\text{void} < 0$), satisfying safety criteria.

---

## 5. Architectural Topology: Comparative Analysis of Arch A, B, and C

The three hexagonal cell architectures investigated in Opt-14 represent distinct geometric allocations of the counter-current coolant channels:

```
    Architecture A               Architecture B               Architecture C
 (Shared Vertex Channels)   (Internal Vertex-Aligned)    (Internal Face-Aligned)

          o                           *                            
       /     \                     /     \                      /    *    \
     /    O    \                 /    O    \                  /      O      \
    o           o               *           *                *               *
     \         /                 \         /                  \      *      /
       \     /                     \     /                      \         /
          o                           *                            
```

1. **Architecture A (Shared Vertex Channels)**:
   - Each channel is positioned at a hexagon vertex, shared by 3 adjacent cells.
   - Effective auxiliary channels per cell: $6 \times \frac{1}{3} = 2$.
   - *Neutronic Impact*: Minimizes moderator displacement. Delivers highest initial $k_\text{eff}$ and longest lifetime, but lower total coolant cross-sectional area and less pronounced axial streaming.
2. **Architecture B (Internal Vertex-Aligned Channels)**:
   - Six independent channels placed inside each cell along the radii pointing to vertices ($\theta = 30^\circ, 90^\circ, \dots$).
   - Effective channels per cell: 6.
   - *Neutronic Impact*: Moderate moderator displacement. Channels are situated at the corners where moderator thickness is greatest.
3. **Architecture C (Internal Face-Aligned Channels, Shifted 30°)**:
   - Six independent channels placed inside each cell along the rays pointing towards the flat faces ($\theta = 0^\circ, 60^\circ, \dots$).
   - Effective channels per cell: 6.
   - *Neutronic Impact*: Channels sit midway between the central pin and the nearest neighboring cell's fuel pin. This arrangement breaks the shortest moderation paths between pins, significantly enhancing axial streaming out of the lattice.
   - *Coolant Void Behavior*: Exhibits the most negative nominal void reactivity ($\Delta\rho_\text{void} \approx -110\text{ to } -150$ pcm in screening), making Architecture C the premier candidate for achieving statistically robust 3-sigma negative void reactivity.
