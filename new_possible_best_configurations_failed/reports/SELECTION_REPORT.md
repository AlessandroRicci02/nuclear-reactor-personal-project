# Geometry selection and recommended candidates

## Decision

**P1_ELLIPTIC_BALANCED is the recommended first simulation candidate.** It offers a controlled departure from TB35: essentially unchanged total moderator mass, exactly conserved reference fuel and fissile inventories, and unchanged total coolant cross-sectional areas. Its estimated solid mass is 578.96 kg, compared with 595.15 kg for the reference under identical area accounting.

This selection is engineering judgment about which experiment is most informative and least disruptive. It is **not** a claim of higher multiplication, longer lifetime, better shutdown margin or lower temperature. Those quantities require the simulations the user will select next.

P2 tests whether more reflector mass can be removed; P3 tests a stronger cooling-surface modification; P4 tests the height/peripheral-moderator tradeoff. They are alternatives, not successive proven improvements.

## What was explored

The deterministic scan contains 36,288 combinations:

| Variable | Values |
|---|---|
| Auxiliary shape/orientation | 6 circular channels at vertices or faces; 6 ellipses with tangential/radial semiaxis ratio 1.5 or 16/9 at vertices or faces; 12 circular channels staggered by 15 degrees |
| Channel-center distance from pin axis | 9.14009, 9.30, 9.50 mm |
| Active height | 44.5, 45.5, 46.5, 48.0 cm |
| Middle-zone fuel cross-sectional multiplier before height correction | 0.96, 1.00, 1.04 |
| Main annular coolant allocation | Uniform or three-zone graded; total area fixed |
| Cylindrical core radius | 19.43568 or 19.53568 cm |
| Peripheral moderator boundary | Full circle, or occupied-cell convex hull expanded by 1.5, 1.0 or 0.5 cm and clipped to the circle |
| Radial BeO thickness outside core circle | 8.8, 9.0, 9.3 cm |
| Each axial BeO thickness | 8.3, 8.8, 9.3 cm |

There are 29,088 combinations passing the study's geometric clearance screen. Of these, 24,138 have modeled solid mass at most 600 kg, and 7,098 are at most 575 kg. These counts are **geometric/inventory screens**, not numbers of viable reactors.

The 0.50 mm minimum ligament, 7 mm outer drum-to-envelope web, and optional 25 kg unallocated mass objective are transparent study choices. They are not fabrication standards or structural qualification limits. The historical reference is retained even though its 0.370 mm inner moderator ligament does not meet this study's preferred clearance.

No weighted score was used to convert moderator kilograms, coolant perimeter or core size into fictitious reactor performance. Low mass can simply indicate that valuable neutron-reflecting or moderating material has been removed. The 552.29 kg minimum in the grid therefore is not the recommended design.

## P1: the first candidate

| Parameter | TB35 reference | P1 |
|---|---:|---:|
| Pin count / pitch | 217 / 22.80 mm | unchanged |
| Active height | 45.50 cm | unchanged |
| Fuel-zone populations, rings 0–2 / 3–5 / 6–8 | 19 / 72 / 126 | 19 / 72 / 126 |
| Fuel radii, inner / middle / outer | 5.850 / 5.850 / 5.850 mm | 5.850 / 5.732 / 5.916 mm |
| Uranium enrichments, by zone | 19.60 / 19.75 / 19.75 wt% | unchanged |
| Main annular gaps, by zone | 0.840 / 0.840 / 0.840 mm | 0.903 / 0.886 / 0.805 mm |
| Auxiliary channels per cell | 6 circles, radius 1.500 mm | 6 tangential ellipses, radial/tangential **semiaxes** 1.225 / 1.837 mm |
| Channel-center radius | 9.140 mm | 9.300 mm |
| Channel azimuths | 0, 60, …, 300 degrees | unchanged |
| Auxiliary area per cell | 42.4115 mm2 | unchanged |
| Total main annular area | 77.3077 cm2 approximately | unchanged; see exact CSV value |
| Radial / each axial BeO thickness | 9.30 / 8.80 cm | 9.00 / 8.30 cm |
| Peripheral YH2 boundary | Full core circle | unchanged |
| Modeled solid mass | 595.15 kg | 578.96 kg |
| Unallocated difference to 600 kg | 4.85 kg | 21.04 kg |

The inner 19-pin fuel volume is held fixed. Four percent of the middle-zone fuel cross section is redistributed to the outer zone, whose cross section increases by 2.286%. Since those two zones have the same uranium enrichment, this conserves both total uranium and U-235, without exceeding the existing enrichment ceiling. The direction is an outward-redistribution hypothesis, not a result from current TB35 pin-power tallies. The opposite zoning direction is included in the scan and isolated controls.

The provisional per-pin main-flow targets, relative to the mean, are 1.08, 1.04 and 0.96508. The main channel areas use these same factors, preserving total area. These are starting targets for an orifice/common-pressure calculation, not an assertion that the passages naturally divide flow this way. Fresh power tallies must determine the eventual allocation.

Why this proposal is worth trying:

1. **No fuel or moderator inventory shortcut.** UN remains 143.3067 kg and U-235 remains 26.7129 kg under the reference composition convention. Total YH2 stays approximately 146.077 kg.
2. **Greater geometric clearance.** The smallest screened moderator ligament rises from 0.370 to 0.742 mm. Six holes are retained, limiting the increase in surface that must eventually be sealed.
3. **Moderate hydraulic change.** The auxiliary wetted perimeter rises 3.09%. For equal flow and viscosity, fully developed laminar auxiliary resistance rises 8.33%. These are geometry-derived ratios, not total-core pressure-drop predictions.
4. **A useful thermal hypothesis.** The sampled middle-cell 95th-percentile distance from moderator points to the nearest auxiliary surface falls from about 3.346 to 3.016 mm. This is a conduction-distance indicator only; it excludes the peripheral moderator and cannot predict temperature.
5. **A recoverable mass tradeoff.** Most of the 16.19 kg mass reduction comes from the BeO allocation. If later transport finds an unacceptable penalty, the radial and axial reflector changes can be reversed independently using the supplied controls.

Potential disadvantages are lower reflector worth, changed drum worth, increased outer-zone power, altered local moderation paths, and more difficult fabrication/coating of noncircular passages. Total moderator mass alone does not fix moderator feedback or reactivity. A negative result would still identify which of these changes should be retained.

The `A15_P1_ORIGINAL_REFLECTORS` control restores both original reflector thicknesses while keeping P1's channels, fuel zones and annular allocation. It provides a direct fallback and helps distinguish a reflector penalty from the effects of the cell redesign. Given the reference's limited reported end-of-life reserve, reflector trimming must remain conditional on transport/depletion results.

## The three alternatives

| Candidate | Change relative to P1 | Reason to retain it | Main unresolved tradeoff |
|---|---|---|---|
| **P2_ELLIPTIC_LIGHT** | Radial BeO thickness 8.80 cm | 572.53 kg modeled solids; 27.47 kg unallocated | Less reflection and an outer drum web of 7.2 mm instead of 9.2 mm; structural/drum consequences unknown |
| **P3_TWELVE_CHANNEL** | 12 circles of radius 1.06066 mm, center radius 9.30 mm, phase 15 degrees | 41.42% more auxiliary perimeter; sampled cell cooling distances improve further | Twice the ideal laminar auxiliary resistance; twice the hole count and more encapsulation surface |
| **P4_TALL_CONTOURED** | Height 46.50 cm; peripheral YH2 follows the occupied-cell hull plus 1.50 cm | Tests lower average linear power and less thick peripheral hydride while conserving fuel | Total YH2 falls to 141.59 kg; leakage, feedback and drum effectiveness may change unfavorably |

P4 fuel radii are 5.787 / 5.670 / 5.853 mm. Extending the fuel length by 1 cm reduces average power per unit length by 2.15% at unchanged total power and pin count. This is a geometric relation, not a predicted peak-temperature decrease. The area outside the contoured moderator but inside the original core cylinder becomes **BeO**, not vacuum. The cylindrical core boundary and drum centers remain unchanged. The contour is represented by explicit halfspaces in its JSON specification.

P3 has 578.82 kg modeled solids, approximately the same as P1. Doubling auxiliary resistance does not double total-core pressure drop: the effect depends on the full branch/manifold/orifice network. Likewise, a larger perimeter does not establish a larger heat-transfer coefficient for all flow regimes.

## Findings that changed the selection

**Actual channel orientation:** the TB35 source uses `HexLattice.orientation = 'y'` and auxiliary-channel azimuths 0+60k degrees. For that lattice the hexagonal face normals are 30+60k degrees. The existing code therefore places the channels toward the vertices, despite some schematic/report labels describing face alignment. Coordinate checks use the installed OpenMC API. P1 keeps the actual reference azimuths. Face-oriented alternatives are separate sensitivity cases, not silently substituted into the reference.

**Axial duct overlap:** the reference main reflector duct radius is 0.860 cm. Auxiliary centers are 0.914009 cm from the pin axis with radius 0.150 cm, so these openings overlap. The old analytic mass formula subtracts their areas twice. Counting their Boolean union once increases the calculated BeO mass by 1.22822 kg, giving 595.15359 kg total; this agrees with the archived stochastic mass within its uncertainty. The original CSG also contains overlapping coolant cells in these regions. Any future input deck should represent a single union or mutually exclusive cells. The new accounting does not add physical material.

**Peripheral mass:** all calculations include the moderator outside the 217 pin-cell polygons. Its baseline mass is 43.44 kg. A small reduction of the bounding cylinder is limited by outer-cell clearance; it cannot remove most of this mass. Contouring changes its distribution explicitly instead.

**Thermal margin:** for a provisional parallel split with 17.5% of the flow in auxiliary passages, depositing 1.5% of 400 kW into that branch raises its bulk outlet to 910.69 K. Only 4.31 K remain to a 915 K moderator target before any film or solid-conduction temperature rise. If the heat fraction is 3%, the bulk outlet alone is 929.37 K. The unknown heat transferred from the hot annulus is therefore consequential. The supplied energy table illustrates this sensitivity; no peak moderator temperature is claimed.

The parallel split is taken as an explicit screening hypothesis from the recent project scripts. Earlier documents also describe serial counter-current flow. The geometry alone does not resolve that inconsistency. Manifold routing and the flow topology must be settled before a coupled thermal calculation; the branch-energy numbers apply only to the stated parallel hypothesis.

## How to decide which reactor is actually better

Future selection should first establish feasibility, then compare non-dominated designs. Mandatory project constraints are 15 full-power years, at least 35 K moderator margin to the provisional 950 K limit (therefore **915 K**, not 920 K), all-drums and worst-stuck-drum margins above 2500/1500 pcm, negative total isothermal and power feedback, pressure drop below 35 kPa, and mass below the agreed 600 kg accounting boundary. These requirements are not established by this folder.

After those checks, compare end-of-life reactivity reserve, worst-case moderator temperature, credible allowance for real hardware, pumping demand and local power peaking. Include uncertainty and off-design conditions; a tiny nominal pass is less attractive than a reproducible margin. More BOL excess reactivity is useful only if shutdown authority and feedback remain acceptable.

The next calculations should compare the same repaired reference and P1 using consistent nuclear data, a suitable depletion chain, common energy deposition accounting, and a thermal model containing the moderator and hot/cold paths. All distinct stuck-drum positions deserve evaluation because ten-drum and sixfold lattice symmetries do not make every position equivalent. The thin central alumina ring must not be treated as demonstrated complete hydrogen encapsulation.

No such simulation has been performed here. The recommended order for the user's next choice is **reference + P1**, then **P2 for mass**, **P3 for cooling**, and **P4 for height/peripheral moderation**.
