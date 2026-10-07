# Reproducible methods and interpretation

## Conserved quantities

All candidates have 217 pins at the reference pitch and the original density/composition prescriptions. Radial zones are rings 0–2 (19 pins), 3–5 (72), and 6–8 (126). With height H, reference H0=45.5 cm and r0=0.585 cm, radii are:

```
r_inner = r0 sqrt(H0/H)
r_middle = r0 sqrt((H0/H) f)
r_outer = r0 sqrt((H0/H) (198 - 72 f)/126)
```

This preserves inner fuel volume and combined middle/outer fuel volume separately. Enrichments remain 0.196 and 0.1975 uranium weight fraction. The uranium fraction 238/252 is intentionally retained from TB35 for an identical comparison; this study does not silently revise its approximate UN composition prescription.

Every pin has the same auxiliary area Aaux = 6 pi (0.150 cm)^2. For N ellipses with tangential-to-radial semiaxis ratio q:

```
a_radial = sqrt(Aaux/(N pi q))
b_tangential = q a_radial
```

Main annular area is conserved in the aggregate over 217 pins. For graded allocation the relative per-pin areas are 1.08, 1.04, and (217-19*1.08-72*1.04)/126. Each annular outer radius follows from its area and the zone-dependent clad radius. The clad and central liner thicknesses remain 0.048 and 0.010 cm.

Fixed cross-sectional area does not mean fixed coolant *volume* when height changes. That distinction is retained in all inventory calculations. No moderator is removed just by exchanging six circular auxiliary holes for equal-area ellipses; local scattering paths still change.

## Geometry and mass

Pin coordinates use the actual OpenMC orientation=y lattice basis. The occupied-cell convex hull is calculated from all pin-cell vertices. A contoured moderator envelope is the intersection of the core circle and the hull offset outward by the specified distance. Its area follows from deterministic polar quadrature. Zero-offset area is independently checked against the convex-hull polygon area.

Fuel, clad, liner and auxiliary openings lie wholly inside their pin cells for selected candidates. The enclosing core cylinder contains every occupied hex cell. Minimum ellipse-to-face clearance uses exact support functions. The auxiliary-to-auxiliary clearance uses a conservative enclosing-circle bound; the selected candidates have positive clearance under that bound. The reported smallest ligament is the minimum of those values and the channel-to-central-liner gap, over all three zones. It is a nominal geometric clearance, before tolerances or encapsulation additions.

Masses are density times deterministic CSG-equivalent volumes. Moderator is partitioned into pin cells and periphery once. BeO includes the radial reflector, drum filler, axial caps and the BeO replacing any trimmed peripheral hydride. B4C occupies the specified ten annular sectors, extending over the active height. There is no double-counting of the total-moderator row.

The axial main opening and auxiliary openings are integrated as a union. For each ellipse, the circle–ellipse overlap is integrated in its local radial coordinate; rotational symmetry makes the overlap independent of azimuth. Auxiliary holes are mutually disjoint. The six-circle reference overlap has an independent analytic circle-lens check. This recovers the discrepancy between the historical analytic and stochastic BeO inventories.

The selected designs retain overlapping *opening shapes* in the axial caps for comparison with the reference, but the union is one coolant volume. Detailed end manifolds, dividers and hot/cold routing are unresolved and will have to be represented for a physical flow topology.

## Geometry-based hydraulic and cooling indicators

For fully developed laminar flow through an elliptical tube of semiaxes a,b, the flow relation is:

```
Q = pi a^3 b^3 DeltaP / [4 mu L (a^2+b^2)]
```

For N identical equal-area passages in parallel, relative to the six-circle reference at the same total auxiliary flow and viscosity:

```
R_aux / R_aux,reference = (H/H0) (N/6) (q + 1/q)/2
```

This is a geometry indicator. It excludes developing-flow effects, roughness, gas property variation, compression, turns, orifices and plenums. No absolute core pressure drop is reported. The perimeter is integrated directly along the ellipse. Diameter-based heat-transfer correlations are not assumed to validate elliptical passages.

Cooling-distance indicators use a deterministic grid of moderator points inside a hex cell and a densely sampled channel boundary. The distance to the nearest auxiliary surface is obtained with a spatial tree. Reported maximum and 95th-percentile values are approximate geometric descriptors; neither is a thermal solution or a rigorous continuum upper bound. Validation doubles both area-grid and boundary resolution for selected middle-zone cells. Peripheral moderator is excluded from this proxy and needs its own heat-transfer calculation.

The branch-energy table conserves total heat and flow for a stated parallel-branch hypothesis. It shows how auxiliary coolant bulk outlet depends on the heat fraction and flow fraction. It does not calculate heat deposition or heat leakage; these must come from transport and conjugate heat transfer. 892/1110 K are comparison boundary targets, not a validated Brayton cycle.

## Selection and uncertainty

Mass, exact geometric areas and conservation identities are computed. Cooling distance is numerically sampled. Flow weights, screening clearances and the shortlist are engineering proposals. Multiplication, power distribution, feedback, lifetime, thermal peaks and electrical output are uncomputed.

The same source material densities are used for every design. Mass precision in CSV/JSON files supports reproduction and conservation checks; decimal places do not describe manufacturing or material-property uncertainty. No Monte Carlo or nuclear-data uncertainty was propagated into nonexistent reactor-performance estimates.

575 kg is an optional study objective leaving 25 kg unallocated below a 600 kg conceptual-solid limit. This is not a demonstrated hardware mass budget. P1 exceeds this preference and remains the recommended first experiment because it retains more reflector than P2; the preference is not presented as a hard physical constraint.

Full encapsulation of auxiliary walls, outer faces and end faces, structural canning, vessel, supports, drum housings/drives and manifolds are not invented or claimed to be present. Their inclusion can alter both mass and neutron economy. The historical 950 K moderator cap is treated as a project criterion, not a universal pressure-independent dissociation temperature.

Before later reactor simulations, retain the audit requirements from the earlier review: reconcile thermal-scattering treatment, eliminate CSG overlapping cells, use meaningful depletion data/time-step convergence, establish power tallies and energy balance, and enforce a feasible power-conversion operating point. These are future work, not silently performed in the present study.
