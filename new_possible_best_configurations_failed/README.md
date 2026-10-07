# Proposed successors to GEOM_OPT15_TB35

Prepared 5 October 2026. **Recommendation: P1_ELLIPTIC_BALANCED is the first candidate to simulate.** It preserves the reference fuel and moderator inventories, makes a modest change to channel geometry, introduces three fuel/flow zones, and reduces reflector mass. It is a proposal for testing, not a demonstrated better reactor.

This folder explores only the three requested directions: channel geometry at fixed coolant area, radial fuel/coolant zoning, and active-height/peripheral-moderator/reflector allocation. All candidates retain 217 pins, the 2.280 cm pitch, UN/N-14 fuel, the reference uranium enrichments, YH2, Nb-1Zr, Al2O3, He, BeO, and ten existing-size B4C/BeO drum cross sections. No inverted moderator-rod topology or nitrogen-isotope change is included.

**36,288 combinations were evaluated by deterministic geometry and inventory calculations. No neutron transport, depletion, CFD, MOOSE/Cardinal thermal solve, or Brayton simulation was run.** Reactor performance fields are explicitly null in the configuration files.

| Configuration | Modeled solid mass, kg | Total YH2, kg | Minimum screened ligament, mm | Purpose |
|---|---:|---:|---:|---|
| TB35_REFERENCE | 595.15 | 146.08 | 0.370 | Same physical reference, consistent union-area accounting |
| **P1_ELLIPTIC_BALANCED** | **578.96** | **146.08** | **0.742** | First choice: modest geometry changes and 16.19 kg modeled mass saving |
| P2_ELLIPTIC_LIGHT | 572.53 | 146.08 | 0.742 | Further reflector reduction; less neutron-reflection and drum-web margin |
| P3_TWELVE_CHANNEL | 578.82 | 146.08 | 0.906 | More auxiliary surface and shorter cooling distances; more hydraulic/sealing cost |
| P4_TALL_CONTOURED | 584.50 | 141.59 | 0.798 | Taller core plus shaped peripheral moderator; greater neutronic uncertainty |

Masses cover the same conceptual solid-material boundary as the TB35 model. Vessel, supports, complete moderator encapsulation, manifolds, actuators and shield are not included. The remaining difference to 600 kg is **unallocated mass**, not proof that these items fit. The historical stochastic reference was 595.02 ± 0.37 kg; the present 595.15 kg is an independent deterministic reconciliation, not a new Monte Carlo result.

Start with [the selection report](reports/SELECTION_REPORT.md), then [the calculation methods and limits](reports/METHODS_AND_LIMITS.md). [Sources](reports/SOURCES.md) distinguish literature principles from this study's own design choices.

The main deliverables are:

- [P1 geometry specification](configurations/P1_ELLIPTIC_BALANCED/geometry.json) and [217-pin map](configurations/P1_ELLIPTIC_BALANCED/pins.csv): explicit dimensions, material densities, zones, channel orientation, reflector boundaries and provisional flow targets.
- The equivalent files in each of the other candidate directories and the reference directory.
- [Shortlist comparison](data/shortlist.csv), [full parameter scan](data/design_space.csv), and [screening counts](data/screening_summary.json).
- [Fifteen isolated-change controls plus P1 with restored reflectors](data/isolated_changes.csv) and their [specifications](configurations/isolated_changes.json), to separate the effects of combined changes in later work.
- [Sampled cooling-distance indicators](data/cooling_distance_proxies.csv) and [branch energy-balance sensitivities](data/branch_energy_balance_assumptions.csv); neither is a predicted moderator temperature.
- [Geometry/accounting verification](data/validation.json) and [SHA-256 hashes of source files](data/input_sha256.json).
- [Cell comparison](figures/unit_cells.png), [core zones and contoured periphery](figures/core_zones_and_periphery.png), and [tradeoff plots](figures/design_tradeoffs.png). SVG versions are also provided.

Reproduce from this directory using the existing OpenMC Python environment:

```bash
/home/aleri/miniforge3/envs/openmc-env/bin/python scripts/explore.py
/home/aleri/miniforge3/envs/openmc-env/bin/python scripts/validate.py
```

The first command regenerates only files in this folder; the second checks the calculations. OpenMC is imported by validation only to check the lattice coordinate convention. The JSON files specify proposed geometries; they are not OpenMC XML input decks or MOOSE meshes.

The next decision is which candidates and simulation campaign to authorize. P1 and the consistent TB35 reference form the most informative initial pair. P2, P3 and P4 answer distinct follow-up questions.
