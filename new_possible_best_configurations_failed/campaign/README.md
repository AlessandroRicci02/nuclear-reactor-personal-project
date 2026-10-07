# P1 simulation campaign

This directory contains the campaign authorized in the user's attached request.
It is separate from the preceding analytic design-space screening and from all
historical results. Completion is established by result files, not by this plan.

## Scope and conventions

- P1 geometry is read from `../configurations/P1_ELLIPTIC_BALANCED/geometry.json`.
- Solid mass includes all material actually represented by the geometry.
  It is **not** total launch mass: pressure vessel, drives, manifolds, supports,
  missing hydrogen barriers, shielding, and power conversion are unmodeled.
- Neutronics: independent source convergence checks, control range and critical
  angle, fuel/moderator/reflector/clad temperature partial derivatives,
  moderator/coolant density derivatives, fixed-geometry isothermal derivative.
  Temperature derivatives hold other quantities fixed unless explicitly stated.
- Production targets: eigenvalue standard deviation at most 20 pcm; feedback
  cases at most 15 pcm per endpoint. These are Monte Carlo sampling errors,
  excluding nuclear-data, geometry and material-property uncertainty.
- Depletion: 400 kW recoverable fission power, thirty spatial fuel compositions
  (three radial zones by ten axial slices), full ENDF/B-VIII.0 depletion chain,
  internal integration steps shorter than the 2.5-year reporting intervals.
  Drums are searched within 0–180 degrees, with 0 pointing outwards. A state
  outside this control range cannot be reported as critical operation.
- Thermal calculations must include heating of moderator and heat transfer from
  the hot central channel to the moderator. The hydraulic boundary conditions
  are assumptions, not a validated closed Brayton cycle.
- Native Cardinal/OpenMC interface tests and file-based reactor thermal
  calculations are identified separately. No claim of a coupled solution is
  made merely because the executables launch.

## Existing installations reused

Standalone OpenMC: `/home/aleri/miniforge3/envs/openmc-env/bin/openmc`, version
0.16.0, commit `617d35a5063c57796b43428bc401e627d2011046`.

Cardinal: `/home/aleri/projects/cardinal/cardinal-opt`, repository commit
`9fb1dffcc095dd6ac9c87c6b16f44635b79cfc46`. It links to its bundled MOOSE framework,
commit `e6a7269a5c8e884ea82a7b64ac687043d219022c`, and to libMesh/PETSc in
`/home/aleri/miniforge3/envs/moose`. Its own OpenMC shared library is version
0.16.1-dev14, commit `150f9ac044031a840c555cfdb2bb2c7f63f23d7c`, built with MPI,
OpenMP and libMesh enabled. This is distinct from standalone OpenMC 0.16.0.

The official `openmc_statepoint` test has executed successfully in
`validation/cardinal_openmc_statepoint`. Changes to the test are recorded there.
This verifies the native interface, not yet reactor thermal feedback.

No solver installation was moved, copied, reinstalled or rebuilt.
Up to four threads for standalone OpenMC, two for small Cardinal tests; only one
production Monte Carlo calculation at a time. Audit evidence: `audit/`.

## Nuclear-data correction

The existing `/home/aleri/nndc_hdf5` and `/home/aleri/openmc-data/nndc_hdf5` have
only 294 K neutron tables for key nuclides, and room-temperature thermal
scattering tables. The WMP overlay does not supply temperature-dependent bound
scattering for YH2. The historical moderator uses `c_H_in_ZrH` as a substitute.
The existing `chain_nep.xml` contains only twelve nuclides.

These concrete limitations require additional data for the requested
temperature coefficients and depletion. Official OpenMC ENDF/B-VIII.0 HDF5
data and its full thermal-branching depletion chain were downloaded **inside
this campaign**; the existing libraries and global environment are unchanged.
Sources are recorded in `data/data_sources.json` and published on
[OpenMC's data page](https://openmc.org/data/).

Use actual hydrogen-in-YH2 and yttrium-in-YH2 scattering, with stochastic
temperature interpolation. Geometry-only comparisons require rerunning TB35
with the same data and corrections; historical and new values are not directly
interchangeable.
