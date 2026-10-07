# Resumed production campaign — 6 October 2026

The controller is `../scripts/p1_campaign.py`. It runs one production Monte
Carlo process at a time, with four OpenMP threads and one BLAS thread. Small
Cardinal thermal solves use one thread. Existing solver installations and
historical results are unchanged.

`status.json` gives the controller PID, active stage, child PID, completed
stages and any errors. `controller_job.json` identifies the initial detached
controller. `logs/campaign.log` and the stage logs contain execution output.
`../comparison.txt` is regenerated every 15 seconds from completed result files.
The report explicitly retains PENDING entries until solver results exist.

## Scheduled work

1. P1 all-out/all-in and common-data TB35 all-out keff, standard error <=20 pcm.
2. Independent P1 seed and longer inactive period for a source-convergence check.
3. P1 and TB35 drum searches. Independent final checks must put the entire
   nominal 95% reactivity interval within +/-60 pcm.
4. Neutron/photon heating calculations at those nominal-temperature drum angles.
5. Heating-driven MOOSE calculations executed by the existing `cardinal-opt`:
   XY moderator/liner meshes with the actual elliptical holes, axial main and
   auxiliary coolant energy balances, and nonlinear RZ fuel/clad FEM.
6. Independent finite-difference endpoints for fuel, moderator, reflector,
   cladding and coolant temperature, moderator and coolant density, and
   isothermal temperature feedback. Endpoint standard errors <=15 pcm.
7. CE/CM full-chain depletion with 30 fuel compositions and 400 kW recoverable
   fission power. Startup poison buildup has shorter steps; subsequent steps
   are <=0.25 equivalent full-power years, reported at 2.5-year intervals.
   Every transport stage adjusts drums within 0–180 degrees if needed. The
   final accepted 95% reactivity interval must fit within +/-75 pcm.

The depletion integration stops when the available drums cannot make the
current composition critical. The last completed and first failed interval
provide a lifetime bracket. It does not call an interpolated estimate an exact
lifespan. The queue continues until a control limit is reached rather than
assigning a lifetime from a predetermined reporting horizon. No critical trajectory is fabricated if the fresh reactor is
outside the available control range.

## Thermal scope

This is a reduced thermal-hydraulic model driven by OpenMC heating files,
executed using Cardinal's MOOSE heat-transfer module. The previous native
Cardinal/OpenMC interface test remains separate. These results are not a
native full-core 3-D CFD solution or a converged two-way reactor calculation.

The model solves heat exchange from the hot annulus through the liner and
moderator to auxiliary coolant, including direct solid heating. Flow split is
prescribed (17.5% auxiliary), not derived from unresolved plenum/orifice losses.
It computes channel friction losses, not total reactor pressure drop.
Peripheral moderator heating is assigned to outer-ring pin cells as an
explicit energy-network closure; peripheral solid temperatures are unresolved.
Reflector/absorber cooling is not modeled. Reported cell temperatures must not
be interpreted as the global moderator peak.

YH2 conductivity 10 W/(m K), Al2O3 conductivity 6 W/(m K), hydraulic-diameter
film correlations and ideal-gas helium are explicit assumptions. Sensitivity
cases use YH2 conductivity 5 and 20 W/(m K) and film coefficients scaled by
0.7 and 1.3. These are scenario variations, not statistical confidence bounds.
Both cell and peak-pin mesh refinement are checked; a >2 K change fails the
stage for review. Pin and global coolant energy balances are checked.

Neutron/photon thermal sources are normalized to 400 kW of summed material
heating; depletion uses 400 kW of recoverable fission energy. This distinction
is recorded in each heating summary, including heat outside the active model.

## Restart and inspection

From `new_possible_best_configurations`:

```bash
/home/aleri/miniforge3/envs/openmc-env/bin/python -u scripts/p1_campaign.py
```

The lock prevents two controllers. Completed stages and completed depletion
intervals are reused. An interrupted interval is recalculated from the previous
material checkpoint. A partially accumulated standalone transport run is
restarted, while completed precision-qualified points are reused.

The detached process survives closing the chat/terminal. It cannot compute
while the laptop sleeps, WSL stops, or Windows shuts down. Restart the command
above after WSL restarts. No startup service or system settings were changed.

## Verification and sources

- `audit/resume_2026-10-06.json`: current executable/version checks.
- `audit/environment.json`: original linked libraries and Cardinal build audit.
- `validation/cardinal_openmc_statepoint/`: successful official native interface test.
- `validation/thermal_synthetic/`: deliberately synthetic numerical verification,
  excluded from reactor performance reporting.
- `validation/depletion_api/`: deliberately low-statistics CE/CM/rotation/export
  verification, excluded from criticality and lifetime reporting.
- Reproduce numerical checks using `scripts/p1_validate_runtime.py depletion`,
  `scripts/p1_validate_runtime.py photons`, and `scripts/p1_thermal.py --validate`
  with the existing `openmc-env/bin/python` (from the parent of this directory).
  The full-chain API check uses about 6 GiB RAM; run it alone when repeating it.
- [OpenMC depletion methods](https://docs.openmc.org/en/stable/usersguide/depletion.html).
- [MOOSE convective heat integral](https://mooseframework.inl.gov/moose/source/postprocessors/ConvectiveHeatTransferSideIntegral.html).
- [Cardinal native solid/fluid coupling](https://cardinal.cels.anl.gov/tutorials/openmc_fluid.html),
  documenting the distinction from this file-based reduced calculation.

Monte Carlo standard errors exclude nuclear-data covariance, geometry,
manufacturing, property and model uncertainties. Historical TB35 feedback signs,
moderator temperature and depletion lifespan are reassessed from their actual
files in the comparison; their old PASS labels are not copied as validation.
