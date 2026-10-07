# TB35 OpenMC reference campaign

This is a new, TB35-only campaign. `configuration/TB35_REFERENCE` is a copied
reference specification; the historical TB35 directory is never altered.

The campaign runs a 5,000,000-sample OpenMC stochastic volume calculation, a
BOL all-drums-out eigenvalue calculation, and the three requested temperature
coefficients: fuel Doppler, moderator at fixed density, and an isothermal total
at fixed geometry/densities. Each eigenvalue endpoint targets a one-standard-
deviation keff uncertainty of 10 pcm. The serial controller uses four OpenMP
threads and writes `results.txt` after each completed stage.

The model uses the official ENDF/B-VIII.0 library and actual YH2 scattering
tables through links to the data preserved in `new_possible_best_configurations_failed`.
This avoids a duplicate multi-gigabyte nuclear-data library.

`results.txt` distinguishes the fixed-density moderator coefficient from a
physical thermo-mechanical expansion feedback. The latter requires verified
YH2 density-versus-temperature data and a specified canning constraint.
