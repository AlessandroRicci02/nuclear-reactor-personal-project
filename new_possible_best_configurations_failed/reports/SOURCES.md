# Sources and attribution

Sources were reviewed on 5 October 2026. The channel dimensions, grid, mass screens and candidate choices in this folder are **new study proposals**, not published validated designs.

1. **A. Peakman and R. Gregg (2020), “The challenges of gas-cooled reactor technology for space propulsion and the development of the JANUS space reactor concept.”** [Author manuscript](https://arxiv.org/html/2004.12116v1). Supports considering radial power shaping, active-core aspect ratio and hydraulic/system penalties together. Its coated fuel, coolant and enrichment differ from TB35; no performance numbers are transferred to these candidates.

2. **B. J. Ade et al. (2022), “Reactor Physics Considerations for Use of Yttrium Hydride Moderator,” Nuclear Science and Engineering 196, 1539–1558.** [ORNL publication record](https://www.ornl.gov/publication/reactor-physics-considerations-use-yttrium-hydride-moderator); [paper DOI](https://doi.org/10.1080/00295639.2022.2035157). Establishes that hydride geometry can materially affect temperature feedback. It motivates explicit feedback checks after redistributing moderator; it does not establish the coefficients of these UN-pin configurations.

3. **OpenMC `Quadric` API.** [Official documentation](https://docs.openmc.org/en/stable/pythonapi/generated/openmc.Quadric.html). Provides a future CSG representation for rotated elliptical cylindrical holes. The installed OpenMC 0.16.0 `HexLattice.get_local_coordinates` method was used directly to verify pin coordinates. No geometry export or transport calculation was run.

4. **Local source model and evidence:** `best_configurations/GEOM_OPT15_TB35/scripts/build_model_tb35.py`, its mass/depletion/thermal scripts and mass CSV files; `PROJECT_STATUS.md`; `ITERATION_LOG.md`; `data/chain_nep.xml`. Their hashes are saved in `data/input_sha256.json`. The baseline source code determines dimensions and materials when schematic labels disagree. Earlier physics claims are not accepted as measured properties of the new candidates.

Ellipse area, support-function clearances, union-area integration, conservation identities and the laminar elliptical-duct relation are explicit analytical calculations documented in METHODS_AND_LIMITS.md. They are not substitutes for neutronic or thermal-hydraulic simulation.
