"""Supplemental YH2-density sensitivity for interpreting moderator feedback."""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import subprocess

import openmc

import build_model as b
from run_bol import ROOT


RUNS = ROOT / "campaign" / "density_sensitivity"
BASE = dict(pitch=2.30, core_h=47.0, refl_rad_th=8.8,
            refl_ax_th=8.5, drum_angles=[0.0] * 10,
            particles=18000, batches=180, inactive=50,
            n_axial_slices=6)


def run(multiplier, seed):
    label = "plus2pct" if multiplier > 1 else "minus2pct"
    directory = RUNS / label
    directory.mkdir(parents=True, exist_ok=True)
    result_path = directory / "result.json"
    if result_path.exists():
        return json.loads(result_path.read_text())
    model = b.build_openmc_model_tb35(str(directory), **BASE,
                                      mod_density_mult=multiplier, seed=seed)
    settings = model["settings"]
    settings.keff_trigger = {"type": "std_dev", "threshold": 0.00030}
    settings.trigger_active = True
    settings.trigger_max_batches = 650
    settings.trigger_batch_interval = 20
    settings.export_to_xml(directory / "settings.xml")
    environment = os.environ.copy()
    environment["OPENMC_CROSS_SECTIONS"] = b.CROSS_SECTIONS_PATH
    environment["OMP_NUM_THREADS"] = "8"
    with (directory / "run.log").open("w") as log:
        completed = subprocess.run([b.OPENMC_EXEC, "-s", "8"], cwd=directory,
                                   stdout=log, stderr=subprocess.STDOUT,
                                   env=environment, check=False)
    if completed.returncode:
        raise RuntimeError(f"OpenMC {label} failed; see {directory/'run.log'}")
    statepoint = max(directory.glob("statepoint.*.h5"),
                     key=lambda p: int(p.name.split(".")[1]))
    with openmc.StatePoint(statepoint) as sp:
        k = float(sp.keff.nominal_value)
        sigma = float(sp.keff.std_dev)
        result = {"case": label, "moderator_density_multiplier": multiplier,
                  "seed": seed, "k_eff": k, "sigma_k": sigma,
                  "reactivity_pcm": 1e5 * (1-1/k),
                  "sigma_reactivity_pcm": 1e5*sigma/k**2,
                  "statepoint": str(statepoint)}
    result_path.write_text(json.dumps(result, indent=2)+"\n")
    print(label, k, "+/-", sigma, flush=True)
    return result


def main():
    minus = run(.98, 2107)
    plus = run(1.02, 2108)
    per_fraction = (plus["reactivity_pcm"]-minus["reactivity_pcm"])/.04
    sigma_per_fraction = math.hypot(plus["sigma_reactivity_pcm"],
                                    minus["sigma_reactivity_pcm"])/.04
    nominal_beta = 3.6e-5
    beta_range = [2.4e-5, 4.95e-5]
    result = {"definition": "central difference of reactivity at fixed moderator temperature 875 K, YH2 density changed by +/-2%; geometry fixed",
              "reactivity_sensitivity_pcm_per_fractional_density": per_fraction,
              "sampling_sigma_pcm_per_fractional_density": sigma_per_fraction,
              "reactivity_sensitivity_pcm_per_percent_density": per_fraction/100,
              "assumed_volumetric_expansion_1_per_K": nominal_beta,
              "exploratory_beta_range_1_per_K": beta_range,
              "estimated_density_feedback_pcm_per_K": -nominal_beta*per_fraction,
              "sampling_sigma_density_feedback_pcm_per_K": nominal_beta*sigma_per_fraction,
              "estimated_feedback_range_beta_only_pcm_per_K":
                [-max(beta_range)*per_fraction, -min(beta_range)*per_fraction],
              "caveat": "An auxiliary sensitivity, not a qualified physical expansion model: YH2 composition and geometry expansion are not modeled. The beta range is a rough 3x conversion of ORNL-reported linear CTE 8-16.5e-6/K for YHx (x 1.60-1.95), not a measurement for this YH2 component.",
              "source": "https://info.ornl.gov/sites/publications/Files/Pub160401.pdf"}
    (ROOT / "density_sensitivity_result.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
