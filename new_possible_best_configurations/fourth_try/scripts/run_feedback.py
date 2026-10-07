"""Restartable BOL temperature-coefficient campaign for fourth_try (ZrH1.2 moderator), drums OUT."""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time

import numpy as np
import openmc

import build_model as b
from run_bol import ROOT


RUNS = ROOT / "campaign" / "feedback"
RUNS.mkdir(parents=True, exist_ok=True)
BASE = dict(pitch=2.30, core_h=47.0, refl_rad_th=8.8,
            refl_ax_th=8.5, drum_angles=[0.0] * 10,
            particles=18000, batches=180, inactive=50,
            n_axial_slices=6)
TARGET_SIGMA_K = 0.00030
MAX_BATCHES = 650
DELTA_K = 100.0


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def run_case(kind, shift, seed):
    name = "reference" if kind == "reference" else f"{kind}_{'plus' if shift > 0 else 'minus'}"
    directory = RUNS / name
    directory.mkdir(parents=True, exist_ok=True)
    request = {"kind": kind, "shift_K": shift, "seed": seed,
               "base": BASE, "target_sigma_k": TARGET_SIGMA_K,
               "max_batches": MAX_BATCHES, "cross_sections": b.CROSS_SECTIONS_PATH,
               "builder_sha256": hashlib.sha256(Path(b.__file__).read_bytes()).hexdigest()}
    signature = hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest()
    result_path = directory / "result.json"
    if result_path.exists():
        result = json.loads(result_path.read_text())
        if result["request_signature"] != signature:
            raise RuntimeError(f"Stale result for {name}; move it aside before rerunning")
        return result

    options = dict(BASE, seed=seed)
    if kind in ("fuel", "total"):
        options["T_fuel_offset"] = shift
    if kind in ("moderator", "total"):
        options["T_mod_offset"] = shift
    if kind == "total":
        options["T_clad_offset"] = shift
        options["T_liner_offset"] = shift
        options["T_refl_offset"] = shift
        midpoint = (np.arange(BASE["n_axial_slices"]) + 0.5) / BASE["n_axial_slices"]
        profile = b.T_INLET_NOM + midpoint * (b.T_OUTLET_NOM - b.T_INLET_NOM)
        options["T_cool_grid"] = np.tile(profile + shift, (b.N_PINS, 1))

    model = b.build_openmc_model_tb35(str(directory), **options)
    settings = model["settings"]
    settings.temperature = {"multipole": False, "method": "interpolation", "tolerance": 10.0}
    settings.keff_trigger = {"type": "std_dev", "threshold": TARGET_SIGMA_K}
    settings.trigger_active = True
    settings.trigger_max_batches = MAX_BATCHES
    settings.trigger_batch_interval = 20
    settings.export_to_xml(directory / "settings.xml")
    write_json(directory / "request.json", request)
    write_json(RUNS / "status.json", {"status": "running", "case": name,
                                      "started_unix": time.time()})
    environment = os.environ.copy()
    environment["OPENMC_CROSS_SECTIONS"] = b.CROSS_SECTIONS_PATH
    environment["OMP_NUM_THREADS"] = "8"
    with (directory / "run.log").open("w") as log:
        completed = subprocess.run([b.OPENMC_EXEC, "-s", "8"], cwd=directory,
                                   stdout=log, stderr=subprocess.STDOUT,
                                   env=environment, check=False)
    if completed.returncode:
        raise RuntimeError(f"OpenMC {name} failed; see {directory/'run.log'}")
    statepoint = max(directory.glob("statepoint.*.h5"),
                     key=lambda p: int(p.name.split(".")[1]))
    with openmc.StatePoint(statepoint) as sp:
        k = float(sp.keff.nominal_value)
        sigma = float(sp.keff.std_dev)
        result = {"case": name, "request_signature": signature,
                  "k_eff": k, "sigma_k": sigma,
                  "reactivity_pcm": 1e5 * (1 - 1/k),
                  "sigma_reactivity_pcm": 1e5 * sigma / k**2,
                  "particles_per_batch": int(sp.n_particles),
                  "active_batches": int(sp.current_batch - sp.n_inactive),
                  "inactive_batches": int(sp.n_inactive),
                  "target_achieved": sigma <= TARGET_SIGMA_K,
                  "statepoint": str(statepoint)}
    write_json(result_path, result)
    print(name, k, "+/-", sigma, "active batches", result["active_batches"],
          "target", result["target_achieved"], flush=True)
    return result


def main():
    cases = [("reference", 0.0, 2100),
             ("fuel", -DELTA_K, 2101), ("fuel", DELTA_K, 2102),
             ("moderator", -DELTA_K, 2103), ("moderator", DELTA_K, 2104),
             ("total", -DELTA_K, 2105), ("total", DELTA_K, 2106)]
    results = {}
    for kind, shift, seed in cases:
        result = run_case(kind, shift, seed)
        results[result["case"]] = result
    coefficients = []
    for kind in ("fuel", "moderator", "total"):
        minus, plus = results[f"{kind}_minus"], results[f"{kind}_plus"]
        alpha = (plus["reactivity_pcm"] - minus["reactivity_pcm"]) / (2 * DELTA_K)
        sigma = math.hypot(plus["sigma_reactivity_pcm"],
                           minus["sigma_reactivity_pcm"]) / (2 * DELTA_K)
        coefficients.append({"kind": kind, "alpha_pcm_per_K": alpha,
                             "sigma_pcm_per_K": sigma,
                             "confidence_95_normal_pcm_per_K": [alpha - 1.96*sigma,
                                                                  alpha + 1.96*sigma],
                             "delta_K": DELTA_K,
                             "minus_case": minus["case"], "plus_case": plus["case"]})
    summary = {"definition": "central difference of rho=(k-1)/k at BOL, all drums OUT; static geometry and solid density; total varies every material temperature and helium density at fixed pressure",
               "moderator_material": "ZrHx (x=1.2)",
               "moderator_thermal_scattering": ["c_H_in_ZrH", "c_Zr_in_ZrH"],
               "reference": results["reference"], "coefficients": coefficients,
               "cross_sections": b.CROSS_SECTIONS_PATH}
    write_json(ROOT / "feedback_results.json", summary)
    write_json(RUNS / "status.json", {"status": "complete", "completed_unix": time.time()})
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        write_json(RUNS / "status.json", {"status": "failed", "error": str(error),
                                          "failed_unix": time.time()})
        raise

