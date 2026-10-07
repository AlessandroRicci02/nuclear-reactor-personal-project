"""Conditional temperature+density feedback for C3 using assumed YH2 expansion."""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import subprocess

import numpy as np
import openmc

import build_model as b
from run_bol import ROOT
from run_feedback import BASE, DELTA_K, MAX_BATCHES, TARGET_SIGMA_K, write_json


RUNS = ROOT / "campaign" / "coupled_feedback"
RUNS.mkdir(parents=True, exist_ok=True)
BETA = 3.6e-5  # assumed volumetric YH2 expansion per K; conditional scenario


def run_case(kind, shift, seed):
    name = f"{kind}_{'plus' if shift > 0 else 'minus'}"
    directory = RUNS / name
    directory.mkdir(parents=True, exist_ok=True)
    request = {"kind": kind, "shift_K": shift, "beta_YH2_per_K": BETA,
               "seed": seed, "base": BASE, "cross_sections": b.CROSS_SECTIONS_PATH,
               "builder_sha256": hashlib.sha256(Path(b.__file__).read_bytes()).hexdigest()}
    signature = hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest()
    path = directory / "result.json"
    if path.exists():
        result = json.loads(path.read_text())
        if result["request_signature"] != signature:
            raise RuntimeError(f"Stale result for {name}")
        return result
    options = dict(BASE, seed=seed, T_mod_offset=shift,
                   mod_density_mult=1-BETA*shift)
    if kind == "total":
        options.update(T_fuel_offset=shift, T_clad_offset=shift,
                       T_liner_offset=shift, T_refl_offset=shift)
        midpoint = (np.arange(BASE["n_axial_slices"]) + .5) / BASE["n_axial_slices"]
        coolant = b.T_INLET_NOM + midpoint * (b.T_OUTLET_NOM-b.T_INLET_NOM)
        options["T_cool_grid"] = np.tile(coolant+shift, (b.N_PINS, 1))
    model = b.build_openmc_model_tb35(str(directory), **options)
    settings = model["settings"]
    settings.keff_trigger = {"type": "std_dev", "threshold": TARGET_SIGMA_K}
    settings.trigger_active = True
    settings.trigger_max_batches = MAX_BATCHES
    settings.trigger_batch_interval = 20
    settings.export_to_xml(directory / "settings.xml")
    write_json(directory / "request.json", request)
    with (directory / "run.log").open("w") as log:
        environment = os.environ.copy()
        environment["OPENMC_CROSS_SECTIONS"] = b.CROSS_SECTIONS_PATH
        environment["OMP_NUM_THREADS"] = "8"
        completed = subprocess.run([b.OPENMC_EXEC, "-s", "8"], cwd=directory,
                                   env=environment, stdout=log,
                                   stderr=subprocess.STDOUT, check=False)
    if completed.returncode:
        raise RuntimeError(f"OpenMC {name} failed; see {directory/'run.log'}")
    statepoint = max(directory.glob("statepoint.*.h5"),
                     key=lambda p: int(p.name.split(".")[1]))
    with openmc.StatePoint(statepoint) as sp:
        k, sigma = float(sp.keff.nominal_value), float(sp.keff.std_dev)
        result = {"case": name, "request_signature": signature,
                  "k_eff": k, "sigma_k": sigma,
                  "reactivity_pcm": 1e5*(1-1/k),
                  "sigma_reactivity_pcm": 1e5*sigma/k**2,
                  "active_batches": int(sp.current_batch-sp.n_inactive),
                  "statepoint": str(statepoint)}
    write_json(path, result)
    print(name, k, "+/-", sigma, flush=True)
    return result


def main():
    coefficients = []
    for index, kind in enumerate(("moderator", "total")):
        minus = run_case(kind, -DELTA_K, 2201+2*index)
        plus = run_case(kind, DELTA_K, 2202+2*index)
        alpha = (plus["reactivity_pcm"]-minus["reactivity_pcm"])/(2*DELTA_K)
        sigma = math.hypot(plus["sigma_reactivity_pcm"],
                           minus["sigma_reactivity_pcm"])/(2*DELTA_K)
        coefficients.append({"kind": kind,
                             "alpha_pcm_per_K": alpha,
                             "sampling_sigma_pcm_per_K": sigma,
                             "confidence_95_normal_pcm_per_K":
                               [alpha-1.96*sigma, alpha+1.96*sigma],
                             "minus_case": minus["case"], "plus_case": plus["case"]})
    summary = {"definition": "central temperature difference at BOL, all drums OUT; YH2 density adjusted as rho(T)=rho_nominal*(1-beta*(T-875K)); geometry fixed; total also varies all material temperatures and helium density at fixed pressure",
               "assumed_volumetric_expansion_1_per_K": BETA,
               "temperature_shift_K": DELTA_K,
               "caveat": "Conditional scenario, not a qualified thermoelastic or hydrogen-retention model; no geometric expansion is included",
               "coefficients": coefficients}
    write_json(ROOT / "coupled_feedback_results.json", summary)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
