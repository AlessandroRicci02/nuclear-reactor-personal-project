"""Persistent OpenMC search for a YH2 NEP neutronic candidate.

The production loop has no candidate-count stop: it exits successfully only
after OUT>1, IN<1 and all three conservative coupled coefficients pass their
one-sided 95% statistical tests. Results are checkpointed case by case.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import random
import subprocess
import sys
import time

import numpy as np
import openmc


ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parents[1]
BUILDER_DIR = PROJECT / "new_possible_best_configurations" / "third_try" / "scripts"
sys.path.insert(0, str(BUILDER_DIR))
import build_model as b  # noqa: E402
from summarize import main as summarize_progress  # noqa: E402


BETA_MIN = 2.4e-5  # /K, weakest assumed volumetric expansion
BETA_MAX = 4.95e-5 # /K, strongest handbook-based screening bound
BASE = dict(r_fuel=.585, clad_th=.048, chan_gap=.084, al2o3_mod_th=.010,
            pitch=2.30, core_h=47.0, refl_ax_th=8.5, refl_rad_th=8.8,
            drum_r=3.95, abs_th=1.50, abs_angle=148.0,
            drum_clearance=.180, r_add=.150,
            reflector_kind="Be", matrix_kind="YH2", peripheral_kind="YH2")
PILOT = dict(particles=4000, batches=100, inactive=30,
             n_axial_slices=2, delta_K=150.0, threads=8)
FINAL = dict(particles=20000, batches=200, inactive=50,
             n_axial_slices=6, delta_K=100.0, threads=8,
             target_sigma_k=.00022, max_batches=1600)


CATALOG = [
    ("matrix_compact_be", dict(matrix_kind="YH2", peripheral_kind="Be", pitch=2.16,
                               core_h=57, refl_rad_th=11, refl_ax_th=10, r_add=.10)),
    ("matrix_compact_beo", dict(matrix_kind="YH2", peripheral_kind="Be", pitch=2.16,
                                core_h=57, refl_rad_th=11, refl_ax_th=10,
                                r_add=.10, reflector_kind="BeO")),
    ("matrix_mid_be", dict(matrix_kind="YH2", peripheral_kind="Be", pitch=2.24,
                           core_h=54, refl_rad_th=10, refl_ax_th=9.5, r_add=.13)),
    ("matrix_mid_beo", dict(matrix_kind="YH2", peripheral_kind="Be", pitch=2.24,
                            core_h=54, refl_rad_th=10, refl_ax_th=9.5,
                            r_add=.13, reflector_kind="BeO")),
    ("matrix_wide_be", dict(matrix_kind="YH2", peripheral_kind="Be", pitch=2.38,
                            core_h=52, refl_rad_th=10, refl_ax_th=9.5, r_add=.16)),
    ("matrix_wide_beo", dict(matrix_kind="YH2", peripheral_kind="Be", pitch=2.38,
                             core_h=52, refl_rad_th=10, refl_ax_th=9.5,
                             r_add=.16, reflector_kind="BeO")),
    ("full_compact_be", dict(pitch=2.16, core_h=54, refl_rad_th=10,
                             refl_ax_th=9, r_add=.10)),
    ("full_compact_beo", dict(pitch=2.16, core_h=54, refl_rad_th=10,
                              refl_ax_th=9, r_add=.10, reflector_kind="BeO")),
    ("full_mid_be", dict(pitch=2.24, core_h=50, refl_rad_th=9.5,
                         refl_ax_th=9, r_add=.13)),
    ("full_mid_beo", dict(pitch=2.24, core_h=50, refl_rad_th=9.5,
                          refl_ax_th=9, r_add=.13, reflector_kind="BeO")),
    ("full_open_be", dict(pitch=2.45, core_h=47, refl_rad_th=9.5,
                          refl_ax_th=9, r_add=.16)),
    ("full_open_beo", dict(pitch=2.45, core_h=47, refl_rad_th=9.5,
                           refl_ax_th=9, r_add=.16, reflector_kind="BeO")),
    ("fuelrich_matrix_be", dict(matrix_kind="YH2", peripheral_kind="Be",
                                pitch=2.30, r_fuel=.60, core_h=58,
                                refl_rad_th=11.5, refl_ax_th=10.5, r_add=.12)),
    ("fuelrich_full_beo", dict(pitch=2.30, r_fuel=.60, core_h=53,
                               refl_rad_th=10.5, refl_ax_th=10,
                               r_add=.12, reflector_kind="BeO")),
    ("peripheral_be", dict(matrix_kind="Be", peripheral_kind="YH2",
                           pitch=2.32, core_h=68, refl_rad_th=12,
                           refl_ax_th=11, r_add=.15, abs_angle=170.0)),
    ("peripheral_beo", dict(matrix_kind="Be", peripheral_kind="YH2",
                            pitch=2.32, core_h=68, refl_rad_th=12,
                            refl_ax_th=11, r_add=.15, abs_angle=170.0,
                            reflector_kind="BeO")),
    ("long_sparse_be", dict(matrix_kind="YH2", peripheral_kind="Be",
                            pitch=2.12, core_h=70, refl_rad_th=12,
                            refl_ax_th=11, r_add=.08, abs_angle=170.0)),
    ("long_sparse_beo", dict(matrix_kind="YH2", peripheral_kind="Be",
                             pitch=2.12, core_h=70, refl_rad_th=12,
                             refl_ax_th=11, r_add=.08, abs_angle=170.0,
                             reflector_kind="BeO")),
]

# The initial matrix-only cores have ample OUT excess reactivity but the
# baseline drums cannot shut several of them down.  Sweep absorber coverage
# explicitly before falling back to the broad deterministic parameter search.
TUNED = [
    ("matrix_compact_strong_be", dict(matrix_kind="YH2", peripheral_kind="Be",
        pitch=2.16, core_h=57, refl_rad_th=11, refl_ax_th=10, r_add=.10,
        abs_angle=172.0, abs_th=1.8)),
    ("matrix_compact_strong_beo", dict(matrix_kind="YH2", peripheral_kind="Be",
        pitch=2.16, core_h=57, refl_rad_th=11, refl_ax_th=10, r_add=.10,
        abs_angle=172.0, abs_th=1.8, reflector_kind="BeO")),
    ("matrix_mid_strong_be", dict(matrix_kind="YH2", peripheral_kind="Be",
        pitch=2.24, core_h=54, refl_rad_th=10, refl_ax_th=9.5, r_add=.13,
        abs_angle=172.0, abs_th=1.8)),
    ("matrix_mid_strong_beo", dict(matrix_kind="YH2", peripheral_kind="Be",
        pitch=2.24, core_h=54, refl_rad_th=10, refl_ax_th=9.5, r_add=.13,
        abs_angle=172.0, abs_th=1.8, reflector_kind="BeO")),
    ("matrix_short_strong_be", dict(matrix_kind="YH2", peripheral_kind="Be",
        pitch=2.16, core_h=47, refl_rad_th=10, refl_ax_th=9.5, r_add=.10,
        abs_angle=172.0, abs_th=1.8)),
    ("matrix_short_strong_beo", dict(matrix_kind="YH2", peripheral_kind="Be",
        pitch=2.16, core_h=47, refl_rad_th=10, refl_ax_th=9.5, r_add=.10,
        abs_angle=172.0, abs_th=1.8, reflector_kind="BeO")),
]


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def catalog_entry(index):
    if index < len(CATALOG):
        name, updates = CATALOG[index]
        return name, BASE | updates
    if index < len(CATALOG) + len(TUNED):
        name, updates = TUNED[index-len(CATALOG)]
        return name, BASE | updates
    rng = random.Random(20261007 + index)
    family = rng.choice(["full", "matrix", "peripheral"])
    pitch = rng.choice([2.12, 2.16, 2.20, 2.24, 2.30, 2.36, 2.42, 2.48, 2.54])
    aux = {2.12:.08, 2.16:.10, 2.20:.11, 2.24:.12, 2.30:.14,
           2.36:.15, 2.42:.16, 2.48:.17, 2.54:.18}[pitch]
    kind = dict(full=("YH2", "YH2"), matrix=("YH2", "Be"),
                peripheral=("Be", "YH2"))[family]
    spec = BASE | dict(pitch=pitch, core_h=rng.choice([47, 52, 57, 62, 68, 74]),
                       refl_rad_th=rng.choice([9.5, 10.5, 11.5, 12.5]),
                       refl_ax_th=rng.choice([9, 10, 11, 12]),
                       r_fuel=rng.choice([.55, .57, .585, .60]),
                       r_add=aux, reflector_kind=rng.choice(["Be", "BeO"]),
                       matrix_kind=kind[0], peripheral_kind=kind[1],
                       abs_angle=rng.choice([148.0, 160.0, 172.0]),
                       abs_th=rng.choice([1.5, 1.8]))
    return f"adaptive_{index:04d}_{family}", spec


def validate_geometry(spec):
    d = .910*(spec["pitch"]/2.270)
    main_outer = sum(spec[k] for k in ("r_fuel", "clad_th", "chan_gap", "al2o3_mod_th"))
    if d-spec["r_add"] <= main_outer+.015:
        return "auxiliary He channel overlaps or approaches main pin layers"
    if d+spec["r_add"] >= spec["pitch"]/2-.025:
        return "auxiliary He channel approaches lattice boundary"
    if spec["refl_rad_th"] <= spec["drum_clearance"]+2*spec["drum_r"]+.20:
        return "control drums exceed radial reflector"
    if spec["matrix_kind"] == "Be" and spec["peripheral_kind"] == "Be":
        return "no YH2 moderator"
    return None


def case_options(spec, case, fidelity):
    options = spec | {k:v for k,v in fidelity.items()
                      if k in ("particles", "batches", "inactive", "n_axial_slices")}
    options["drum_angles"] = [180.0 if case == "in" else 0.0]*10
    if case in ("out", "in"):
        return options
    sign = -1 if case.endswith("minus") else 1
    shift = sign*fidelity["delta_K"]
    if case.startswith("fuel") or case.startswith("total"):
        options["T_fuel_offset"] = shift
    if case.startswith("moderator") or case.startswith("total"):
        options["T_mod_offset"] = shift
        options["mod_density_mult"] = 1-fidelity.get("beta", BETA_MIN)*shift
    if case.startswith("total"):
        options.update(T_clad_offset=shift, T_liner_offset=shift,
                       T_refl_offset=shift)
        s = fidelity["n_axial_slices"]
        fraction = (np.arange(s)+.5)/s
        coolant = b.T_INLET_NOM + fraction*(b.T_OUTLET_NOM-b.T_INLET_NOM)
        options["T_cool_grid"] = np.tile(coolant+shift, (b.N_PINS, 1))
    return options


def run_case(name, spec, case, stage, fidelity, sequence):
    directory = ROOT / "campaign" / name / stage / case
    directory.mkdir(parents=True, exist_ok=True)
    request = dict(candidate=name, spec=spec, case=case, stage=stage,
                   fidelity=fidelity, beta_min=BETA_MIN,
                   builder_sha256=hashlib.sha256(Path(b.__file__).read_bytes()).hexdigest(),
                   xs=b.CROSS_SECTIONS_PATH)
    if "beta" in fidelity:
        request["beta"] = fidelity["beta"]
    signature = hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest()
    path = directory / "result.json"
    if path.exists():
        cached = json.loads(path.read_text())
        if cached["request_sha256"] == signature:
            return cached
        raise RuntimeError(f"stale cache in {path}")
    options = case_options(spec, case, fidelity)
    options["seed"] = (100_000 + sequence*100 +
        {"out":0, "in":1, "moderator_minus":2, "moderator_plus":3,
         "total_minus":4, "total_plus":5, "fuel_minus":6, "fuel_plus":7}[case])
    model = b.build_openmc_model_tb35(str(directory), **options)
    if "target_sigma_k" in fidelity:
        settings = model["settings"]
        settings.keff_trigger = {"type":"std_dev",
                                 "threshold":fidelity["target_sigma_k"]}
        settings.trigger_active = True
        settings.trigger_max_batches = fidelity["max_batches"]
        settings.trigger_batch_interval = 20
        settings.export_to_xml(directory / "settings.xml")
    atomic_json(directory / "request.json", request)
    atomic_json(ROOT / "status.json", dict(status="running", candidate=name,
                                           stage=stage, case=case, index=sequence,
                                           updated_unix=time.time()))
    environment = os.environ.copy()
    environment["OPENMC_CROSS_SECTIONS"] = b.CROSS_SECTIONS_PATH
    environment["OMP_NUM_THREADS"] = str(fidelity["threads"])
    with (directory / "run.log").open("w") as log:
        result = subprocess.run([b.OPENMC_EXEC, "-s", str(fidelity["threads"])],
                                cwd=directory, env=environment,
                                stdout=log, stderr=subprocess.STDOUT,
                                check=False)
    if result.returncode:
        raise RuntimeError(f"OpenMC failed for {name}/{stage}/{case}; see {directory/'run.log'}")
    statepoint = max(directory.glob("statepoint.*.h5"),
                     key=lambda p:int(p.name.split(".")[1]))
    with openmc.StatePoint(statepoint) as sp:
        k = float(sp.keff.nominal_value)
        sigma = float(sp.keff.std_dev)
        out = dict(candidate=name, stage=stage, case=case,
                   request_sha256=signature, k_eff=k, sigma_k=sigma,
                   rho_pcm=1e5*(1-1/k), sigma_rho_pcm=1e5*sigma/k**2,
                   batches=int(sp.current_batch), particles=int(sp.n_particles),
                   statepoint=str(statepoint))
    atomic_json(path, out)
    print(name, stage, case, f"k={k:.6f} +/- {sigma:.6f}", flush=True)
    return out


def alpha(minus, plus, delta):
    value = (plus["rho_pcm"]-minus["rho_pcm"])/(2*delta)
    sigma = math.hypot(plus["sigma_rho_pcm"],minus["sigma_rho_pcm"])/(2*delta)
    return dict(alpha_pcm_K=value, sigma_pcm_K=sigma,
                upper_95_pcm_K=value+1.645*sigma,
                lower_95_pcm_K=value-1.645*sigma)


def evaluate(name, spec, stage, fidelity, sequence):
    results = {}
    out = run_case(name,spec,"out",stage,fidelity,sequence)
    results["out"] = out
    if out["k_eff"]+1.96*out["sigma_k"] <= 1:
        return dict(pass_all=False, reason="OUT clearly subcritical", results=results)
    inside = run_case(name,spec,"in",stage,fidelity,sequence)
    results["in"] = inside
    if inside["k_eff"]-1.96*inside["sigma_k"] >= 1:
        return dict(pass_all=False, reason="IN clearly supercritical", results=results)
    coefficients = {}
    for kind in ("moderator", "total", "fuel"):
        low = run_case(name,spec,kind+"_minus",stage,fidelity,sequence)
        high = run_case(name,spec,kind+"_plus",stage,fidelity,sequence)
        results[low["case"]] = low
        results[high["case"]] = high
        coeff = alpha(low, high, fidelity["delta_K"])
        coefficients[kind] = coeff
        if coeff["lower_95_pcm_K"] >= 0:
            return dict(pass_all=False, reason=f"{kind} clearly nonnegative",
                        results=results, coefficients=coefficients)
        if stage == "confirmation" and coeff["upper_95_pcm_K"] >= 0:
            return dict(pass_all=False,
                        reason=f"{kind} not negative at 95% confidence",
                        results=results, coefficients=coefficients)
    if stage == "pilot":
        return dict(pass_all=False, reason="pilot promising: proceed to confirmation",
                    results=results, coefficients=coefficients, promote=True)
    conditions = dict(out=out["k_eff"]-1.96*out["sigma_k"]>1,
                      inside=inside["k_eff"]+1.96*inside["sigma_k"]<1,
                      **{kind:coefficients[kind]["upper_95_pcm_K"]<0
                         for kind in ("fuel", "moderator", "total")})
    return dict(pass_all=all(conditions.values()),
                reason="all five criteria pass" if all(conditions.values()) else
                       "confirmation confidence/margin insufficient",
                conditions=conditions, results=results, coefficients=coefficients)


def evaluate_beta_max(name, spec, sequence):
    """Check the other end of the plausible YHx expansion range for finalists."""
    fidelity = FINAL | {"beta": BETA_MAX}
    results, coefficients = {}, {}
    for kind in ("moderator", "total"):
        low = run_case(name,spec,kind+"_minus","beta_max",fidelity,sequence)
        high = run_case(name,spec,kind+"_plus","beta_max",fidelity,sequence)
        results[low["case"]] = low
        results[high["case"]] = high
        coefficients[kind] = alpha(low,high,fidelity["delta_K"])
        if coefficients[kind]["upper_95_pcm_K"] >= 0:
            return dict(pass_all=False,reason=f"{kind} not negative at beta_max",
                        results=results,coefficients=coefficients)
    return dict(pass_all=True,reason="moderator and total negative at beta_max",
                results=results,coefficients=coefficients)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    atomic_json(ROOT / "candidate_catalog.json",
                [{"name":name,"parameters":BASE|updates}
                 for name,updates in CATALOG + TUNED])
    if args.dry_run:
        print("catalogued",len(CATALOG)+len(TUNED),"distinct configurations")
        for name,updates in CATALOG + TUNED:
            print(name,validate_geometry(BASE|updates) or "geometry checks pass")
        return
    solution = ROOT / "solution.json"
    if solution.exists():
        print("Qualified solution already recorded:",solution)
        return
    index = 0
    while True:
        name, spec = catalog_entry(index)
        reason = validate_geometry(spec)
        folder = ROOT / "campaign" / name
        folder.mkdir(parents=True,exist_ok=True)
        completed = folder / "assessment.json"
        if completed.exists():
            prior = json.loads(completed.read_text())
            if prior.get("spec") == spec and prior.get("status") == "rejected":
                index += 1
                continue
            if prior.get("status") == "qualified":
                atomic_json(solution, prior)
                return
        atomic_json(folder / "candidate.json",dict(index=index,name=name,parameters=spec,
                                                   geometry_check=reason or "passed"))
        if reason:
            assessment = dict(candidate=name,index=index,status="invalid_geometry",
                              reason=reason)
        else:
            pilot = evaluate(name,spec,"pilot",PILOT,index)
            if pilot.get("promote"):
                confirmation = evaluate(name,spec,"confirmation",FINAL,index)
                beta_max = (evaluate_beta_max(name,spec,index)
                            if confirmation["pass_all"] else None)
                assessment = dict(candidate=name,index=index,spec=spec,
                                  pilot=pilot,confirmation=confirmation,
                                  beta_max=beta_max,
                                  status="qualified" if beta_max and beta_max["pass_all"]
                                  else "rejected")
            else:
                assessment = dict(candidate=name,index=index,spec=spec,
                                  pilot=pilot,status="rejected")
        atomic_json(folder / "assessment.json",assessment)
        with (ROOT / "search_journal.jsonl").open("a") as journal:
            journal.write(json.dumps(dict(candidate=name,index=index,status=assessment["status"],
                                          reason=(assessment.get("beta_max") or
                                                  assessment.get("confirmation") or
                                                  assessment.get("pilot") or assessment).get("reason")))+"\n")
        summarize_progress()
        print(name,assessment["status"],flush=True)
        if assessment["status"] == "qualified":
            atomic_json(solution,assessment)
            atomic_json(ROOT / "status.json",dict(status="qualified",candidate=name,
                                                   index=index,completed_unix=time.time()))
            summarize_progress()
            return
        index += 1


if __name__ == "__main__":
    main()
