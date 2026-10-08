"""Criticality-controlled 400 kWth OpenMC depletion, sampled every 2.5 y.

The result is conditional on the fixed BOL temperature field. Each interval
uses the critical drum angle found at its beginning, with the full ENDF/B-VIII.0
depletion chain. The angle is changed again at the next 2.5-year state.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

import openmc
import openmc.deplete


ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[1] / "new_possible_best_configurations/third_try/scripts"))
import build_model as b  # noqa: E402

SPEC = json.loads((ROOT / "solution.json").read_text())["spec"]
CHAIN = ROOT.parents[1] / "new_possible_best_configurations_failed/campaign/data/chain_endfb80_thermal.xml"
POWER_W = 400_000.0
STEP_Y = 2.5
STEP_D = STEP_Y * 365.25
THREADS = 8
RUN = HERE / "depletion"
RUN.mkdir(parents=True, exist_ok=True)


def atomic(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def model_options(f_in, f_out, theta, seed, particles, batches, inactive):
    return SPEC | dict(f_in_mat=f_in, f_out_mat=f_out,
                       drum_angles=[theta] * 10, n_axial_slices=6,
                       seed=seed, particles=particles, batches=batches,
                       inactive=inactive)


def initial_fuel():
    mats = openmc.Materials.from_xml(
        ROOT / "campaign/matrix_short_strong_be/confirmation/out/materials.xml")
    return next(m for m in mats if m.id == 1), next(m for m in mats if m.id == 2)


def burned_fuel(burn_dir):
    marker = burn_dir / "burn_complete.json"
    if not marker.exists():
        return None
    data = json.loads(marker.read_text())
    result = openmc.deplete.Results(data["depletion_results"])
    mats = result.export_to_materials(-1, path=burn_dir / "burned_materials.xml")
    return next(m for m in mats if m.id == 1), next(m for m in mats if m.id == 2)


def run_k(step, label, f_in, f_out, theta, seed, target=.00034):
    directory = step / label
    directory.mkdir(parents=True, exist_ok=True)
    result_path = directory / "result.json"
    if result_path.exists():
        return json.loads(result_path.read_text())
    opts = model_options(f_in, f_out, theta, seed, 16000, 200, 50)
    built = b.build_openmc_model_tb35(str(directory), **opts)
    settings = built["settings"]
    settings.keff_trigger = {"type": "std_dev", "threshold": target}
    settings.trigger_active = True
    settings.trigger_max_batches = 1200
    settings.trigger_batch_interval = 20
    settings.export_to_xml(directory / "settings.xml")
    atomic(RUN / "status.json", dict(status="running", stage="transport",
                                    step=str(step.name), case=label,
                                    angle_deg=theta, started_unix=time.time()))
    env = os.environ.copy()
    env.update(OPENMC_CROSS_SECTIONS=b.CROSS_SECTIONS_PATH,
               OMP_NUM_THREADS=str(THREADS), OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1")
    with (directory / "run.log").open("w") as log:
        cp = subprocess.run([b.OPENMC_EXEC, "-s", str(THREADS)], cwd=directory,
                            env=env, stdout=log, stderr=subprocess.STDOUT)
    if cp.returncode:
        raise RuntimeError(f"OpenMC failed: {directory/'run.log'}")
    statepoint = max(directory.glob("statepoint.*.h5"),
                     key=lambda p: int(p.name.split(".")[1]))
    with openmc.StatePoint(statepoint) as sp:
        outcome = dict(angle_deg=theta, k_eff=float(sp.keff.n),
                       sigma_k=float(sp.keff.s), batches=int(sp.current_batch),
                       statepoint=str(statepoint))
    atomic(result_path, outcome)
    print(step.name, label, theta, outcome["k_eff"], "+/-",
          outcome["sigma_k"], flush=True)
    return outcome


def critical_angle(step, f_in, f_out, k0, k180, i):
    if k0["k_eff"] <= 1 or k180["k_eff"] >= 1:
        raise RuntimeError("Critical angle is not bracketed by OUT and IN")
    low = (0.0, k0["k_eff"])
    high = (180.0, k180["k_eff"])
    evaluations = []
    for j in range(12):
        guess = low[0] + (1-low[1])*(high[0]-low[0])/(high[1]-low[1])
        if not low[0]+0.01 < guess < high[0]-0.01:
            guess = (low[0]+high[0])/2
        guess = round(guess, 3)
        point = run_k(step, f"critical_{j:02d}_{guess:.3f}deg", f_in, f_out,
                      guess, seed=21000+i*100+j)
        evaluations.append(point)
        if abs(point["k_eff"]-1) <= .00050 and point["sigma_k"] <= .00034:
            return point, evaluations
        if point["k_eff"] > 1:
            low = (guess, point["k_eff"])
        else:
            high = (guess, point["k_eff"])
    best = min(evaluations, key=lambda p: abs(p["k_eff"]-1))
    if abs(best["k_eff"]-1) <= .00075:
        return best, evaluations
    raise RuntimeError(f"Critical-angle solve failed at {step}: {best}")


def nuclide_mass_kg(mats, nuclide):
    value = 0.0
    for m in mats:
        if nuclide in m.get_nuclides():
            value += m.get_mass(nuclide=nuclide)
    return value/1000


def burn(step, f_in, f_out, angle, i):
    directory = step / "burn_to_next"
    directory.mkdir(parents=True, exist_ok=True)
    recycled = burned_fuel(directory)
    if recycled:
        return recycled
    atomic(RUN / "status.json", dict(status="running", stage="depletion",
                                    step=str(step.name), angle_deg=angle,
                                    started_unix=time.time()))
    opts = model_options(f_in, f_out, angle, seed=31000+i,
                         particles=24000, batches=120, inactive=40)
    built = b.build_openmc_model_tb35(str(directory), **opts)
    model = openmc.Model(materials=built["materials"],
                         geometry=built["geometry"], settings=built["settings"])
    os.environ["OPENMC_CROSS_SECTIONS"] = b.CROSS_SECTIONS_PATH
    os.environ["OMP_NUM_THREADS"] = str(THREADS)
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    old_cwd = Path.cwd()
    os.chdir(directory)
    try:
        operator = openmc.deplete.CoupledOperator(
            model, str(CHAIN), normalization_mode="fission-q")
        integrator = openmc.deplete.CECMIntegrator(
            operator, [STEP_D], power=POWER_W, timestep_units="d")
        integrator.integrate()
    finally:
        os.chdir(old_cwd)
    path = directory / "depletion_results.h5"
    if not path.exists():
        raise RuntimeError(f"Depletion did not write {path}")
    atomic(directory / "burn_complete.json",
           dict(depletion_results=str(path), power_W=POWER_W,
                interval_years=STEP_Y, angle_deg=angle,
                chain=str(CHAIN), completed_unix=time.time()))
    return burned_fuel(directory)


def save_rows(rows):
    if not rows:
        return
    with (HERE / "depletion_2p5yr.csv").open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main(max_years):
    assert CHAIN.exists(), CHAIN
    f_in, f_out = initial_fuel()
    rows = []
    max_steps = int(max_years/STEP_Y)
    for i in range(max_steps+1):
        year = i*STEP_Y
        step = RUN / f"step_{i:02d}_{year:.1f}yr"
        step.mkdir(parents=True, exist_ok=True)
        if i:
            previous = RUN / f"step_{i-1:02d}_{(i-1)*STEP_Y:.1f}yr/burn_to_next"
            f_in, f_out = burned_fuel(previous) or (f_in, f_out)
        k0 = run_k(step, "out_0deg", f_in, f_out, 0.0,
                   seed=11000+i*100)
        k180 = run_k(step, "in_180deg", f_in, f_out, 180.0,
                     seed=11001+i*100)
        if k0["k_eff"] <= 1:
            status = "exhausted_or_uncertain"
            critical = None
        elif k180["k_eff"] >= 1:
            status = "cannot_shutdown_or_uncertain"
            critical = None
        else:
            critical, _ = critical_angle(step, f_in, f_out, k0, k180, i)
            status = "critical_at_sample"
        row = dict(years=year, status=status,
                   k_out=k0["k_eff"], sigma_k_out=k0["sigma_k"],
                   k_in=k180["k_eff"], sigma_k_in=k180["sigma_k"],
                   critical_angle_deg=critical["angle_deg"] if critical else "",
                   k_critical=critical["k_eff"] if critical else "",
                   sigma_k_critical=critical["sigma_k"] if critical else "",
                   u235_kg=nuclide_mass_kg((f_in,f_out), "U235"),
                   pu239_kg=nuclide_mass_kg((f_in,f_out), "Pu239"),
                   xe135_g=1000*nuclide_mass_kg((f_in,f_out), "Xe135"))
        rows.append(row)
        save_rows(rows)
        atomic(RUN / "status.json", dict(status=status, year=year,
                                         result=row, updated_unix=time.time()))
        if not critical:
            break
        if i < max_steps:
            f_in, f_out = burn(step, f_in, f_out, critical["angle_deg"], i)
    if len(rows)>1 and rows[-1]["k_out"]<=1<rows[-2]["k_out"]:
        a, z = rows[-2], rows[-1]
        crossing = a["years"] + (1-a["k_out"])*(z["years"]-a["years"])/(z["k_out"]-a["k_out"])
        lifetime = dict(status="bracketed", interval_years=[a["years"],z["years"]],
                        linear_interpolation_years=crossing,
                        note="Interpolation is conditional on fixed thermal model and power; not a continuous-control burnup solution.")
    else:
        lifetime = dict(status="not_bracketed", last_sample_years=rows[-1]["years"])
    atomic(HERE / "lifespan_result.json", lifetime)
    atomic(RUN / "status.json", dict(status="finished", lifetime=lifetime,
                                     completed_unix=time.time()))
    print(json.dumps(lifetime, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-years", type=float, default=50.0)
    args = parser.parse_args()
    try:
        main(args.max_years)
    except Exception as exc:
        atomic(RUN / "status.json", dict(status="failed", error=str(exc),
                                         failed_unix=time.time()))
        raise
