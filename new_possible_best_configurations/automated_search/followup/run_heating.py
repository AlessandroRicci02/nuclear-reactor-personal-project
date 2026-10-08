"""Neutron/photon heat deposition map for the qualified YH2 matrix core."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import openmc


ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[1] / "new_possible_best_configurations/third_try/scripts"))
import build_model as b  # noqa: E402

SPEC = json.loads((ROOT / "solution.json").read_text())["spec"]
RUN = HERE / "heating_transport"
N_SLICES = 6
N_PINS = 217
POWER_W = 400_000.0


def main(theta):
    RUN.mkdir(parents=True, exist_ok=True)
    result_path = HERE / "heating_result.json"
    if result_path.exists():
        old = json.loads(result_path.read_text())
        if abs(old["theta_deg"] - theta) > 1e-6:
            raise RuntimeError("Existing heating result has another drum angle")
        print(json.dumps(old,indent=2)); return
    built = b.build_openmc_model_tb35(str(RUN), **(SPEC | dict(
        n_axial_slices=N_SLICES, particles=12000, batches=180,
        inactive=50, seed=51001, drum_angles=[theta]*10)))
    built["settings"].photon_transport = True
    # OpenMC 0.16 does not support windowed multipole with photon transport.
    # The temperature-interpolated cross sections remain enabled.
    built["settings"].temperature["multipole"] = False
    built["settings"].sourcepoint = {"write": False}
    built["settings"].export_to_xml(RUN / "settings.xml")
    geometry = built["geometry"]
    cells = list(geometry.get_all_cells().values())
    fuel_cells = [built["fuel_cells_map"][(p,s)]
                  for p in range(N_PINS) for s in range(N_SLICES)]
    mod_cells = sorted((c for c in cells if isinstance(c.fill, openmc.Material)
                        and c.fill.id == 6), key=lambda c: c.id)
    if len(mod_cells) != N_PINS:
        raise RuntimeError(f"Expected {N_PINS} YH2 cells, found {len(mod_cells)}")
    mesh = openmc.RegularMesh()
    cr = built["core_r"]
    h = SPEC["core_h"]
    mesh.lower_left = (-cr, -cr, -h/2)
    mesh.upper_right = (cr, cr, h/2)
    mesh.dimension = (1, 1, N_SLICES)
    all_mat = openmc.Tally(name="all_material_heating")
    all_mat.filters = [openmc.MaterialFilter(list(built["materials"]))]
    all_mat.scores = ["heating"]
    fuel = openmc.Tally(name="fuel_pin_slice_heating")
    fuel.filters = [openmc.CellFilter(fuel_cells)]
    fuel.scores = ["heating"]
    mod = openmc.Tally(name="yh2_pin_slice_heating")
    mod.filters = [openmc.CellFilter(mod_cells), openmc.MeshFilter(mesh)]
    mod.scores = ["heating"]
    openmc.Tallies([all_mat, fuel, mod]).export_to_xml(RUN / "tallies.xml")
    (RUN / "manifest.json").write_text(json.dumps(dict(theta_deg=theta,
        power_W=POWER_W, particles=12000, batches=180, inactive=50,
        photon_transport=True, spec=SPEC),indent=2)+"\n")
    env=os.environ.copy()
    env.update(OPENMC_CROSS_SECTIONS=b.CROSS_SECTIONS_PATH,
               OMP_NUM_THREADS="2",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
    with (RUN / "run.log").open("w") as log:
        cp=subprocess.run([b.OPENMC_EXEC,"-s","2"],cwd=RUN,env=env,
                          stdout=log,stderr=subprocess.STDOUT)
    if cp.returncode:
        raise RuntimeError(f"OpenMC heating run failed; see {RUN/'run.log'}")
    spath=max(RUN.glob("statepoint.*.h5"),key=lambda p:int(p.name.split(".")[1]))
    with openmc.StatePoint(spath) as sp:
        mat_t=sp.get_tally(name="all_material_heating")
        fuel_t=sp.get_tally(name="fuel_pin_slice_heating")
        mod_t=sp.get_tally(name="yh2_pin_slice_heating")
        mat_raw=np.asarray(mat_t.mean).ravel()
        total_raw=float(mat_raw.sum())
        scale=POWER_W/total_raw
        mat_names=[m.name for m in built["materials"]]
        mat_power={name:float(val*scale) for name,val in zip(mat_names,mat_raw)}
        fuel_map=np.asarray(fuel_t.mean).reshape(N_PINS,N_SLICES)*scale
        mod_map=np.asarray(mod_t.mean).reshape(N_PINS,N_SLICES)*scale
        fuel_sd=np.asarray(fuel_t.std_dev).reshape(N_PINS,N_SLICES)*scale
        mod_sd=np.asarray(mod_t.std_dev).reshape(N_PINS,N_SLICES)*scale
        k=float(sp.keff.n);sk=float(sp.keff.s)
    np.savez(HERE/"heating_maps_W.npz",fuel_W=fuel_map,
             moderator_W=mod_map,sigma_fuel_W=fuel_sd,
             sigma_moderator_W=mod_sd)
    result=dict(theta_deg=theta,power_W=POWER_W,statepoint=str(spath),
                k_eff=k,sigma_k=sk,normalization_W_per_eV_source=scale,
                material_power_W=mat_power,fuel_power_W=float(fuel_map.sum()),
                moderator_power_W=float(mod_map.sum()),
                other_power_W=float(POWER_W-fuel_map.sum()-mod_map.sum()),
                moderator_relative_sigma_max=float(np.max(mod_sd/np.maximum(mod_map,1e-12))),
                fuel_relative_sigma_max=float(np.max(fuel_sd/np.maximum(fuel_map,1e-12))))
    result_path.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--theta",type=float,required=True)
    args=parser.parse_args()
    main(args.theta)
