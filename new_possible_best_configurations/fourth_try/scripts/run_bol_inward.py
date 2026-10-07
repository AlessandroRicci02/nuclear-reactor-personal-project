"""BOL check with all absorber drums facing the core (drums IN, 180 deg) for fourth_try."""
from __future__ import annotations

import json
import os
from pathlib import Path

import openmc

import build_model as b
from run_bol import ROOT, PARAMETERS


def main():
    run = ROOT / "campaign" / "c4_bol_inward"
    run.mkdir(parents=True, exist_ok=True)
    parameters = PARAMETERS | {"drum_angles": [180.0] * 10, "seed": 1011}
    model = b.build_openmc_model_tb35(str(run), **parameters)
    model["settings"].temperature = {"multipole": False, "method": "interpolation",
                                     "tolerance": 10.0}
    model["settings"].export_to_xml(run / "settings.xml")
    os.environ["OPENMC_CROSS_SECTIONS"] = b.CROSS_SECTIONS_PATH
    os.environ["OMP_NUM_THREADS"] = "8"
    openmc.run(cwd=str(run), threads=8, openmc_exec=b.OPENMC_EXEC, output=False)
    statepoint_path = max(run.glob("statepoint.*.h5"), key=lambda p: int(p.name.split(".")[1]))
    with openmc.StatePoint(statepoint_path) as statepoint:
        keff = statepoint.keff
        k_val = float(keff.nominal_value)
        k_sig = float(keff.std_dev)
        rho_pcm = 1e5 * (k_val - 1.0) / k_val
        sdm_pcm = 1e5 * (1.0 - k_val) / k_val
        result = {
            "case": "fourth_try_C4_BOL_all_drums_IN",
            "k_eff": k_val,
            "k_eff_sigma": k_sig,
            "reactivity_pcm": rho_pcm,
            "shutdown_margin_pcm": sdm_pcm,
            "confidence_95_normal": [k_val - 1.96 * k_sig, k_val + 1.96 * k_sig],
            "parameters": parameters,
            "cross_sections": b.CROSS_SECTIONS_PATH,
        }
    (ROOT / "bol_inward_result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

