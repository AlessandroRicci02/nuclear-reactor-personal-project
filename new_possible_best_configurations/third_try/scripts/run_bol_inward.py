"""Reproduce the C3 BOL check with all absorber drums facing the core."""
from __future__ import annotations

import json
import os

import openmc

import build_model as b
from run_bol import ROOT, PARAMETERS


def main():
    run = ROOT / "campaign" / "c3_bol_inward"
    run.mkdir(parents=True, exist_ok=True)
    parameters = PARAMETERS | {"drum_angles": [180.0] * 10, "seed": 1011}
    model = b.build_openmc_model_tb35(str(run), **parameters)
    model["settings"].temperature = {"multipole": True, "method": "nearest",
                                     "tolerance": 1500.0}
    model["settings"].export_to_xml(run / "settings.xml")
    os.environ["OPENMC_CROSS_SECTIONS"] = b.CROSS_SECTIONS_PATH
    openmc.run(cwd=str(run), threads=8, openmc_exec=b.OPENMC_EXEC, output=False)
    with openmc.StatePoint(run / "statepoint.170.h5") as statepoint:
        keff = statepoint.keff
        result = {
            "case": "third_try_C3_BOL_all_drums_IN",
            "k_eff": float(keff.nominal_value),
            "k_eff_sigma": float(keff.std_dev),
            "confidence_95_normal": [float(keff.nominal_value - 1.96 * keff.std_dev),
                                     float(keff.nominal_value + 1.96 * keff.std_dev)],
            "parameters": parameters,
            "cross_sections": b.CROSS_SECTIONS_PATH,
        }
    (ROOT / "bol_inward_result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
