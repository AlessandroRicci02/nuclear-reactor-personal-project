"""Independent OpenMC stochastic-volume check of the C3 modeled mass."""
from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
import shutil

import openmc

import build_model as b
from run_bol import ROOT


SOURCE = ROOT / "campaign" / "c3_bol"
RUN = ROOT / "campaign" / "mass_volume"
SAMPLES = 5_000_000


def main():
    RUN.mkdir(parents=True, exist_ok=True)
    for name in ("geometry.xml", "materials.xml", "settings.xml"):
        shutil.copy2(SOURCE / name, RUN / name)
    settings = openmc.Settings.from_xml(RUN / "settings.xml")
    radius = 19.45 * (2.30 / 2.270) + 8.8
    half_height = 47.0/2 + 8.5
    lower = [-radius, -radius, -half_height]
    upper = [radius, radius, half_height]
    materials = openmc.Materials.from_xml(RUN / "materials.xml")
    settings.volume_calculations = [openmc.VolumeCalculation(list(materials),
                                     SAMPLES, lower, upper)]
    settings.export_to_xml(RUN / "settings.xml")
    os.environ["OPENMC_CROSS_SECTIONS"] = b.CROSS_SECTIONS_PATH
    if not (RUN / "volume_1.h5").exists():
        openmc.calculate_volumes(cwd=str(RUN), threads=2,
                                 openmc_exec=b.OPENMC_EXEC, output=False)
    volumes = openmc.VolumeCalculation.from_hdf5(RUN / "volume_1.h5")
    box = (2*radius)**2 * 2*half_height
    rows = []
    for material in materials:
        volume = volumes.volumes[material.id]
        density = material.get_mass_density()
        rows.append({"material": material.name,
                     "volume_cm3": float(volume.nominal_value),
                     "volume_sigma_cm3": float(volume.std_dev),
                     "density_g_cm3": float(density),
                     "mass_kg": float(volume.nominal_value*density/1000)})
    with (ROOT / "mass_breakdown.csv").open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    solids = [row for row in rows if not row["material"].startswith("He Coolant")]
    modeled = sum(row["mass_kg"] for row in solids)
    helium = sum(row["mass_kg"] for row in rows if row not in solids)
    mean_density = sum(row["volume_cm3"] / box * row["density_g_cm3"] for row in solids)
    mean_sq_density = sum(row["volume_cm3"] / box * row["density_g_cm3"]**2 for row in solids)
    sigma = box / 1000 * math.sqrt(max(0.0, mean_sq_density - mean_density**2) / SAMPLES)
    b.PITCH = 2.30
    b.H_CORE = 47.0
    b.REFL_RAD_TH = 8.8
    b.REFL_AX_TH = 8.5
    analytic = b.compute_mass_breakdown()["M_tot_kg"]
    result = {"method": "OpenMC stochastic material-volume calculation",
              "samples": SAMPLES, "modeled_dry_mass_kg": modeled,
              "sampling_sigma_kg": sigma,
              "analytic_dry_mass_kg": analytic,
              "stochastic_minus_analytic_kg": modeled-analytic,
              "helium_inventory_kg": helium,
              "hardware_allowance_kg_not_modeled": 65.0,
              "contingency_allowance_kg_not_modeled": 15.0,
              "screening_total_with_allowances_kg": modeled + 80.0,
              "excluded": "unmodeled shielding, power conversion, radiators, external piping and vehicle structure"}
    (ROOT / "mass_result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
