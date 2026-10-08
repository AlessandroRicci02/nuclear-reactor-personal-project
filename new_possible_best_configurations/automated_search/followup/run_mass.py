"""Stochastic material-volume mass of the qualified OpenMC geometry."""
from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
import shutil

import openmc


ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
SOURCE = ROOT / "campaign/matrix_short_strong_be/confirmation/out"
RUN = HERE / "mass_volume"
SAMPLES = 5_000_000


def main():
    RUN.mkdir(parents=True, exist_ok=True)
    for name in ("geometry.xml", "materials.xml", "settings.xml"):
        shutil.copy2(SOURCE / name, RUN / name)
    settings = openmc.Settings.from_xml(RUN / "settings.xml")
    spec = json.loads((ROOT / "solution.json").read_text())["spec"]
    core_r = 19.45 * (spec["pitch"] / 2.270)
    radius = core_r + spec["refl_rad_th"]
    half_height = spec["core_h"]/2 + spec["refl_ax_th"]
    lower = [-radius, -radius, -half_height]
    upper = [radius, radius, half_height]
    materials = openmc.Materials.from_xml(RUN / "materials.xml")
    settings.volume_calculations = [openmc.VolumeCalculation(list(materials),
                                     SAMPLES, lower, upper)]
    settings.export_to_xml(RUN / "settings.xml")
    os.environ["OPENMC_CROSS_SECTIONS"] = str(
        Path("/home/aleri/OpenMC_projects/first_test/new_possible_best_configurations_failed/campaign/data/endfb-viii.0-hdf5/cross_sections.xml"))
    if not (RUN / "volume_1.h5").exists():
        openmc.calculate_volumes(cwd=str(RUN), threads=2,
                                 openmc_exec="/home/aleri/miniforge3/envs/openmc-env/bin/openmc",
                                 output=False)
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
    with (HERE / "mass_breakdown.csv").open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    solids = [row for row in rows if not row["material"].startswith("He Coolant")]
    modeled = sum(row["mass_kg"] for row in solids)
    helium = sum(row["mass_kg"] for row in rows if row not in solids)
    mean_density = sum(row["volume_cm3"] / box * row["density_g_cm3"] for row in solids)
    mean_sq_density = sum(row["volume_cm3"] / box * row["density_g_cm3"]**2 for row in solids)
    sigma = box / 1000 * math.sqrt(max(0.0, mean_sq_density - mean_density**2) / SAMPLES)
    result = {"method": "OpenMC stochastic material-volume calculation",
              "samples": SAMPLES, "modeled_dry_mass_kg": modeled,
              "sampling_sigma_kg": sigma, "helium_inventory_kg": helium,
              "core_plus_helium_kg": modeled+helium,
              "optional_unmodeled_hardware_allowance_kg": 65.0,
              "optional_contingency_allowance_kg": 15.0,
              "screening_mass_with_allowances_kg": modeled+80.0,
              "excluded": "unmodeled shielding, controls/drive hardware, power conversion, radiators, external piping and vehicle structure"}
    (HERE / "mass_result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
