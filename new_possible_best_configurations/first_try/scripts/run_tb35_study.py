#!/usr/bin/env python3
"""TB35-only OpenMC mass and BOL feedback campaign.

All calculations are independent, serial and restartable.  The results report
is rewritten atomically after every completed calculation so an interrupted
campaign never presents a partial value as final.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time

import numpy as np
import openmc

from tb35_model import build, CAMPAIGN, DEFAULT_XS, NOMINAL, OPENMC

ROOT = Path(__file__).resolve().parents[1]
RUNS = CAMPAIGN / 'neutronics'
MASS = CAMPAIGN / 'mass'
STATUS = CAMPAIGN / 'status.json'
RESULTS = ROOT / 'results.txt'
TARGET_SIGMA = 1.0e-4  # 10 pcm in k-effective at one standard deviation
MASS_SAMPLES = 5_000_000


def atomic_json(path: Path, value: dict | list) -> None:
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def read_json(path: Path):
    return json.loads(path.read_text()) if path.exists() else None


def update_status(stage: str, **more) -> None:
    value = read_json(STATUS) or {}
    value.update(stage=stage, updated_unix=time.time(), **more)
    atomic_json(STATUS, value)
    write_report()


def run_point(name: str, *, temperatures=None, seed=1):
    directory = RUNS / name
    directory.mkdir(parents=True, exist_ok=True)
    request = dict(
        configuration='TB35_REFERENCE', theta_deg=0.0,
        target_sigma=TARGET_SIGMA, particles=20_000, initial_batches=180,
        inactive_batches=60, seed=seed, temperatures=temperatures or NOMINAL,
        cross_sections=str(DEFAULT_XS),
    )
    result_path = directory / 'result.json'
    signature = hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest()
    if result_path.exists():
        result = read_json(result_path)
        if result['request_sha256'] != signature:
            raise RuntimeError(f'Existing result has a different request: {directory}')
        return result

    model, _ = build(directory, theta=0.0, temperatures=temperatures,
                     particles=20_000, batches=180, inactive=60, seed=seed,
                     target_sigma=TARGET_SIGMA)
    model.settings.trigger_max_batches = 20_000
    model.settings.export_to_xml(directory / 'settings.xml')
    env = os.environ.copy()
    env.update(OPENMC_CROSS_SECTIONS=str(DEFAULT_XS), OMP_NUM_THREADS='4',
               OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    update_status(name, current_run=name, current_directory=str(directory))
    started = time.time()
    with (directory / 'run.log').open('w') as log:
        completed = subprocess.run([OPENMC, '-s', '4'], cwd=directory, env=env,
                                   stdout=log, stderr=subprocess.STDOUT)
    if completed.returncode:
        raise RuntimeError(f'OpenMC failed for {name}; see {directory}/run.log')
    statepoint = max(directory.glob('statepoint.*.h5'),
                     key=lambda path: int(path.name.split('.')[1]))
    with openmc.StatePoint(statepoint) as point:
        keff, sigma = point.keff.n, point.keff.s
        entropy = np.asarray(point.entropy)
        result = dict(
            name=name, keff=keff, sigma=sigma, sigma_pcm_k=1e5*sigma,
            reactivity_pcm=1e5*(1.0 - 1.0/keff),
            reactivity_sigma_pcm=1e5*sigma/keff**2,
            batches=int(point.current_batch), requested_batches=int(point.n_batches),
            inactive=int(point.n_inactive), particles=int(point.n_particles),
            wall_seconds=time.time()-started, statepoint=str(statepoint),
            request=request, request_sha256=signature,
            target_achieved=bool(sigma <= TARGET_SIGMA),
            entropy_inactive_last20_mean=float(np.mean(entropy[max(0, point.n_inactive-20):point.n_inactive])),
            entropy_active_first20_mean=float(np.mean(entropy[point.n_inactive:point.n_inactive+20])),
        )
    if not result['target_achieved']:
        raise RuntimeError(f'{name} missed its 10 pcm target: {sigma:.6g}')
    atomic_json(result_path, result)
    return result


def coefficient(name: str, minus: dict, plus: dict, span_K: float) -> dict:
    value = (plus['reactivity_pcm'] - minus['reactivity_pcm']) / span_K
    sigma = math.hypot(plus['reactivity_sigma_pcm'], minus['reactivity_sigma_pcm']) / span_K
    return dict(
        name=name, value_pcm_per_K=value, sigma_pcm_per_K=sigma,
        span_K=span_K,
        positive_at_95pct=bool(value - 1.96*sigma > 0),
        negative_at_95pct=bool(value + 1.96*sigma < 0),
        definition='Central finite difference at all-drums-out, fixed geometry and fixed nonperturbed material densities.',
        minus_statepoint=minus['statepoint'], plus_statepoint=plus['statepoint'],
    )


def run_mass() -> dict:
    result_path = MASS / 'result.json'
    if result_path.exists():
        return read_json(result_path)
    MASS.mkdir(parents=True, exist_ok=True)
    model, info = build(MASS, theta=0.0, particles=1000, batches=10, inactive=5)
    cfg = info['config']
    radius = cfg['metrics']['outer_radius_cm']
    height = cfg['metrics']['total_height_cm']
    model.settings.run_mode = 'volume'
    model.settings.volume_calculations = [openmc.VolumeCalculation(
        list(model.materials), MASS_SAMPLES, [-radius, -radius, -height/2],
        [radius, radius, height/2])]
    model.settings.export_to_xml(MASS / 'settings.xml')
    env = os.environ.copy()
    env.update(OPENMC_CROSS_SECTIONS=str(DEFAULT_XS), OMP_NUM_THREADS='4',
               OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    update_status('mass', current_run='mass', current_directory=str(MASS))
    with (MASS / 'run.log').open('w') as log:
        completed = subprocess.run([OPENMC, '-s', '4'], cwd=MASS, env=env,
                                   stdout=log, stderr=subprocess.STDOUT)
    if completed.returncode:
        raise RuntimeError(f'OpenMC volume calculation failed; see {MASS}/run.log')
    calculation = openmc.VolumeCalculation.from_hdf5(MASS / 'volume_1.h5')
    grouped = {}
    for material in model.materials:
        group = 'UN fuel' if material.name.startswith('UN_zone') else material.name
        entry = grouped.setdefault(group, dict(volume_cm3=0.0, density_g_cm3=material.get_mass_density()))
        entry['volume_cm3'] += calculation.volumes[material.id].n
    domain = (2*radius)**2 * height
    rows = []
    for name, entry in grouped.items():
        probability = entry['volume_cm3'] / domain
        mass = entry['volume_cm3'] * entry['density_g_cm3'] / 1000
        sigma = domain * entry['density_g_cm3'] / 1000 * math.sqrt(probability*(1-probability)/MASS_SAMPLES)
        rows.append(dict(material=name, volume_cm3=entry['volume_cm3'], density_g_cm3=entry['density_g_cm3'], mass_kg=mass, sigma_kg=sigma))
    solid_rows = [row for row in rows if row['material'] != 'Helium']
    mass = sum(row['mass_kg'] for row in solid_rows)
    rho_bar = sum(row['volume_cm3']/domain*row['density_g_cm3'] for row in solid_rows)
    rho2_bar = sum(row['volume_cm3']/domain*row['density_g_cm3']**2 for row in solid_rows)
    sigma = domain/1000 * math.sqrt((rho2_bar-rho_bar**2)/MASS_SAMPLES)
    with (MASS / 'breakdown.csv').open('w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=rows[0])
        writer.writeheader(); writer.writerows(rows)
    result = dict(
        method='OpenMC stochastic volume calculation', samples=MASS_SAMPLES,
        modeled_solid_mass_kg=mass, sigma_kg=sigma,
        analytic_geometry_mass_kg=cfg['metrics']['modeled_solids_kg'],
        analytic_minus_stochastic_sigma=(cfg['metrics']['modeled_solids_kg']-mass)/sigma,
        helium_mass_kg=grouped['Helium']['volume_cm3']*grouped['Helium']['density_g_cm3']/1000,
        excludes='Pressure vessel, drum drives/shrouds, end manifolds, supports, additional hydrogen encapsulation, shielding, conversion system.',
    )
    atomic_json(result_path, result)
    return result


def write_report() -> None:
    mass = read_json(MASS / 'result.json')
    bol = read_json(RUNS / 'bol_0deg' / 'result.json')
    feedback = read_json(CAMPAIGN / 'feedback_coefficients.json') or []
    state = read_json(STATUS) or {}
    lines = [
        'TB35_REFERENCE — OpenMC BOL mass and temperature-feedback campaign',
        'Configuration: configuration/TB35_REFERENCE (copied into this workspace).',
        'Nuclear data: official ENDF/B-VIII.0 HDF5, temperature interpolation, c_H_in_YH2 and c_Y_in_YH2 thermal scattering.',
        f"Campaign status: {state.get('stage', 'not started')}",
        'Reported keff uncertainties are one-standard-deviation Monte Carlo sampling errors only; they exclude nuclear-data, geometry and material-property uncertainty.',
        '', 'MASS',
    ]
    if mass:
        lines.extend([
            f"Modeled solid mass: {mass['modeled_solid_mass_kg']:.3f} +/- {mass['sigma_kg']:.3f} kg (OpenMC {mass['samples']:,}-sample volume calculation).",
            f"Analytic geometry mass: {mass['analytic_geometry_mass_kg']:.3f} kg; difference: {mass['analytic_minus_stochastic_sigma']:.2f} sigma.",
            f"Helium inventory: {mass['helium_mass_kg']:.5f} kg.",
            'Excludes: '+mass['excludes'],
        ])
    else: lines.append('PENDING')
    lines.extend(['', 'BOL KEFF, ALL DRUMS OUT (0 degrees)'])
    if bol:
        lines.append(f"keff = {bol['keff']:.7f} +/- {bol['sigma']:.7f} ({bol['sigma_pcm_k']:.2f} pcm in k); rho = {bol['reactivity_pcm']:+.1f} +/- {bol['reactivity_sigma_pcm']:.1f} pcm.")
        lines.append(f"Histories: {bol['particles']:,} particles/generation; {bol['inactive']} inactive; {bol['batches']} total batches. Target 10 pcm: achieved.")
    else: lines.append('PENDING')
    lines.extend(['', 'THERMAL COEFFICIENTS'])
    if feedback:
        for item in feedback:
            sign = 'positive at 95%' if item['positive_at_95pct'] else 'negative at 95%' if item['negative_at_95pct'] else 'sign unresolved at 95%'
            lines.append(f"{item['name']}: {item['value_pcm_per_K']:+.5f} +/- {item['sigma_pcm_per_K']:.5f} pcm/K; {sign}.")
        lines.extend([
            '',
            'Definitions: fuel is a +/-200 K Doppler finite difference; moderator is a +/-100 K finite difference at fixed moderator density; total is a +/-100 K isothermal finite difference at fixed geometry, fixed material densities and fixed helium pressure.',
            'The total coefficient is therefore a transport-temperature coefficient, not a mechanically coupled thermal-expansion coefficient. It must not be used to claim the sign of physical moderator expansion feedback without a temperature-dependent YH2 density/geometry model.',
        ])
    else: lines.append('PENDING')
    temporary = RESULTS.with_suffix('.txt.tmp')
    temporary.write_text('\n'.join(lines)+'\n')
    temporary.replace(RESULTS)


def main() -> None:
    update_status('starting', pid=os.getpid())
    run_mass(); update_status('mass complete')
    run_point('bol_0deg', seed=41001); update_status('BOL keff complete')
    cases = [
        ('fuel', 200.0, 42001),
        ('moderator', 100.0, 43001),
    ]
    output = []
    for material, delta, seed in cases:
        minus = run_point(f'{material}_minus', temperatures={material: NOMINAL[material]-delta}, seed=seed)
        update_status(f'{material} minus complete')
        plus = run_point(f'{material}_plus', temperatures={material: NOMINAL[material]+delta}, seed=seed+1)
        output.append(coefficient(f'{material}_temperature_fixed_density', minus, plus, 2*delta))
        atomic_json(CAMPAIGN / 'feedback_coefficients.json', output); write_report()
        update_status(f'{material} coefficient complete')
    iso_minus = run_point('total_isothermal_minus', temperatures={key:value-100. for key,value in NOMINAL.items()}, seed=44001)
    update_status('total isothermal minus complete')
    iso_plus = run_point('total_isothermal_plus', temperatures={key:value+100. for key,value in NOMINAL.items()}, seed=44002)
    output.append(coefficient('total_isothermal_fixed_geometry', iso_minus, iso_plus, 200.0))
    atomic_json(CAMPAIGN / 'feedback_coefficients.json', output)
    update_status('complete', completed_unix=time.time(), current_run=None, current_directory=None)


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        update_status('failed', error=str(error))
        raise
