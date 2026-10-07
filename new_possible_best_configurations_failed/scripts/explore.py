#!/usr/bin/env python3
"""Deterministic geometry/inventory screening, NOT reactor-performance simulation.

Run with /home/aleri/miniforge3/envs/openmc-env/bin/python scripts/explore.py
All lengths in cm unless a field explicitly states otherwise. No OpenMC/MOOSE run.
"""
from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.spatial import ConvexHull, cKDTree

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parent
PITCH = 2.280
H0 = 45.50
RF0 = 0.585
RCORE0 = 19.45 * PITCH / 2.270
COUNTS = np.array([19, 72, 126])
ENR = np.array([0.196, 0.1975, 0.1975])
CLAD = 0.048
GAP0 = 0.084
LINER = 0.010
R_AUX0 = 0.150
A_AUX0 = 6 * math.pi * R_AUX0**2
A_ANN0 = math.pi * ((RF0 + CLAD + GAP0)**2 - (RF0 + CLAD)**2)
A_HEX = math.sqrt(3) / 2 * PITCH**2
R_DUCT = 0.860
R_DRUM = 3.95
T_ABS = 1.50
A_ABS = (148 / 360) * math.pi * (R_DRUM**2 - (R_DRUM - T_ABS)**2)
DRUM_CLEARANCE = 0.180
RHO = {'fuel': 13.50, 'clad': 8.57, 'liner': 3.95, 'moderator': 4.30,
       'reflector': 2.90, 'absorber': 2.52}
SHAPES = {
    # Actual OpenMC orientation=y: faces normal at 30+60k degrees.
    'C6_VERTEX': {'count': 6, 'ratio': 1.0, 'phase_deg': 0.0},
    'C6_FACE': {'count': 6, 'ratio': 1.0, 'phase_deg': 30.0},
    'E6_VERTEX_150': {'count': 6, 'ratio': 1.5, 'phase_deg': 0.0},
    'E6_VERTEX_178': {'count': 6, 'ratio': 16 / 9, 'phase_deg': 0.0},
    'E6_FACE_150': {'count': 6, 'ratio': 1.5, 'phase_deg': 30.0},
    'E6_FACE_178': {'count': 6, 'ratio': 16 / 9, 'phase_deg': 30.0},
    'C12_STAGGERED': {'count': 12, 'ratio': 1.0, 'phase_deg': 15.0},
}


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def write_csv(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(records[0]))
        w.writeheader()
        w.writerows(records)


def pin_centres():
    """Coordinates from OpenMC HexLattice orientation=y, not old plot metadata."""
    out = []
    for a in range(-8, 9):
        for b in range(-8, 9):
            ring = max(abs(a), abs(b), abs(a + b))
            if ring <= 8:
                out.append(dict(ring=ring, lattice_a=a, lattice_b=b,
                                x_cm=math.sqrt(3)/2*PITCH*a,
                                y_cm=PITCH*(b + 0.5*a),
                                zone=0 if ring <= 2 else (1 if ring <= 5 else 2)))
    return sorted(out, key=lambda x: (x['ring'], math.atan2(x['y_cm'], x['x_cm'])))


PINS = pin_centres()
ANGLES = np.arange(6)*math.pi/3
HEX_VERTICES = (PITCH/math.sqrt(3))*np.column_stack((np.cos(ANGLES), np.sin(ANGLES)))
HEX_NORMALS = np.column_stack((np.cos(ANGLES+math.pi/6), np.sin(ANGLES+math.pi/6)))
ALL_VERTICES = np.vstack([HEX_VERTICES + [p['x_cm'], p['y_cm']] for p in PINS])
HULL = ConvexHull(ALL_VERTICES)


@lru_cache(None)
def moderator_envelope_area(core_r, offset):
    """Intersection of the core circle and outward-offset occupied-cell hull.

    offset=-1 means the original full circular YH2 envelope. For offset>=0,
    the complement inside the same core cylinder is solid BeO, never void.
    Areas by deterministic polar quadrature; error is integration error only.
    """
    if offset < 0:
        return math.pi*core_r**2, 0.0
    normals = HULL.equations[:, :2]
    distances = -HULL.equations[:, 2] + offset

    def f(theta):
        dots = normals @ np.array([math.cos(theta), math.sin(theta)])
        mask = dots > 1e-14
        radius = min(core_r, float(np.min(distances[mask]/dots[mask])))
        return radius**2 / 2

    # Sixfold symmetry, split for robust quadrature across envelope transitions.
    pieces = [quad(f, a, b, epsabs=1e-9, epsrel=1e-10, limit=150)
              for a, b in zip(np.linspace(0, math.pi/3, 49)[:-1],
                              np.linspace(0, math.pi/3, 49)[1:])]
    return 6*sum(p[0] for p in pieces), 6*sum(p[1] for p in pieces)


def channel_parameters(shape):
    s = SHAPES[shape]
    # rad=tang for circles; for ellipses tang/rad=ratio. Exact fixed area.
    rad = math.sqrt(A_AUX0/(s['count']*math.pi*s['ratio']))
    tang = rad*s['ratio']
    return s['count'], rad, tang, s['phase_deg']


@lru_cache(None)
def channel_geometry(shape, centre):
    n, a, b, phase = channel_parameters(shape)
    angles = np.deg2rad(phase) + np.arange(n)*2*math.pi/n
    centres = centre*np.column_stack((np.cos(angles), np.sin(angles)))
    radial = np.column_stack((np.cos(angles), np.sin(angles)))
    tangent = np.column_stack((-np.sin(angles), np.cos(angles)))
    projections = centres @ HEX_NORMALS.T
    supports = np.sqrt((a*radial @ HEX_NORMALS.T)**2 + (b*tangent @ HEX_NORMALS.T)**2)
    to_face = float(np.min(PITCH/2 - projections - supports))
    perimeter = quad(lambda t: math.sqrt(a*a*math.sin(t)**2+b*b*math.cos(t)**2),
                     0, 2*math.pi, epsabs=1e-12)[0]
    pair_bound = min(np.linalg.norm(centres[i]-centres[j])-2*max(a,b)
                     for i in range(n) for j in range(i))
    # Union of axial central duct and auxiliary holes: do not double count overlap.
    lo, hi = max(-R_DUCT, centre-a), min(R_DUCT, centre+a)
    overlap = 0.0
    if hi > lo:
        def strip(x):
            ey = b*math.sqrt(max(0, 1-((x-centre)/a)**2))
            cy = math.sqrt(max(0, R_DUCT**2-x*x))
            return 2*min(ey,cy)
        overlap = quad(strip, lo, hi, epsabs=1e-11, limit=150)[0]
    axial_union = math.pi*R_DUCT**2 + A_AUX0 - n*overlap
    # Hagen-Poiseuille exact ratio for equal-area ellipse bundle; same mu and L.
    resistance_ratio = (n/6)*(a*a+b*b)/(2*a*b)
    return dict(count=n, a_cm=a, b_cm=b, phase_deg=phase,
                centre_radius_cm=centre, perimeter_cm=n*perimeter,
                perimeter_ratio=n*perimeter/(12*math.pi*R_AUX0),
                min_face_ligament_cm=to_face,
                min_aux_pair_ligament_bound_cm=float(pair_bound),
                axial_union_area_cm2=axial_union,
                axial_overlap_area_cm2=n*overlap,
                laminar_aux_resistance_ratio=resistance_ratio)


def zone_geometry(height, middle_area_factor, allocation):
    """Keep inner19 and remaining198 fuel volumes separately invariant.

    Consequently UN, HM, and U235 inventories all remain invariant at fixed
    enrichment. Annular area is conserved over all 217 pins, not per zone.
    """
    factors = np.array([1.0, middle_area_factor,
                        (198-72*middle_area_factor)/126])
    radii = RF0*np.sqrt(H0/height*factors)
    if allocation == 'graded':
        weights = np.array([1.08, 1.04, (217-19*1.08-72*1.04)/126])
    else:
        weights = np.ones(3)
    r_clad = radii+CLAD
    ann_areas = A_ANN0*weights
    r_ann = np.sqrt(r_clad**2+ann_areas/math.pi)
    return radii, r_clad, r_ann, r_ann+LINER, ann_areas, weights


@lru_cache(None)
def fuel_and_cell(height, middle_area_factor, allocation):
    rf, rc, ra, rm, areas, weights = zone_geometry(height, middle_area_factor, allocation)
    fuel = math.pi*rf**2*height
    clad = math.pi*(rc**2-rf**2)*height
    liner = math.pi*(rm**2-ra**2)*height
    cellmod = (A_HEX-math.pi*rm**2-A_AUX0)*height
    return rf, rc, ra, rm, areas, weights, fuel, clad, liner, cellmod


def evaluate(c):
    h, r = c['height_cm'], c['core_radius_cm']
    ch = channel_geometry(c['shape'], c['channel_centre_cm'])
    rf, rc, ra, rm, ann, weights, vf, vc, vl, vm = fuel_and_cell(
        h, c['middle_fuel_area_factor'], c['annulus_allocation'])
    a_env, err = moderator_envelope_area(r, c['peripheral_hull_offset_cm'])
    a_periph = a_env-217*A_HEX
    v_abs = 10*A_ABS*h
    r_out = r+c['radial_reflector_cm']
    v_beo_active = (math.pi*r_out**2-a_env)*h-v_abs
    v_beo_ax = 2*c['axial_reflector_cm']*(math.pi*r_out**2-217*ch['axial_union_area_cm2'])
    masses = dict(fuel_kg=float(COUNTS@vf)*RHO['fuel']/1000,
                  cladding_kg=float(COUNTS@vc)*RHO['clad']/1000,
                  liner_kg=float(COUNTS@vl)*RHO['liner']/1000,
                  moderator_cells_kg=float(COUNTS@vm)*RHO['moderator']/1000,
                  moderator_periphery_kg=a_periph*h*RHO['moderator']/1000,
                  beo_kg=(v_beo_active+v_beo_ax)*RHO['reflector']/1000,
                  b4c_kg=v_abs*RHO['absorber']/1000)
    mass = sum(masses.values())
    hm = masses['fuel_kg']*238/252  # Match reference's stated weight fractions.
    fissile = float(COUNTS@(vf*ENR))*RHO['fuel']/1000*238/252
    inner_gap = c['channel_centre_cm']-ch['a_cm']-float(max(rm))
    minlig = min(inner_gap, ch['min_face_ligament_cm'],
                 ch['min_aux_pair_ligament_bound_cm'])
    outer_web = c['radial_reflector_cm']-DRUM_CLEARANCE-2*R_DRUM
    geometry_ok = (minlig >= 0.05-1e-10 and outer_web >= 0.7-1e-10
                   and r >= np.linalg.norm(ALL_VERTICES, axis=1).max()+0.03)
    reasons = []
    if minlig < 0.05-1e-10: reasons.append('ligament_below_study_0.50mm')
    if outer_web < 0.7-1e-10: reasons.append('outer_drum_web_below_study_7mm')
    if not r >= np.linalg.norm(ALL_VERTICES, axis=1).max()+0.03:
        reasons.append('occupied_cells_too_close_to_core_boundary')
    if mass > 600: reasons.append('modeled_solids_above_600kg')
    if mass > 575: reasons.append('less_than_study_25kg_hardware_allowance')
    return dict(**c, **masses, modeled_solids_kg=mass, hm_kg=hm, u235_kg=fissile,
                moderator_total_kg=masses['moderator_cells_kg']+masses['moderator_periphery_kg'],
                nominal_unallocated_mass_kg=600-mass,
                geometry_screen_pass=bool(geometry_ok),
                study_575kg_screen_pass=bool(geometry_ok and mass <=575),
                reasons=';'.join(reasons),
                min_ligament_mm=minlig*10, inner_ligament_mm=inner_gap*10,
                cell_face_ligament_mm=ch['min_face_ligament_cm']*10,
                outer_drum_web_mm=outer_web*10,
                auxiliary_perimeter_ratio=ch['perimeter_ratio'],
                auxiliary_resistance_ratio_at_same_mu=ch['laminar_aux_resistance_ratio']*h/H0,
                auxiliary_area_cm2_per_cell=A_AUX0,
                annular_area_cm2_all_pins=float(COUNTS@ann),
                fuel_radius_inner_mm=rf[0]*10, fuel_radius_middle_mm=rf[1]*10,
                fuel_radius_outer_mm=rf[2]*10,
                outer_radius_cm=r_out, total_height_cm=h+2*c['axial_reflector_cm'],
                moderator_envelope_area_cm2=a_env, area_quadrature_error_cm2=err,
                axial_coolant_overlap_union_correction_cm2_per_pin=ch['axial_overlap_area_cm2'])


def config(name, **kwargs):
    d = dict(id=name, shape='C6_VERTEX', channel_centre_cm=.910*PITCH/2.270,
             height_cm=H0, middle_fuel_area_factor=1.0, annulus_allocation='uniform',
             core_radius_cm=RCORE0, peripheral_hull_offset_cm=-1.0,
             radial_reflector_cm=9.30, axial_reflector_cm=8.80)
    d.update(kwargs)
    return d


def explicit_spec(c):
    e = evaluate(c)
    ch = channel_geometry(c['shape'], c['channel_centre_cm'])
    rf, rc, ra, rm, areas, weights = zone_geometry(c['height_cm'],
        c['middle_fuel_area_factor'], c['annulus_allocation'])
    zones = []
    for i in range(3):
        zones.append(dict(zone=i, rings=[[0,2],[3,5],[6,8]][i], pins=int(COUNTS[i]),
                          uranium_u235_weight_fraction=ENR[i], fuel_radius_cm=rf[i],
                          clad_outer_radius_cm=rc[i], annulus_outer_radius_cm=ra[i],
                          liner_outer_radius_cm=rm[i], annulus_gap_mm=10*(ra[i]-rc[i]),
                          annulus_area_cm2_per_pin=areas[i],
                          provisional_main_flow_per_pin_relative_to_mean=weights[i],
                          auxiliary_flow_per_pin_relative_to_mean=1.0))
    return dict(schema_version=1, status='PROPOSED_GEOMETRY_NOT_TRANSPORT_OR_THERMALLY_VALIDATED',
                parameters=c, metrics=e, zones=zones, auxiliary_channels=ch,
                geometry=dict(length_unit='cm', pin_count=217, pitch_cm=PITCH,
                    lattice_orientation='OpenMC HexLattice orientation=y',
                    unit_cell_face_normal_degrees=[30+60*k for k in range(6)],
                    axial_main_duct_radius_cm=R_DUCT,
                    axial_coolant_rule='Boolean union of central duct and auxiliary holes; overlap counted once',
                    peripheral_rule='full circle if offset<0; otherwise intersection of circle and offset occupied-cell convex hull',
                    peripheral_complement_material='BeO',
                    hull_halfspaces_nx_ny_d=[[*x[:2],-x[2]] for x in HULL.equations],
                    central_liner_thickness_cm=LINER, fuel_clad_thickness_cm=CLAD,
                    drum_count=10, drum_radius_cm=R_DRUM, absorber_thickness_cm=T_ABS,
                    absorber_sector_degrees=148., b10_atom_fraction=.9,
                    drum_centres_radius_cm=c['core_radius_cm']+DRUM_CLEARANCE+R_DRUM,
                    drum_azimuth_degrees=[36*k for k in range(10)]),
                materials=dict(densities_g_cm3=RHO, fuel='UN with reference N14 weight prescription',
                    moderator='YH2.0', cladding='Nb-1Zr', liner='Al2O3',
                    reflector='BeO', absorber='B4C', coolant='He4'),
                boundaries=dict(power_W=400000, pressure_Pa=2500000,
                    inlet_K=892, target_mixed_outlet_K=1110,
                    cp_J_kgK=5193.16, total_mass_flow_kg_s=400000/(5193.16*218),
                    interpretation='Common screening boundary conditions, not an established Brayton operating point'),
                flow=dict(topology='Provisional parallel fuel-annulus and moderator-auxiliary branches',
                    auxiliary_fraction=.175,
                    pressure_balance='Common inlet/outlet pressure; orifice coefficients and plenum routing unresolved',
                    note='Flow weights are design targets, not calculated natural flow. Reassess after pin heating tallies.'),
                unresolved=['neutron multiplication and all feedbacks', 'pin and moderator heat deposition',
                    'full hydrogen encapsulation including auxiliary walls and end faces',
                    'pressure vessel, supports, drum shrouds/drives and end manifolds',
                    'common-pressure flow solution and coupled temperature fields',
                    'validated depletion and feasible Brayton cycle'],
                performance=dict(keff=None, lifetime_years=None, peak_moderator_K=None,
                    core_pressure_drop_kPa=None, net_electrical_kW=None))


def coolant_distance(c, zone, grid_n=141, boundary_n=720):
    """Sampled geometry-only distance to cold auxiliary surfaces, not temperature.

    Includes per-cell moderator only. Does not assess peripheral moderator.
    Finite grid and boundary sampling: report as approximate, not upper bound.
    """
    _, _, _, rm, _, _ = zone_geometry(c['height_cm'],c['middle_fuel_area_factor'],c['annulus_allocation'])
    n,a,b,phase = channel_parameters(c['shape'])
    t=np.arange(boundary_n)*2*math.pi/boundary_n
    clouds=[]
    pts=np.column_stack([x.ravel() for x in np.meshgrid(np.linspace(-PITCH/math.sqrt(3),PITCH/math.sqrt(3),grid_n),
                                                        np.linspace(-PITCH/2,PITCH/2,grid_n))])
    mask=(pts@HEX_NORMALS.T <= PITCH/2+1e-12).all(axis=1)
    mask &= (np.linalg.norm(pts,axis=1)>=rm[zone])
    for phi in np.deg2rad(phase)+np.arange(n)*2*math.pi/n:
        u=np.array([math.cos(phi),math.sin(phi)]); v=np.array([-math.sin(phi),math.cos(phi)])
        ctr=c['channel_centre_cm']*u
        clouds.append(ctr+a*np.cos(t)[:,None]*u+b*np.sin(t)[:,None]*v)
        rel=pts-ctr
        mask &= ((rel@u/a)**2+(rel@v/b)**2>=1)
    d=cKDTree(np.vstack(clouds)).query(pts[mask])[0]*10
    return dict(zone=zone, sample_count=len(d), max_distance_to_aux_mm=float(max(d)),
                p95_distance_to_aux_mm=float(np.percentile(d,95)),
                area_grid_n=grid_n, boundary_samples_per_channel=boundary_n)


def generate_scan():
    records=[]
    grid=itertools.product(SHAPES, [.9140088105726871,.93,.95], [44.5,45.5,46.5,48.0],
        [.96,1.0,1.04], ['uniform','graded'], [RCORE0-.10,RCORE0], [-1.,1.5,1.,.5],
        [8.8,9.0,9.3], [8.3,8.8,9.3])
    for i,(shape,centre,h,f,alloc,r,off,rad,ax) in enumerate(grid):
        c=config(f'SCAN_{i:05d}',shape=shape,channel_centre_cm=centre,height_cm=h,
            middle_fuel_area_factor=f,annulus_allocation=alloc,core_radius_cm=r,
            peripheral_hull_offset_cm=off,radial_reflector_cm=rad,axial_reflector_cm=ax)
        records.append(evaluate(c))
    write_csv(ROOT/'data/design_space.csv',records)
    summary=dict(designs_evaluated=len(records),
                 geometry_pass=sum(r['geometry_screen_pass'] for r in records),
                 modeled_mass_le600_and_geometry_pass=sum(r['geometry_screen_pass'] and r['modeled_solids_kg']<=600 for r in records),
                 study575_and_geometry_pass=sum(r['study_575kg_screen_pass'] for r in records),
                 min_modeled_mass_kg=min(r['modeled_solids_kg'] for r in records),
                 max_modeled_mass_kg=max(r['modeled_solids_kg'] for r in records),
                 fuel_inventory_spread_kg=max(r['fuel_kg'] for r in records)-min(r['fuel_kg'] for r in records),
                 u235_inventory_spread_kg=max(r['u235_kg'] for r in records)-min(r['u235_kg'] for r in records),
                 ranking='No neutronic/thermal ranking or weighted pseudo-performance score applied')
    write_json(ROOT/'data/screening_summary.json',summary)
    return records,summary


def main():
    records,summary=generate_scan()
    build_deliverables(records, summary)
    print(json.dumps(summary,indent=2))


def selected_designs():
    common=dict(shape='E6_VERTEX_150', channel_centre_cm=.93,
                middle_fuel_area_factor=.96, annulus_allocation='graded',
                radial_reflector_cm=9.0, axial_reflector_cm=8.3)
    return [config('TB35_REFERENCE'),
            config('P1_ELLIPTIC_BALANCED', **common),
            config('P2_ELLIPTIC_LIGHT', **(common | {'radial_reflector_cm':8.8})),
            config('P3_TWELVE_CHANNEL', **(common | {'shape':'C12_STAGGERED'})),
            config('P4_TALL_CONTOURED', **(common | {'height_cm':46.5,
                                                   'peripheral_hull_offset_cm':1.5}))]


def build_deliverables(records, summary):
    selected=selected_designs()
    shortlist=[]
    distances=[]
    for c in selected:
        spec=explicit_spec(c)
        distances.extend([dict(id=c['id'],**coolant_distance(c,z,161,1440)) for z in range(3)])
        write_json(ROOT/'configurations'/c['id']/'geometry.json',spec)
        pins=[]
        for i,p in enumerate(PINS):
            z=spec['zones'][p['zone']]
            pins.append(dict(pin_id=i,**p,fuel_radius_cm=z['fuel_radius_cm'],
                clad_outer_radius_cm=z['clad_outer_radius_cm'],
                annulus_outer_radius_cm=z['annulus_outer_radius_cm'],
                liner_outer_radius_cm=z['liner_outer_radius_cm'],
                uranium_u235_weight_fraction=z['uranium_u235_weight_fraction'],
                main_flow_weight=z['provisional_main_flow_per_pin_relative_to_mean']))
        write_csv(ROOT/'configurations'/c['id']/'pins.csv',pins)
        shortlist.append(spec['metrics'])
    write_csv(ROOT/'data/shortlist.csv',shortlist)
    write_csv(ROOT/'data/cooling_distance_proxies.csv',distances)

    # Changes isolated against the reference; all computed under identical accounting.
    controls=[config('A00_REFERENCE'),
              config('A01_PLACEMENT_ONLY',channel_centre_cm=.93),
              config('A02_ELLIPSE_ONLY',shape='E6_VERTEX_150'),
              config('A03_FACE_ROTATION_ONLY',shape='C6_FACE'),
              config('A04_TWELVE_ONLY',shape='C12_STAGGERED'),
              config('A05_FUEL_ZONING_ONLY',middle_fuel_area_factor=.96),
              config('A06_REVERSE_FUEL_ZONING_ONLY',middle_fuel_area_factor=1.04),
              config('A07_ANNULUS_ALLOCATION_ONLY',annulus_allocation='graded'),
              config('A08_HEIGHT_SHORT_ONLY',height_cm=44.5),
              config('A09_HEIGHT_TALL_ONLY',height_cm=46.5),
              config('A10_HEIGHT_48_ONLY',height_cm=48.),
              config('A11_RADIAL_REFLECTOR_ONLY',radial_reflector_cm=9.),
              config('A12_AXIAL_REFLECTOR_ONLY',axial_reflector_cm=8.3),
              config('A13_PERIPHERY_ONLY',peripheral_hull_offset_cm=1.5),
              config('A14_CORE_RADIUS_ONLY',core_radius_cm=RCORE0-.10),
              config('A15_P1_ORIGINAL_REFLECTORS',shape='E6_VERTEX_150',channel_centre_cm=.93,
                     middle_fuel_area_factor=.96,annulus_allocation='graded')]
    write_json(ROOT/'configurations/isolated_changes.json',[explicit_spec(c) for c in controls])
    write_csv(ROOT/'data/isolated_changes.csv',[evaluate(c) for c in controls])

    energy=[]
    for fraction, heat_fraction in itertools.product([.10,.175,.25],[.015,.03,.06]):
        cold=892+218*heat_fraction/fraction
        hot=892+218*(1-heat_fraction)/(1-fraction)
        energy.append(dict(auxiliary_flow_fraction=fraction,
            assumed_fraction_of_total_heat_entering_auxiliary_branch=heat_fraction,
            auxiliary_bulk_outlet_K=cold,main_bulk_outlet_K=hot,
            mixed_outlet_K=fraction*cold+(1-fraction)*hot,
            remaining_K_to_915_before_any_film_or_conduction_drop=915-cold))
    write_csv(ROOT/'data/branch_energy_balance_assumptions.csv',energy)

    # Preserve a compact, falsifiable evidence trail without modifying original files.
    inputs=['best_configurations/GEOM_OPT15_TB35/scripts/build_model_tb35.py',
            'best_configurations/GEOM_OPT15_TB35/scripts/compute_precise_mass_tb35.py',
            'best_configurations/GEOM_OPT15_TB35/scripts/run_thermodynamics_and_moose_cardinal_tb35.py',
            'best_configurations/GEOM_OPT15_TB35/scripts/run_criticality_depletion_tb35.py',
            'best_configurations/GEOM_OPT15_TB35/data/tb35_precise_mass_breakdown.csv',
            'best_configurations/GEOM_OPT15_TB35/data/tb35_stochastic_volume_mass.csv',
            'data/chain_nep.xml','PROJECT_STATUS.md','ITERATION_LOG.md']
    manifest={s:hashlib.sha256((PROJECT/s).read_bytes()).hexdigest() for s in inputs}
    write_json(ROOT/'data/input_sha256.json',manifest)
    generate_figures(selected, records)


def generate_figures(selected, records):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Ellipse, Polygon, Wedge, Patch
    plt.rcParams.update({'font.size':10,'figure.dpi':120})
    colors={'fuel':'#be9d61','clad':'#8b98a2','liner':'#e18438','mod':'#e4ecd5',
            'cool':'#a4d8ed','beo':'#ccd5e2'}
    fig,axs=plt.subplots(1,3,figsize=(13,4.6))
    for ax,c in zip(axs,[selected[0],selected[1],selected[3]]):
        spec=explicit_spec(c);z=spec['zones'][1]
        ax.add_patch(Polygon(HEX_VERTICES*10,fc=colors['mod'],ec='#435340',lw=1.5))
        for r,col in [(z['liner_outer_radius_cm'],'liner'),(z['annulus_outer_radius_cm'],'cool'),
                      (z['clad_outer_radius_cm'],'clad'),(z['fuel_radius_cm'],'fuel')]:
            ax.add_patch(Circle((0,0),r*10,fc=colors[col],ec='#49525d',lw=.6))
        n,a,b,phase=channel_parameters(c['shape'])
        for phi in phase+np.arange(n)*360/n:
            angle=math.radians(phi)
            ax.add_patch(Ellipse((10*c['channel_centre_cm']*math.cos(angle),
                                 10*c['channel_centre_cm']*math.sin(angle)),20*a,20*b,
                                angle=phi,fc=colors['cool'],ec='#2685a5',lw=1))
        ax.set_aspect('equal');ax.set_xlim(-14,14);ax.set_ylim(-13,13)
        ax.set_xlabel('x [mm]');ax.set_ylabel('y [mm]')
        ax.set_title(c['id'].replace('_','\n',1)+'\nMiddle-zone cell',fontsize=10)
    fig.suptitle('Equal auxiliary coolant area; actual OpenMC lattice orientation',fontsize=13)
    fig.legend(handles=[Patch(fc=colors[k],label=label) for k,label in
               [('fuel','UN'),('clad','Nb-1Zr'),('cool','He'),('liner','Al2O3'),('mod','YH2')]],
               loc='lower center',ncol=5,frameon=False)
    fig.tight_layout(rect=(0,.06,1,1))
    for ext in ['png','svg']:fig.savefig(ROOT/'figures'/f'unit_cells.{ext}',bbox_inches='tight')
    plt.close(fig)

    c=selected[-1];spec=explicit_spec(c)
    fig,ax=plt.subplots(figsize=(8,8))
    r=c['core_radius_cm'];ro=r+c['radial_reflector_cm']
    ax.add_patch(Circle((0,0),ro,fc=colors['beo'],ec='#34455d',lw=1.5))
    normals=HULL.equations[:,:2];d=-HULL.equations[:,2]+c['peripheral_hull_offset_cm']
    theta=np.linspace(0,2*math.pi,1441)
    dirs=np.column_stack((np.cos(theta),np.sin(theta)))
    dots=dirs@normals.T
    rho=np.minimum(r,np.min(np.divide(d[None,:],dots,out=np.full_like(dots,np.inf),where=dots>1e-14),axis=1))
    ax.add_patch(Polygon(rho[:,None]*dirs,fc=colors['mod'],ec='#557244',lw=1.1))
    zonecols=['#a5653a','#c19450','#dfba78']
    for p in PINS:
        ax.add_patch(Circle((p['x_cm'],p['y_cm']),spec['zones'][p['zone']]['liner_outer_radius_cm'],fc=colors['cool'],lw=0))
        ax.add_patch(Circle((p['x_cm'],p['y_cm']),spec['zones'][p['zone']]['fuel_radius_cm'],fc=zonecols[p['zone']],lw=0))
    rd=r+DRUM_CLEARANCE+R_DRUM
    for phi in np.arange(10)*36:
        x,y=rd*np.cos(np.deg2rad(phi)),rd*np.sin(np.deg2rad(phi))
        ax.add_patch(Circle((x,y),R_DRUM,fc='#edf1f6',ec='#657486',lw=.8))
        ax.add_patch(Wedge((x,y),R_DRUM,phi-74,phi+74,width=T_ABS,fc='#495663'))
    ax.add_patch(Circle((0,0),r,fill=False,ec='#63754e',ls='--',lw=.6))
    ax.set_aspect('equal');ax.set_xlim(-ro-1,ro+1);ax.set_ylim(-ro-1,ro+1)
    ax.set_xlabel('x [cm]');ax.set_ylabel('y [cm]')
    ax.set_title('P4: three fuel zones and contoured peripheral moderator\nAuxiliary holes omitted at this scale; drums shown outward')
    ax.legend(handles=[Patch(fc=zonecols[i],label=f'Fuel rings {label}')
                       for i,label in enumerate(['0–2','3–5','6–8'])]
              +[Patch(fc=colors['mod'],label='YH2'),Patch(fc=colors['beo'],label='BeO')],
              loc='lower center',ncol=3,fontsize=8)
    fig.tight_layout()
    for ext in ['png','svg']:fig.savefig(ROOT/'figures'/f'core_zones_and_periphery.{ext}',bbox_inches='tight')
    plt.close(fig)

    fig,axs=plt.subplots(1,2,figsize=(13,4.8))
    valid=[x for x in records if x['geometry_screen_pass']]
    for shape,color in [('C6_VERTEX','#8692a3'),('E6_VERTEX_150','#2a9d8f'),('C12_STAGGERED','#b75d69')]:
        pts=[x for x in valid if x['shape']==shape]
        axs[0].scatter([x['modeled_solids_kg'] for x in pts],[x['moderator_total_kg'] for x in pts],
                       s=5,alpha=.16,c=color,label=shape)
    for c in selected:
        e=evaluate(c)
        axs[0].scatter(e['modeled_solids_kg'],e['moderator_total_kg'],c='black',s=35)
        axs[0].annotate(c['id'].split('_')[0],(e['modeled_solids_kg'],e['moderator_total_kg']),xytext=(3,5),textcoords='offset points')
    axs[0].axvline(600,color='red',ls='--',lw=1);axs[0].axvline(575,color='gray',ls=':',lw=1)
    axs[0].set_xlabel('Modeled solid mass [kg]');axs[0].set_ylabel('Total YH2 inventory [kg]')
    axs[0].set_title('Mass and moderation tradeoff — no predicted reactivity')
    axs[0].legend(fontsize=8)
    rows=[evaluate(c) for c in selected]
    x=np.arange(len(rows));w=.36
    axs[1].bar(x-w/2,[e['auxiliary_perimeter_ratio'] for e in rows],w,label='Wetted perimeter / reference',color='#2a9d8f')
    axs[1].bar(x+w/2,[e['auxiliary_resistance_ratio_at_same_mu'] for e in rows],w,label='Laminar auxiliary resistance / reference',color='#e5ad65')
    axs[1].set_xticks(x,[c['id'].split('_')[0] for c in selected]);axs[1].set_ylabel('Ratio')
    axs[1].set_title('More coolant surface comes with a hydraulic cost')
    axs[1].legend(fontsize=8);axs[1].grid(axis='y',alpha=.2)
    fig.tight_layout()
    for ext in ['png','svg']:fig.savefig(ROOT/'figures'/f'design_tradeoffs.{ext}',bbox_inches='tight')
    plt.close(fig)


if __name__=='__main__':
    main()
