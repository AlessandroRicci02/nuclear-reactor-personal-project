#!/usr/bin/env python3
"""Independent invariants and geometric accounting checks; no reactor solves."""
import json
import math
from pathlib import Path

import numpy as np

from explore import (ROOT, H0, RF0, PITCH, R_AUX0, R_DUCT, A_AUX0, A_ANN0,
                     A_HEX, RCORE0, COUNTS, PINS, HULL, RHO, SHAPES,
                     ALL_VERTICES, selected_designs, config, evaluate,
                     explicit_spec, channel_geometry, channel_parameters,
                     moderator_envelope_area, coolant_distance)


def main():
    checks=[]
    def check(name, condition, detail):
        if not condition:
            raise AssertionError(f'{name}: {detail}')
        checks.append(dict(check=name,status='PASS',detail=detail))

    check('pin_count_and_zones',len(PINS)==217 and
          [sum(p['zone']==z for p in PINS) for z in range(3)]==[19,72,126],
          '217 pins; radial groups 19/72/126')
    centres=np.array([[p['x_cm'],p['y_cm']] for p in PINS])
    sep=np.linalg.norm(centres[:,None,:]-centres[None,:,:],axis=2)
    sep[sep==0]=np.inf
    check('nearest_neighbor_pitch',abs(sep.min()-PITCH)<1e-12,str(sep.min()))

    # Actual OpenMC API, no XML generation or neutron transport.
    import openmc
    lat=openmc.HexLattice();lat.orientation='y';lat.center=(0.,0.);lat.pitch=[PITCH]
    error=max(np.linalg.norm(lat.get_local_coordinates((p['x_cm'],p['y_cm'],0.),
                               (p['lattice_a'],p['lattice_b'],0))) for p in PINS)
    check('coordinates_match_openmc_lattice',error<1e-12,str(error))

    # Independent exact two-circle intersection for the historical baseline.
    d=.910*PITCH/2.270;r=R_AUX0;R=R_DUCT
    lens=(r*r*math.acos((d*d+r*r-R*R)/(2*d*r))
          +R*R*math.acos((d*d+R*R-r*r)/(2*d*R))
          -.5*math.sqrt((-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R)))
    ch=channel_geometry('C6_VERTEX',d)
    check('axial_duct_union_matches_exact_lens',
          abs(ch['axial_overlap_area_cm2']-6*lens)<1e-9,
          f'Overlap counted once: {6*lens:.10f} cm2/cell')
    base=evaluate(config('base'))
    correction=217*2*8.8*6*lens*2.90/1000
    check('reference_mass_correction',abs(base['modeled_solids_kg']-(593.9253679639015+correction))<1e-7,
          f'{base["modeled_solids_kg"]:.9f} kg vs historical analytic 593.925367964 kg')
    check('agreement_with_archived_stochastic_total',
          abs(base['modeled_solids_kg']-595.0156784996518)<2*.367,
          'Within archived two-sigma total mass interval; this is not a new Monte Carlo calculation')
    # Convex hull area has an independent polygon formula provided by Qhull.
    hull_area,_=moderator_envelope_area(100.,0.)
    check('polar_envelope_area_vs_convex_hull',abs(hull_area-HULL.volume)<1e-6,
          f'Polar {hull_area:.9f}; Qhull {HULL.volume:.9f} cm2')
    for shape in SHAPES:
        n,a,b,_=channel_parameters(shape)
        check(f'fixed_auxiliary_area_{shape}',abs(n*math.pi*a*b-A_AUX0)<1e-12,
              f'{n*math.pi*a*b:.12f} cm2 per cell')
    for c in selected_designs():
        e=evaluate(c);s=explicit_spec(c)
        fuel_volume=sum(z['pins']*math.pi*z['fuel_radius_cm']**2*c['height_cm'] for z in s['zones'])
        u235=sum(z['pins']*math.pi*z['fuel_radius_cm']**2*c['height_cm']*13.5/1000*238/252*
                  z['uranium_u235_weight_fraction'] for z in s['zones'])
        ann_area=sum(z['pins']*z['annulus_area_cm2_per_pin'] for z in s['zones'])
        check(c['id']+'_inventory',abs(fuel_volume-217*math.pi*RF0**2*H0)<1e-8 and
              abs(u235-base['u235_kg'])<1e-10 and abs(ann_area-217*A_ANN0)<1e-12,
              'Fuel volume, U235 mass, and total main annular area conserved')
        check(c['id']+'_flow_normalization',abs(sum(z['pins']*
              z['provisional_main_flow_per_pin_relative_to_mean'] for z in s['zones'])-217)<1e-10,
              'Per-pin target flow weights normalize to 217')
        check(c['id']+'_containment',c['core_radius_cm']>np.linalg.norm(ALL_VERTICES,axis=1).max(),
              'All occupied hex cells retained; no cropped fuel or channels')
        # Full envelope volume closure: solids plus main/auxiliary coolant.
        solid_volume=(e['fuel_kg']/13.5+e['cladding_kg']/8.57+e['liner_kg']/3.95+
                      e['moderator_total_kg']/4.3+e['beo_kg']/2.9+e['b4c_kg']/2.52)*1000
        coolant_volume=((217*A_AUX0+ann_area)*c['height_cm'] + 217*2*c['axial_reflector_cm']*
                        channel_geometry(c['shape'],c['channel_centre_cm'])['axial_union_area_cm2'])
        total_volume=math.pi*e['outer_radius_cm']**2*e['total_height_cm']
        check(c['id']+'_volume_closure',abs(solid_volume+coolant_volume-total_volume)<1e-7,
              f'Residual {solid_volume+coolant_volume-total_volume:.3g} cm3')
        if c['id']!='TB35_REFERENCE':
            check(c['id']+'_study_geometry',e['geometry_screen_pass'] and e['modeled_solids_kg']<600,
                  '0.50 mm study ligament and 7 mm outer drum web criteria; modeled solids <600 kg')
        low=coolant_distance(c,1,161,1440)
        high=coolant_distance(c,1,321,2880)
        check(c['id']+'_distance_sampling',
              abs(low['max_distance_to_aux_mm']-high['max_distance_to_aux_mm'])<.05 and
              abs(low['p95_distance_to_aux_mm']-high['p95_distance_to_aux_mm'])<.05,
              f'Max {low["max_distance_to_aux_mm"]:.5f}/{high["max_distance_to_aux_mm"]:.5f} mm; '
              f'p95 {low["p95_distance_to_aux_mm"]:.5f}/{high["p95_distance_to_aux_mm"]:.5f} mm')
        check(c['id']+'_uncomputed_physics_null',all(x is None for x in s['performance'].values()),
              'No invented reactor performance values')
    result=dict(status='PASS_GEOMETRY_AND_ACCOUNTING_ONLY',checks=checks,
                count=len(checks),openmc_used_for='Coordinate API check only; no transport or volume run',
                excludes='Neutronics, depletion, thermal hydraulics, structural qualification, hydrogen retention')
    (ROOT/'data/validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'{len(checks)} geometry/accounting checks passed. No reactor simulation executed.')


if __name__=='__main__':
    main()
