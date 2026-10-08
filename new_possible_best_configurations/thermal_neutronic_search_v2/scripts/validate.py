"""Physics/geometry checks; synthetic thermal loads are not candidate results."""
import argparse
import json
import math
import os
import subprocess
import numpy as np
from catalog import ROOT,catalog,layers,pins,main as write_catalog
from model import build,OPENMC
from thermal import scenario,cardinal_check,radial


def main(transport_check=False):
    write_catalog();cc=catalog();out={'candidate_count':len(cc),'checks':[]}
    assert len({c['id'] for c in cc})==30
    for c in cc:
        s=c['spec'];ll=layers(s);pp=pins(s)
        assert s['pitch_cm']>2*ll[-1]['r_outer_cm']
        assert s['core_radius_cm']>max(math.hypot(p['x_cm'],p['y_cm'])+ll[-1]['r_outer_cm'] for p in pp)
        assert s['drum_centres_radius_cm']-s['drum_radius_cm']>s['core_radius_cm']
        assert s['drum_centres_radius_cm']+s['drum_radius_cm']<s['core_radius_cm']+s['radial_reflector_cm']
        assert 2*s['drum_centres_radius_cm']*math.sin(math.pi/s['drum_count'])>2*s['drum_radius_cm']
        assert 850<=s['inlet_K']<=900 and 1100<=s['outlet_K']<=1150
        for temp in (-100.,100.):
            factor=1+3.9e-5*temp;r=s['fuel_radius_cm']*math.sqrt(factor)
            assert r<ll[1]['r_outer_cm']
            assert abs((math.pi*r*r*13.5/factor)-(math.pi*s['fuel_radius_cm']**2*13.5))<1e-10
    out['checks'].append('All 30: pin/drum separation, containment, cycle bounds, fuel mass conservation')
    q=np.array([100.,200.]);t,qi,qo=radial(q,np.array([850.,860.]),np.array([880.,900.]),.010,.012,10.,np.array([800.,1000.]),np.array([500.,700.]))
    assert np.allclose(qi+qo,q,rtol=1e-9)
    out['checks'].append('Two-sided cylindrical conduction: heat removal equals source')
    for i in (0,6):
        c=cc[i];s=c['spec'];n=len(pins(s));z=8
        maps={x['name']:np.zeros((n,z)) for x in layers(s) if x['material']!='void'}
        maps['fuel'][:]=380000/n/z;maps['moderator'][:]=10000/n/z
        maps['reflector_allocated']=np.full((n,z),5000/n/z)
        maps['reflector_conservative_allocated']=np.full((n,z),15000/n/z)
        errors={k:np.zeros_like(v) for k,v in maps.items()}
        for worst in (False,True):
            result,data=scenario(s,maps,errors,dict(material_W={'reflector':5000.}),worst=worst)
            assert abs(result['mixed_outlet_K']-s['outlet_K'])<1e-7
            assert result['hydraulic_closure_valid']
            if not worst:
                check=cardinal_check(ROOT/'validation/cardinal',c['id'],result,data)
                assert check['max_moderator_K']<=result['max_moderator_K']+1
            out['checks'].append(dict(candidate=c['id'],synthetic_load=True,worst=worst,
                source_W=400000,analytic_peak_K=result['max_moderator_K'],mixed_outlet_K=result['mixed_outlet_K']))
    c=cc[0];s=c['spec']
    for theta in (0.,180.):
        d=ROOT/'validation'/f'final_angle{theta:.0f}'
        m,info=build(d,s,theta=theta,particles=100,batches=12,inactive=3)
        center=s['drum_centres_radius_cm'];offset=s['drum_radius_cm']-s['absorber_thickness_cm']/2
        for label,x in [('outward',center+offset),('inward',center-offset)]:
            path=m.geometry.find((x,0.,0.));name=path[-1].fill.name
            assert (name=='B4C')==(label==('outward' if theta==0 else 'inward'))
        for material in info['materials']:
            if material.name.startswith('UN_'):
                ad=material.get_nuclide_atom_densities()
                ratio=(ad['U235']+ad['U238'])/(ad['N14']+ad['N15'])
                assert abs(ratio-1)<1e-6
        if transport_check and theta==0:
            env=os.environ.copy();env['OMP_NUM_THREADS']='2'
            with (d/'geometry_debug.log').open('w') as log:
                cp=subprocess.run([OPENMC,'-s','2','--geometry-debug'],cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT)
            assert cp.returncode==0,'OpenMC geometry-debug failed'
    out['checks'].append('Absorber orientation at OUT/IN and UN atomic stoichiometry verified')
    out['geometry_debug_run']=transport_check
    (ROOT/'validation/validation_report.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--transport-check',action='store_true')
    main(p.parse_args().transport_check)
