"""Restartable S2-R1 OpenMC temperature finite differences and mass audit."""
from __future__ import annotations
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time
import openmc
from s2_model import build, ROOT, XS, TEMP

CAM=ROOT/'campaign';RUNS=CAM/'feedback';RUNS.mkdir(parents=True,exist_ok=True)
OPENMC='/home/aleri/miniforge3/envs/openmc-env/bin/openmc'
TARGET=.00045  # sampling uncertainty in k, 1 sigma; refine if signs unresolved
PARTICLES=10000
PAIRS=[('fuel',200.),('moderator',100.),('total',100.)]

def atomic(path,obj):
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2)+'\n');tmp.replace(path)

def run_case(name,change,seed):
    d=RUNS/name;d.mkdir(exist_ok=True)
    req=dict(name=name,change=change,seed=seed,target_sigma=TARGET,particles=PARTICLES,
             xs=str(XS),config_sha256=hashlib.sha256((ROOT/'configuration.json').read_bytes()).hexdigest())
    sig=hashlib.sha256(json.dumps(req,sort_keys=True).encode()).hexdigest()
    rpath=d/'result.json'
    if rpath.exists():
        r=json.loads(rpath.read_text())
        if r['request_sha256']!=sig:raise RuntimeError('Result request mismatch: '+name)
        return r
    model,_=build(d,temperatures={k:v+change[k] for k,v in TEMP.items() if k in change},
                  particles=PARTICLES,batches=180,inactive=60,seed=seed,slices=6)
    model.settings.keff_trigger={'type':'std_dev','threshold':TARGET}
    model.settings.trigger_active=True;model.settings.trigger_max_batches=1800
    model.settings.trigger_batch_interval=20
    model.settings.export_to_xml(d/'settings.xml')
    threads=int(os.environ.get('S2_THREADS','4'))
    env=os.environ.copy();env.update(OPENMC_CROSS_SECTIONS=str(XS),OMP_NUM_THREADS='4',
                                     OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    env['OMP_NUM_THREADS']=str(threads)
    atomic(CAM/'status.json',dict(stage='running',case=name,started_unix=time.time()))
    with (d/'run.log').open('w') as log:
        cp=subprocess.run([OPENMC,'-s',str(threads)],cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT)
    if cp.returncode:raise RuntimeError(f'{name} failed; see {d}/run.log')
    s=max(d.glob('statepoint.*.h5'),key=lambda p:int(p.name.split('.')[1]))
    with openmc.StatePoint(s) as sp:
        k=float(sp.keff.n);sigma=float(sp.keff.s)
        r=dict(request=req,request_sha256=sig,keff=k,sigma_k=sigma,
               reactivity_pcm=1e5*(1-1/k),sigma_reactivity_pcm=1e5*sigma/k**2,
               active_batches=int(sp.current_batch-sp.n_inactive),inactive_batches=int(sp.n_inactive),
               particles=int(sp.n_particles),statepoint=str(s),target_achieved=sigma<=TARGET)
    atomic(rpath,r);return r

def mass():
    d=CAM/'mass';v=openmc.VolumeCalculation.from_hdf5(d/'volume_1.h5')
    mats=openmc.Materials.from_xml(d/'materials.xml')
    rows=[]
    box=(2*28.5)**2*65
    for m in mats:
        vol=v.volumes[m.id]
        rho=m.get_mass_density()
        rows.append(dict(material=m.name,volume_cm3=vol.n,sigma_volume_cm3=vol.s,
                         density_g_cm3=rho,mass_kg=vol.n*rho/1000))
    with (d/'mass_breakdown.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    solid=[r for r in rows if not r['material'].startswith('He_')]
    mass=sum(x['mass_kg'] for x in solid)
    rhobar=sum(x['volume_cm3']/box*x['density_g_cm3'] for x in solid)
    rho2bar=sum(x['volume_cm3']/box*x['density_g_cm3']**2 for x in solid)
    sigma=box/1000*math.sqrt((rho2bar-rhobar**2)/2_000_000)
    ref=json.loads((ROOT/'screening.json').read_text())['solid_mass_kg']
    out=dict(method='OpenMC stochastic material volumes; 2,000,000 samples',
             solid_mass_kg=mass,sampling_sigma_kg=sigma,analytic_mass_kg=ref,
             difference_kg=mass-ref,helium_mass_kg=sum(x['mass_kg'] for x in rows if x['material'].startswith('He_')),
             prescribed_hardware_allowance_kg=65.,prescribed_contingency_kg=15.,
             conditional_total_kg=mass+80.,conditional_total_sampling_sigma_kg=sigma,
             excluded='Shield, conversion equipment, radiators, external piping and vehicle structure; actual hardware dimensions unknown')
    atomic(d/'result.json',out);return out

def main():
    mr=mass();print('Mass:',mr,flush=True)
    coeffs=[]
    for j,(which,delta) in enumerate(PAIRS):
        pair=[]
        for sign in (-1,1):
            name=which+('_minus' if sign<0 else '_plus')
            changes={('moderator' if which=='moderator' else 'fuel'):sign*delta} if which!='total' else {k:sign*delta for k in TEMP}
            # Changes values are deltas; run_case applies them to nominal.
            r=run_case(name,changes,7100+2*j+(sign+1)//2)
            print(name,r['keff'],r['sigma_k'],r['active_batches'],flush=True)
            pair.append(r)
        lo,hi=pair;value=(hi['reactivity_pcm']-lo['reactivity_pcm'])/(2*delta)
        sigma=math.hypot(lo['sigma_reactivity_pcm'],hi['sigma_reactivity_pcm'])/(2*delta)
        coeffs.append(dict(kind=which,alpha_pcm_K=value,sampling_sigma_pcm_K=sigma,
             low_95_pcm_K=value-1.96*sigma,high_95_pcm_K=value+1.96*sigma,
             minus=lo['request']['name'],plus=hi['request']['name'],delta_K=delta,
             definition='Eigenvalue reactivity central difference; all drums 0 deg; fixed geometry and solid densities; helium fixed pressure; total all material temperatures varied.'))
        atomic(CAM/'feedback_coefficients.json',coeffs)
    atomic(CAM/'status.json',dict(stage='complete',completed_unix=time.time()))

if __name__=='__main__':
    try:main()
    except Exception as exc:
        atomic(CAM/'status.json',dict(stage='failed',error=str(exc),failed_unix=time.time()))
        raise
