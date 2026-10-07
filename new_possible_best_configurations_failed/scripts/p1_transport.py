#!/usr/bin/env python3
"""Resumable, logged standalone OpenMC calculations; no historical writes."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import numpy as np
import openmc
from p1_model import build, CAMPAIGN, OPENMC, DEFAULT_XS, NOMINAL


def run_point(name, sigma=2e-4, **kwargs):
    d=CAMPAIGN/'neutronics'/name;d.mkdir(parents=True,exist_ok=True)
    request=dict(sigma=sigma,**kwargs)
    fuel_file=kwargs.pop('fuel_file',None)
    if fuel_file:
        request['fuel_sha256']=hashlib.sha256(Path(fuel_file).read_bytes()).hexdigest()
        kwargs['fuel_materials']={m.id:m for m in openmc.Materials.from_xml(fuel_file) if m.depletable}
    signature=hashlib.sha256(json.dumps(request,sort_keys=True).encode()).hexdigest()
    resultpath=d/'result.json'
    if resultpath.exists():
        r=json.loads(resultpath.read_text())
        if r['request_sha256']!=signature:raise RuntimeError(f'Different request already exists: {d}')
        if not r.get('target_achieved',False):raise RuntimeError(f'Cached point missed target: {d}')
        with openmc.StatePoint(r['statepoint']) as sp:
            r['batches']=int(sp.current_batch)
            r['requested_batches']=int(sp.n_batches)
        resultpath.write_text(json.dumps(r,indent=2)+'\n')
        print('CACHE',name,r['keff'],r['sigma'],flush=True);return r
    args=dict(particles=20000,batches=180,inactive=60,seed=173,target_sigma=sigma)
    args.update(kwargs)
    model,info=build(d,**args)
    if sigma:
        model.settings.trigger_max_batches=20000
        model.settings.export_to_xml(d/'settings.xml')
    command=[OPENMC,'-s','4'];env=os.environ.copy();env['OPENMC_CROSS_SECTIONS']=str(kwargs.get('xs',DEFAULT_XS));env['OMP_NUM_THREADS']='4'
    start=time.time()
    status={'name':name,'directory':str(d),'started_unix':start,'command':command,'request':request}
    (CAMPAIGN/'current_transport.json').write_text(json.dumps(status,indent=2)+'\n')
    print('RUN',name,'target sigma',sigma,flush=True)
    with (d/'run.log').open('w') as log:
        proc=subprocess.Popen(command,cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT)
        status['pid']=proc.pid;(CAMPAIGN/'current_transport.json').write_text(json.dumps(status,indent=2)+'\n')
        ret=proc.wait()
    if ret:raise RuntimeError(f'OpenMC exit {ret}; see {d}/run.log')
    spath=max(d.glob('statepoint.*.h5'),key=lambda p:int(p.name.split('.')[1]))
    with openmc.StatePoint(spath) as sp:
        k=sp.keff.n;sd=sp.keff.s;ent=np.asarray(sp.entropy)
        r={'keff':k,'sigma':sd,'sigma_pcm_k':1e5*sd,'reactivity_pcm':1e5*(1-1/k),
           'reactivity_sigma_pcm':1e5*sd/k**2,'batches':int(sp.current_batch),'requested_batches':int(sp.n_batches),'inactive':int(sp.n_inactive),
           'particles':int(sp.n_particles),'wall_seconds':time.time()-start,'statepoint':str(spath),
           'request_sha256':signature,'request':request,'target_achieved':sigma is None or sd<=sigma,
           'entropy_inactive_last20_mean':float(np.mean(ent[max(0,sp.n_inactive-20):sp.n_inactive])),
           'entropy_active_first20_mean':float(np.mean(ent[sp.n_inactive:sp.n_inactive+20]))}
        np.savetxt(d/'entropy.csv',np.column_stack([np.arange(1,len(ent)+1),ent]),delimiter=',',header='generation,entropy',comments='')
    resultpath.write_text(json.dumps(r,indent=2)+'\n')
    status.update(completed_unix=time.time(),returncode=ret,result=str(resultpath))
    (CAMPAIGN/'current_transport.json').write_text(json.dumps(status,indent=2)+'\n')
    print('DONE',name,k,sd,'seconds',r['wall_seconds'],flush=True)
    if not r['target_achieved']:raise RuntimeError(f'Statistical target not achieved in {name}')
    return r


def feedback(theta=0.):
    # Independent random seeds: uncertainty propagation assumes independent runs.
    out=[]
    for i,(name,key,delta) in enumerate([('fuel_temperature','fuel',200.),('moderator_temperature','moderator',100.),('reflector_temperature','reflector',100.),('clad_temperature','clad',100.),('coolant_temperature_fixed_density','coolant',100.)]):
        results=[]
        for j,sign in enumerate([-1,1]):
            T=NOMINAL[key]+sign*delta
            opts={'temperatures':{key:T}}
            if key=='coolant':opts['coolant_density']=T/NOMINAL['coolant']
            results.append(run_point(f'feedback/{name}_{"minus" if sign<0 else "plus"}',sigma=1.5e-4,theta=theta,seed=6001+100*i+2*j,**opts))
        out.append(coefficient(name,results,2*delta,'pcm/K'))
        save_feedback(out)
    for i,(name,arg,delta) in enumerate([('moderator_density','moderator_density',.02),('coolant_density','coolant_density',.5)]):
        r=[run_point(f'feedback/{name}_{"minus" if s<0 else "plus"}',sigma=1.5e-4,theta=theta,seed=7101+100*i+2*j,**{arg:1+s*delta}) for j,s in enumerate([-1,1])]
        out.append(coefficient(name,r,200*delta,'pcm/percent_density'));save_feedback(out)
    r=[run_point(f'feedback/isothermal_{"minus" if s<0 else "plus"}',sigma=1.5e-4,theta=theta,seed=7501+2*j,temperatures={k:v+s*100 for k,v in NOMINAL.items()}) for j,s in enumerate([-1,1])]
    out.append(coefficient('isothermal_fixed_geometry_fixed_helium_pressure',r,200,'pcm/K'));save_feedback(out)
    return out


def coefficient(name,results,dx,units):
    a,b=results;value=(b['reactivity_pcm']-a['reactivity_pcm'])/dx
    sigma=np.hypot(a['reactivity_sigma_pcm'],b['reactivity_sigma_pcm'])/dx
    return dict(name=name,value=value,sigma=float(sigma),units=units,span=dx,
                negative_at_95pct=bool(value+1.96*sigma<0),positive_at_95pct=bool(value-1.96*sigma>0),
                reference='Nominal prescribed temperatures; fixed drum angle recorded in endpoint manifests; geometry and nonperturbed densities held fixed',
                minus_statepoint=a['statepoint'],plus_statepoint=b['statepoint'])


def save_feedback(out):
    (CAMPAIGN/'data/feedback_coefficients.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('task',choices=['pilot','bol','feedback']);a=ap.parse_args()
    if a.task=='pilot':
        for ref,name in [(False,'p1'),(True,'tb35')]:
            run_point(name+'_pilot_out',sigma=None,reference=ref,theta=0.,particles=10000,batches=160,inactive=60,seed=1103)
    elif a.task=='bol':
        for theta in [0.,180.]:run_point(f'p1_bol_{theta:.0f}deg',theta=theta,seed=1201+int(theta))
        run_point('tb35_corrected_bol_0deg',reference=True,theta=0.,seed=1601)
    else:feedback()
