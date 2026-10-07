"""Statistically verified drum searches in isolated campaign directories."""
import json
import math
from p1_model import CAMPAIGN
from p1_transport import run_point


def critical(name, reference=False, fuel_file=None, endpoints=None):
    opts={'reference':reference}
    if fuel_file: opts['fuel_file']=str(fuel_file)
    points=[]
    for j,angle in enumerate((0.,180.)):
        r=endpoints[j] if endpoints else run_point(name+f'/endpoint_{angle:.0f}',theta=angle,seed=21001+2*j,**opts)
        points.append((angle,r))
    out,ins=points[0][1],points[1][1]
    result={'points':points,'theta_deg':None,'critical':False}
    if out['keff']+1.96*out['sigma']<1:
        result['status']='SUBCRITICAL_AT_ALL_DRUM_ANGLES'
    elif ins['keff']-1.96*ins['sigma']>1:
        result['status']='SUPERCRITICAL_AT_ALL_DRUM_ANGLES'
    else:
        lo,hi=points
        for i in range(24):
            # Safeguarded interpolation; every proposed angle gets a transport solve.
            den=lo[1]['reactivity_pcm']-hi[1]['reactivity_pcm']
            f=lo[1]['reactivity_pcm']/den if den>0 else .5
            theta=lo[0]+max(.15,min(.85,f))*(hi[0]-lo[0])
            r=run_point(name+f'/search_{i:02d}',sigma=1.5e-4,theta=theta,seed=21101+37*i,**opts)
            points.append((theta,r))
            if abs(r['reactivity_pcm'])+1.96*r['reactivity_sigma_pcm']<=60.:
                # Fresh seed protects against acceptance driven by the search noise.
                check=run_point(name+f'/verify_{i:02d}',sigma=1.5e-4,theta=theta,seed=41001+41*i,**opts)
                points.append((theta,check))
                if abs(check['reactivity_pcm'])+1.96*check['reactivity_sigma_pcm']<=60.:
                    result.update(theta_deg=theta,critical=True,status='VERIFIED_WITHIN_60_PCM_AT_95_PERCENT',verification=check)
                    break
                r=check
            if r['keff']>1:lo=(theta,r)
            else:hi=(theta,r)
        else:raise RuntimeError(f'Critical search did not converge: {name}')
    result['points']=points
    path=CAMPAIGN/'neutronics'/name/'control.json';path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2)+'\n')
    return result
