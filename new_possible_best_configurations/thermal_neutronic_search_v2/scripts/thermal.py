"""Counterflow serial He energy balances, shield radiation, YH2 radial bounds.

An operating-state screen, not full-core CFD. Cardinal verifies the hottest
RZ sleeve. Acceptance uses the conservative radial slice bound, not the lower
axially smoothed Cardinal peak. The uncertainty envelope is explicitly bounded.
"""
from pathlib import Path
import csv
import json
import math
import os
import subprocess
import numpy as np
import openmc
from catalog import layers,pins
from model import SLICES

CARDINAL='/home/aleri/projects/cardinal/cardinal-opt'
CP=5193.16;SIGMA=5.670374419e-8


def extract(directory,info,spec):
    d=Path(directory);n=len(info['pins']);s=SLICES
    state=max(d.glob('statepoint.*.h5'),key=lambda p:int(p.name.split('.')[1]))
    with openmc.StatePoint(state,autolink=False) as sp:
        t=sp.get_tally(name='all_material_heating')
        values=np.asarray(t.mean).ravel();scale=spec['power_W']/float(values.sum())
        names=[m.name for m in info['materials']]
        material={name:float(v*scale) for name,v in zip(names,values)}
        maps={};sigma={}
        for name in info['scored']:
            t=sp.get_tally(name=name+'_heating')
            maps[name]=np.asarray(t.mean).reshape(n,s)*scale
            sigma[name]=np.asarray(t.std_dev).reshape(n,s)*scale
        for name in info.get('scored_end',{}):
            t=sp.get_tally(name='end_'+name+'_heating')
            values=np.asarray(t.mean).reshape(n,2)*scale;std=np.asarray(t.std_dev).reshape(n,2)*scale
            maps[name][:,0]+=values[:,0];maps[name][:,-1]+=values[:,1]
            sigma[name][:,0]+=std[:,0];sigma[name][:,-1]+=std[:,1]
        rt=sp.get_tally(name='reflector_spatial_heating');mesh=rt.find_filter(openmc.MeshFilter).mesh
        lower=np.asarray(mesh.lower_left);width=(np.asarray(mesh.upper_right)-lower)/np.asarray(mesh.dimension)
        centers=np.asarray([lower+(np.asarray(index)-.5)*width for index in mesh.indices])
        rv=np.asarray(rt.mean).ravel()*scale;rs=np.asarray(rt.std_dev).ravel()*scale
        xy=np.asarray([[p['x_cm'],p['y_cm']] for p in info['pins']])
        dist=((centers[:,None,:2]-xy[None,:,:])**2).sum(axis=2)
        nearest=np.argsort(dist,axis=1)[:,:3]
        zs=np.clip(((centers[:,2]+spec['height_cm']/2)/spec['height_cm']*s).astype(int),0,s-1)
        nominal_ref=np.zeros((n,s));worst_ref=np.zeros((n,s));ref_sigma=np.zeros((n,s))
        for i,z in enumerate(zs):
            nominal_ref[nearest[i,0],z]+=rv[i];ref_sigma[nearest[i,0],z]+=rs[i]
            # Give each of the three nearest sleeves the full upper bin power
            # for the conservative routing envelope; no ideal cold reflector.
            for p in nearest[i]:worst_ref[p,z]+=rv[i]+3*rs[i]
        if abs(rv.sum()-material['reflector'])>.01*spec['power_W']:
            raise RuntimeError('Reflector spatial/material heating closure failed')
        maps['reflector_allocated']=nominal_ref;sigma['reflector_allocated']=ref_sigma
        maps['reflector_conservative_allocated']=worst_ref;sigma['reflector_conservative_allocated']=np.zeros_like(worst_ref)
    # Sixfold rotation is an exact symmetry of hex pins + 12 identical drums.
    # Keep the mean standard deviation as a conservative uncertainty estimate;
    # do not assume independent normalized heat bins without their covariance.
    groups={}
    for pin in info['pins']:
        q,r=pin['q'],pin['r'];orbit=[]
        for _ in range(6):orbit.append((q,r));q,r=-r,q+r
        groups.setdefault(min(orbit),[]).append(pin['id'])
    for group in groups.values():
        for name in maps:
            maps[name][group]=maps[name][group].mean(axis=0)
            sigma[name][group]=sigma[name][group].mean(axis=0)
    np.savez(d/'heating_W.npz',**maps,**{'sigma_'+k:v for k,v in sigma.items()})
    result=dict(material_W=material,normalization_W_per_eV=scale,statepoint=str(state),
        note='400 kW deposited-heat normalization; photon leakage is not added to coolant power.')
    (d/'deposition.json').write_text(json.dumps(result,indent=2)+'\n')
    return maps,sigma,result


def props(T,pressure):
    mu=3.674e-7*np.asarray(T)**.70
    k=2.682e-3*(1+1.123e-8*pressure)*np.asarray(T)**(.71*(1-2e-9*pressure))
    return mu,k


def film(massflow,ri,ro,T,pressure,mult=1.):
    area=math.pi*(ro**2-ri**2);dh=2*(ro-ri);mu,k=props(T,pressure)
    re=massflow*dh/(area*mu);pr=mu*CP/k
    # Use a conservative laminar closure through the transition regime.
    f=np.where(re<3000.,64/np.maximum(re,1.),(.79*np.log(np.maximum(re,3000.))-1.64)**-2)
    nut=(f/8)*(re-1000)*pr/(1+12.7*np.sqrt(f/8)*(pr**(2/3)-1))
    nu=np.where(re<4000.,4.0,np.maximum(4.,nut))
    return nu*k/dh*mult,re,f,area,dh


def radial(q,ti,to,ri,ro,k,hi,ho):
    """Exact uniform-source cylindrical sleeve with two Robin boundaries."""
    if np.all(np.asarray(ho)==0):
        rise=q/(2*math.pi*ri*hi)+q*(2*ro**2*math.log(ro/ri)-(ro**2-ri**2))/(4*k*math.pi*(ro**2-ri**2))
        return ti+rise,q,np.zeros_like(q)
    qvol=q/(math.pi*(ro**2-ri**2))
    a11=k/ri-hi*math.log(ri);a12=-hi
    a21=-k/ro-ho*math.log(ro);a22=-ho
    b1=qvol*ri/2-hi*(qvol*ri**2/(4*k)+ti)
    b2=-qvol*ro/2-ho*(qvol*ro**2/(4*k)+to)
    det=a11*a22-a12*a21
    c1=(b1*a22-a12*b2)/det;c2=(a11*b2-b1*a21)/det
    temp=lambda r:-qvol*r*r/(4*k)+c1*np.log(r)+c2
    rpeak=np.sqrt(np.maximum(0,2*k*c1/np.maximum(qvol,1e-30)))
    rpeak=np.clip(rpeak,ri,ro)
    peak=np.maximum(np.maximum(temp(ri),temp(ro)),temp(rpeak))
    qi=2*math.pi*ri*hi*(temp(ri)-ti);qo=2*math.pi*ro*ho*(temp(ro)-to)
    return peak,qi,qo


def scenario(spec,maps,errors,deposition,*,worst=False):
    ll={x['name']:x for x in layers(spec)};n=len(pins(spec));s=SLICES
    h=spec['height_cm']/100;dz=h/s
    ri=ll['moderator']['r_inner_cm']/100;ro=ll['moderator']['r_outer_cm']/100
    rf=spec['fuel_radius_cm']/100
    power=spec['power_W'];totalflow=power/CP/(spec['outlet_K']-spec['inlet_K'])
    fu=maps['fuel'];pinw=fu.sum(axis=1);pinw/=pinw.sum()
    flatten=.30 if worst else .10
    flows=totalflow*((1-flatten)*pinw+flatten/n);cap=flows*CP
    shape=fu/np.maximum(fu.sum(axis=1)[:,None],1e-30)
    qmod=np.zeros_like(fu)
    coldnames=['moderator','cold_wall','Mo_shield','inner_can','inner_barrier',
               'outer_barrier','outer_can','outer_jacket','cold_He','outer_cold_He','expansion_clearance']
    for name in coldnames:
        if name in maps:qmod+=maps[name]+(3*errors[name] if worst else 0)
    # Route ALL reflector gamma heating to the cold side, including axial and
    # external reflector. This is a conservative energy routing assumption.
    qmod+=maps['reflector_conservative_allocated'] if worst else maps['reflector_allocated']
    if worst:qmod*=1.05
    if qmod.sum()>=power:raise ValueError('Cold-side deposition exceeds total power')
    qhot=(power-qmod.sum())*pinw[:,None]*shape
    hotri=ll['hot_He']['r_inner_cm']/100;hotro=ll['hot_He']['r_outer_cm']/100
    coldri=ll['cold_He']['r_inner_cm']/100;coldro=ll['cold_He']['r_outer_cm']/100
    hi,re,fr,area,dh=film(flows,coldri,coldro,900.,spec['pressure_Pa'],.70 if worst else 1.)
    resistance=1/(2*math.pi*coldro*hi)
    resistance+=math.log(ri/coldro)/(2*math.pi*6.) # conservative for metal+barrier
    heff=1/(2*math.pi*ri*resistance)
    outer=spec['outer_cold_channel'];fm=np.ones(n)
    ho=np.zeros(n)
    if outer:
        or_i=ll['outer_cold_He']['r_inner_cm']/100;or_o=ll['outer_cold_He']['r_outer_cm']/100
        ain=area;ao=math.pi*(or_o**2-or_i**2)
        fm=np.full(n,float(np.mean(ain/(ain+ao))))
        if worst:fm=np.clip(fm+.15,.25,.85)
        hi,re,fr,area,dh=film(flows*fm,coldri,coldro,900.,spec['pressure_Pa'],.70 if worst else 1.)
        heff=1/(2*math.pi*ri*(1/(2*math.pi*coldro*hi)+math.log(ri/coldro)/(2*math.pi*6.)))
        hout,_,_,_,_=film(flows*(1-fm),or_i,or_o,900.,spec['pressure_Pa'],.70 if worst else 1.)
        # Explicit outer clearance/canning limits cooling: gas gap conduction
        # is included, so the outer cold channel is not credited perfectly.
        kgas=float(props(900.,spec['pressure_Pa'])[1])
        rgap=ll['expansion_clearance']['r_outer_cm']/100
        resistance=math.log(rgap/ro)/(2*math.pi*kgas)+math.log(or_i/rgap)/(2*math.pi*6.)+1/(2*math.pi*or_i*hout)
        ho=1/(2*math.pi*ro*resistance)
    km=5. if worst else 10.;eps=.30 if worst else .15
    bridge=(35. if worst else 12.)/n/s
    mainwall=ll['hot_wall']['r_outer_cm']/100
    shieldi=ll['vacuum_inner']['r_outer_cm']/100;shieldo=ll['Mo_shield']['r_outer_cm']/100
    coldwall=ll['vacuum_outer']['r_outer_cm']/100
    denom=1/eps+(mainwall/shieldi)*(1/eps-1)+(mainwall/shieldo)/eps+(mainwall/coldwall)*(1/eps-1)
    hh,_,_,_,_=film(flows,hotri,hotro,1100.,spec['pressure_Pa'],.7 if worst else 1.)
    # Basis flux partitions of the radial cell, evaluated per pin for films.
    _,fi,_=radial(np.ones(n),np.zeros(n),np.zeros(n),ri,ro,km,heff,ho)
    _,gneg,_=radial(np.zeros(n),np.ones(n),np.zeros(n),ri,ro,km,heff,ho)
    conductance=-gneg*dz
    leak=np.full((n,s),20000/n/s);ci=cap*fm;co=cap*(1-fm)
    def network(currentleak):
        ei=np.zeros((n,s+1));eo=np.zeros_like(ei)
        ei[:,s]=spec['inlet_K'];eo[:,s]=spec['inlet_K']
        cm=np.zeros((n,s));om=np.zeros_like(cm)
        for j in reversed(range(s)):
            q=qmod[:,j];li=currentleak[:,j]
            if outer:
                qi=(fi*q-conductance*(ei[:,j+1]-eo[:,j+1]+li/(2*ci)-q/(2*co)))/(1+conductance/(2*ci)+conductance/(2*co))
                ei[:,j]=ei[:,j+1]+(qi+li)/ci;eo[:,j]=eo[:,j+1]+(q-qi)/co
            else:ei[:,j]=ei[:,j+1]+(q+li)/ci;eo[:,j]=ei[:,j]
            cm[:,j]=(ei[:,j]+ei[:,j+1])/2;om[:,j]=(eo[:,j]+eo[:,j+1])/2
        hotedges=np.zeros_like(ei);hotedges[:,0]=fm*ei[:,0]+(1-fm)*eo[:,0]
        for j in range(s):hotedges[:,j+1]=hotedges[:,j]+(qhot[:,j]-currentleak[:,j])/cap
        hm=(hotedges[:,:-1]+hotedges[:,1:])/2
        # Gamma-heated hot jacket can exceed bulk gas; bounded wall offset.
        qwall=maps['hot_wall']+(3*errors['hot_wall'] if worst else 0)
        wall=hm+qwall/(2*math.pi*mainwall*dz*hh[:,None])+(50 if worst else 10)
        return ei,eo,cm,om,hotedges,hm,wall
    for it in range(100):
        ei,eo,cm,om,he,hm,wall=network(leak)
        new=np.maximum(0.,2*math.pi*mainwall*dz*SIGMA*(wall**4-cm**4)/denom)+bridge*np.maximum(0,wall-cm)
        # Include both axial reflector jackets and their exposed hot-return
        # bore, rather than charging only the active-core insulation area.
        ah=spec['axial_reflector_cm']/100
        new[:,0]+=np.maximum(0.,2*math.pi*mainwall*ah*SIGMA*(wall[:,0]**4-cm[:,0]**4)/denom)
        new[:,-1]+=np.maximum(0.,2*math.pi*mainwall*ah*SIGMA*((he[:,-1]+(50 if worst else 10))**4-spec['inlet_K']**4)/denom)
        if np.max(abs(new-leak))<.005:leak=new;break
        leak=.5*leak+.5*new
    else:raise RuntimeError('Counterflow/insulation energy balance did not converge')
    ei,eo,cm,om,he,hm,wall=network(leak)
    peaks=[]
    for j in range(s):
        # End-of-slice coolant maxima give a conservative bound on the local
        # Robin boundary instead of reducing the peak with midpoint averaging.
        ti=np.maximum(ei[:,j],ei[:,j+1]);to=np.maximum(eo[:,j],eo[:,j+1])
        peak,_,_=radial(qmod[:,j]/dz,ti,to,ri,ro,km,heff,ho)
        peaks.append(peak)
    peaks=np.array(peaks).T;where=np.unravel_index(peaks.argmax(),peaks.shape)
    mixed=float(np.dot(flows,he[:,-1])/totalflow)
    if abs(mixed-spec['outlet_K'])>.05:raise RuntimeError('Thermal energy closure error')
    qfuel=fu/dz
    rclad=ll['clad']['r_outer_cm']/100;rgap=ll['fuel_gap']['r_outer_cm']/100
    kg=float(props(1200,spec['pressure_Pa'])[1])
    rr=1/(2*math.pi*rclad*hh[:,None])+math.log(rclad/rgap)/(2*math.pi*35.)+math.log(rgap/rf)/(2*math.pi*kg)
    fuelT=hm+qfuel*rr+qfuel/(8*math.pi*(10 if worst else 15))
    field=dict(fuel=fuelT.tolist(),moderator=peaks.tolist(),hot_he=hm.tolist(),cold_he=cm.tolist(),outer_cold_he=om.tolist(),
        clad=float(np.mean(hm+qfuel/(2*math.pi*rclad*hh[:,None]))),hot_wall=float(np.mean(wall)),
        cold_wall=float(np.mean(cm)+5),foil=float(np.mean((wall+cm)/2)),
        can=float(np.mean(peaks)),barrier=float(np.mean(peaks)),reflector=float(np.mean(cm)+10),absorber=float(np.mean(cm)+10))
    def hydraulics(mflow,r1,r2,temp):
        _,_,friction,a,dhh=film(mflow,r1,r2,temp,spec['pressure_Pa'])
        density=spec['pressure_Pa']*.004002602/(8.314462618*temp)
        velocity=mflow/(density*a)
        mach=velocity/np.sqrt((5/3)*8.314462618/.004002602*temp)
        drop=(friction*(h+2*spec['axial_reflector_cm']/100)/dhh+4)*density*velocity**2/2
        return mach,drop
    hotmach,hotdrop=hydraulics(flows,hotri,hotro,np.maximum(he.max(axis=1),850.))
    coldmach,colddrop=hydraulics(flows*fm,coldri,coldro,np.maximum(ei.max(axis=1),850.))
    if outer:
        omach,odrop=hydraulics(flows*(1-fm),or_i,or_o,np.maximum(eo.max(axis=1),850.))
        coldmach=np.maximum(coldmach,omach);colddrop=np.maximum(colddrop,odrop)
    maxmach=float(np.maximum(hotmach,coldmach).max())
    dpfrac=float((hotdrop+colddrop).max()/spec['pressure_Pa'])
    summary=dict(worst=worst,max_moderator_K=float(peaks.max()),pin=int(where[0]),slice=int(where[1]),
        mixed_inlet_K=spec['inlet_K'],mixed_outlet_K=mixed,max_hot_outlet_K=float(he[:,-1].max()),
        cold_heat_W=float(qmod.sum()+leak.sum()),insulation_leak_W=float(leak.sum()),
        cold_return_max_K=float(max(ei[:,0].max(),eo[:,0].max())),max_fuel_K=float(fuelT.max()),
        k_YH2_W_mK=km,emissivity=eps,bridge_conductance_total_W_K=bridge*n*s,
        source_uncertainty='3 sigma per scored cold bin + 5% deposition allowance' if worst else 'mean tally',
        flow_flattening_fraction=flatten,modeled_helium_mass_flow_kg_s=totalflow,
        max_channel_Mach=maxmach,max_pressure_drop_fraction=dpfrac,
        hydraulic_closure_valid=maxmach<.30 and dpfrac<.15,
        field=field)
    return summary,dict(qmod=qmod,cold_edges=ei,outer_edges=eo,hi=heff,ho=ho,ri=ri,ro=ro,h=h,k=km)


def cardinal_check(directory,label,summary,data):
    d=Path(directory)/label;d.mkdir(parents=True,exist_ok=True)
    p=summary['pin'];h=data['h'];edges=np.linspace(0,h,SLICES+1)
    qvol=data['qmod'][p]/(math.pi*(data['ro']**2-data['ri']**2)*(h/SLICES))
    fmt=lambda a:' '.join(f'{float(x):.12g}' for x in a)
    inp=f'''[Mesh]
  coord_type = RZ
  [cell]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 160
    xmin = {data['ri']}
    xmax = {data['ro']}
    ymin = 0
    ymax = {h}
  []
[]
[Variables/T]
  initial_condition = 900
[]
[Functions]
  [source]
    type = PiecewiseConstant
    axis = y
    direction = left
    x = '{fmt(edges)}'
    y = '{fmt(np.r_[qvol,qvol[-1]])}'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '{fmt(edges)}'
    y = '{fmt(data['cold_edges'][p])}'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '{fmt(edges)}'
    y = '{fmt(data['outer_edges'][p])}'
  []
[]
[Kernels]
  [conduction]
    type = HeatConduction
    variable = T
  []
  [heat]
    type = BodyForce
    variable = T
    function = source
  []
[]
[BCs]
  [inner]
    type = ConvectiveFluxFunction
    variable = T
    boundary = left
    coefficient = {data['hi'][p]}
    T_infinity = inner_bulk
  []
'''
    if data['ho'][p]>0:
        inp+=f'''  [outer]
    type = ConvectiveFluxFunction
    variable = T
    boundary = right
    coefficient = {data['ho'][p]}
    T_infinity = outer_bulk
  []
'''
    inp+=f'''[]
[Materials/thermal]
  type = GenericConstantMaterial
  prop_names = thermal_conductivity
  prop_values = {data['k']}
[]
[Postprocessors/Tmax]
  type = NodalExtremeValue
  variable = T
  value_type = max
[]
[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-10
[]
[Outputs]
  csv = true
  exodus = true
[]
'''
    (d/'case.i').write_text(inp)
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    with (d/'run.log').open('w') as log:
        cp=subprocess.run([CARDINAL,'-i','case.i','--n-threads=1'],cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT)
    if cp.returncode:raise RuntimeError(f'Cardinal failure: {d}/run.log')
    with (d/'case_out.csv').open() as f:r=list(csv.DictReader(f))[-1]
    return dict(max_moderator_K=float(r['Tmax']),input=str(d/'case.i'),bound_K=summary['max_moderator_K'])


def evaluate(directory,spec,maps,errors,deposition):
    nominal,nd=scenario(spec,maps,errors,deposition)
    worst,wd=scenario(spec,maps,errors,deposition,worst=True)
    checks=[cardinal_check(directory,'nominal',nominal,nd),cardinal_check(directory,'conservative',worst,wd)]
    upper=max(nominal['max_moderator_K'],worst['max_moderator_K'],*(x['max_moderator_K'] for x in checks))+5.
    return dict(nominal=nominal,conservative=worst,cardinal_checks=checks,
        acceptance_upper_K=upper,pass_thermal=upper<950. and nominal['hydraulic_closure_valid'] and worst['hydraulic_closure_valid'],
        uncertainty_scope='Monte Carlo source envelope + bounded thermal closures + 5 K numerical allowance. Not nuclear-data covariance or full-system uncertainty.',
        limitations=['Prescribed/orificed serial flow, no manifold CFD','Intact vacuum shield only; bridges bounded 12--35 W/K',
                    'Reflector heat routed conservatively to sleeves; no full-core conduction','Containment/material compatibility and structural lifetime unqualified'])
