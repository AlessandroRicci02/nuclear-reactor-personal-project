"""OpenMC photon-heating -> serial-flow network -> Cardinal YH2 RZ solves.

Prescribed material conductivities, emissivity, and helium film closures are
uncertain model inputs. The output explicitly distinguishes simulated heating
from conditional thermal predictions. Not an operating reactor while k<1.
"""
from __future__ import annotations
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import numpy as np
import openmc
from s2_model import ROOT

CAM=ROOT/'campaign';HD=CAM/'heating';TD=CAM/'thermal';TD.mkdir(exist_ok=True)
CARDINAL='/home/aleri/projects/cardinal/cardinal-opt'
P=json.loads((ROOT/'configuration.json').read_text())['parameters']
N=217;S=6;CP=P['helium_cp_J_kgK'];FLOW=P['power_W']/CP/(P['mixed_outlet_target_K']-P['inlet_K'])
MFLOW=FLOW/N;CAP=MFLOW*CP
SIGMA=5.670374419e-8
RI=.00927;RO=.01087;RCOOL=.00897;RMID=.00747;RSI=.0077575;RSO=.0077825;RCOLD=.00807
H=P['active_height_cm']/100;DZ=H/S;HISO=(P['active_height_cm']+2*P['axial_reflector_cm'])/100

def write(path,x):Path(path).write_text(json.dumps(x,indent=2)+'\n')

def get_heating(sp,name):
    t=sp.get_tally(name=name+'_heating');return np.asarray(t.mean).ravel(),np.asarray(t.std_dev).ravel()

def extraction():
    state=max(HD.glob('statepoint.*.h5'),key=lambda x:int(x.name.split('.')[1]))
    with openmc.StatePoint(state) as sp:
        mt=sp.get_tally(name='all_material_heating')
        names=[x.name for x in openmc.Materials.from_xml(HD/'materials.xml')]
        values=np.asarray(mt.mean).ravel()
        total=float(values.sum())
        scale=P['power_W']/total
        materials={name:float(value*scale) for name,value in zip(names,values)}
        arrays={};errors={}
        for name in ('fuel','moderator','hot_He','cold_He','core_Be')+[f'layer_{i}' for i in (3,5,7,9,10,12,13)]:
            x,sd=get_heating(sp,name)
            arrays[name]=(x*scale).reshape((N,S) if name in ('fuel','moderator') else (1,) if name=='core_Be' else (N,))
            errors[name]=(sd*scale).reshape(arrays[name].shape)
        out=dict(statepoint=str(state),keff=float(sp.keff.n),sigma_k=float(sp.keff.s),
            raw_heating_eV_per_source=total,normalization_W_per_eV_per_source=scale,
            material_W=materials,scored_W={k:float(v.sum()) for k,v in arrays.items()},
            moderator_relative_sampling_sigma_max=float(np.max(errors['moderator']/np.maximum(arrays['moderator'],1e-12))))
    np.savez(TD/'heating_W.npz',**arrays,**{'sigma_'+k:v for k,v in errors.items()})
    write(TD/'deposition_summary.json',out)
    return arrays,out

def radiation(Th,Tc,eps):
    d=1/eps+(RMID/RSI)*(1/eps-1)+(RMID/RSO)/eps+(RMID/RCOLD)*(1/eps-1)
    return 2*math.pi*RMID*DZ*SIGMA*(Th**4-Tc**4)/d

def he_props(T):
    mu=3.674e-7*T**.70
    k=2.682e-3*(1+1.123e-8*P['pressure_Pa'])*T**(.71*(1-2e-9*P['pressure_Pa']))
    return mu,k

def h_cold(T=850.,nu=5.):
    mu,k=he_props(T);area=math.pi*(.00897**2-.00837**2);dh=.0012
    re=MFLOW*dh/(area*mu)
    h=nu*k/dh
    # Convective surface at r=8.97 mm, conductive Nb and Al2O3 before YH2.
    resistance=1/(2*math.pi*RCOOL*h)+math.log(.00917/RCOOL)/(2*math.pi*40.)+math.log(RI/.00917)/(2*math.pi*6.)
    heff=1/(2*math.pi*RI*resistance)
    return dict(Re=re,k_He_W_mK=k,h_fluid_W_m2K=h,h_effective_at_YH2_W_m2K=heff,
                viscosity_Pa_s=mu,assumed_Nusselt=nu)

def flow_case(arrays,deposition,eps=.15,bridge_W_K=10.,gas_W_K=2.,be_mode='core',
              kmod=10.,nu=5.,hot_wall_offset_K=0.,max_iter=40):
    a=arrays
    f=a['fuel'];fp=f.sum(axis=1);weights=fp/fp.sum()
    shape=np.divide(f,fp[:,None],out=np.full_like(f,1/S),where=fp[:,None]>0)
    mat=deposition['material_W']
    # Actual OpenMC all-material heating is normalized to 400 kW. Assume all
    # unassigned energy ultimately reaches the hot branch. This keeps total
    # heat exactly conserved but is not a reflector cooling calculation.
    direct=(a['moderator'].copy()
            +shape*(a['layer_7']+a['layer_9']+a['layer_10']+a['layer_12']+
                    a['layer_13']+a['layer_5']+a['cold_He'])[:,None])
    qbe=float(a['core_Be'][0]) if be_mode=='core' else mat['Be_reflector']
    direct+=qbe*weights[:,None]*shape
    qdirect=float(direct.sum())
    if qdirect>P['power_W']:raise RuntimeError('Cold direct load exceeds total deposited power')
    hot_source=(P['power_W']-qdirect)*weights[:,None]*shape
    wallpower=a['layer_3'][:,None]*shape
    leak=np.full((N,S),30_000./N/S)
    coldmid=np.full((N,S),825.)
    hotmid=np.full((N,S),1000.)
    bridge=(bridge_W_K+gas_W_K)/N/S
    # Weakly coupled fixed point for serial counterflow and gray-body exchange.
    for it in range(max_iter):
        cold=direct+leak
        cold_edges=np.zeros((N,S+1));cold_edges[:,S]=P['inlet_K']
        for j in reversed(range(S)):
            cold_edges[:,j]=cold_edges[:,j+1]+cold[:,j]/CAP
        cnew=(cold_edges[:,:-1]+cold_edges[:,1:])/2
        hot=hot_source-leak
        hot_edges=np.zeros((N,S+1));hot_edges[:,0]=cold_edges[:,0]
        for j in range(S):hot_edges[:,j+1]=hot_edges[:,j]+hot[:,j]/CAP
        hnew=(hot_edges[:,:-1]+hot_edges[:,1:])/2
        Th=hnew+wallpower/(2*math.pi*RMID*DZ*1100.)+hot_wall_offset_K
        # He and wall temperatures are not equal. Taking Tc=cold bulk makes
        # the radiative leak larger than with a hotter cold wall.
        newleak=np.maximum(0.,radiation(Th,cnew,eps))+bridge*np.maximum(0.,Th-cnew)
        if np.max(abs(newleak-leak))<.01:
            leak=newleak;break
        leak=.6*leak+.4*newleak
    else:raise RuntimeError('Thermal counterflow iteration did not converge')
    # Recompute consistently at final leak.
    cold=direct+leak
    cold_edges=np.zeros((N,S+1));cold_edges[:,S]=P['inlet_K']
    for j in reversed(range(S)):cold_edges[:,j]=cold_edges[:,j+1]+cold[:,j]/CAP
    cnew=(cold_edges[:,:-1]+cold_edges[:,1:])/2
    hot=hot_source-leak
    hot_edges=np.zeros((N,S+1));hot_edges[:,0]=cold_edges[:,0]
    for j in range(S):hot_edges[:,j+1]=hot_edges[:,j]+hot[:,j]/CAP
    hnew=(hot_edges[:,:-1]+hot_edges[:,1:])/2
    Th=hnew+wallpower/(2*math.pi*RMID*DZ*1100.)+hot_wall_offset_K
    qmod=direct # Strongly conservative: all direct cold-side and allocated Be
    # deposited in YH2. Radiation/bridges heat cold He only.
    hinfo=h_cold(float(np.max(cnew)),nu)
    h_eff=hinfo['h_effective_at_YH2_W_m2K']
    yvolume=math.pi*(RO**2-RI**2)*DZ
    qvol=qmod/yvolume
    # Exact radial steady solution for uniform heat in each axial slice,
    # adiabatic outside and inner convection. Axial conduction would flatten
    # local peaks, so this is a per-slice peak upper screen under given sources.
    # qvol*(...) has K. This expression is intentionally evaluated directly.
    conduction=qvol*(2*RO**2*math.log(RO/RI)-(RO**2-RI**2))/(4*kmod)
    film=qmod/(2*math.pi*RI*DZ*h_eff)
    peak=cnew+film+conduction
    idx=np.unravel_index(np.argmax(peak),peak.shape)
    mixed=float(np.mean(hot_edges[:,S])) # equal per-pin mass flows
    summary=dict(heating_case=be_mode,eps=eps,bridge_W_K=bridge_W_K,gas_W_K=gas_W_K,
                 hot_wall_offset_K=hot_wall_offset_K,
                 k_YH2_W_mK=kmod,film=hinfo,cold_direct_W=qdirect,
                 cold_leak_W=float(leak.sum()),cold_total_W=float(cold.sum()),
                 cold_mixed_outlet_K=float(np.mean(cold_edges[:,0])),
                 hot_mixed_outlet_K=mixed,max_hot_outlet_K=float(hot_edges[:,S].max()),
                 min_hot_outlet_K=float(hot_edges[:,S].min()),
                 max_hot_wall_estimate_K=float(Th.max()),
                 max_cold_bulk_K=float(cnew.max()),
                 max_YH2_analytic_K=float(peak.max()),
                 max_YH2_pin_id=int(idx[0]+1),max_YH2_slice=int(idx[1]),
                 max_YH2_film_rise_K=float(film[idx]),
                 max_YH2_conduction_rise_K=float(conduction[idx]),
                 total_heat_balance_error_W=float((hot_edges[:,S].mean()-P['inlet_K'])*FLOW*CP-P['power_W']),
                 radiation_model='one Mo shield gray-cylinder formula with prescribed emissivity; no vacuum-loss event',
                 worst_pin_axial_bulk_edges_K=cold_edges[idx[0]].tolist(),
                 worst_pin_direct_W_by_slice=direct[idx[0]].tolist(),
                 worst_pin_total_cold_W_by_slice=cold[idx[0]].tolist(),
                 model_limitations=['400 kW imposed heat even though BOL k<1','All unassigned material heat routed to hot He',
                   'Be core heat allocated by fuel power; no full-core conduction',
                   'Uniform per-pin mass flow; no manifold/orifice calculation',
                   'Film coefficient uses assumed Nusselt; no CFD',
                   'Radiation uses coolant-midpoint temperatures as wall estimate',
                   'Vacuum integrity, contacts, gamma tally systematics and material properties unqualified'])
    np.savez(TD/f'fields_{be_mode}_{eps}_{kmod}_{nu}_{bridge_W_K}_{hot_wall_offset_K}.npz',cold_bulk=cnew,hot_bulk=hnew,
             cold_edges=cold_edges,hot_edges=hot_edges,moderator_source_W=qmod,
             moderator_peak_slice_K=peak,cold_leak_W=leak)
    return summary,(idx, qmod[idx[0]],cold_edges[idx[0]])

def cardinal_pin(label,qby,cedges,h_eff,kmod=10.,nx=20,ny=120):
    d=TD/'cardinal'/label;d.mkdir(parents=True,exist_ok=True)
    edges=np.linspace(0,H,S+1)
    # z is bottom-to-top; cold fluid flows top-to-bottom.
    vals=lambda a:' '.join(f'{v:.13g}' for v in a)
    qvol=qby/(math.pi*(RO**2-RI**2)*DZ)
    inp=f'''[Mesh]
  coord_type = RZ
  [moderator]
    type = GeneratedMeshGenerator
    dim = 2
    nx = {nx}
    ny = {ny}
    xmin = {RI}
    xmax = {RO}
    ymin = 0
    ymax = {H}
  []
[]
[Variables/T]
  initial_condition = 850
[]
[Functions]
  [heat]
    type = PiecewiseConstant
    axis = y
    direction = left
    x = '{vals(edges)}'
    y = '{vals(np.r_[qvol,qvol[-1]])}'
  []
  [bulk]
    type = PiecewiseLinear
    axis = y
    x = '{vals(edges)}'
    y = '{vals(cedges)}'
  []
[]
[Kernels]
  [conduction]
    type = HeatConduction
    variable = T
  []
  [source]
    type = BodyForce
    variable = T
    function = heat
  []
[]
[BCs]
  [inner]
    type = ConvectiveFluxFunction
    variable = T
    boundary = left
    coefficient = {h_eff}
    T_infinity = bulk
  []
[]
[Materials]
  [conductivity]
    type = GenericConstantMaterial
    prop_names = thermal_conductivity
    prop_values = {kmod}
  []
  [bulk_temperature]
    type = GenericFunctionMaterial
    prop_names = bulk_temperature
    prop_values = bulk
  []
[]
[Postprocessors]
  [T_YH2_max_K]
    type = NodalExtremeValue
    variable = T
    value_type = max
  []
  [Q_removed_W]
    type = ConvectiveHeatTransferSideIntegral
    T_solid = T
    T_fluid = bulk_temperature
    htc = {h_eff}
    boundary = left
  []
[]
[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-9
[]
[Outputs]
  csv = true
[]
'''
    (d/'moderator.i').write_text(inp)
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    with (d/'run.log').open('w') as log:
        p=subprocess.run([CARDINAL,'-i','moderator.i','--n-threads=1'],cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT)
    if p.returncode:raise RuntimeError(f'Cardinal failed: {d}/run.log')
    with (d/'moderator_out.csv').open() as f:r={k:float(v) for k,v in list(csv.DictReader(f))[-1].items()}
    r.update(expected_heat_W=float(sum(qby)),relative_energy_error=(r['Q_removed_W']-float(sum(qby)))/float(sum(qby)),
             k_YH2_W_mK=kmod,nx=nx,ny=ny)
    if abs(r['relative_energy_error'])>.02:raise RuntimeError(f'Cardinal energy balance failed: {r}')
    write(d/'result.json',r);return r

def main():
    arrays,dep=extraction()
    cases=[]
    for label,opts in [('baseline',{}),('oxidized',dict(eps=.3)),('hotter_wall',dict(hot_wall_offset_K=100.)),
                       ('high_bridge',dict(bridge_W_K=40)),('low_film',dict(nu=3.5)),
                       ('low_k_YH2',dict(kmod=5.)),('all_Be_to_cold',dict(be_mode='all'))]:
        r,details=flow_case(arrays,dep,**opts)
        r['label']=label;cases.append(r)
        print(label,r['max_YH2_analytic_K'],r['hot_mixed_outlet_K'],flush=True)
        if label in ('baseline','oxidized','all_Be_to_cold'):
            idx,qby,cedges=details
            r['cardinal_peak_pin']=cardinal_pin(label,qby,cedges,r['film']['h_effective_at_YH2_W_m2K'])
            if label=='baseline':
                r['cardinal_peak_pin_refined']=cardinal_pin('baseline_refined',qby,cedges,r['film']['h_effective_at_YH2_W_m2K'],nx=40,ny=240)
    write(TD/'thermal_results.json',dict(status='complete',cases=cases,
          disclaimer='Conditional temperature at imposed 400 kW; candidate is subcritical at all-drums-out BOL.'))

if __name__=='__main__':main()
