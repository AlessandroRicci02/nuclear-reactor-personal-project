"""OpenMC heating -> Cardinal/MOOSE cell FEM -> axial helium balances.

Reduced thermal-hydraulic model, NOT a native 3-D conjugate CFD calculation.
Actual elliptical passages and Al2O3 liner are meshed in XY. Linear material
properties permit exact superposition of three FEM response fields per zone.
Axial fluid balances include heat exchanged from the hot annulus through the
liner/moderator to auxiliary coolant. Full-core radial conduction, manifolds,
and reflector cooling remain outside this reduced model.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import numpy as np
from scipy.spatial import Delaunay
from scipy.optimize import brentq
import openmc
from p1_model import CAMPAIGN, CARDINAL, load_config

CP=5193.16
TIN=892.
PRESSURE=2.5e6
K_MOD=10. # W/m/K, explicit constant-property assumption, sensitivity below
K_LINER=6.


def write_json(path,value):
    Path(path).write_text(json.dumps(value,indent=2)+'\n')


def polygon_contains(xy,poly):
    # Convex outer hexagon, counterclockwise vertices.
    good=np.ones(len(xy),bool)
    for a,b in zip(poly,np.roll(poly,-1,axis=0)):
        good &= cross2(b-a,xy-a)>=-1e-13
    return good


def cross2(a,b):
    return a[...,0]*b[...,1]-a[...,1]*b[...,0]


def cell_mesh(path,cfg,zone,spacing=.0004):
    """Recover every geometric edge by Delaunay midpoint insertion; audit areas."""
    z=cfg['zones'][zone];aux=cfg['auxiliary_channels'];pitch=cfg['geometry']['pitch_cm']*.01
    rn=z['annulus_outer_radius_cm']*.01;rl=z['liner_outer_radius_cm']*.01
    angle=np.arange(6)*math.pi/3
    hexagon=pitch/math.sqrt(3)*np.column_stack((np.cos(angle),np.sin(angle)))
    pts=[];segments=[];curves=[]
    def curve(poly,tag):
        offset=len(pts);pts.extend(poly.tolist());curves.append(np.asarray(poly))
        for i in range(len(poly)):segments.append((offset+i,offset+(i+1)%len(poly),tag))
    hexpts=[]
    for a,b in zip(hexagon,np.roll(hexagon,-1,axis=0)):
        for f in np.arange(math.ceil(np.linalg.norm(b-a)/spacing))/math.ceil(np.linalg.norm(b-a)/spacing):hexpts.append(a+f*(b-a))
    curve(np.array(hexpts),12)
    n=max(int(120*.0004/spacing),int(math.ceil(2*math.pi*rl/spacing))*2)
    th=np.arange(n)*2*math.pi/n
    curve(rn*np.column_stack((np.cos(th),np.sin(th))),10)
    curve(rl*np.column_stack((np.cos(th),np.sin(th))),0)
    holes=[]
    for j in range(aux['count']):
        phi=math.radians(aux['phase_deg'])+j*2*math.pi/aux['count'];c,s=math.cos(phi),math.sin(phi)
        rot=np.array([[c,-s],[s,c]])
        center=aux['centre_radius_cm']*.01*np.array([c,s])
        a,b=aux['a_cm']*.01,aux['b_cm']*.01
        ne=max(72,int(72*.0004/spacing));t=np.arange(ne)*2*math.pi/ne
        poly=np.column_stack((a*np.cos(t),b*np.sin(t)))@rot.T+center
        curve(poly,11);holes.append((center,rot,a,b))
    def inside(xy):
        ok=polygon_contains(xy,hexagon)&(np.linalg.norm(xy,axis=1)>=rn-1e-12)
        for center,rot,a,b in holes:
            p=(xy-center)@rot;ok &= (p[:,0]/a)**2+(p[:,1]/b)**2>=1-1e-10
        return ok
    v=np.arange(-pitch,pitch+spacing/2,spacing)
    x,y=np.meshgrid(v,v);grid=np.column_stack((x.ravel(),y.ravel()))
    # Do not put interior points in the very thin liner ring.
    grid=grid[inside(grid)&(np.linalg.norm(grid,axis=1)>rl+spacing*.15)]
    # Keep grid nodes away from constrained segments, especially exact matches
    # on horizontal hex faces at the refined spacing.
    keep_grid=np.ones(len(grid),bool)
    for va,vb in zip(hexagon,np.roll(hexagon,-1,axis=0)):
        keep_grid &= cross2(vb-va,grid-va)/np.linalg.norm(vb-va)>spacing*.1
    for center,rot,aa,bb in holes:
        local=(grid-center)@rot
        keep_grid &= (local[:,0]/aa)**2+(local[:,1]/bb)**2>1+spacing*.1/min(aa,bb)
    grid=grid[keep_grid]
    pts.extend(grid.tolist())
    for iteration in range(20):
        xy=np.asarray(pts);dt=Delaunay(xy)
        edges={tuple(sorted((int(a),int(b)))) for tri in dt.simplices for a,b in zip(tri,np.roll(tri,-1))}
        missing=[s for s in segments if tuple(sorted(s[:2])) not in edges]
        if not missing:break
        bad=set(missing);new=[]
        for a,b,tag in segments:
            if (a,b,tag) in bad:
                mid=len(pts);pts.append(((xy[a]+xy[b])/2).tolist());new.extend([(a,mid,tag),(mid,b,tag)])
            else:new.append((a,b,tag))
        segments=new
    else:raise RuntimeError('Failed to recover thermal mesh boundary edges')
    tris=dt.simplices;cent=xy[tris].mean(axis=1)
    # Use polygonal boundaries for domain membership to retain boundary slivers.
    from matplotlib.path import Path as PolyPath
    keep=PolyPath(curves[0]).contains_points(cent,radius=1e-12)
    for poly in [curves[1]]+curves[3:]:keep &= ~PolyPath(poly).contains_points(cent,radius=-1e-13)
    tris=tris[keep];cent=cent[keep]
    blocks=np.where(np.linalg.norm(cent,axis=1)<rl,2,1)
    areas=np.abs(cross2(xy[tris[:,1]]-xy[tris[:,0]],xy[tris[:,2]]-xy[tris[:,0]]))/2
    boundary=[s for s in segments if s[2]]
    used=set(tris.ravel().tolist());assert all(a in used and b in used for a,b,_ in boundary)
    edge_count={}
    for tri in tris:
        for a,b in zip(tri,np.roll(tri,-1)):
            key=tuple(sorted((int(a),int(b))));edge_count[key]=edge_count.get(key,0)+1
    assert {e for e,n in edge_count.items() if n==1}=={tuple(sorted(s[:2])) for s in boundary}
    with Path(path).open('w') as f:
        f.write('$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n5\n2 1 "moderator"\n2 2 "liner"\n1 10 "main"\n1 11 "aux"\n1 12 "outer"\n$EndPhysicalNames\n')
        f.write(f'$Nodes\n{len(xy)}\n')
        for i,(x,y) in enumerate(xy,1):f.write(f'{i} {x:.15g} {y:.15g} 0\n')
        f.write(f'$EndNodes\n$Elements\n{len(boundary)+len(tris)}\n')
        eid=0
        for a,b,tag in boundary:eid+=1;f.write(f'{eid} 1 2 {tag} {tag} {a+1} {b+1}\n')
        for tri,block in zip(tris,blocks):eid+=1;f.write(f'{eid} 2 2 {block} {block} '+ ' '.join(str(n+1) for n in tri)+'\n')
        f.write('$EndElements\n')
    exact_mod=math.sqrt(3)/2*pitch**2-math.pi*rl**2-aux['count']*math.pi*aux['a_cm']*aux['b_cm']*1e-4
    exact_lin=math.pi*(rl**2-rn**2)
    report=dict(nodes=len(xy),elements=len(tris),moderator_area_m2=float(areas[blocks==1].sum()),liner_area_m2=float(areas[blocks==2].sum()),exact_moderator_area_m2=exact_mod,exact_liner_area_m2=exact_lin,spacing_m=spacing,all_boundary_edges_verified=True)
    assert abs(report['moderator_area_m2']/exact_mod-1)<.005
    assert abs(report['liner_area_m2']/exact_lin-1)<.005
    write_json(Path(path).with_suffix('.json'),report)
    return report


def helium(T):
    mu=3.674e-7*T**.70
    k=2.682e-3*(1+1.123e-8*PRESSURE)*T**(.71*(1-2e-9*PRESSURE))
    rho=PRESSURE*.004002602/(8.314462618*T)
    return mu,k,rho


def channel_coefficients(cfg,z,mdot_main,mdot_aux):
    mu,k,rho=helium(1001.)
    ri=z['clad_outer_radius_cm']*.01;ro=z['annulus_outer_radius_cm']*.01
    area=math.pi*(ro*ro-ri*ri);dh=2*(ro-ri)
    aux=cfg['auxiliary_channels'];aa=aux['count']*math.pi*aux['a_cm']*aux['b_cm']*1e-4
    da=4*aa/(aux['perimeter_cm']*.01)
    re_m=mdot_main*dh/(area*mu);re_a=mdot_aux*da/(aa*mu)
    def nu(re,lam):
        if re<=2300:return lam
        pr=CP*mu/k;f=(.79*math.log(max(re,3000))-1.64)**-2
        turbulent=(f/8)*(max(re,3000)-1000)*pr/(1+12.7*math.sqrt(f/8)*(pr**(2/3)-1))
        w=np.clip((re-2300)/1700,0,1)
        return (1-w)*lam+w*max(lam,turbulent)
    return dict(h_main=nu(re_m,5.385)*k/dh,h_aux=nu(re_a,3.66)*k/da,Re_main=re_m,Re_aux=re_a,main_area=area,aux_area=aa,dh_main=dh,dh_aux=da)


def execute_cardinal(directory,input_name):
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    command=[CARDINAL,'-i',input_name,'--n-threads=1']
    with (directory/(input_name+'.log')).open('w') as f:
        p=subprocess.run(command,cwd=directory,env=env,stdout=f,stderr=subprocess.STDOUT)
    if p.returncode:raise RuntimeError(f'Cardinal exit {p.returncode}: {directory/input_name}.log')
    with (directory/(Path(input_name).stem+'_out.csv')).open() as f:row=list(csv.DictReader(f))[-1]
    return {k:float(v) for k,v in row.items()}


def basis(directory,cfg,zone,mdot_main,mdot_aux,spacing=.0004,kmod=K_MOD,hmult=1.):
    directory.mkdir(parents=True,exist_ok=True)
    meta=cell_mesh(directory/'cell.msh',cfg,zone,spacing)
    c=channel_coefficients(cfg,cfg['zones'][zone],mdot_main,mdot_aux)
    hm,ha=c['h_main']*hmult,c['h_aux']*hmult
    responses=[];nodal=[]
    for label,qm,ql,tm in [('mod',1.,0.,0.),('liner',0.,1.,0.),('exchange',0.,0.,1.)]:
        text=f'''[Mesh]
  type = FileMesh
  file = cell.msh
[]
[Variables/T]
  initial_condition = 0
[]
[Kernels]
  [diffusion]
    type = HeatConduction
    variable = T
  []
  [mod_source]
    type = BodyForce
    variable = T
    value = {qm/meta['moderator_area_m2']:.16g}
    block = '1'
  []
  [liner_source]
    type = BodyForce
    variable = T
    value = {ql/meta['liner_area_m2']:.16g}
    block = '2'
  []
[]
[Materials]
  [moderator]
    type = GenericConstantMaterial
    prop_names = thermal_conductivity
    prop_values = {kmod}
    block = 1
  []
  [liner]
    type = GenericConstantMaterial
    prop_names = thermal_conductivity
    prop_values = {K_LINER}
    block = 2
  []
[]
[BCs]
  [main]
    type = ConvectiveFluxFunction
    variable = T
    boundary = 10
    coefficient = {hm:.16g}
    T_infinity = {tm}
  []
  [aux]
    type = ConvectiveFluxFunction
    variable = T
    boundary = 11
    coefficient = {ha:.16g}
    T_infinity = 0
  []
[]
[Postprocessors]
  [Q_main]
    type = ConvectiveHeatTransferSideIntegral
    T_solid = T
    boundary = 10
    T_fluid = {tm}
    htc = {hm:.16g}
  []
  [Q_aux]
    type = ConvectiveHeatTransferSideIntegral
    T_solid = T
    boundary = 11
    T_fluid = 0
    htc = {ha:.16g}
  []
  [T_mod_avg]
    type = ElementAverageValue
    variable = T
    block = 1
  []
  [T_liner_avg]
    type = ElementAverageValue
    variable = T
    block = 2
  []
[]
[VectorPostprocessors/nodes]
  type = NodalValueSampler
  variable = T
  sort_by = id
[]
[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10
  l_tol = 1e-10
[]
[Outputs]
  csv = true
[]
'''
        (directory/f'{label}.i').write_text(text)
        r=execute_cardinal(directory,f'{label}.i')
        residual=r['Q_main']+r['Q_aux']-qm-ql
        if abs(residual)>1e-6:raise RuntimeError(f'FEM heat balance failed {residual}')
        samples=sorted(directory.glob(f'{label}_out_nodes_*.csv'))
        data=np.genfromtxt(samples[-1],delimiter=',',names=True)
        nodal.append(data['T']);responses.append(r)
    write_json(directory/'basis.json',dict(coefficients=c,responses=responses,kmod=kmod,k_liner=K_LINER,htc_multiplier=hmult,mesh=meta))
    return dict(coefficients=c,responses=responses,nodal=np.asarray(nodal),kmod=kmod,htc_multiplier=hmult)


def heating(statepoint,directory):
    with openmc.StatePoint(statepoint) as sp:
        material=sp.get_tally(name='material_heating')
        total=float(material.mean.sum())
        if total<=0:raise RuntimeError('Nonpositive heating tally')
        scale=400000./total
        fields={name:sp.get_tally(name=name+'_heating').mean.ravel().reshape(217,10)*scale for name in ('fuel','clad','liner','moderator')}
        errors={name:sp.get_tally(name=name+'_heating').std_dev.ravel().reshape(217,10)*scale for name in fields}
        material_ids=material.filters[0].bins
        matpower={int(mid):float(q)*scale for mid,q in zip(material_ids,material.mean.ravel())}
        peripheral=max(0.,matpower[3]-fields['moderator'].sum())
        summary=dict(statepoint=str(statepoint),normalization='Sum of neutron+photon deposited material heating scaled to 400 kW; not identical to recoverable-fission depletion normalization',raw_deposited_eV_per_source=total,scale_W_per_eV=scale,pin_block_heat_W={k:float(v.sum()) for k,v in fields.items()},peripheral_moderator_heat_W=peripheral,unmodeled_reflector_absorber_coolant_heat_W=400000-sum(float(v.sum()) for v in fields.values())-peripheral,max_cell_relative_sampling_sigma={k:float(np.max(errors[k]/np.maximum(v,1e-20))) for k,v in fields.items()})
    directory.mkdir(parents=True,exist_ok=True)
    write_json(directory/'heating_summary.json',summary)
    np.savez(directory/'heating_arrays.npz',**fields,**{'sigma_'+k:v for k,v in errors.items()})
    return fields,summary


def solve_flows(cfg,pins,fields,bases,peripheral_W):
    H=cfg['parameters']['height_cm']*.01;dz=H/10
    total_m=cfg['boundaries']['total_mass_flow_kg_s'];fraction=cfg['flow']['auxiliary_fraction']
    weights=np.array([cfg['zones'][p['zone']]['provisional_main_flow_per_pin_relative_to_mean'] for p in pins])
    mm=total_m*(1-fraction)*weights/weights.sum();ma=np.full(217,total_m*fraction/217)
    qmod=fields['moderator'].copy()
    # Peripheral heat collected by the outermost pin ring, explicit network closure.
    outer=[p['pin_id'] for p in pins if p['ring']==8]
    qmod[outer,:]+=peripheral_W/(len(outer)*10)
    records=[];grid={k:np.zeros((217,10)) for k in ('fuel','clad','liner','moderator','coolant')}
    for p in pins:
        pi=p['pin_id'];z=cfg['zones'][p['zone']];b=bases[p['zone']];c=b['coefficients'];rr=b['responses']
        tmi=tai=TIN;dp_m=dp_a=0.
        for s in range(10):
            pm=qmod[pi,s]/dz;pl=fields['liner'][pi,s]/dz;pf=(fields['fuel'][pi,s]+fields['clad'][pi,s])/dz
            # Q_main=a*Pmod+b*Pliner+g*(Tm-Ta), W per metre.
            A=rr[0]['Q_main']*pm+rr[1]['Q_main']*pl;g=rr[2]['Q_main']
            cm=mm[pi]*CP;ca=ma[pi]*CP
            mat=np.array([[1-dz*g/(2*cm),dz*g/(2*cm)],[dz*g/(2*ca),1-dz*g/(2*ca)]])
            rhs=np.array([tmi+dz*(pf+A)/(2*cm),tai+dz*(pm+pl-A)/(2*ca)])
            tm,ta=np.linalg.solve(mat,rhs)
            q_main=A+g*(tm-ta);q_aux=pm+pl-q_main
            tmo=tmi+dz*(pf+q_main)/cm;tao=tai+dz*q_aux/ca
            Tnodes=ta+b['nodal'][0]*pm+b['nodal'][1]*pl+b['nodal'][2]*(tm-ta)
            rf=z['fuel_radius_cm']*.01;rc=z['clad_outer_radius_cm']*.01
            tc=tm+pf/(2*math.pi*rc*c['h_main']*b['htc_multiplier'])
            tf_s=tc+pf*math.log(rc/rf)/(2*math.pi*(41.8+.0125*tc))
            tf_peak=(tf_s**1.41+1.41/1.37*(fields['fuel'][pi,s]/dz)/(4*math.pi))**(1/1.41)
            grid['fuel'][pi,s]=(tf_s+tf_peak)/2;grid['clad'][pi,s]=(tc+tf_s)/2
            for name in ('moderator','liner'):
                key='T_mod_avg' if name=='moderator' else 'T_liner_avg'
                grid[name][pi,s]=ta+rr[0][key]*pm+rr[1][key]*pl+rr[2][key]*(tm-ta)
            grid['coolant'][pi,s]=(tm*c['main_area']+ta*c['aux_area'])/(c['main_area']+c['aux_area'])
            # Fully developed laminar exact shape conductance, turbulent Darcy blend.
            def pressure_drop(mdot,T,area,dh,kind):
                mu,_,rho=helium(T);re=mdot*dh/(area*mu);G=mdot/area
                if kind=='main':
                    ri=z['clad_outer_radius_cm']*.01;ro=z['annulus_outer_radius_cm']*.01
                    den=ro**4-ri**4-(ro**2-ri**2)**2/math.log(ro/ri)
                    lam=8*mu*dz*mdot/(rho*math.pi*den)
                else:
                    ax=cfg['auxiliary_channels'];aa=ax['a_cm']*.01;bb=ax['b_cm']*.01
                    lam=4*mu*dz*(aa*aa+bb*bb)*mdot/(rho*ax['count']*math.pi*aa**3*bb**3)
                turb=(.79*math.log(max(re,3000))-1.64)**-2*dz/dh*G*G/(2*rho)
                w=float(np.clip((re-2300)/1700,0,1));return (1-w)*lam+w*turb
            dp_m+=pressure_drop(mm[pi],tm,c['main_area'],c['dh_main'],'main')
            dp_a+=pressure_drop(ma[pi],ta,c['aux_area'],c['dh_aux'],'aux')
            records.append(dict(pin_id=pi,zone=p['zone'],ring=p['ring'],slice=s,T_main_bulk_K=tm,T_aux_bulk_K=ta,T_main_out_K=tmo,T_aux_out_K=tao,T_fuel_peak_1D_K=tf_peak,T_clad_outer_K=tc,T_cell_solid_peak_K=float(Tnodes.max()),T_mod_avg_K=grid['moderator'][pi,s],Q_to_main_W=q_main*dz,Q_to_aux_W=q_aux*dz,Q_fuel_clad_W=pf*dz,main_friction_cumulative_Pa=dp_m,aux_friction_cumulative_Pa=dp_a,mdot_main_kg_s=mm[pi],mdot_aux_kg_s=ma[pi]))
            tmi,tai=tmo,tao
    heat_removed=sum((r['mdot_main_kg_s']*CP*(r['T_main_out_K']-TIN)+r['mdot_aux_kg_s']*CP*(r['T_aux_out_K']-TIN)) for r in records if r['slice']==9)
    expected=sum(float(v.sum()) for v in fields.values())+peripheral_W
    assert abs(heat_removed-expected)<1e-6*expected
    out=dict(peak_fuel_1D_K=max(r['T_fuel_peak_1D_K'] for r in records),peak_clad_outer_1D_K=max(r['T_clad_outer_K'] for r in records),peak_pin_cell_solid_K=max(r['T_cell_solid_peak_K'] for r in records),max_main_friction_Pa=max(r['main_friction_cumulative_Pa'] for r in records),max_aux_friction_Pa=max(r['aux_friction_cumulative_Pa'] for r in records),heat_removed_W=heat_removed,energy_residual_W=heat_removed-expected,mixed_active_coolant_outlet_K=TIN+heat_removed/(total_m*CP))
    return records,grid,out


def pin_fem(directory,cfg,records,fields,pin_id,refinement=1):
    """Actual RZ nonlinear Cardinal solve, with cell-tally heat and solved coolant."""
    directory.mkdir(parents=True,exist_ok=True)
    rows=[r for r in records if r['pin_id']==pin_id];z=cfg['zones'][rows[0]['zone']]
    H=cfg['parameters']['height_cm']*.01;rf=z['fuel_radius_cm']*.01;rc=z['clad_outer_radius_cm']*.01
    edges=np.linspace(0,H,11);mid=(edges[:-1]+edges[1:])/2
    # Piecewise constant exactly preserves slice power and aligns with mesh planes.
    q=fields['fuel'][pin_id]/(math.pi*rf*rf*H/10)
    qc=fields['clad'][pin_id]/(math.pi*(rc*rc-rf*rf)*H/10)
    main=np.r_[TIN,[r['T_main_out_K'] for r in rows]]
    c=channel_coefficients(cfg,z,rows[0]['mdot_main_kg_s'],rows[0]['mdot_aux_kg_s'])
    def values(x):return ' '.join(f'{v:.14g}' for v in x)
    inp=f'''[Mesh]
  coord_type = RZ
  [fuel]
    type = GeneratedMeshGenerator
    dim = 2
    nx = {24*refinement}
    ny = {60*refinement}
    xmin = 0
    xmax = {rf}
    ymin = 0
    ymax = {H}
    boundary_name_prefix = fuel
  []
  [fuel_id]
    type = SubdomainIDGenerator
    input = fuel
    subdomain_id = 1
  []
  [clad]
    type = GeneratedMeshGenerator
    dim = 2
    nx = {8*refinement}
    ny = {60*refinement}
    xmin = {rf}
    xmax = {rc}
    ymin = 0
    ymax = {H}
    boundary_name_prefix = clad
  []
  [clad_id]
    type = SubdomainIDGenerator
    input = clad
    subdomain_id = 2
  []
  [stitch]
    type = StitchMeshGenerator
    inputs = 'fuel_id clad_id'
    stitch_boundaries_pairs = 'fuel_right clad_left'
  []
[]
[Variables/T]
  initial_condition = 1000
[]
[Functions]
  [qfuel]
    type = PiecewiseConstant
    axis = y
    direction = left
    x = '{values(edges)}'
    y = '{values(np.r_[q,q[-1]])}'
  []
  [qclad]
    type = PiecewiseConstant
    axis = y
    direction = left
    x = '{values(edges)}'
    y = '{values(np.r_[qc,qc[-1]])}'
  []
  [bulk]
    type = PiecewiseLinear
    axis = y
    x = '{values(edges)}'
    y = '{values(main)}'
  []
[]
[Kernels]
  [conduction]
    type = HeatConduction
    variable = T
  []
  [fuel_source]
    type = BodyForce
    variable = T
    block = 1
    function = qfuel
  []
  [clad_source]
    type = BodyForce
    variable = T
    block = 2
    function = qclad
  []
[]
[BCs/convection]
  type = ConvectiveFluxFunction
  variable = T
  boundary = clad_right
  coefficient = {c['h_main']}
  T_infinity = bulk
[]
[Materials]
  [bulk_temperature]
    type = GenericFunctionMaterial
    prop_names = bulk_temperature
    prop_values = bulk
    block = 2
  []
  [fuel]
    type = ParsedMaterial
    property_name = thermal_conductivity
    coupled_variables = T
    expression = '1.37*T^0.41'
    block = 1
  []
  [clad]
    type = ParsedMaterial
    property_name = thermal_conductivity
    coupled_variables = T
    expression = '41.8+0.0125*T'
    block = 2
  []
[]
[Postprocessors]
  [T_fuel_max]
    type = NodalExtremeValue
    variable = T
    value_type = max
    block = 1
  []
  [T_clad_max]
    type = NodalExtremeValue
    variable = T
    value_type = max
    block = 2
  []
  [Q_removed]
    type = ConvectiveHeatTransferSideIntegral
    T_solid = T
    T_fluid = bulk_temperature
    htc = {c['h_main']}
    boundary = clad_right
  []
[]
[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-8
[]
[Outputs]
  csv = true
[]
'''
    (directory/'pin.i').write_text(inp)
    r=execute_cardinal(directory,'pin.i');expected=float(fields['fuel'][pin_id].sum()+fields['clad'][pin_id].sum())
    r.update(pin_id=pin_id,expected_power_W=expected,relative_energy_residual=(r['Q_removed']-expected)/expected)
    if abs(r['relative_energy_residual'])>.005:raise RuntimeError(f'Pin FEM energy mismatch: {r}')
    write_json(directory/'result.json',r);return r


def run(statepoint,reference=False,validation=False):
    cfg,pins=load_config(reference);name='tb35' if reference else 'p1'
    d=CAMPAIGN/('validation/thermal_synthetic' if validation else 'thermal')/name;d.mkdir(parents=True,exist_ok=True)
    if validation:
        # This case tests numerical plumbing only; it is never a reactor result.
        fields={k:np.full((217,10),v) for k,v in dict(fuel=170.,clad=1.,liner=.2,moderator=3.).items()};heat={'peripheral_moderator_heat_W':1000.}
    else:fields,heat=heating(statepoint,d)
    flow=cfg['boundaries']['total_mass_flow_kg_s'];frac=cfg['flow']['auxiliary_fraction'];bases=[]
    norm=sum(z['pins']*z['provisional_main_flow_per_pin_relative_to_mean'] for z in cfg['zones'])
    for z in cfg['zones']:
        bases.append(basis(d/f"zone_{z['zone']}",cfg,z['zone'],flow*(1-frac)*z['provisional_main_flow_per_pin_relative_to_mean']/norm,flow*frac/217))
    records,grid,result=solve_flows(cfg,pins,fields,bases,heat['peripheral_moderator_heat_W'])
    with (d/'channels.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    write_json(d/'temperature_grid.json',{k:v.tolist() for k,v in grid.items()})
    peak=max(records,key=lambda r:r['T_fuel_peak_1D_K'])['pin_id']
    selected=sorted(set([peak,0]+[next(p['pin_id'] for p in pins if p['zone']==z) for z in range(3)]))
    result['pin_fem']=[pin_fem(d/f'pin_{pi}',cfg,records,fields,pi) for pi in selected]
    refined=pin_fem(d/f'pin_{peak}_refined',cfg,records,fields,peak,2)
    coarse=next(r for r in result['pin_fem'] if r['pin_id']==peak)
    result['pin_mesh_change_peak_K']=refined['T_fuel_max']-coarse['T_fuel_max']
    if abs(result['pin_mesh_change_peak_K'])>2.:raise RuntimeError('Pin thermal mesh refinement changed peak by >2 K')
    worst=max(records,key=lambda r:r['T_cell_solid_peak_K'])['zone']
    z=cfg['zones'][worst];mm=flow*(1-frac)*z['provisional_main_flow_per_pin_relative_to_mean']/norm;ma=flow*frac/217
    fine_basis=basis(d/f'zone_{worst}_refined',cfg,worst,mm,ma,spacing=.0002)
    fine_bases=list(bases);fine_bases[worst]=fine_basis
    _,_,fine_summary=solve_flows(cfg,pins,fields,fine_bases,heat['peripheral_moderator_heat_W'])
    result['cell_mesh_change_peak_K']=fine_summary['peak_pin_cell_solid_K']-result['peak_pin_cell_solid_K']
    if abs(result['cell_mesh_change_peak_K'])>2.:raise RuntimeError('Cell thermal mesh refinement changed peak by >2 K')
    if not validation:
        result['property_sensitivities']=[]
        for label,km,hm in [('kmod_5',5.,1.),('kmod_20',20.,1.),('htc_70pct',10.,.7),('htc_130pct',10.,1.3)]:
            varied=[]
            for z in cfg['zones']:
                varied.append(basis(d/f'sensitivity_{label}'/f"zone_{z['zone']}",cfg,z['zone'],flow*(1-frac)*z['provisional_main_flow_per_pin_relative_to_mean']/norm,ma,kmod=km,hmult=hm))
            _,_,variation=solve_flows(cfg,pins,fields,varied,heat['peripheral_moderator_heat_W'])
            result['property_sensitivities'].append(dict(case=label,kmod_W_mK=km,htc_multiplier=hm,**variation))
    result['validation_only']=validation
    result['model']='XY moderator/liner FEM response fields + 1D axial helium energy/friction + selected RZ fuel/clad FEM; file-based one-way OpenMC heating transfer'
    result['limitations']=['Constant assumed YH2 conductivity 10 W/m/K and Al2O3 6 W/m/K; not validated irradiated properties','Prescribed parallel flow split 17.5% auxiliary; no plenum/orifice pressure-balance solution','Hydraulic-diameter heat-transfer closures; no resolved ellipse CFD','Friction pressure losses only; no manifolds/spacers/acceleration losses','Peripheral heat allocated to outer-ring cells; peripheral solid temperatures and full-core peak are unresolved','Reflector/absorber heat and cooling excluded from active thermal model','No two-way reactor temperature feedback; cell profiles at prescribed nominal neutronic temperatures','Monte Carlo heating statistics recorded separately; thermal systematic uncertainties are not a statistical confidence interval']
    write_json(d/'summary.json',result)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--statepoint');p.add_argument('--reference',action='store_true');p.add_argument('--validate',action='store_true');a=p.parse_args()
    if not a.statepoint and not a.validate:p.error('--statepoint or --validate required')
    print(json.dumps(run(a.statepoint,a.reference,a.validate),indent=2))
