"""Explicit UN/He/vacuum/Mo/canned-YH2 CSG, with mass-conserving feedback."""
from pathlib import Path
import math
import numpy as np
import openmc
from catalog import layers,pins

XS=Path('/home/aleri/OpenMC_projects/first_test/new_possible_best_configurations_failed/campaign/data/endfb-viii.0-hdf5/cross_sections.xml')
OPENMC='/home/aleri/miniforge3/envs/openmc-env/bin/openmc'
SLICES=8


def initial_field(spec):
    n=len(pins(spec));s=SLICES
    cold=spec['inlet_K']+15.
    return dict(fuel=np.full((n,s),1350.).tolist(),moderator=np.full((n,s),cold+20.).tolist(),
        hot_he=np.tile(np.linspace(cold+30,spec['outlet_K']-10,s),(n,1)).tolist(),
        cold_he=np.full((n,s),cold).tolist(),outer_cold_he=np.full((n,s),cold).tolist(),clad=1250.,hot_wall=1100.,
        cold_wall=cold+5,foil=1000.,can=cold+15.,barrier=cold+15.,
        reflector=cold+10.,absorber=cold+10.)


def build(directory,spec,*,field=None,theta=0.,particles=5000,batches=140,
          inactive=40,seed=1001,heating=False,target=None,max_batches=2000,
          perturb_kind=None,delta_K=0.,beta_mod=0.,beta_fuel=0.,volume=False):
    d=Path(directory);d.mkdir(parents=True,exist_ok=True)
    openmc.reset_auto_ids();openmc.config['cross_sections']=str(XS)
    pp=pins(spec);n=len(pp);ll=layers(spec);idx={x['name']:j for j,x in enumerate(ll)}
    T=initial_field(spec) if field is None else {k:np.asarray(v).copy() if isinstance(v,list) else v for k,v in field.items()}
    for key in ('fuel','moderator','hot_he','cold_he','outer_cold_he'):T[key]=np.asarray(T[key],dtype=float)
    if perturb_kind:
        keys=['fuel'] if perturb_kind=='fuel' else ['moderator'] if perturb_kind=='moderator' else list(T)
        for key in keys:T[key]=T[key]+delta_K
    vf=1+beta_fuel*delta_K if perturb_kind in ('fuel','total') else 1.
    vm=1+beta_mod*delta_K if perturb_kind in ('moderator','total') else 1.
    radii=[x['r_outer_cm'] for x in ll]
    radii[0]*=math.sqrt(vf)
    mi=ll[idx['moderator']]['r_inner_cm'];mo=ll[idx['moderator']]['r_outer_cm']
    radii[idx['moderator']]=math.sqrt(mi**2+(mo**2-mi**2)*vm)
    if not all(a<b for a,b in zip(radii,radii[1:])):
        raise ValueError('Thermal expansion exceeds explicit radial clearance')
    mats=[]
    def mat(name,rho,items,key,sabs=()):
        m=openmc.Material(name=name);m.set_density('g/cm3',rho)
        for nu,q,unit in items:
            if nu in ('Nb','Zr','Al','O','Mo','C'):m.add_element(nu,q,unit)
            else:m.add_nuclide(nu,q,unit)
        for sab in sabs:m.add_s_alpha_beta(sab)
        m.temperature=float(np.mean(T[key]));mats.append(m);return m
    fuel={}
    for zone,enr in [('inner',spec['enrichment_inner']),('outer',spec['enrichment_outer'])]:
        f15=spec['nitrogen15_atom_fraction'];mn=(1-f15)*14.003074+f15*15.000109
        mu=1/(enr/235.0439299+(1-enr)/238.0507884)
        fuel[zone]=mat('UN_'+zone,spec['fuel_density_g_cm3']/vf,[('U235',enr*mu/(mu+mn),'wo'),
            ('U238',(1-enr)*mu/(mu+mn),'wo'),('N14',(1-f15)*14.003074/(mu+mn),'wo'),
            ('N15',f15*15.000109/(mu+mn),'wo')],'fuel')
    mod=mat('YH2',spec['moderator_density_g_cm3']/vm,[('Y89',1.,'ao'),('H1',spec['moderator_x'],'ao')],
            'moderator',('c_H_in_YH2','c_Y_in_YH2'))
    if spec['moderator_x']!=2.:raise ValueError('Only stoichiometric YH2 TSL available in validated local library')
    refl=mat('reflector',1.85 if spec['reflector']=='Be' else 2.9,
        [('Be9',1.,'ao')] if spec['reflector']=='Be' else [('Be9',1.,'ao'),('O16',1.,'ao')],
        'reflector',('c_Be',) if spec['reflector']=='Be' else ('c_Be_in_BeO','c_O_in_BeO'))
    af=spec['b10_atom_fraction']
    absorber=mat('B4C',2.52,[('B10',4*af,'ao'),('B11',4*(1-af),'ao'),('C',1.,'ao')],'absorber')
    lm={};gas_layers={}
    for layer in ll:
        name=layer['name'];material=layer['material']
        if name=='fuel':continue
        if material=='void':lm[name]=None;continue
        if material=='YH2':lm[name]=mod;continue
        if material.startswith('He'):
            key='hot_he' if material=='He_hot' else 'outer_cold_he' if name=='outer_cold_He' else 'cold_he'
            rho=spec['pressure_Pa']*.004002602/(8.314462618*float(np.mean(T[key])))/1000
            lm[name]=mat(name,rho,[('He4',1.,'ao')],key)
            gas_layers[name]=[]
            for s in range(SLICES):
                ts=float(np.mean(T[key][:,s]))
                rh=spec['pressure_Pa']*.004002602/(8.314462618*ts)/1000
                gm=mat(f'{name}_z{s}',rh,[('He4',1.,'ao')],key)
                gm.temperature=ts;gas_layers[name].append(gm)
        elif material=='Nb1Zr':
            key='clad' if name=='clad' else 'hot_wall' if name=='hot_wall' else 'cold_wall' if name=='cold_wall' else 'can'
            lm[name]=mat(name,8.57,[('Nb',.99,'wo'),('Zr',.01,'wo')],key)
        elif material=='Mo':lm[name]=mat(name,10.28,[('Mo',1.,'ao')],'foil')
        else:lm[name]=mat(name,3.95,[('Al',2.,'ao'),('O',3.,'ao')],'barrier')
    h=spec['height_cm'];ah=spec['axial_reflector_cm'];cr=spec['core_radius_cm']
    zb=openmc.ZPlane(z0=-h/2);zt=openmc.ZPlane(z0=h/2)
    zlo=openmc.ZPlane(z0=-h/2-ah,boundary_type='vacuum');zhi=openmc.ZPlane(z0=h/2+ah,boundary_type='vacuum')
    active=+zb&-zt;end=(+zlo&-zb)|(+zt&-zhi)
    rr=[openmc.ZCylinder(r=x) for x in radii]
    split=[openmc.ZPlane(z0=-h/2+j*h/SLICES) for j in range(SLICES+1)]
    cells=[];pinregions=[];scored={x['name']:[] for x in ll if x['material']!='void'};scored_end={}
    def shell(i):return -rr[i] if i==0 else +rr[i-1]&-rr[i]
    for pin in pp:
        p=pin['id'];local=[];zone='inner' if pin['ring']<=2 else 'outer'
        for j,layer in enumerate(ll):
            name=layer['name'];fill=fuel[zone] if name=='fuel' else lm.get(name)
            for s in range(SLICES):
                actualfill=gas_layers[name][s] if name in gas_layers else fill
                c=openmc.Cell(name=f'{name}_p{p}_s{s}',fill=actualfill,
                    region=shell(j)&+split[s]&-split[s+1])
                if fill:
                    key='fuel' if name=='fuel' else 'moderator' if name=='moderator' else 'hot_he' if layer['material']=='He_hot' else 'outer_cold_he' if name=='outer_cold_He' else 'cold_he' if layer['material']=='He_cold' else None
                    if key:c.temperature=float(T[key][p,s])
                    scored[name].append(c)
                local.append(c)
        local.append(openmc.Cell(fill=refl,region=+rr[-1]&active))
        # The hot bore, vacuum shield and cold jacket continue through BOTH
        # axial reflectors. Replacing these layers by Be would short insulation.
        passage=[]
        endnames=['hot_He','hot_wall','vacuum_inner','Mo_shield','vacuum_outer',
            'cold_wall','cold_He','inner_can','inner_barrier','outer_barrier',
            'outer_can','outer_cold_He','outer_jacket']
        for name in endnames:
            if name not in idx:continue
            j=idx[name];rad=-rr[j] if name=='hot_He' else shell(j)
            for eslab in (+zlo&-zb,+zt&-zhi):
                reg=rad&eslab;passage.append(reg)
                ec=openmc.Cell(fill=lm.get(name),region=reg);local.append(ec)
                if lm.get(name):scored_end.setdefault(name,[]).append(ec)
        local.append(openmc.Cell(fill=refl,region=end&~openmc.Union(passage)))
        u=openmc.Universe(cells=local)
        cylinder=openmc.ZCylinder(x0=pin['x_cm'],y0=pin['y_cm'],r=radii[-1])
        reg=-cylinder&+zlo&-zhi
        c=openmc.Cell(fill=u,region=reg);c.translation=(pin['x_cm'],pin['y_cm'],0.)
        cells.append(c);pinregions.append(-cylinder)
    core=openmc.ZCylinder(r=cr);outer=openmc.ZCylinder(r=cr+spec['radial_reflector_cm'],boundary_type='vacuum')
    cells.append(openmc.Cell(fill=refl,region=-core&~openmc.Union(pinregions)&+zlo&-zhi))
    dr=spec['drum_radius_cm'];di=dr-spec['absorber_thickness_cm'];half=math.radians(spec['absorber_sector_deg']/2)
    dc=openmc.ZCylinder(r=dr);din=openmc.ZCylinder(r=di)
    a=openmc.Plane(a=math.sin(half),b=-math.cos(half));bb=openmc.Plane(a=math.sin(half),b=math.cos(half))
    ar=-dc&+din&+a&+bb
    du=openmc.Universe(cells=[openmc.Cell(fill=absorber,region=ar),openmc.Cell(fill=refl,region=-dc&~ar)])
    drumregions=[]
    for i in range(spec['drum_count']):
        az=360*i/spec['drum_count'];ph=math.radians(az)
        x=spec['drum_centres_radius_cm']*math.cos(ph);y=spec['drum_centres_radius_cm']*math.sin(ph)
        cylinder=openmc.ZCylinder(x0=x,y0=y,r=dr);reg=-cylinder&active
        c=openmc.Cell(fill=du,region=reg);c.translation=(x,y,0.);c.rotation=(0.,0.,az+theta)
        cells.append(c);drumregions.append(reg)
    cells.append(openmc.Cell(fill=refl,region=+core&-outer&active&~openmc.Union(drumregions)))
    cells.append(openmc.Cell(fill=refl,region=+core&-outer&end))
    settings=openmc.Settings();settings.particles=particles;settings.batches=batches;settings.inactive=inactive;settings.seed=seed
    settings.source=openmc.IndependentSource(space=openmc.stats.Box((-cr,-cr,-h/2),(cr,cr,h/2)),constraints={'fissionable':True})
    settings.temperature={'method':'interpolation','tolerance':10.,'multipole':not heating}
    settings.photon_transport=heating;settings.sourcepoint={'write':False}
    settings.output={'tallies':heating}
    if target:
        settings.keff_trigger={'type':'std_dev','threshold':target};settings.trigger_active=True
        settings.trigger_max_batches=max_batches;settings.trigger_batch_interval=20
    tallies=[]
    if heating:
        for name,subset in scored.items():
            t=openmc.Tally(name=name+'_heating');t.filters=[openmc.CellFilter(subset)];t.scores=['heating'];tallies.append(t)
        for name,subset in scored_end.items():
            t=openmc.Tally(name='end_'+name+'_heating');t.filters=[openmc.CellFilter(subset)];t.scores=['heating'];tallies.append(t)
        t=openmc.Tally(name='all_material_heating');t.filters=[openmc.MaterialFilter(mats)];t.scores=['heating'];tallies.append(t)
        mesh=openmc.RegularMesh();ro=cr+spec['radial_reflector_cm']
        mesh.lower_left=(-ro,-ro,-h/2-ah);mesh.upper_right=(ro,ro,h/2+ah);mesh.dimension=(24,24,10)
        t=openmc.Tally(name='reflector_spatial_heating')
        t.filters=[openmc.MaterialFilter([refl]),openmc.MeshFilter(mesh)];t.scores=['heating'];tallies.append(t)
    materials=openmc.Materials(mats);materials.cross_sections=str(XS)
    model=openmc.Model(openmc.Geometry(cells),materials,settings,openmc.Tallies(tallies))
    if volume:
        settings.run_mode='volume';ro=cr+spec['radial_reflector_cm']
        settings.volume_calculations=[openmc.VolumeCalculation(mats,3000000,[-ro,-ro,-h/2-ah],[ro,ro,h/2+ah])]
    model.export_to_xml(d)
    return model,dict(scored=scored,scored_end=scored_end,materials=mats,layers=ll,pins=pp,field=T)
