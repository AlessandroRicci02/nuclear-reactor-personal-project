"""Geometria concettuale e bilanci analitici; NON esegue simulazioni nucleari/TH."""
from pathlib import Path
import csv
import json
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, Rectangle, Patch

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / 'figures'
FIG.mkdir(exist_ok=True)
P = dict(name='S2_YH_SLEEVES_BE_COLD_FIRST', status='CONCEPT_NOT_VALIDATED',
         power_W=400000, pitch_cm=2.30, rings=8, active_height_cm=50.0,
         core_radius_cm=19.70, radial_reflector_cm=8.80, axial_reflector_cm=8.50,
         drum_count=10, drum_radius_cm=3.95, drum_centres_radius_cm=23.83,
         absorber_thickness_cm=1.50, absorber_sector_deg=148.0,
         b10_atom_fraction=0.90, pressure_Pa=2500000,
         inlet_K=700.0, mixed_outlet_target_K=880.0, helium_cp_J_kgK=5193.16,
         hardware_allowance_kg=60.0, contingency_kg=15.0,
         nominal_moderator_peak_target_K=900.0, moderator_limit_K=950.0,
         life_goal_full_power_years=18.0, minimum_life_full_power_years=15.0)
# r_outer [cm], density [g/cm3], material, color. All annuli are disjoint.
LAYERS = [
 ('Combustibile UN', .585, 13.5, 'UN', '#bd5d45'),
 ('Guaina combustibile', .633, 8.57, 'Nb-1Zr', '#525e70'),
 ('Elio: secondo passaggio', .717, 0., 'He', '#ecaf65'),
 ('Parete separatrice', .737, 8.57, 'Nb-1Zr', '#525e70'),
 ('Elio: primo passaggio', .837, 0., 'He', '#64bcdf'),
 ('Contenitore interno YH2', .857, 8.57, 'Nb-1Zr', '#525e70'),
 ('Barriera interna proposta', .867, 3.95, 'Al2O3', '#eee1bd'),
 ('Moderatore YH2', 1.087, 4.30, 'YH2', '#738dbe'),
 ('Barriera esterna proposta', 1.097, 3.95, 'Al2O3', '#eee1bd'),
 ('Contenitore esterno YH2', 1.117, 8.57, 'Nb-1Zr', '#525e70'),
]
pins=[]
for q in range(-P['rings'], P['rings']+1):
    for r in range(-P['rings'], P['rings']+1):
        ring=max(abs(q),abs(r),abs(q+r))
        if ring<=P['rings']:
            pins.append(dict(id=len(pins)+1,q=q,r=r,ring=ring,
                x_cm=P['pitch_cm']*(q+r/2),y_cm=P['pitch_cm']*math.sqrt(3)*r/2,
                u235_weight_fraction=.196 if ring<=2 else .1975))
n=len(pins); assert n==217
h=P['active_height_cm']; az=P['axial_reflector_cm']; rc=P['core_radius_cm']
ro=rc+P['radial_reflector_cm']; ht=h+2*az; pi=math.pi
min_distance=min(math.hypot(a['x_cm']-b['x_cm'],a['y_cm']-b['y_cm'])
                 for i,a in enumerate(pins) for b in pins[i+1:])
assert min_distance>2*LAYERS[-1][1]
assert max(math.hypot(p['x_cm'],p['y_cm']) for p in pins)+LAYERS[-1][1]<rc
dc=P['drum_centres_radius_cm']; dr=P['drum_radius_cm']
assert dc-dr>rc and dc+dr<ro and 2*dc*math.sin(pi/P['drum_count'])>2*dr
rows=[]
def add(name,volume,density):
    assert volume>=0
    rows.append(dict(component=name,volume_cm3=volume,density_g_cm3=density,mass_kg=volume*density/1000))
r0=0
for name,r,density,_,_ in LAYERS:
    if density: add(name,n*pi*(r*r-r0*r0)*h,density)
    r0=r
# The moderator and fuel stop at active ends. Through each end reflector:
# hot bore r<.717, divider .717-.737, cold annulus .737-.837,
# duct wall .837-.857; Be occupies the remainder outside .857.
add('Pareti condotti nei due riflettori assiali',n*pi*((.737**2-.717**2)+(.857**2-.837**2))*2*az,8.57)
vb=P['drum_count']*P['absorber_sector_deg']/360*pi*(dr**2-(dr-P['absorber_thickness_cm'])**2)*h
add('B4C settori tamburi',vb,2.52)
add('Be interstiziale e periferico nel nocciolo',(pi*rc**2-n*pi*1.117**2)*h,1.85)
add('Be riflettori assiali interni al raggio nocciolo',(pi*rc**2-n*pi*.857**2)*2*az,1.85)
add('Be riflettore radiale e corpi tamburi',pi*(ro**2-rc**2)*ht-vb,1.85)
solid=sum(r['mass_kg'] for r in rows)
gasvol=n*pi*((.717**2-.633**2)+(.837**2-.737**2))*h+n*pi*(.717**2+.837**2-.737**2)*2*az
closure=sum(r['volume_cm3'] for r in rows)+gasvol-pi*ro**2*ht
assert abs(closure)<1e-8
fuel=rows[0]['mass_kg']
# Explicit, approximate historical mass convention, not a depletion composition.
hm=fuel*238/252
enr=sum(p['u235_weight_fraction'] for p in pins)/n
flow=P['power_W']/P['helium_cp_J_kgK']/(P['mixed_outlet_target_K']-P['inlet_K'])
metrics=dict(solid_mass_kg=solid,concept_mass_with_allowances_kg=solid+75,
    remaining_to_600_kg=600-solid-75,maximum_hardware_plus_contingency_kg=600-solid,
    fuel_UN_kg=fuel,approx_HM_kg=hm,approx_U235_kg=hm*enr,
    moderator_kg=next(r['mass_kg'] for r in rows if r['component']=='Moderatore YH2'),
    outer_radius_cm=ro,total_reflector_height_cm=ht,helium_flow_kg_s=flow,
    active_hot_flow_area_cm2=n*pi*(.717**2-.633**2),
    active_cold_flow_area_cm2=n*pi*(.837**2-.737**2),
    cell_clearance_mm=10*(min_distance-2*1.117),
    outer_cell_clearance_mm=10*(rc-8*P['pitch_cm']-1.117),
    drum_inner_web_mm=10*(dc-dr-rc),drum_outer_web_mm=10*(ro-dc-dr),
    mass_volume_closure_error_cm3=closure,
    modeled_gas_volume_cm3=gasvol,
    energy_15y_MWd=.4*365.25*15,burnup_15y_MWd_kgHM=.4*365.25*15/hm,
    energy_18y_MWd=.4*365.25*18,burnup_18y_MWd_kgHM=.4*365.25*18/hm)
spec=dict(parameters=P,radial_layers=[dict(name=a,r_outer_cm=b,density_g_cm3=c,material=d) for a,b,c,d,_ in LAYERS],
          metrics=metrics,notes=['Dimensions nominal; no thermal expansion or tolerances.',
          'Nb-1Zr/Al2O3 hydrogen containment is a proposal, not qualified.',
          '60 kg hardware + 15 kg contingency are allowances, not modeled parts.',
          'No keff, feedback, lifetime or peak temperature calculated.'])
(ROOT/'configuration.json').write_text(json.dumps(spec,indent=2,ensure_ascii=False)+'\n')
for filename,data in [('pins.csv',pins),('mass_breakdown.csv',rows)]:
    with (ROOT/filename).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0])); w.writeheader();w.writerows(data)
(ROOT/'screening.json').write_text(json.dumps(metrics,indent=2)+'\n')
plt.rcParams.update({'font.size':10,'axes.titlesize':13,'figure.facecolor':'white'})
def save(fig,name):
    fig.savefig(FIG/(name+'.png'),dpi=180,bbox_inches='tight')
    fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight')
    plt.close(fig)

fig,ax=plt.subplots(figsize=(10,9))
ax.add_patch(Circle((0,0),ro,facecolor='#b7cdb9',edgecolor='#334b3f'))
ax.add_patch(Circle((0,0),rc,fill=False,linestyle='--',edgecolor='#334b3f'))
for p in pins:
    for _,r,_,_,color in reversed(LAYERS): ax.add_patch(Circle((p['x_cm'],p['y_cm']),r,facecolor=color,linewidth=0))
for i in range(P['drum_count']):
    angle=i*360/P['drum_count']; xy=(dc*math.cos(math.radians(angle)),dc*math.sin(math.radians(angle)))
    ax.add_patch(Circle(xy,dr,fill=False,edgecolor='#304b38',linewidth=1))
    ax.add_patch(Wedge(xy,dr,angle-74,angle+74,width=1.5,facecolor='#403348'))
ax.set(xlim=(-ro-2,ro+2),ylim=(-ro-2,ro+2),aspect='equal',xlabel='x [cm]',ylabel='y [cm]',
       title='S2 — Sezione attiva proposta, tamburi a 0° (OUT)\n217 celle; YH₂ confinato in manicotti; matrice e riflettore Be')
ax.legend(handles=[Patch(color=c,label=t) for t,c in [('UN','#bd5d45'),('YH₂','#738dbe'),('Be','#b7cdb9'),('B₄C','#403348'),('Elio freddo','#64bcdf')]],loc='upper right',fontsize=8)
fig.text(.5,.02,'Geometria concettuale: criticità e coefficienti NON calcolati. Vessel e attuatori non disegnati.',ha='center',fontsize=9)
save(fig,'01_core_xy')

fig,(ax,txt)=plt.subplots(1,2,figsize=(12,6),gridspec_kw={'width_ratios':[1,1.3]})
for _,r,_,_,color in reversed(LAYERS): ax.add_patch(Circle((0,0),r*10,facecolor=color,edgecolor='white',linewidth=.3))
ax.set(xlim=(-12,12),ylim=(-12,12),aspect='equal',xlabel='x [mm]',ylabel='y [mm]',title='Cella coassiale — quote nominali')
txt.axis('off'); textlines=['RAGGI ESTERNI (dal centro):','']+[f'{r*10:5.2f} mm   {name}' for name,r,_,_,_ in LAYERS]
textlines+=['',f'Passo: {P["pitch_cm"]*10:.2f} mm; altezza attiva: {h*10:.0f} mm',
            'Gap combustibile/guaina e tolleranze: da progettare.',
            'Barriere H e chiusure terminali: da qualificare.',
            'Elio freddo a contatto con la parete interna del moderatore.']
txt.text(0,1,'\n'.join(textlines),va='top',fontsize=10,linespacing=1.6)
fig.tight_layout();save(fig,'02_unit_cell')

fig,(ax,txt)=plt.subplots(1,2,figsize=(12,7),gridspec_kw={'width_ratios':[1,1.1]})
# A single-cell longitudinal section, radial width exaggerated by aspect ratio.
ax.add_patch(Rectangle((-1.15,-ht/2),2.3,ht,facecolor='#b7cdb9'))
for _,r,_,_,color in reversed(LAYERS): ax.add_patch(Rectangle((-r,-h/2),2*r,h,facecolor=color))
for z in (-ht/2,h/2):
    for left,right,color in [(-.857,.857,'#525e70'),(-.837,.837,'#64bcdf'),(-.737,.737,'#525e70'),(-.717,.717,'#ecaf65')]:
        ax.add_patch(Rectangle((left,z),right-left,az,facecolor=color))
for x in (-.787,.787): ax.annotate('',xy=(x,-h/2+1),xytext=(x,h/2-1),arrowprops=dict(arrowstyle='->',color='#07557b',lw=2))
for x in (-.675,.675): ax.annotate('',xy=(x,h/2-1),xytext=(x,-h/2+1),arrowprops=dict(arrowstyle='->',color='#943918',lw=2))
ax.set(xlim=(-1.2,1.2),ylim=(-ht/2-3,ht/2+3),xlabel='r firmato [cm] — scala radiale dilatata',ylabel='z [cm]',title='Sezione longitudinale di una cella')
txt.axis('off'); txt.text(0,1,
    'PERCORSO SERIALE PROPOSTO\n\n'
    '1. Ingresso superiore: He a 700 K, 2.5 MPa.\n'
    '2. Primo passaggio ↓ nell’anello freddo,\n   adiacente al contenitore del moderatore.\n'
    '3. Collettore inferiore: inversione del flusso.\n'
    '4. Secondo passaggio ↑ lungo la guaina UN.\n'
    '5. Uscita superiore separata: miscela a 880 K.\n\n'
    f'Portata totale di bilancio: {flow:.4f} kg/s.\n'
    'Stessa portata nei due passaggi; nessuno split\n'
    'arbitrario tra combustibile e moderatore.\n\n'
    'Collettori, setti di separazione e orifizi\n'
    'sono requisiti di progetto, qui non dimensionati.\n\n'
    'Tra anelli caldo/freddo: parete metallica,\n'
    'NON isolamento ideale. Lo scambio interno\n'
    'deve entrare nel bilancio coniugato.\n\n'
    '900 K: obiettivo picco YH₂, NON risultato.\n'
    '880 K: target medio He, NON limite locale.',va='top',linespacing=1.5)
fig.tight_layout();save(fig,'03_axial_and_flow')

fig,(ax,ax2)=plt.subplots(1,2,figsize=(12,5))
groups=[('UN',fuel),('YH₂',metrics['moderator_kg']),('Be',sum(r['mass_kg'] for r in rows if r['component'].startswith('Be '))),
        ('Altri solidi',solid-fuel-metrics['moderator_kg']-sum(r['mass_kg'] for r in rows if r['component'].startswith('Be '))),
        ('Hardware*',60),('Riserva*',15)]
bottom=0
for (name,value),color in zip(groups,['#bd5d45','#738dbe','#b7cdb9','#777777','#ceb783','#eee1bd']):
    ax.bar('S2',value,bottom=bottom,label=f'{name}: {value:.1f} kg',color=color);bottom+=value
ax.bar('TB35\nsoli solidi',595.154,color='#b0b0b0')
ax.axhline(600,color='firebrick',linestyle='--');ax.set(ylabel='Massa [kg]',ylim=(0,650),title='Bilancio geometrico + assegnazioni*')
ax.legend(fontsize=8,loc='upper left',bbox_to_anchor=(0,-.1))
loads=np.array([20,60,100,160,240])
temps=700+loads*1000/flow/P['helium_cp_J_kgK']
ax2.plot(loads,temps,'o-',color='#167da5',label='Bulk uscita primo passaggio')
ax2.axhline(900,color='#947328',ls='--',label='Obiettivo picco YH₂')
ax2.axhline(950,color='firebrick',ls='--',label='Limite picco YH₂')
ax2.set(xlabel='Calore netto al primo passaggio [kW], ipotizzato',ylabel='Temperatura [K]',ylim=(680,970),title='Solo bilancio energetico: non predice il picco')
ax2.legend(fontsize=8)
fig.text(.5,-.18,'*Hardware e riserva non dimensionati. Film, conduzione, potenze locali e scambio fra passaggi ancora da calcolare.',ha='center',fontsize=9)
fig.tight_layout();save(fig,'04_mass_and_energy_screen')
print(json.dumps(metrics,indent=2))
