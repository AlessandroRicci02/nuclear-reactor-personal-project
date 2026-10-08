from pathlib import Path
import json
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Wedge
from catalog import layers,pins


def figures(directory,candidate):
    d=Path(directory)/'figures';d.mkdir(parents=True,exist_ok=True)
    s=candidate['spec'];ll=layers(s)
    colors={'UN':'#af543c','Nb1Zr':'#777c83','He_hot':'#edaa50','He_cold':'#74c3d6',
        'YH2':'#738cb9','Al2O3':'#e6dab2','Mo':'#a59564','void':'#ffffff'}
    fig,ax=plt.subplots(figsize=(8,8))
    for layer in reversed(ll):
        ax.add_patch(Circle((0,0),layer['r_outer_cm'],color=colors[layer['material']],ec='black',lw=.3))
    r=ll[-1]['r_outer_cm'];ax.set(xlim=(-r*1.15,r*1.15),ylim=(-r*1.15,r*1.15),aspect='equal',
        xlabel='cm',ylabel='cm',title=candidate['id']+' — pin section')
    ax.text(1.02,.99,'\n'.join(f"{i+1}. {x['name']}" for i,x in enumerate(ll)),
        transform=ax.transAxes,va='top',fontsize=8)
    fig.savefig(d/'pin_section.png',dpi=180,bbox_inches='tight');plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,9));cr=s['core_radius_cm'];ro=cr+s['radial_reflector_cm']
    ax.add_patch(Circle((0,0),ro,color='#d6e2cb'));ax.add_patch(Circle((0,0),cr,fill=False,ec='black'))
    for p in pins(s):
        ax.add_patch(Circle((p['x_cm'],p['y_cm']),r,color='#738cb9',ec='black',lw=.15))
        ax.add_patch(Circle((p['x_cm'],p['y_cm']),s['fuel_radius_cm'],color='#af543c'))
    for i in range(s['drum_count']):
        a=360*i/s['drum_count'];ar=math.radians(a);dc=s['drum_centres_radius_cm']
        xy=(dc*math.cos(ar),dc*math.sin(ar));dr=s['drum_radius_cm'];half=s['absorber_sector_deg']/2
        ax.add_patch(Circle(xy,dr,fill=False,ec='black',lw=.6))
        ax.add_patch(Wedge(xy,dr,a-half,a+half,width=s['absorber_thickness_cm'],color='#343b46'))
    ax.set(xlim=(-ro*1.05,ro*1.05),ylim=(-ro*1.05,ro*1.05),aspect='equal',xlabel='cm',ylabel='cm',
        title=candidate['id']+' — absorber OUT (0°)')
    fig.savefig(d/'core_section.png',dpi=180,bbox_inches='tight');plt.close(fig)


def write_found(directory,result):
    d=Path(directory);c=result['candidate'];s=c['spec'];th=result['thermal'];metrics=c['metrics']
    lines=[f"CANDIDATO CHE SUPERA I FILTRI BOL: {c['id']}",
        "Superamento di questa campagna di simulazioni, non qualifica regolatoria.",
        f"Moderatore YH{s['moderator_x']:.1f}; combustibile UN; riflettore {s['reflector']}.",
        f"Azoto del fuel: {100*s['nitrogen15_atom_fraction']:.3f}% atomico N15, restante N14.",
        f"{metrics['pin_count']} pin, altezza {s['height_cm']} cm, passo {s['pitch_cm']} cm.",
        f"He: {s['inlet_K']} K ingresso / {s['outlet_K']} K uscita; potenza {s['power_W']/1000} kWth.",
        f"Massa geometrica preliminare superiore: {metrics['preliminary_solid_mass_upper_kg']:.1f} kg; schermatura e hardware esterno esclusi.",
        f"Inventario U235 {metrics['initial_U235_kg']:.2f} kg; burnup richiesto in 10 anni {metrics['burnup_10yr_MWd_kgHM']:.2f} MWd/kgHM. Lifespan non simulato.",
        f"k OUT = {result['out']['keff']:.6f} +/- {result['out']['sigma_k']:.6f}",
        f"k IN = {result['inward']['keff']:.6f} +/- {result['inward']['sigma_k']:.6f}",
        f"Angolo critico {result['critical']['theta']:.3f} gradi: k = {result['critical']['keff']:.6f} +/- {result['critical']['sigma_k']:.6f}",
        f"Tmax moderatore nominale {th['nominal']['max_moderator_K']:.2f} K; limite superiore di accettazione {th['acceptance_upper_K']:.2f} K.",
        "Coefficienti in pcm/K (reattivita' 1-1/k; differenze centrali a tamburi bloccati all'angolo critico):"]
    for a in result['feedback']+result['confirmation']:
        lines.append(f"{a['kind']}: {a['alpha_pcm_K']:.4f} +/- {a['sigma_pcm_K']:.4f}, superiore 4sigma {a['upper_bound_pcm_K']:.4f}; betaY={a['beta_mod_per_K']}, betaUN={a['beta_fuel_per_K']}, delta={a['delta_K']} K")
    lines += ["Espansione radiale conservativa della massa entro giochi espliciti; canning e altezza fissi.",
        "Incertezze nucleari di libreria e qualifica materiali/meccanica non incluse. Vedere protocollo e risultati termici completi.",
        "Le simulazioni sono state fermate al primo candidato che supera tutti i filtri."]
    (d/'report.txt').write_text('\n'.join(lines)+'\n');figures(d,c)
