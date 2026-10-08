# Reattore Nucleare Spaziale NEP (400 kWth / 100+ kWe)
## Geometria e Architettura Nucleare Definitiva: `C01_matrix_217_yh20_880_1120`

Benvenuti nell'archivio consolidato e definitivo del reattore nucleare spaziale progettato per missioni di **Propulsione Elettrica Nucleare (Nuclear Electric Propulsion - NEP)** nello spazio profondo.

Questo repository raccoglie tutti i dati progettuali, i risultati delle simulazioni Monte Carlo ad alta fedeltà (**OpenMC**), i benchmark multifisici e termoidraulici agli elementi finiti (**MOOSE / Cardinal / Griffin**), l'analisi di burnup/depletion a lungo termine con controllo attivo di criticità, e l'intera galleria di 13 figure ad alta risoluzione (300 DPI).

---

## 1. Sintesi Prestazionale e Rispetto dei Requisiti di Progetto

Tutti i vincoli imposti per la missione NEP sono pienamente verificati e superati:

| Parametro / Vincolo di Progetto | Valore Ottenuto / Calcolato | Specifiche di Progetto | Margine di Sicurezza | Esito |
| :--- | :--- | :--- | :--- | :---: |
| **Potenza Termica Nocciolo** | **400.0 kWth** | 400.0 kWth | Nominale | **PASS** |
| **Potenza Elettrica Netta** | **136.39 kWe** | $\ge 100.0\text{ kWe}$ | **+36.39 kWe (+36.4%)** | **PASS** |
| **Rendimento Termodinamico Globale** | **34.10 %** | $\ge 30.0\text{ \%}$ | **+4.10 %** | **PASS** |
| **Moderatore Solido Impiegato** | **$\text{YH}_{2.0}$ ($\rho = 4.30\text{ g/cm}^3$)** | $\text{YH}_x$ ($x$ a scelta) | Idruro a max densità H | **PASS** |
| **Temperatura Ingresso Gas (Elio)** | **$880.0\text{ K}$** | $850 - 900\text{ K}$ | Perfettamente centrato | **PASS** |
| **Temperatura Uscita Gas (Elio)** | **$1120.0\text{ K}$** | $1100 - 1150\text{ K}$ | Perfettamente centrato | **PASS** |
| **Portata Massica Refrigerante (He)** | **$0.32093\text{ kg/s}$** | Circuito a $2.5\text{ MPa}$ | Ottimizzata per $\Delta T = 240\text{ K}$ | **PASS** |
| **$k_{\text{eff}}$ BOL Drums OUT ($0^\circ$, $2\sigma$)** | **$1.05671$** ($1.05945 \pm 0.00137$) | $> 1.00000$ | **$+5611.8\text{ pcm}$ excess** | **PASS** |
| **$k_{\text{eff}}$ BOL Drums IN ($180^\circ$, $2\sigma$)** | **$0.99811$** ($0.99551 \pm 0.00130$) | $< 1.00000$ | **$+451.3\text{ pcm}$ SDM** | **PASS** |
| **Coeff. Temp. Combustibile ($\alpha_{\text{fuel}}$, $2\sigma$)** | **$-1.2198\text{ pcm/K}$** ($-2.062 \pm 0.421$) | $< 0.0000\text{ pcm/K}$ | Doppler intrinsecamente stabile | **PASS** |
| **Coeff. Temp. Moderatore ($\alpha_{\text{mod}}$, $2\sigma$)** | **$-0.2518\text{ pcm/K}$** ($-1.382 \pm 0.565$) | $< 0.0000\text{ pcm/K}$ | Spettro termico stabilizzante | **PASS** |
| **Coeff. Temp. Totale Isotermico ($\alpha_{\text{tot}}$, $2\sigma$)** | **$-1.0476\text{ pcm/K}$** ($-2.465 \pm 0.709$) | $< 0.0000\text{ pcm/K}$ | Stabilità intrinseca garantita | **PASS** |
| **Temperatura Picco Moderatore $\text{YH}_{2.0}$ (FEM)** | **$935.31\text{ K}$** | $< 950.00\text{ K}$ | **$+14.69\text{ K}$ buffer termico** | **PASS** |
| **Temperatura Picco Guaina Nb-1Zr (FEM)** | **$1226.82\text{ K}$** | $< 1350.00\text{ K}$ | **$+123.18\text{ K}$ creep buffer** | **PASS** |
| **Temperatura Picco Combustibile UN (FEM)** | **$1239.61\text{ K}$** | $< 2073.00\text{ K}$ | **$+833.39\text{ K}$ dissoc. margin**| **PASS** |
| **Caduta di Pressione Core ($\Delta P$)** | **$14.67\text{ kPa}$** | $< 35.00\text{ kPa}$ | **$+20.33\text{ kPa}$ margine pompaggio**| **PASS** |
| **Massa Totale Reattore (Core Assembly)** | **$636.48\text{ kg}$** | $< 650.0\text{ kg}$ | **$-13.52\text{ kg}$ sotto il target** | **PASS** |
| **Vita Operativa Verificata (Burnup)** | **$25.00\text{ anni}$** | $> 10.0\text{ anni}$ | **$+15.0\text{ anni}$ oltre il requisito** | **PASS** |
| **Vita Operativa Estrapolata (Zero-crossing)**| **$44.80\text{ anni}$** | $> 10.0\text{ anni}$ | **$+34.8\text{ anni}$ operatività** | **PASS** |

---

## 2. Architettura del Nocciolo e Innovazioni Tecnologiche

### 2.1 Cella Combustibile a 5 Strati Concentrici
1. **Pastiglia di combustibile UN** ($r = 0.585\text{ cm}$): Nitruro di Uranio ad altissima conducibilità termica ($k \approx 20\text{ W/m-K}$) e densità uranica ($14.3\text{ g/cm}^3$). Arricchimento fissile pari al $19.75\text{ wt\%}$ U-235 (HALEU non proliferante).
2. **Guaina in Nb-1Zr** ($t = 0.048\text{ cm}$): Lega refrattaria biocompatibile con l'elio e ad altissima resistenza al creep termico fino a oltre $1350\text{ K}$.
3. **Canale anulare He caldo** ($t = 0.084\text{ cm}$): Refrigerante elio puro a $2.5\text{ MPa}$ con flusso ascendente verso il plenum superiore ($T_{\text{out}} = 1120\text{ K}$).
4. **Barriera ceramica di allumina $\text{Al}_2\text{O}_3$** ($t = 100\ \mu\text{m} = 0.010\text{ cm}$): Rivestimento denso a bassa permeabilità che blocca totalmente la migrazione e diffusione dell'idrogeno dal moderatore verso il circuito di elio caldo.
5. **Matrice solida continua di $\text{YH}_{2.0}$**: Reticolo esagonale a passo $2.28\text{ cm}$ con **6 canali ausiliari longitudinali di raffreddamento** ($r = 1.65\text{ mm}$ a $r_{\text{offset}} = 0.760\text{ cm}$ dal centro pin) alimentati con elio fresco a $880\text{ K}$.

### 2.2 Riflettore e Sistema di Controllo a Tamburi
- **Riflettore neutronico**: Ossido di Berillio (BeO, $\rho = 3.01\text{ g/cm}^3$) con spessore radiale di $9.30\text{ cm}$ e spessori assiali di $8.80\text{ cm}$ superiore e inferiore.
- **Tamburi di Controllo**: 10 cilindri rotanti in BeO ($r = 3.95\text{ cm}$) posizionati a raggio $22.23\text{ cm}$, ciascuno provvisto di un settore assorbitore in Carburo di Boro $\text{B}_4\text{C}$ (arricchito al $90\%\ ^{10}\text{B}$) di spessore $1.80\text{ cm}$ e arco angolare di $155.0^\circ$.

---

## 3. Struttura dell'Archivio e Indice dei File

```
definitive_final_geometry/
├── README.md                                    # Questa guida completa
├── report.txt                                   # Relazione tecnica esaustiva
├── results.json                                 # Risultati numerici strutturati
│
├── figures/                                     # Galleria di 13 figure ad alta risoluzione (300 DPI)
│   ├── 01_reactor_radial_xy.png                 # Sezione trasversale del reattore (X-Y)
│   ├── 02_reactor_axial_xz.png                  # Sezione longitudinale assiale del reattore (X-Z)
│   ├── 03_unit_cell_5layer_detailed.png         # Schema ingegneristico cella a 5 strati
│   ├── 04_control_drum_mechanism.png            # Meccanismo tamburo di controllo BeO/B4C
│   ├── 05_openmc_csg_slice_xy.png               # Rendering numerico OpenMC CSG (piano X-Y)
│   ├── 06_openmc_csg_slice_xz.png               # Rendering numerico OpenMC CSG (piano X-Z)
│   ├── 07_depletion_keff_and_excess_reactivity.png # Evoluzione keff e reattività (0-25 anni)
│   ├── 08_depletion_control_drum_critical_angle.png# Angolo critico tamburi vs tempo
│   ├── 09_depletion_isotopic_inventory_evolution.png# Consumo U-235 e produzione Pu-239
│   ├── 10_depletion_burnup_and_fima.png         # Burnup specifico (MWd/kgHM) e FIMA
│   ├── 11_moose_cardinal_temperature_profiles_pins.png# Benchmark termico FEM barre combustibile
│   ├── 12_moose_cardinal_moderator_thermal_field.png # Campo termico FEM nel moderatore YH2
│   └── 13_closed_brayton_cycle_ts_diagram.png   # Diagramma T-s del ciclo Brayton chiuso
│
├── data/                                        # Tabelle numeriche formattate
│   ├── c01_bol_neutronics_and_feedbacks.csv     # Parametri neutronici BOL e 2-sigma bounds
│   ├── c01_criticality_depletion_1yr_steps.csv  # 26 step annuali di burnup e criticità
│   ├── c01_moose_cardinal_pin_benchmarks.csv    # Benchmark FEM barre combustibile
│   ├── c01_brayton_cycle_statepoints.csv        # Punti termodinamici ciclo Brayton
│   └── results.json                             # Risultati macchina di C01
│
├── neutronics_openmc/                           # Script per neutronica OpenMC
│   ├── model.py                                 # Generatore modello CSG OpenMC
│   ├── thermal.py                               # Modulo termoidraulico ausiliario
│   └── generate_high_precision_plots_c01.py     # Generatore rendering 300 DPI figure geometriche
│
├── depletion_analysis/                          # Analisi di depletion a step annuali
│   ├── run_1yr_criticality_depletion.py         # Script OpenMC di burnup con controllo tamburi
│   └── results_depletion.txt                    # Report testuale della depletion
│
├── thermohydraulics_moose_cardinal/             # Modellazione termica FEM MOOSE/Cardinal
│   ├── run_thermodynamics_and_moose_cardinal_c01.py # Script risolutore FEM e ciclo Brayton
│   ├── c01_thermohydraulic_moose_report.txt     # Report testuale termico e FEM
│   └── cardinal_inputs/                         # File di input (.i) e output (.csv) FEM
│
└── configurations/                              # Cronologia e catalogo di ricerca
    ├── configurations.csv                       # Catalogo delle configurazioni esplorate
    ├── configurations.json                      # Parametri costruttivi di tutte le varianti
    └── progress.csv                             # Tabella di marcia dello screening
```

---

## 4. Galleria delle Figure ad Alta Risoluzione

Tutte le figure sono state renderizzate a **300 DPI** per garantire la massima leggibilità di quote, etichette e legende:

| File Figura | Descrizione Fisica e Dettagli |
| :--- | :--- |
| [01_reactor_radial_xy.png](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/01_reactor_radial_xy.png) | Sezione trasversale X-Y del nocciolo: evidenzia i 217 pin esagonali, la matrice YH2.0, i canali di bypass, il riflettore di BeO e i 10 tamburi rotanti. |
| [02_reactor_axial_xz.png](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/02_reactor_axial_xz.png) | Sezione longitudinale X-Z: altezza attiva di $46\text{ cm}$, riflettori assiali di $8.8\text{ cm}$, plenum gas e contenimento in Nb-1Zr. |
| [03_unit_cell_5layer_detailed.png](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/03_unit_cell_5layer_detailed.png) | Dettaglio microscopico della cella elementare: pastiglia UN, guaina Nb-1Zr, canale He, barriera $\text{Al}_2\text{O}_3$ ($100\ \mu\text{m}$) e matrice $\text{YH}_{2.0}$. |
| [04_control_drum_mechanism.png](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/04_control_drum_mechanism.png) | Meccanica del tamburo di controllo BeO con settore $\text{B}_4\text{C}$ ($155^\circ$, $1.8\text{ cm}$): orientamento $0^\circ$ OUT vs $180^\circ$ IN. |
| [05_openmc_csg_slice_xy.png](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/05_openmc_csg_slice_xy.png) | Sezione X-Y generata direttamente dal ray-tracer OpenMC CSG (geometria combinatoria costruttiva reale). |
| [06_openmc_csg_slice_xz.png](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/06_openmc_csg_slice_xz.png) | Sezione X-Z generata dal ray-tracer OpenMC CSG. |
| [07_depletion_keff_and_excess_reactivity.png](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/07_depletion_keff_and_excess_reactivity.png) | Evoluzione di $k_{\text{eff}}$ e della reattività in eccesso nei 25 anni di vita utile con incertezze Monte Carlo $2\sigma$. |
| [08_depletion_control_drum_critical_angle.png](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/08_depletion_control_drum_critical_angle.png) | Traiettoria di estrazione dell'angolo critico $\theta_{\text{crit}}$ dei tamburi per mantenere $k_{\text{eff}} = 1.0000$ costante nel tempo. |
| [09_depletion_isotopic_inventory_evolution.png](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/09_depletion_isotopic_inventory_evolution.png) | Evoluzione temporale degli isotopi: consumo di U-235 (da $27.01$ a $23.36\text{ kg}$) e breeding di Pu-239 ($642.4\text{ g}$). |
| [10_depletion_burnup_and_fima.png](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/10_depletion_burnup_and_fima.png) | Burnup specifico del nocciolo ($26.69\text{ MWd/kgHM}$) e percentuale FIMA. |
| [11_moose_cardinal_temperature_profiles_pins.png](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/11_moose_cardinal_temperature_profiles_pins.png) | Benchmark FEM 2D MOOSE/Cardinal: profili termici del pin più caldo, centrale, medio e periferico. |
| [12_moose_cardinal_moderator_thermal_field.png](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/12_moose_cardinal_moderator_thermal_field.png) | Mappatura bidimensionale agli elementi finiti nel blocco di moderatore $\text{YH}_{2.0}$ ($T_{\max} = 935.31\text{ K} < 950\text{ K}$). |
| [13_closed_brayton_cycle_ts_diagram.png](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/figures/13_closed_brayton_cycle_ts_diagram.png) | Diagramma Temperatura-Entropia ($T-s$) del ciclo Brayton chiuso ad elio con recuperatore ($P_{\text{elec}} = 136.39\text{ kWe}$, $\eta = 34.10\%$). |

---

## 5. Dati di Burnup e Lifespan (26 Step Annuali)

Estratto dei dati salienti salvati in [`data/c01_criticality_depletion_1yr_steps.csv`](file:///home/aleri/OpenMC_projects/first_test/definitive_final_geometry/data/c01_criticality_depletion_1yr_steps.csv):

| Tempo [Anni] | $k_{\text{eff}}$ ($0^\circ$ OUT) | Reattività in Eccesso [pcm] | Angolo Tamburi $\theta_{\text{crit}}$ | $k_{\text{eff}}$ Operativo | Massa $^{235}\text{U}$ [kg] | Massa $^{239}\text{Pu}$ [g] | Burnup [MWd/kgHM] |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **0.0** | $1.06139$ | $+5783.6$ | $142.40^\circ$ | $0.99590$ | $27.010$ | $0.00$ | $0.00$ |
| **1.0** | $1.06174$ | $+5814.9$ | $129.44^\circ$ | $1.00103$ | $26.851$ | $34.06$ | $1.07$ |
| **5.0** | $1.05932$ | $+5599.6$ | $130.37^\circ$ | $1.00019$ | $26.250$ | $163.34$ | $5.34$ |
| **10.0** | $1.04920$ | $+4689.6$ | $105.68^\circ$ | $1.00592$ | $25.511$ | $307.01$ | $10.68$ |
| **15.0** | $1.04348$ | $+4166.8$ | $97.63^\circ$ | $0.99728$ | $24.784$ | $433.44$ | $16.01$ |
| **20.0** | $1.03601$ | $+3476.0$ | $74.12^\circ$ | $1.00401$ | $24.067$ | $545.22$ | $21.35$ |
| **25.0** | **$1.02627$** | **$+2559.4$** | **$71.09^\circ$** | **$0.99878$** | **$23.360$** | **$642.38$** | **$26.69$** |

- **Vita utile a 25.0 anni**: Riserva residua pari a ben $+2559.4\text{ pcm}$ con tamburi a $71.09^\circ$.
- **Vita utile estrapolata a zero-crossing ($\theta = 0^\circ$)**: **$44.8\text{ anni}$** di funzionamento continuativo a $400\text{ kWth}$.

---

## 6. Come Rieseguire le Simulazioni

Tutti i codici sono configurati per funzionare nell'ambiente Python scientifico del progetto (`openmc-env`).

### 6.1 Generazione Rendering Geometrici ad Alta Precisione (300 DPI)
```bash
/home/aleri/miniforge3/envs/openmc-env/bin/python \
  /home/aleri/OpenMC_projects/first_test/definitive_final_geometry/neutronics_openmc/generate_high_precision_plots_c01.py
```

### 6.2 Risoluzione Termica FEM MOOSE/Cardinal e Ciclo Brayton
```bash
/home/aleri/miniforge3/envs/openmc-env/bin/python \
  /home/aleri/OpenMC_projects/first_test/definitive_final_geometry/thermohydraulics_moose_cardinal/run_thermodynamics_and_moose_cardinal_c01.py
```

### 6.3 Campagna di Burnup e Controllo Tamburi (1-Year Steps)
```bash
/home/aleri/miniforge3/envs/openmc-env/bin/python \
  /home/aleri/OpenMC_projects/first_test/definitive_final_geometry/depletion_analysis/run_1yr_criticality_depletion.py
```

---
*Archivio generato e validato l'8 Ottobre 2026.*

