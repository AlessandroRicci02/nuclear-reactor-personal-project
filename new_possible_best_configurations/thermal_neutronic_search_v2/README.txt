CAMPAGNA V2: 30 REATTORI NEP CON MODERATORE YH2 TERMICAMENTE SEPARATO
Creata il 2026-10-08. Proposte concettuali, non progetti qualificati.

REQUISITI UTENTE
  He ingresso 850--900 K, uscita miscelata 1100--1150 K.
  k_eff BOL OUT >1; k_eff BOL IN <1.
  Massimo locale YH2 <950 K.
  Coefficienti del moderatore, combustibile e totale negativi con incertezze.
  Nessun tetto rigido di 650 kg. Probabilita' di vita >10 anni da valutare
  tramite inventario e burnup richiesto; nessuna depletion in questa campagna.
Potenza conservata dai tentativi precedenti: 400 kWth; He puro a 3 MPa.

FILE E AVVIO
  configurations.json: catalogo completo di 30 candidati, parametri e motivi.
  configurations.csv: tabella leggibile dei principali parametri.
  scripts/auto_search.py: controller con checkpoint e lock anti-doppio avvio.
  campaign/<id>/: input, log, statepoint e motivo di scarto di ogni candidato.
  status.json / progress.csv / controller.log: stato e avanzamento.
  found_reactor/<id>/: report, configurazione, risultati e figure del primo
  candidato che supera tutti i filtri. Nessun altro trasporto dopo questo punto.
  exhausted_no_solution: tutte le 30 proposte esaurite senza una soluzione.
  numerical_failure: errore numerico/software; non viene etichettato come
  fallimento fisico. Il controller si ferma per evitare 30 scarti artificiali.

Comando da questa directory:
  /home/aleri/miniforge3/envs/openmc-env/bin/python -u scripts/auto_search.py
Il riavvio riutilizza solo risultati con hash identico di codice e richiesta.

LEZIONI DEL WORKSPACE
Il candidato matrix_short_strong_be e' SCARTATO: il limite inferiore fra le
13 stime termiche precedenti era 1242.6 K, rispetto a 950 K. La sua precedente
"qualifica" era soltanto neutronica. Non e' utilizzabile per questi requisiti.
La depletion vecchia e' fallita dopo il primo intervallo: nessuna durata
verificata e nessuna simulazione vecchia attualmente attiva.
TB35 e C3 hanno mostrato coefficienti del moderatore positivi. Il concetto
second_try aveva separazione termica promettente, ma combustibile/inventario
e geometria insufficienti per la criticita'. Non trasferiamo i suoi risultati.

SCELTA DELLE 30 CONFIGURAZIONI
Cinque famiglie, sei varianti ciascuna:
  1) 271 pin, manicotto intermedio, raffreddamento freddo su due lati;
  2) 331 pin, fuel piu' grande e manicotto sottile, raffreddamento interno;
  3) 217 pin, manicotto spesso e doppio raffreddamento;
  4) 271 pin, nocciolo allungato con raffreddamento interno;
  5) 271 pin con riflettore BeO e doppio raffreddamento.
Variamo altezza, spessore del manicotto, temperature del ciclo e dimensione/
spessore dei tamburi B4C. Le varianti 3 e 6 aggiungono un anello di pin.
Ogni famiglia ha un caso UN con azoto naturale, uno al 95% N15 e quattro al
99.5% N15 atomico. Densita' operative assunte: UN 13.1 e YH2 4.20 g/cm3,
compatibili con una riduzione rispetto ai valori freddi dei vecchi modelli;
da confermare sul materiale reale. Le dimensioni tabulate sono nominali a
caldo, non quote di fabbricazione a temperatura ambiente. Composizione e
stoichiometria U:N=1:1 sono normalizzate con le masse isotopiche reali.
N15 riduce le catture parassite rispetto a N14; non e' una garanzia di criticita'
o feedback negativo. Richiede fabbricazione e approvvigionamento specifici.
Il pitch viene derivato dai veri strati, con 0.8 mm
di separazione minima fra pin. Tamburi 12, settore 175 gradi, B10 95%.
Il catalogo esplicita tutte le dimensioni e l'inventario fissile.
Tutte le proposte usano YH2.0. La libreria locale ENDF/B-VIII.0 contiene
S(alpha,beta) per YH2; non viene sostituita con ZrH ne' riutilizzata tacitamente
per una composizione YHx arbitraria. Una futura campagna x<2 richiede TSL e
proprieta' coerenti o una specifica valutazione dell'errore di tale approssimazione.

ARCHITETTURA TERMICA
L'elio entra dall'estremita' calda esterna e percorre prima i canali freddi
lungo i manicotti YH2; all'estremita' opposta torna nel canale caldo attorno
al fuel e raggiunge la temperatura finale richiesta. Nel doppio raffreddamento
i due canali freddi si ricombinano prima del passaggio sul fuel.
Camera evacuata con foglio Mo fra i passaggi caldo e freddo. La geometria
neutronica include pareti Nb1Zr, foil Mo, barriera Al2O3 proposta e clearance.
La camera isolante e i canali proseguono attraverso entrambi i riflettori
assiali; la deposizione nelle pareti assiali e la loro area radiativa sono
incluse nel carico termico, non sostituite da una continuita' ideale in Be.
Il gioco esterno del manicotto e' riempito di He; la sua resistenza termica
limita esplicitamente il credito del canale esterno. Non e' un contatto ideale.
Distribuzione di portata fra pin con orifizi: 90% proporzionale alla potenza
integrata del pin e 10% uniforme, sensibilita' conservativa 70%/30%.
Gli orifizi/collettori e i dettagli meccanici sono requisiti di progetto ancora
da sviluppare. Non si assegna una temperatura fredda al solido per decreto.

ORDINE DEL CONTROLLER / SCARTO IMMEDIATO
  1. OUT e poi IN BOL preliminari, con raffinamento se il segno e' incerto.
     Fallimento di uno: scarta subito, niente calcoli termici/coefficienti.
  2. Angolo critico intermedio; mappa OpenMC neutroni+fotoni del calore.
  3. Bilanci di energia dei canali seriali, radiazione sul foil Mo, conduzione
     nel manicotto e Robin convettivo. Verifica Cardinal RZ del pin piu' caldo.
     Temperatura superiore non sotto 950 K: scarta subito.
  4. Ripete trasporto/termica e OUT/IN alle temperature risultanti fino alla
     coerenza del campo (<8 K He/YH2, <5% fuel, <3 gradi tamburi), massimo
     tre iterazioni; mancata convergenza e' un risultato irrisolto.
  5. Coefficienti differenziali fuel, moderatore e totale all'angolo critico
     e al campo operativo; interrompe al primo coefficiente che non passa.
  6. Conferma indipendente dei tre coefficienti a passo termico dimezzato.
     Primo successo: ferma il loop e prepara found_reactor/<id>/.

INCERTEZZE / CRITERI DI ACCETTAZIONE
Eigenvalue: k_OUT -4*sigma >1 e k_IN +4*sigma <1. Non basta il valore centrale.
Coefficienti: alpha +4*sigma_alpha <0. Reattivita' rho=1-1/k; propagazione
dei due trasporti indipendenti; semi distinti per ogni caso e conferma.
Raffinamento automatico sigma_k fino a 6e-5 per coefficienti difficili, con
limite 12000 batches x 20000 particelle per i raffinamenti piu' fini. Se il segno resta incerto, il candidato
e' scartato come IRRISOLTO, non dichiarato fisicamente positivo.
Passi principali: +/-50 K moderatore/totale, +/-100 K fuel; conferma +/-25
e +/-50 K. I coefficienti sono feedback locali a tamburi bloccati.
Espansione massa-conservativa: densita' inversa al volume, fuel che si espande
nel suo gioco e manicotto che si espande radialmente nella clearance, con
raggio interno e lunghezza assunti vincolati dal contenitore. Beta volumetrico
YH2 2.4--4.95e-5/K, UN 2.1--3.9e-5/K; moderatore/fuel a tre livelli, totale
ai quattro vertici delle due bande. Queste bande sono ipotesi di screening,
da rimpiazzare con caratterizzazione del lotto e analisi termomeccanica.
Altri solidi: temperatura variabile, geometria fissata. Non e' un modello
completo della dilatazione dell'assemblaggio o del comportamento in transitorio.

Termica: tally medio nominale; scenario conservativo con sorgenti fredde
media+3*sigma per bin e ulteriore 5%. Una mesh fotonica 24x24x10 localizza
il calore del riflettore: nominalmente al manicotto piu' vicino, nello scenario
conservativo la potenza superiore di ogni bin e' assegnata interamente a tutti
e tre i manicotti piu' vicini (sovrastima intenzionale della potenza fredda).
kYH2=5 vs 10 W/mK, film 70%, emissivita' 0.30 vs 0.15,
ponti/dispersioni 35 vs 12 W/K, parete calda +50 vs +10 K, redistribuzione
di portata. La temperatura per-slice usa il massimo bulk ai due estremi e
il massimo radiale analitico; Cardinal risolve anche la conduzione assiale.
Accettazione: massimo fra limiti analitici/Cardinal +5 K <950 K.
Sono bande ingegneristiche esplicite, NON un intervallo probabilistico completo.
La simmetria rotazionale di ordine sei riduce il rumore della mappa, conservando
la media delle sigma originarie come limite prudente senza assumere covarianze.
Incertezze di dati nucleari, degrado per irraggiamento, integrita' della camera
evacuata, compatibilita' canning/barriera, stress e CFD non sono quantificate.
Un successo permette ulteriori studi; non certifica che un regolatore approvi.

MASSA E POSSIBILE DURATA
L'inventario fuel piu' grande e la minore potenza specifica motivano lo studio
di una vita >10 anni. Il consumo energetico fissile e il burnup in 10 anni
sono soltanto controlli di inventario: veleni, margine di reattivita', deriva
dei coefficienti e burnup del B4C richiederanno depletion in una fase successiva.
La massa tabulata e' un limite geometrico preliminare dei solidi modellati;
schermo per equipaggio, attuatori, vessel completo, collettori, turbomacchine,
radiatori e strutture esterne sono esclusi. Nessun candidato passa/fallisce
per massa nella nuova campagna. Non si confonde massa del core con sistema NEP.

LETTERATURA PRIMARIA CONSULTATA E USO
[1] Wright et al., NASA/Sandia, Direct-Drive Gas-Cooled Reactor Power System:
    Concept and Preliminary Testing (2003).
    https://ntrs.nasa.gov/citations/20030061231
    UN a pin e Brayton diretto He/Xe, gas 1144 K e cavita' isolanti per superfici
    sotto 900 K: precedente per separare temperature di gas e strutture.
    Questa campagna He puro + YH2 e' una nostra derivazione, non quel reattore.
[2] Peakman & Gregg (2020), JANUS space reactor concept.
    https://arxiv.org/abs/2004.12116
    Gas-cooled NEP/Brayton e compromessi fra fuel, pressione, flusso e massa.
    Non viene assunto come benchmark neutronicamente equivalente.
[3] Ade et al. (2022), Reactor Physics Considerations for Use of Yttrium Hydride.
    https://doi.org/10.1080/00295639.2022.2035157
    https://impact.ornl.gov/en/publications/reactor-physics-considerations-for-use-of-yttrium-hydride-moderat/
    TCR YH1.85 e rischio di coefficiente positivo; motivo per simulare ogni
    coefficiente e non usare composizione/idruro come garanzia del segno.
[4] Chapman et al. (2021), Thermal neutron scattering measurements and modeling.
    https://www.osti.gov/biblio/1775234
    Scattering di YHx a temperature elevate e limiti/effetti anarmonici nelle
    valutazioni: motivo per tracciare la TSL e dichiarare il bias non quantificato.
[5] ORNL/TCR, Handbook on the Material Properties of Yttrium Hydride (2021).
    https://info.ornl.gov/sites/publications/Files/Pub160401.pdf
    Proprieta' termofisiche, composizione, stabilita' e dati post-irraggiamento
    ancora limitati: motivo per le bande conservative e futura caratterizzazione.
[6] NASA/DOE, Nuclear Power Concepts and Development Strategies for High-Power
    Electric Propulsion Missions to Mars (2021).
    https://ntrs.nasa.gov/citations/20210016968
    Il sistema NEP comprende reattore, schermatura, conversione e rigetto del
    calore: la nostra massa geometrica non rappresenta l'intero sistema.
[7] KTH, Uranium Nitride Fuel — ricerca primaria sul combustibile e fabbricazione.
    https://www.reactor.sci.kth.se/une/about-the-project-1.793757
    N14 assorbitore forte nello spettro termico; UN con N15 richiede riduzione
    dei costi di fabbricazione. Motivo per variare isotopia dell'azoto anziche'
    soltanto aggiungere massa o cambiare tacitamente l'arricchimento dell'uranio.
