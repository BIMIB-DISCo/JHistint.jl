# JHistint.jl - Julia Histopathology Interface

Interfaccia Julia per implementazione delle REST API disponibili sul portale CDSA (Cancer Slide Digital Archive) per il download di immagini istologiche reperibili nel TCGA (The Cancer Genome Atlas). Il Cancer Slide Digital Archive (CDSA) è una piattaforma web per il supporto, la condivisione e l'analisi di dati patologici digitali. Attualmente ospita oltre 23.000 immagini associate ai dati disponibili nel «The Cancer Genome Atlas» Data Portal.  

Link d'accesso al CDSA: [Clicca qui](https://api.digitalslidearchive.org/#collections)    

Link repository contenente i dati mappati nel portale: [Clicca qui](https://cancer.digitalslidearchive.org/#!/CDSA/acc/TCGA-OR-A5J1)

Link guida all'utilizzo delle API: [Clicca qui](https://api.digitalslidearchive.org/api/v1)

## Struttura del Package
* Le folder `case` e `collection` memorizzano in formato `.json` i metadati relativi ai singoli casi e alle collezioni disponibili nel TCGA Data Portal. La folder `collection` è composta come segue:  
  * `collectionlist.jsn` = Memorizza i dati d'accesso (metadati) delle collezioni (Project in TCGA).  
  * `colletion_name.jsn` = Memorizza i dati d'accesso (metadati) relativi alla singola collezione. Il file `.json` viene generato in base alla collezione scelta dall'utente.
* La folder `case` è composta come segue:
  * `collection_name.jsn` = Memorizza tutti i metadati relativi ai casi associati alla collezione selezionata dall'utente.  
* La folder `slides` memorizza le immagini istologiche relative ai singoli casi. Le immagini sono organizzate in base al formato (`.svs`), alla collezione (`TCGA-chol`, `TCGA-esca`, etc.) e al singolo caso da analizzare (`TCGA-2H-A9GF`, `TCGA-2H-A9GG`, etc.). All'interno di ogni folder relativa al caso sono memorizzate le slide compresse in file `.zip`. Il formato delle singole slide è .svs. La denominazione delle folder inerenti ai casi coincide con i valori del campo `Case ID` riportato nel TCGA Data Portal. La denominazione dei file `.zip` situati in ogni folder fa riferimento all'attributo `Sample ID` associato al paziente. La denominazione della slide è data dalla concatenazione degli attributi `Slide ID` e `Slide UUID` consultabili nella sezione inferiore della pagina web dedicata al caso generico `TCGA-XX-YYYY`.

```
Esempio: TCGA-02-0001-01C-01-TS1.zip  
  - 02 = si riferisce al TSS (Tissue Source Site).  
  - 0001 = si riferisce al codice associato al Participant, stringa alfanumerica.  
  - 01 = si riferisce al Sample Type. I valori associati ai campioni aventi tumori sono nell'intervallo 01 - 09. 10 - 19 indica l'intervallo dedicato a campioni normali non malati. 20 - 29 indica campioni attualmente sotto controllo.  
  - C = si riferisce al campo Vial relativo all'ordinamento del campione nella sequenza di campioni. I valori variano tra A - Z.  
  - 01 = si riferisce al campo Portion relativo all'ordinamento delle porzioni analizzate associate ad un campione. Assume valori nell'intervallo 01-99.  
  - TS1 = si riferisce al campo Slide relativo al tipo di immagine. I valori assumbili sono TS (Top Slide), BS (Bottom Slide) e MS (Middle Slide). Il valore alfanumerico indica l'ordinamento della slide.
```

## Collezioni JHistint 
Le collezioni disponibili sono:  
  * TCGA-BRCA = Breast Invasive Carcinoma (Breast)
  * TCGA-OV = Ovarian Serous Cystadenocarcinoma (Ovary)
  * TCGA-LUAD = Lung Adenocarcinoma (Bronchus and Lung)
  * TCGA-UCEC = Uterine Corpus Endometrial Carcinoma (Corpus uteri)
  * TCGA-GBM = Glioblastoma Multiforme (Brain)
  * TCGA-HSNC = Head and Neck Squamous Cell Carcinoma (Larynx, Lip, Tonsil, Gum, Other and unspecified parths of mouth)
  * TCGA-KIRC = Kidney Renal Clear Cell Carcinoma (Kidney)
  * TCGA-LGG = Brain Lower Grade Glioma (Brain)
  * TCGA-LUSC = Lung Squamous Cell Carcinoma (Bronchus and lung)
  * TCGA-TCHA = Thyroid Carcinoma (Thyroid gland)
  * TCGA-PRAD = Prostate Adenocarcinoma (Prostate gland)
  * TCGA-SKCM = Skin Cutaneous Melanoma (Skin)
  * TCGA-COAD = Colon Adenocarcinoma (Colon)
  * TCGA-STAD = Stomach Adenocarcinoma (Stomach)
  * TCGA-BLCA = Bladder Urothelial Carcinoma (Bladder)
  * TCGA-LIHC = Liver Hepatocellular Carcinoma (Liver and intrahepatic bile ducts)  
  * TCGA-CESC = Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma (Cervix uteri)
  * TCGA-KIRP = Kidney Renal Papillary Cell Carcinoma (Kidney)
  * TCGA-SARC = Sarcoma (Various)
  * TCGA-ESCA = Esophageal Carcinoma (Esophagus)
  * TCGA-PAAD = Pancreatic Adenocarcinoma (Pancreas)
  * TCGA-READ = Rectum Adenocarcinoma (Rectum)
  * TCGA-PCPG = Pheochromocytoma and Paraganglioma (Adrenal gland)
  * TCGA-TGCT = Testicular Germ Cell Tumors (Testis)
  * TCGA-THYM = Thymoma (Thymus)
  * TCGA-ACC = Adrenocortical Carcinoma -Adenomas and Adenocarcinomas (Adrenal gland)
  * TCGA-MESO = Mesothelioma (Heart, mediastinum and pleura)
  * TCGA-UVM = Uveal Melanoma (Eye and adnexa)
  * TCGA-KICH = Kidney Chromophobe (Kidney)
  * TCGA-UCS = Uterine Carcinosarcoma (Uterus, NOS)
  * TCGA-CHOL = Cholangiocarcinoma (Liver and intrahepatic bile ducts, Other and unspecified part of biliary track)
  * TCGA-DLBC = Lymphoid Neoplasm Diffuse Large B-cell Lymphoma (Various)
Per il download di una collezione specifica è sufficiente indicare il nome della collezione: `BRCA`, `OV`, `LUAD`.

## Installazione del Package
Il package `JHistint` è disponibile nei Julia Registries, quindi installabile come segue:
```
julia > using Pkg
julia > Pkg.add("JHistint")
julia > using JHistint
```
In alternativa, digitare `]` nel Julia REPL ed eseguire:
```
(@v1.8) pkg > add JHistint
(@v1.8) pkg > using JHistint
``` 
## Funzioni Download Slides (JHistint.jl)

```@docs
download_single_collection(collection_name::AbstractString)
```

```@docs
download_all_collection()
```

## Funzioni di supporto API (functions.jl)
```@docs
download_collection_values(filepath::AbstractString)
```

```@docs
extract_collection_values(filepath::AbstractString)
```

```@docs
download_project_infos(filepath::AbstractString, collection_name::AbstractString)
```

```@docs
extract_project_id(filepath::AbstractString)
```

```@docs
getCasesForProject(filepath_case::AbstractString, project_id::AbstractString)
```

```@docs
download_zip(link::AbstractString, filepath::AbstractString)
```

