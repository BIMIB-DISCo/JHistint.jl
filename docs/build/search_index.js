var documenterSearchIndex = {"docs":
[{"location":"#JHistint.jl-Julia-Histopathology-Interface","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"","category":"section"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"Interfaccia Julia per implementazione delle REST API disponibili sul portale CDSA (Cancer Slide Digital Archive) per il download di immagini istologiche reperibili nel TCGA (The Cancer Genome Atlas). Il Cancer Slide Digital Archive (CDSA) è una piattaforma web per il supporto, la condivisione e l'analisi di dati patologici digitali. Attualmente ospita oltre 23.000 immagini associate ai dati disponibili nel «The Cancer Genome Atlas» Data Portal.  ","category":"page"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"Link GitHub repository: JHistint.jl","category":"page"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"Link d'accesso al CDSA: Clicca qui    ","category":"page"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"Link repository contenente i dati mappati nel portale: Clicca qui","category":"page"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"Link guida all'utilizzo delle API: Clicca qui","category":"page"},{"location":"#Struttura-del-Package","page":"JHistint.jl - Julia Histopathology Interface","title":"Struttura del Package","text":"","category":"section"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"Le folder case e collection memorizzano in formato .json i metadati relativi ai singoli casi e alle collezioni disponibili nel TCGA Data Portal. La folder collection è composta come segue:  \ncollectionlist.jsn = Memorizza i dati d'accesso (metadati) delle collezioni (Project in TCGA).  \ncolletion_name.jsn = Memorizza i dati d'accesso (metadati) relativi alla singola collezione. Il file .json viene generato in base alla collezione scelta dall'utente.\nLa folder case è composta come segue:\ncollection_name.jsn = Memorizza tutti i metadati relativi ai casi associati alla collezione selezionata dall'utente.  \nLa folder slides memorizza le immagini istologiche relative ai singoli casi. Le immagini sono organizzate in base al formato (.svs), alla collezione (TCGA-chol, TCGA-esca, etc.) e al singolo caso da analizzare (TCGA-2H-A9GF, TCGA-2H-A9GG, etc.). All'interno di ogni folder relativa al caso sono memorizzate le slide compresse in file .zip. Il formato delle singole slide è .svs. La denominazione delle folder inerenti ai casi coincide con i valori del campo Case ID riportato nel TCGA Data Portal. La denominazione dei file .zip situati in ogni folder fa riferimento all'attributo Sample ID associato al paziente. La denominazione della slide è data dalla concatenazione degli attributi Slide ID e Slide UUID consultabili nella sezione inferiore della pagina web dedicata al caso generico TCGA-XX-YYYY.","category":"page"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"Esempio: TCGA-02-0001-01C-01-TS1.zip  \n  - 02 = si riferisce al TSS (Tissue Source Site).  \n  - 0001 = si riferisce al codice associato al Participant, stringa alfanumerica.  \n  - 01 = si riferisce al Sample Type. I valori associati ai campioni aventi tumori sono nell'intervallo 01 - 09. 10 - 19 indica l'intervallo dedicato a campioni normali non malati. 20 - 29 indica campioni attualmente sotto controllo.  \n  - C = si riferisce al campo Vial relativo all'ordinamento del campione nella sequenza di campioni. I valori variano tra A - Z.  \n  - 01 = si riferisce al campo Portion relativo all'ordinamento delle porzioni analizzate associate ad un campione. Assume valori nell'intervallo 01-99.  \n  - TS1 = si riferisce al campo Slide relativo al tipo di immagine. I valori assumbili sono TS (Top Slide), BS (Bottom Slide) e MS (Middle Slide). Il valore alfanumerico indica l'ordinamento della slide.","category":"page"},{"location":"#Collezioni-JHistint","page":"JHistint.jl - Julia Histopathology Interface","title":"Collezioni JHistint","text":"","category":"section"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"Le collezioni disponibili sono:  ","category":"page"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"TCGA-BRCA = Breast Invasive Carcinoma (Breast)\nTCGA-OV = Ovarian Serous Cystadenocarcinoma (Ovary)\nTCGA-LUAD = Lung Adenocarcinoma (Bronchus and Lung)\nTCGA-UCEC = Uterine Corpus Endometrial Carcinoma (Corpus uteri)\nTCGA-GBM = Glioblastoma Multiforme (Brain)\nTCGA-HSNC = Head and Neck Squamous Cell Carcinoma (Larynx, Lip, Tonsil, Gum, Other and unspecified parths of mouth)\nTCGA-KIRC = Kidney Renal Clear Cell Carcinoma (Kidney)\nTCGA-LGG = Brain Lower Grade Glioma (Brain)\nTCGA-LUSC = Lung Squamous Cell Carcinoma (Bronchus and lung)\nTCGA-TCHA = Thyroid Carcinoma (Thyroid gland)\nTCGA-PRAD = Prostate Adenocarcinoma (Prostate gland)\nTCGA-SKCM = Skin Cutaneous Melanoma (Skin)\nTCGA-COAD = Colon Adenocarcinoma (Colon)\nTCGA-STAD = Stomach Adenocarcinoma (Stomach)\nTCGA-BLCA = Bladder Urothelial Carcinoma (Bladder)\nTCGA-LIHC = Liver Hepatocellular Carcinoma (Liver and intrahepatic bile ducts)  \nTCGA-CESC = Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma (Cervix uteri)\nTCGA-KIRP = Kidney Renal Papillary Cell Carcinoma (Kidney)\nTCGA-SARC = Sarcoma (Various)\nTCGA-ESCA = Esophageal Carcinoma (Esophagus)\nTCGA-PAAD = Pancreatic Adenocarcinoma (Pancreas)\nTCGA-READ = Rectum Adenocarcinoma (Rectum)\nTCGA-PCPG = Pheochromocytoma and Paraganglioma (Adrenal gland)\nTCGA-TGCT = Testicular Germ Cell Tumors (Testis)\nTCGA-THYM = Thymoma (Thymus)\nTCGA-ACC = Adrenocortical Carcinoma -Adenomas and Adenocarcinomas (Adrenal gland)\nTCGA-MESO = Mesothelioma (Heart, mediastinum and pleura)\nTCGA-UVM = Uveal Melanoma (Eye and adnexa)\nTCGA-KICH = Kidney Chromophobe (Kidney)\nTCGA-UCS = Uterine Carcinosarcoma (Uterus, NOS)\nTCGA-CHOL = Cholangiocarcinoma (Liver and intrahepatic bile ducts, Other and unspecified part of biliary track)\nTCGA-DLBC = Lymphoid Neoplasm Diffuse Large B-cell Lymphoma (Various)","category":"page"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"Per il download di una collezione specifica è sufficiente indicare il nome della collezione: BRCA, OV, LUAD.","category":"page"},{"location":"#Installazione-del-Package","page":"JHistint.jl - Julia Histopathology Interface","title":"Installazione del Package","text":"","category":"section"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"Il package JHistint è disponibile nei Julia Registries, quindi installabile come segue:","category":"page"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"julia > using Pkg\njulia > Pkg.add(\"JHistint\")\njulia > using JHistint","category":"page"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"In alternativa, digitare ] nel Julia REPL ed eseguire:","category":"page"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"(@v1.8) pkg > add JHistint\n(@v1.8) pkg > using JHistint","category":"page"},{"location":"#Funzioni-Download-Slides-(JHistint.jl)","page":"JHistint.jl - Julia Histopathology Interface","title":"Funzioni Download Slides (JHistint.jl)","text":"","category":"section"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"download_single_collection(collection_name::AbstractString)","category":"page"},{"location":"#JHistint.download_single_collection-Tuple{AbstractString}","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.download_single_collection","text":"download_single_collection(collection_name::AbstractString)\n\nFunzione per il download delle slides istologiche associate ad una collezione disponibile nel TCGA.\n\nArgomenti\n\ncollection_name::AbstractString = Collezione di dati TCGA di cui scaricare le slides istologiche.\n\nNote\n\nLa funzione valuta l'argomento collection_name, in caso di collezione non valida considera la configurazione del file Config.toml. Il valore impostato nel package è default.\n\n# Esempi con input validi\njulia> JHistint.download_single_collection(\"acc\")\njulia> JHistint.download_single_collection(\"bLca\")\n\n# Esempi con input non validi\njulia> JHistint.download_single_collection(\"ac\")\njulia> JHistint.download_single_collection(\"\")\n\n\n\n\n\n","category":"method"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"download_all_collection()","category":"page"},{"location":"#JHistint.download_all_collection-Tuple{}","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.download_all_collection","text":"download_all_collection()\n\nFunzione per il download delle slides istologiche associate a tutte le collezioni disponbile nel TCGA.\n\n# Esempi con input validi\njulia> JHistint.download_all_collection()\n\n\n\n\n\n","category":"method"},{"location":"#Funzioni-Segmentazione-Slides-(JHistint.jl)","page":"JHistint.jl - Julia Histopathology Interface","title":"Funzioni Segmentazione Slides (JHistint.jl)","text":"","category":"section"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"slide_cell_segmentation(collection_name::AbstractString)","category":"page"},{"location":"#JHistint.slide_cell_segmentation-Tuple{AbstractString}","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.slide_cell_segmentation","text":"slide_cell_segmentation(collection_name::AbstractString)\n\nFunzione per esecuzione della segmentazione cellulare delle slide istopatologiche presenti nel database JHistint_DB associate al nome della collezione fornita come argomento.\n\nArgomenti\n\ncollection_name::AbstractString = Collezione di slide TCGA di cui effettuare la segmentazione cellulare.\n\nNote\n\nPer ogni slide nel DB viene eseguita la segmentazione delle cellule utilizzando la funzione apply_segmentation e il risultato viene salvato nel DB utilizzando la funzione load_seg_slide. Infine, viene stampato un messaggio di conferma per ogni slide segmentata.\n\n# Esempi con input validi\njulia> JHistint.slide_cell_segmentation(\"acc\")\njulia> JHistint.slide_cell_segmentation(\"bLca\")\n\n# Esempi con input non validi\njulia> JHistint.slide_cell_segmentation(\"ac\")\njulia> JHistint.slide_cell_segmentation(\"\")\n\n\n\n\n\n","category":"method"},{"location":"#Funzioni-di-supporto-per-Segmentazione-(segmentationManager.jl)","page":"JHistint.jl - Julia Histopathology Interface","title":"Funzioni di supporto per Segmentazione (segmentationManager.jl)","text":"","category":"section"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"apply_segmentation(slide_info::Tuple{String, Array{ColorTypes.RGB{FixedPointNumbers.N0f8}, 3}, String})","category":"page"},{"location":"#JHistint.apply_segmentation-Tuple{Tuple{String, Array{RGB{N0f8}, 3}, String}}","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.apply_segmentation","text":"apply_segmentation(slide_info::Tuple{String, Array{RGB{}, 3}, String})\n\nLa funzione esegue la segmentazione di un'immagine istologica.\n\nArgomenti\n\nslide_info::Tuple{String, Array{ColorTypes.RGB{FixedPointNumbers.N0f8}, 3}, String}: Tupla contenente l'id della slide, l'immagine stessa ottenuto dal DB e il percorso del file dell'immagine originale.\n\nValori di ritorno\n\nfilepath_seg: Il percorso all'interno del quale è memorizzata la slide segmentata.\nsegmented_slide: La slide istologica segmentata in formato .tif.\n\nNote\n\nLa funzione utilizza l'algoritmo di segmentazione watershed per segmentare l'immagine in diversi raggruppamenti di pixel. La segmentazione viene eseguita utilizzando una trasformazione delle caratteristiche dell'immagine (feature_transform) e l'etichettatura dei componenti connessi. Viene quindi calcolata la distanza tra le diverse regioni e viene costruito un grafo di adiacenza delle regioni, utilizzando la funzione region_adjacency_graph. Viene inoltre assegnato un colore casuale a ciascuna regione segmentata. Infine, viene salvata un'immagine .tif segmentata e viene restituito il percorso del file della slide segmentata e la slide stessa.\n\n\n\n\n\n","category":"method"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"get_random_color(seed)","category":"page"},{"location":"#JHistint.get_random_color-Tuple{Any}","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.get_random_color","text":"get_random_color(seed)\n\nFunzione per restituire un colore casuale in formato RGB a 8 bit, utilizzando un seme specificato.\n\nArgomenti\n\nseed = Un intero utilizzato per inizializzare il generatore di numeri casuali. Se due chiamate alla funzione utilizzano lo stesso seme, verrà generato lo stesso colore.\n\nValore di ritorno\n\nLa funzione restituisce un colore casuale in formato RGB a 8 bit.\n\n\n\n\n\n","category":"method"},{"location":"#Funzioni-di-supporto-DB-(dbManager.jl)","page":"JHistint.jl - Julia Histopathology Interface","title":"Funzioni di supporto DB (dbManager.jl)","text":"","category":"section"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"insert_record_DB(col_name::AbstractString,\n                        cas_name::AbstractString,\n                        tcga_case_id::AbstractString,\n                        sin_cas_name::AbstractString,\n                        tcga_slide_id::AbstractString,\n                        link_slide::AbstractString,\n                        filepath_zip::AbstractString,\n                        filepath_svs::AbstractString)","category":"page"},{"location":"#JHistint.insert_record_DB-NTuple{8, AbstractString}","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.insert_record_DB","text":"insert_record_DB(col_name::AbstractString,\n                      cas_name::AbstractString,\n                      tcga_case_id::AbstractString,\n                      sin_cas_name::AbstractString,\n                      tcga_slide_id::AbstractString,\n                      link_slide::AbstractString,\n                      filepath_zip::AbstractString,\n                      filepath_svs::AbstractString)\n\nFunzione per la memorizzazione nel database JHistint_DB delle informazioni associate ad ogni slide scaricata dal CDSA.\n\nArgomenti\n\ncol_name::AbstractString = Nome della collezione.\ncas_name::AbstractString = Nome del caso, coincide con il CASE-NAME stampato a video dal package.\ntcga_case_id::AbstractString = ID utilizzato dal TCGA per identificare il caso.\nsin_cas_name::AbstractString = Nome del singolo caso, coincide con lo SLIDE-ID stampato a video dal package.\ntcga_slide_id::AbstractString = ID utilizzato dal TCGA per identificare la slide.\nlink_slide::AbstractString = Link alle API per il download della slide.\nfilepath_zip::AbstractString = Percorso in cui è memorizzato il file .zip.\nfilepath_svs::AbstractString = Percorso in cui è memorizzato il file .svs.\n\nNote\n\nDati disponibili nel database JHistint_DB per ogni slide:\n\ncollection_name TEXT = Nome della collezione.\ncase_name TEXT = Nome del caso.\nTCGA_caseID TEXT = ID utilizzato dal TCGA per identificare il caso.\nslide_ID TEXT = Nome del singolo caso.\nTCGA_slideID TEXT UNIQUE = ID utilizzato dal TCGA per identificare la slide, UNIQUE evita la generazione di duplicati.\nslide_path_folder_zip TEXT = Percorso in cui è memorizzato il file .zip.\nslide_path_folder_svs TEXT = Percorso in cui è memorizzato il file .svs.\nslide_path_api TEXT = Link alle API per il download della slide.\nslide_path_folder_seg TEXT = Percorso in cui è memorizzato il file .tif segmentata.\nslide_svs BLOB = Slide istopatologica (immagine).\nslide_seg BLOB = Slide istopatologica segmentata (immagine).\nslide_info_TSS TEXT = Informazioni sulla slide - Tissue Source Site.\nslide_info_participant_code TEXT = Informazioni sulla slide - Codice associato al Participant, stringa alfanumerica.\nslide_info_sample_type TEXT = Informazioni sulla slide - Sample Type. I valori associati ai campioni aventi tumori sono nell'intervallo 01 - 09. 10 - 19 indica l'intervallo dedicato a campioni normali non malati. 20 - 29 indica campioni attualmente sotto controllo.\nslide_info_vial TEXT = Informazioni sulla slide - Vial. Relativo all'ordinamento del campione nella sequenza di campioni. I valori variano tra A - Z.\nslide_info_portion TEXT = Informazioni sulla slide - Portion. Relativo all'ordinamento delle porzioni analizzate associate ad un campione. Assume valori nell'intervallo 01-99.\nslide_info_type TEXT = Informazioni sulla slide - Tipo di immagine. I valori assumbili sono TS (Top Slide), BS (Bottom Slide) e MS (Middle Slide). Il valore alfanumerico indica l'ordinamento della slide.\n\n\n\n\n\n","category":"method"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"query_extract_slide_svs(collection_name::AbstractString)","category":"page"},{"location":"#JHistint.query_extract_slide_svs-Tuple{AbstractString}","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.query_extract_slide_svs","text":"query_extract_slide_svs(collection_name::AbstractString)\n\nLa funzione interroga il JHistint_DB ed estrae la lista di slide associate al nome della collezione fornita come argomento.\n\nArgomenti\n\ncollection_name::AbstractString: Nome della collezione di slide da ricercare nel JHistint_DB.\n\nValore di ritorno\n\nslide_list: Lista di tuple, ognuna delle quali contiene l'ID della slide, il file .svs della slide e il percorso della cartella contenente il file .svs.\n\n\n\n\n\n","category":"method"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"load_seg_slide(filepath_seg::AbstractString, segmented_slide::Array{ColorTypes.RGB{Float32}, 3}, slide_id::AbstractString)","category":"page"},{"location":"#JHistint.load_seg_slide-Tuple{AbstractString, Array{RGB{Float32}, 3}, AbstractString}","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.load_seg_slide","text":"load_seg_slide(filepath_seg::AbstractString, segmented_slide::Array{RGB{Float32}, 3}, slide_id::AbstractString)\n\nLa funzione aggiorna il JHistint_DB con il percorso del file dell'immagine segmentata e l'immagine segmentata.\n\nArgomento\n\nfilepath_seg::AbstractString: Percorso del file dell'immagine segmentata da aggiungere al DB.\nsegmented_slide::Array{ColorTypes.RGB{Float32}, 3}: Immagine segmentata da aggiungere al DB.\nslide_id::AbstractString: ID della slide da aggiornare con le informazioni dell'immagine segmentata.\n\n\n\n\n\n","category":"method"},{"location":"#Funzioni-di-supporto-API-(apiManager.jl)","page":"JHistint.jl - Julia Histopathology Interface","title":"Funzioni di supporto API (apiManager.jl)","text":"","category":"section"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"download_collection_values(filepath::AbstractString)","category":"page"},{"location":"#JHistint.download_collection_values-Tuple{AbstractString}","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.download_collection_values","text":"download_collection_values(filepath::AbstractString)\n\nFunzione per il download dei dati delle collezioni disponibili nel TCGA.\n\nArgomenti\n\nfilepath::AbstractString = Percorso in cui salvare il file .json ottenuto dalle API disponibili nel CDSA.\n\nNote\n\nL'API richiede la definizione del parentType e del parentId. Il parentId specifica l'identificativo della collezione. La collezione di immagini associate al TCGA è identificata dal codice: 5b9ef8e3e62914002e454c39. L'utilizzo del limit=0 imposta l'assenza di limiti nel file interrogato, garantendo il download del file in modo completo. L'API appartiene alla categoria per gestire le folder memorizzate nel repository. Il file scaricato è .json.\n\n\n\n\n\n","category":"method"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"extract_collection_values(filepath::AbstractString)","category":"page"},{"location":"#JHistint.extract_collection_values-Tuple{AbstractString}","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.extract_collection_values","text":"extract_collection_values(filepath::AbstractString)\n\nFunzione per estrarre i valori delle collezioni di dati dal file .json scaricato dalla funzione download_collection_values.\n\nArgomenti\n\nfilepath::AbstractString = Percorso in cui è memorizzato il file collectionlist.json.\n\nValore di ritorno\n\ncollection_values::Array{String} = Lista di collezioni di dati disponibili nel TCGA.\n\n\n\n\n\n","category":"method"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"download_project_infos(filepath::AbstractString, collection_name::AbstractString)","category":"page"},{"location":"#JHistint.download_project_infos-Tuple{AbstractString, AbstractString}","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.download_project_infos","text":"download_project_infos(filepath::AbstractString, collection_name::AbstractString)\n\nFunzione per il download dei metadati associati alla collezione selezionata all'avvio.\n\nArgomenti\n\nfilepath::AbstractString = Percorso in cui salvare il file .json associato alla collezione. Il file è indicato\n\ncon la dicitura collection_name.json.\n\ncollection_name::AbstractString = Nome della collezione di cui scaricare le slides.\n\nNote\n\nL'API richiede la definizione del parentType, del parentId e del name. L'attributo name identifica il nome della collezione di cui si desidera prelevare i dati (esempio: chol, esca, etc.). L'API appartiene alla categoria per gestire le folder memorizzate nel repository. Il file scaricato è .json.\n\n\n\n\n\n","category":"method"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"extract_project_id(filepath::AbstractString)","category":"page"},{"location":"#JHistint.extract_project_id-Tuple{AbstractString}","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.extract_project_id","text":"extract_project_id(filepath::AbstractString)\n\nFunzione per estrarre il valore del id dai metadati della collezione selezionata all'avvio.\n\nArgomenti\n\nfilepath::AbstractString =  Percorso in cui è memorizzato il file collection_name.json.\n\nValori di ritorno\n\nproject_id =  id della collezione.\n\n\n\n\n\n","category":"method"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"getCasesForProject(filepath_case::AbstractString, project_id::AbstractString)","category":"page"},{"location":"#JHistint.getCasesForProject-Tuple{AbstractString, AbstractString}","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.getCasesForProject","text":"getCasesForProject(filepath_case::AbstractString, project_id::AbstractString)\n\nFunzione per il download dei metadati associati ai casi della collezione selezionata all'avvio.\n\nArgomenti\n\nfilepath::AbstractString =  Percorso in cui salvare il file .json associato ai casi della collezione. Il file è indicato\n\ncon la dicitura collection_name.json.\n\nproject_id::AbstractString =  id della collezione.\n\nValori di ritorno\n\ncasesID_values =  Lista di id di tutti i casi della collezione.\ncasesNAME_values = Lista di name di tutti i casi della collezione.\n\nNote\n\nL'API richiede la definizione del parentType e del parentId. L'attributo parentType è impostato a folder data la struttura del repository. Il parentId è configurato definendo l'identificativo della collezione scelta. Il file scaricato è .json.\n\n\n\n\n\n","category":"method"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"download_zip(link::AbstractString, filepath::AbstractString)","category":"page"},{"location":"#JHistint.download_zip-Tuple{AbstractString, AbstractString}","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.download_zip","text":"download_zip(link::AbstractString, filepath::AbstractString)\n\nFunzione per il download delle slides istologiche in formato .zip associate ai casi della collezione selezionata all'avvio.\n\nArgomenti\n\nlink::AbstractString = URL per l'accesso alla API per il download delle slides.\nfilepath::AbstractString =  Percorso in cui salvare il file .zip.\n\n\n\n\n\n","category":"method"},{"location":"#Funzioni-di-supporto-ZIP-(zipManager.jl)","page":"JHistint.jl - Julia Histopathology Interface","title":"Funzioni di supporto ZIP (zipManager.jl)","text":"","category":"section"},{"location":"","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.jl - Julia Histopathology Interface","text":"extract_slide(filepath_zip::AbstractString)","category":"page"},{"location":"#JHistint.extract_slide-Tuple{AbstractString}","page":"JHistint.jl - Julia Histopathology Interface","title":"JHistint.extract_slide","text":"extract_slide(filepath_zip::AbstractString)\n\nFunzione per estrarre il contenuto dei file .zip scaricati dal CDSA.\n\nArgomenti\n\nfilepath_zip::AbstractString = Percorso in cui è salvato il file .zip relativo al singolo caso.\n\n\n\n\n\n","category":"method"}]
}
