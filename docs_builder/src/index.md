# JHistint.jl - Julia Histopathology Interface

Julia interface for implementing the REST APIs available on the Cancer Slide Digital Archive (CDSA) portal for downloading histological images available in The Cancer Genome Atlas (TCGA). The Cancer Slide Digital Archive (CDSA) is a web platform for support, sharing, and analysis of digital pathological data. Currently, it hosts over 23,000 images associated with the data available on "The Cancer Genome Atlas" Data Portal. The library includes functions for managing image-processing algorithms for cell segmentation, constructing the adjacency matrix, and interfacing with the J-Space.jl package.

Link GitHub repository: [JHistint.jl](https://github.com/niccolo99mandelli/JHistint.jl)

Link GitHub repository: [J-Space.jl](https://github.com/niccolo99mandelli/J-Space.jl)

CDSA Portal: [Click Here](https://api.digitalslidearchive.org/#collections)    

Repository containing the data mapped in the portal: [Click Here](https://cancer.digitalslidearchive.org/#!/CDSA/acc/TCGA-OR-A5J1)

Guide to using the APIs: [Click Here](https://api.digitalslidearchive.org/api/v1)

## Package Structure
* The `case` and `collection` folders store metadata in `.json` format for individual cases and collections available on the TCGA Data Portal. The `collection` folder is structured as follows:
  * `collectionlist.json` = Stores access data (metadata) for collections (Projects in TCGA).
  * `collection_name.json` = Stores access data (metadata) for a single collection. The `.json` file is generated based on the collection chosen by the user.
* The `case` folder is structured as follows:
  * `collection_name.json` = Stores all metadata related to cases associated with the collection selected by the user.
* The `slides` folder stores histological images related to individual cases. The images are organized based on collection (`TCGA-chol`, `TCGA-esca`, etc.), and the individual case being analyzed (`TCGA-2H-A9GF`, `TCGA-2H-A9GG`, etc.). Within each folder related to the case, the slides are stored in compressed `.zip` files. The format of each individual slide is `.tif`. The folder names related to the cases correspond to the values of the `Case ID` field listed in the TCGA Data Portal. The names of the `.zip` files located in each folder refer to the `Sample ID` attribute associated with the patient. The slide name is given by concatenating the `Slide ID` and `Slide UUID` attributes that can be found in the lower section of the web page dedicated to the generic case `TCGA-XX-YYYY`.

```
Example: TCGA-02-0001-01C-01-TS1.zip  
  - 02 = refers to the TSS (Tissue Source Site).  
  - 0001 = refers to the code associated with the Participant, an alphanumeric string.  
  - 01 = refers to the Sample Type. The values associated with tumor samples are in the range 01-09. 10-19 indicates the range for non-diseased normal samples. 20-29 indicates samples currently under control.  
  - C = refers to the Vial field related to the ordering of the sample in the sample sequence. Values range from A-Z.  
  - 01 = refers to the Portion field related to the ordering of the analyzed portions associated with a sample. It takes values in the range 01-99.  
  - TS1 = refers to the Slide field related to the type of image. The values that can be assumed are TS (Top Slide), BS (Bottom Slide), and MS (Middle Slide). The alphanumeric value indicates the slide ordering.
```

## JHistint Collections
The available collections are:
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
To download a specific collection, just indicate the name of the collection: `BRCA`, `OV`, `LUAD`.

## Package Installation
The `JHistint` package is available in the Julia Registries and can be installed as follows:
```
julia > using Pkg
julia > Pkg.add("JHistint")
julia > using JHistint
```
Otherwise, type `]` in the Julia REPL and execute:
```
(@v1.8) pkg > add JHistint
(@v1.8) pkg > using JHistint
```

## Download Slides main functions (JHistint.jl)
```@docs
download_single_collection(collection_name::AbstractString)
```

```@docs
download_all_collection()
```
## Cell Segmentation Slides main functions (JHistint.jl)
```@docs
slide_cell_segmentation_without_download(collection_name::AbstractString)
```

```@docs
slide_cell_segmentation_with_download(collection_name::AbstractString)
```
## Support Functions for Cell Segmentation (segmentationManager.jl)
```@docs
apply_segmentation_without_download(slide_info::Tuple{String, Vector{UInt8}, String})
```

```@docs
apply_segmentation_with_download(slide_info::Tuple{String, Vector{UInt8}, String})
```

```@docs
get_random_color(seed)
```

## Support Functions for Graph (graphManager.jl)
```@docs
save_adjacency_matrix(matrix::Matrix{Int64}, filepath_matrix::AbstractString)
```

```@docs
weighted_graph_to_adjacency_matrix(G::SimpleWeightedGraph{Int64, Float64}, n::Int64)
```

## DB Support Functions (dbManager.jl)
```@docs
insert_record_DB(col_name::AbstractString,
                        cas_name::AbstractString,
                        tcga_case_id::AbstractString,
                        sin_cas_name::AbstractString,
                        tcga_slide_id::AbstractString,
                        link_slide::AbstractString,
                        filepath_zip::AbstractString,
                        filepath_svs::AbstractString)
```

```@docs
query_extract_slide_svs(collection_name::AbstractString)
```

```@docs
load_seg_slide(filepath_seg::AbstractString, filepath_matrix::AbstractString, matrix::Matrix{Int64}, slide_id::AbstractString)
```

## API Support Functions (apiManager.jl)
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

## ZIP Support Functions (zipManager.jl)
```@docs
extract_slide(filepath_zip::AbstractString)
```
