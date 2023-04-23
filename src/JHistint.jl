module JHistint

# Packages
using HTTP
using JSON
using ZipFile
using TOML
using SQLite
using DataFrames
using Images
using ImageSegmentation
using ImageMagick
using ImageView
using FileIO
using Random
using IndirectArrays
using Graphs
using LightGraphs
using SimpleWeightedGraphs

# Exported Functions
export download_single_collection
export download_all_collection
export slide_cell_segmentation_without_download
export slide_cell_segmentation_with_donwload

# Included Files
include("apiManager.jl")
include("dbManager.jl")
include("zipManager.jl")
include("segmentationManager.jl")

# Line of code to input the name of collection from command line. Not usable in Package
# collection_name = select_collection_name(collection_list)

"""
    download_single_collection(collection_name::AbstractString)

Funzione per il download delle slides istologiche associate ad una collezione disponibile nel TCGA.

# Argomenti
- `collection_name::AbstractString` = Collezione di dati TCGA di cui scaricare le slides istologiche.

# Note
La funzione valuta l'argomento `collection_name`, in caso di collezione non valida considera
la configurazione del file `Config.toml`. Il valore impostato nel package è `default`.
```julia
# Esempi con input validi
julia> JHistint.download_single_collection("acc")
julia> JHistint.download_single_collection("bLca")
```
```julia
# Esempi con input non validi
julia> JHistint.download_single_collection("ac")
julia> JHistint.download_single_collection("")
```
"""
function download_single_collection(collection_name::AbstractString)
    # Check the value of the parameter
    filepath_collection = joinpath(@__DIR__, "..", "collection", "collectionlist.jsn")
    download_collection_values(filepath_collection)
    collection_list = extract_collection_values(filepath_collection)

    if lowercase(collection_name) in collection_list

        collection_name = lowercase(collection_name)
        # Project Management (TCGA-OR-A5J1, TCGA-OR-A5J2, etc.)
        filepath_collection = joinpath(@__DIR__, "..", "collection", "$(collection_name).jsn")
        download_project_infos(filepath_collection, collection_name)
        project_id = extract_project_id(filepath_collection)
        filepath_case = joinpath(@__DIR__, "..", "case", "$(collection_name).jsn")
        casesID_values, casesNAME_values = getCasesForProject(filepath_case, project_id)

        # Slides Management
        if isdir(joinpath(@__DIR__, "..", "slides", "$collection_name"))
            println("Update data ...")
        else
            mkdir(joinpath(@__DIR__, "..", "slides", "$collection_name"))
        end
        for (i, j) in zip(casesID_values, casesNAME_values)
            single_casesID_values, single_casesNAME_values = getCasesForProject(filepath_case, i)
            for (x, y) in zip(single_casesID_values, single_casesNAME_values)
                if isdir(joinpath(@__DIR__, "..", "slides", "$(collection_name)", "$j"))
                else
                    mkdir(joinpath(@__DIR__, "..", "slides", "$(collection_name)", "$j"))
                end
                filepath_slides = joinpath(@__DIR__, "..", "slides", "$(collection_name)", "$j", "$(y).zip")
                link_slides = "https://api.digitalslidearchive.org/api/v1/folder/$x/download"
                download_zip(link_slides, filepath_slides)
                filepath_svs = extract_slide(filepath_slides)
                insert_record_DB(collection_name, j, i, y, x, link_slides, filepath_slides, filepath_svs)
                println("DOWNLOAD Slide complete: CASE NAME = $j - SLIDE ID = $y")
            end
        end
    else
        println("ERROR : $collection_name - collection not avaiable. Using collection_name field in Config.toml ...")

        # Line of code for definition of "collection_name" from Config.toml
        filepath_config = joinpath(@__DIR__, "..", "Config.toml")
        config = TOML.parsefile(filepath_config)
        collection_name = config["collection_name"]
        if collection_name != "default" && lowercase(collection_name) in collection_list

            collection_name = lowercase(collection_name)
            # Project Management (TCGA-OR-A5J1, TCGA-OR-A5J2, etc.)
            filepath_collection = joinpath(@__DIR__, "..", "collection", "$(collection_name).jsn")
            download_project_infos(filepath_collection, collection_name)
            project_id = extract_project_id(filepath_collection)
            filepath_case = joinpath(@__DIR__, "..", "case", "$(collection_name).jsn")
            casesID_values, casesNAME_values = getCasesForProject(filepath_case, project_id)

            # Slides Management
            if isdir(joinpath(@__DIR__, "..", "slides", "$collection_name"))
                println("Update data ...")
            else
                mkdir(joinpath(@__DIR__, "..", "slides", "$collection_name"))
            end
            for (i, j) in zip(casesID_values, casesNAME_values)
                single_casesID_values, single_casesNAME_values = getCasesForProject(filepath_case, i)
                for (x, y) in zip(single_casesID_values, single_casesNAME_values)
                    if isdir(joinpath(@__DIR__, "..", "slides", "$(collection_name)", "$j"))
                    else
                        mkdir(joinpath(@__DIR__, "..", "slides", "$(collection_name)", "$j"))
                    end
                    filepath_slides = joinpath(@__DIR__, "..", "slides", "$(collection_name)", "$j", "$(y).zip")
                    link_slides = "https://api.digitalslidearchive.org/api/v1/folder/$x/download"

                    download_zip(link_slides, filepath_slides)
                    filepath_svs = extract_slide(filepath_slides)
                    insert_record_DB(collection_name, j, i, y, x, link_slides, filepath_slides, filepath_svs)
                    println("DOWNLOAD Slide complete: CASE NAME = $j - SLIDE ID = $y")
                end
            end
        else
            println("ERROR : $collection_name - collection not avaiable. Change collection_name field in Config.toml.")
        end
    end
end

"""
    download_all_collection()

Funzione per il download delle slides istologiche associate a tutte le collezioni disponbile nel TCGA.

```julia
# Esempi con input validi
julia> JHistint.download_all_collection()
```
"""
function download_all_collection()
    # Collection Management (acc, blca, etc.)
    filepath_collection = joinpath(@__DIR__, "..", "collection", "collectionlist.jsn")
    download_collection_values(filepath_collection)
    collection_list=extract_collection_values(filepath_collection)

    for collection_name in collection_list
        # Project Management (TCGA-OR-A5J1, TCGA-OR-A5J2, etc.)
        filepath_collection = joinpath(@__DIR__, "..", "collection", "$(collection_name).jsn")
        download_project_infos(filepath_collection, collection_name)
        project_id = extract_project_id(filepath_collection)
        filepath_case = joinpath(@__DIR__, "..", "case", "$(collection_name).jsn")
        casesID_values, casesNAME_values = getCasesForProject(filepath_case, project_id)

        # Slides Management
        if isdir(joinpath(@__DIR__, "..", "slides", "$collection_name"))
            println("Update data ...")
        else
            mkdir(joinpath(@__DIR__, "..", "slides", "$collection_name"))
        end
        for (i, j) in zip(casesID_values, casesNAME_values)
            single_casesID_values, single_casesNAME_values = getCasesForProject(filepath_case, i)
            for (x, y) in zip(single_casesID_values, single_casesNAME_values)
                if isdir(joinpath(@__DIR__, "..", "slides", "$(collection_name)", "$j"))
                else
                    mkdir(joinpath(@__DIR__, "..", "slides", "$(collection_name)", "$j"))
                end
                filepath_slides = joinpath(@__DIR__, "..", "slides", "$(collection_name)", "$j", "$(y).zip")
                link_slides = "https://api.digitalslidearchive.org/api/v1/folder/$x/download"

                download_zip(link_slides, filepath_slides)
                filepath_svs = extract_slide(filepath_slides)
                insert_record_DB(collection_name, j, i, y, x, link_slides, filepath_slides, filepath_svs)
                println("DOWNLOAD Slide complete: CASE NAME = $j - SLIDE ID = $y")
            end
        end
    end
end

"""
    slide_cell_segmentation_with_download(collection_name::AbstractString)

Funzione per esecuzione della segmentazione cellulare delle slide istopatologiche presenti nel database `JHistint_DB` associate al nome della collezione fornita come argomento.
La funzione effettua il download della slide segmentata. Il risultato viene posto nella directory in cui è presente la slide.
Dopo la generazione della slide segmentata la funzione procede con la costruzione e il salvataggio del rispettivo grafo e della matrice di adiacenza.

# Argomenti
- `collection_name::AbstractString` = Collezione di slide TCGA di cui effettuare la segmentazione cellulare.

# Note
Per ogni slide nel DB viene eseguita la segmentazione delle cellule utilizzando la funzione `apply_segmentation_with_download` e il percorso in cui è salvato il risultato viene memorizzato nel DB utilizzando la funzione `load_seg_slide`.
Il processo di segmentazione è definito in 6 step:
- LOAD SLIDE ... (slide_id)
- APPLY SEGMENTATION ... (slide_id)
- BUILD SEGMENTED SLIDE ... (slide_id)
- BUILD GRAPH ... (slide_id)
- BUILD & SAVE ADJACENCY MATRIX ... (slide_id)
- SAVE SEGMENTED SLIDE ... (slide_id)
La matrice di adiacenza viene memorizzata nello stesso percorso della immagine originale e di quella segmentata in formato testuale.
Infine, viene stampato un messaggio di conferma per ogni slide segmentata.
```julia
# Esempi con input validi
julia> JHistint.slide_cell_segmentation("acc")
julia> JHistint.slide_cell_segmentation("bLca")
```
```julia
# Esempi con input non validi
julia> JHistint.slide_cell_segmentation("ac")
julia> JHistint.slide_cell_segmentation("")
```
"""
function slide_cell_segmentation_with_download(collection_name::AbstractString)
    # Check the value of the parameter
    filepath_collection = joinpath(@__DIR__, "..", "collection", "collectionlist.jsn")
    download_collection_values(filepath_collection)
    collection_list = extract_collection_values(filepath_collection)

    if lowercase(collection_name) in collection_list
        collection_name = lowercase(collection_name)
        slide_list = query_extract_slide_svs(collection_name)
        if isempty(slide_list)
            println("ERROR : No match for $collection_name - collection in DB. Download the collection before.")
        else
            for record in slide_list
                filepath_seg, filepath_matrix, matrix = apply_segmentation_with_download(record)
                load_seg_slide(filepath_seg, filepath_matrix, matrix, record[1])
                slide_id = record[1]
                println("SEGMENTATION SLIDE, BUILD GRAPH & MATRIX complete for SLIDE ID = $slide_id")
            end
        end
    else
        println("ERROR : $collection_name - collection not avaiable. Retry with a new collection.")
    end
end

function slide_cell_segmentation_without_download(collection_name::AbstractString)
    # Check the value of the parameter
    filepath_collection = joinpath(@__DIR__, "..", "collection", "collectionlist.jsn")
    download_collection_values(filepath_collection)
    collection_list = extract_collection_values(filepath_collection)

    if lowercase(collection_name) in collection_list
        collection_name = lowercase(collection_name)
        slide_list = query_extract_slide_svs(collection_name)
        if isempty(slide_list)
            println("ERROR : No match for $collection_name - collection in DB. Download the collection before.")
        else
            for record in slide_list
                filepath_matrix, matrix = apply_segmentation_without_download(record)
                load_seg_slide("not saved", filepath_matrix, matrix, record[1])
                slide_id = record[1]
                println("SEGMENTATION SLIDE, BUILD GRAPH & MATRIX complete for SLIDE ID = $slide_id")
            end
        end
    else
        println("ERROR : $collection_name - collection not avaiable. Retry with a new collection.")
    end
end
end
