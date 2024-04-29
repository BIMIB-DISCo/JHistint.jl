### -*- Mode: Julia -*-

### JHistint -- Julia Histopathology Interface.
### JHistint.jl

module JHistint

### Packages
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
using CSV
using Luxor
using Karnak
using MetaGraphs
using Plots
using VoronoiCells
using GeometryBasics

### Exported Functions
export download_single_collection
export download_all_collection
export slide_cell_segmentation_without_download
export slide_cell_segmentation_with_download
export start_segmentation_SOPHYSM_tessellation
export start_segmentation_SOPHYSM_graph

### Included Files
include("DirectoryManager.jl")

include("apiManager.jl")
include("dbManager.jl")
include("zipManager.jl")
include("segmentationManager.jl")
include("graphManager.jl")
include("noiseManager.jl")
include("tessellationManager.jl")

### Main Functions
"""
    download_single_collection(collection_name::AbstractString, path_to_save::AbstractString)

Function for downloading histological slides in SOPYHSM_app associated
with a collection available in TCGA.

# Arguments
- `collection_name::AbstractString` = Collection of TCGA data to download the
histological slides.
- `path_to_save::AbstractString` = Local folder path for saving
histological slides.

# Notes
The function evaluates the `collection_name` argument, and in case of an
invalid collection, considers the configuration in the
`Config.toml` file. The value set in the package is `default`.
```julia
# Examples with valid input
julia> JHistint.download_single_collection("acc", "C:\\...")
julia> JHistint.download_single_collection("bLca", "C:\\...")
```
```julia
# Examples with invalid input
julia> JHistint.download_single_collection("ac", "C:\\...")
julia> JHistint.download_single_collection("", "C:\\...")
```
"""
function download_single_collection(collection_name::AbstractString, path_to_save::AbstractString)
    DirectoryManager.set_environment()
    # Check the value of the parameter
    filepath_collection_list = joinpath(DirectoryManager.CONFIG_DIR, "collections", "collectionlist.json")
    download_collection_values(filepath_collection_list)

    # List with all possible collection values
    collection_list = extract_collection_values(filepath_collection_list)
    collection_name = lowercase(collection_name)

    if collection_name in collection_list
        # Project Management (TCGA-OR-A5J1, TCGA-OR-A5J2, etc.)
        filepath_collection = joinpath(DirectoryManager.CONFIG_DIR, "collections", "$(collection_name).json")
        download_project_infos(filepath_collection, collection_name)
        project_id = extract_project_id(filepath_collection)
        filepath_case = joinpath(DirectoryManager.CONFIG_DIR, "cases", "$(collection_name).json")
        casesID_values, casesNAME_values = getCasesForProject(filepath_case, project_id)

        # Slides Management
        if isdir(joinpath(path_to_save, "$collection_name"))
            println("Update data ...")
        else
            mkdir(joinpath(path_to_save, "$collection_name"))
        end
        for (i, j) in zip(casesID_values, casesNAME_values)
            single_casesID_values, single_casesNAME_values = getCasesForProject(filepath_case, i)
            for (x, y) in zip(single_casesID_values, single_casesNAME_values)
                if !isdir(joinpath(path_to_save, "$(collection_name)", "$j"))
                    mkdir(joinpath(path_to_save, "$(collection_name)", "$j"))
                end
                filepath_slides = joinpath(path_to_save, "$(collection_name)", "$j", "$(y).zip")
                link_slides = "https://api.digitalslidearchive.org/api/v1/folder/$x/download"
                download_zip(link_slides, filepath_slides)
                filepath_svs = extract_slide(filepath_slides)
                insert_record_DB_SOPHYSM(collection_name,
                                            j, i, y, x,
                                            link_slides,
                                            filepath_slides,
                                            filepath_svs,
                                            path_to_save)
                println("DOWNLOAD Slide complete: CASE NAME = $j - SLIDE ID = $y")
            end
        end

    else
        return "error"
    end
end

"""
    download_all_collection()

Function for downloading histological slides associated with all
collections available in TCGA.

```julia
# Examples with valid input
julia> JHistint.download_all_collection()
```
"""
function download_all_collection()
    # Collection Management (acc, blca, etc.)
    filepath_collection_list = joinpath(@__DIR__, "..", "collection", "collectionlist.json")
    download_collection_values(filepath_collection_list)
    collection_list=extract_collection_values(filepath_collection_list)

    for collection_name in collection_list
        # Project Management (TCGA-OR-A5J1, TCGA-OR-A5J2, etc.)
        filepath_collection_list = joinpath(@__DIR__, "..", "collection", "$(collection_name).json")
        download_project_infos(filepath_collection_list, collection_name)
        project_id = extract_project_id(filepath_collection_list)
        filepath_case = joinpath(@__DIR__, "..", "case", "$(collection_name).json")
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
                insert_record_DB(collection_name,
                                    j, i, y, x,
                                    link_slides,
                                    filepath_slides,
                                    filepath_svs)
                println("DOWNLOAD Slide complete: CASE NAME = $j - SLIDE ID = $y")
            end
        end
    end
end

"""
    download_all_collection_SOPHYSM(path_to_save::AbstractString)

Function for downloading histological slides associated with all
collections available in TCGA.

# Arguments
- `path_to_save::AbstractString` = Local folder path for saving
histological slides.

```julia
# Examples with valid input
julia> JHistint.download_all_collection_SOPHYSM("C:\\...")
```
"""
function download_all_collection_SOPHYSM(path_to_save::AbstractString)
    # Collection Management (acc, blca, etc.)
    filepath_collection_list = joinpath(@__DIR__, "..", "collection", "collectionlist.json")
    download_collection_values(filepath_collection_list)
    collection_list=extract_collection_values(filepath_collection_list)

    for collection_name in collection_list
        # Project Management (TCGA-OR-A5J1, TCGA-OR-A5J2, etc.)
        filepath_collection_list = joinpath(@__DIR__, "..", "collection", "$(collection_name).json")
        download_project_infos(filepath_collection_list, collection_name)
        project_id = extract_project_id(filepath_collection_list)
        filepath_case = joinpath(@__DIR__, "..", "case", "$(collection_name).json")
        casesID_values, casesNAME_values = getCasesForProject(filepath_case, project_id)

        # Slides Management
        if isdir(joinpath(path_to_save, "$collection_name"))
            println("Update data ...")
        else
            mkdir(joinpath(path_to_save, "$collection_name"))
        end
        for (i, j) in zip(casesID_values, casesNAME_values)
            single_casesID_values, single_casesNAME_values = getCasesForProject(filepath_case, i)
            for (x, y) in zip(single_casesID_values, single_casesNAME_values)
                if isdir(joinpath(path_to_save, "$(collection_name)", "$j"))
                else
                    mkdir(joinpath(path_to_save, "$(collection_name)", "$j"))
                end
                filepath_slides = joinpath(path_to_save, "$(collection_name)", "$j", "$(y).zip")
                link_slides = "https://api.digitalslidearchive.org/api/v1/folder/$x/download"

                download_zip(link_slides, filepath_slides)
                filepath_svs = extract_slide(filepath_slides)
                insert_record_DB_SOPHYSM(collection_name,
                                            j, i, y, x,
                                            link_slides,
                                            filepath_slides,
                                            filepath_svs,
                                            path_to_save)
                println("DOWNLOAD Slide complete: CASE NAME = $j - SLIDE ID = $y")
            end
        end
    end
end

"""
    slide_cell_segmentation_without_download(collection_name::AbstractString)

Function for performing cell segmentation on histopathological slides present
in the `JHistint_DB` database associated with the collection name provided as
an argument. After generating the segmented slide, the function
proceeds with constructing and saving the corresponding graph and adjacency matrix.

# Arguments
- `collection_name::AbstractString` = Collection of TCGA data to download
the histological slides.

# Notes
The function utilizes the `JHistint_DB` database for performing cell
segmentation on the histopathological slides associated with the provided
collection name. It generates a segmented slide and constructs a corresponding
graph and adjacency matrix. The output files are saved in a user-defined
directory. The function may take a considerable amount of time to complete,
depending on the size of the slides and the complexity of the segmentation algorithm.
For each slide in the database, cell segmentation is performed using the
`apply_segmentation_without_download` function, and the path where the result
is saved is stored in the database using the `load_seg_slide` function.
The segmentation process is defined in 4 steps:
- LOAD SLIDE ... (slide_id)
- APPLY SEGMENTATION ... (slide_id)
- BUILD GRAPH ... (slide_id)
- BUILD & SAVE ADJACENCY MATRIX ... (slide_id)
- J-SPACE features ... (slide_id)
The adjacency matrix is saved in the same directory as the original image in
text format.
Finally, a confirmation message is printed for each segmented slide. Unlike
the `slide_cell_segmentation_with_download` function,
this function does not involve the creation and download of the segmented image.
```julia
# Examples with valid input
julia> JHistint.slide_cell_segmentation_without_download("acc")
julia> JHistint.slide_cell_segmentation_without_download("bLca")
```
```julia
# Examples with invalid input
julia> JHistint.slide_cell_segmentation_without_download("ac")
julia> JHistint.slide_cell_segmentation_without_download("")
```
"""
function slide_cell_segmentation_without_download(collection_name::AbstractString)
    # Check the value of the parameter
    filepath_collection_list = joinpath(@__DIR__, "..", "collection", "collectionlist.json")
    download_collection_values(filepath_collection_list)
    collection_list = extract_collection_values(filepath_collection_list)

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
                println("")
                # J-SPACE Interface
                filepath_file_JSPACE = replace(filepath_matrix, ".txt" => "_Files_JSpace")
                if isdir(filepath_file_JSPACE)
                    # do nothing
                else
                    mkdir(filepath_file_JSPACE) # create directory for saving J_Space files
                end

                filepath_plot_JSPACE = replace(filepath_matrix, ".txt" => "_Plots_JSpace")
                if isdir(filepath_plot_JSPACE)
                    # do nothing
                else
                    mkdir(filepath_plot_JSPACE) # create directory for saving J_Space plots
                end

                filepath_reference_JSPACE = replace(filepath_matrix, ".txt" => "_reference.fasta")
                filepath_dataframe_labels = replace(filepath_matrix, r"....$" => "_dataframe_labels.csv")
                filepath_dataframe_edges = replace(filepath_matrix, r"....$" => "_dataframe_edges.csv")
                Start_J_Space(filepath_reference_JSPACE,
                                filepath_matrix,
                                filepath_file_JSPACE,
                                filepath_plot_JSPACE,
                                slide_id,
                                filepath_dataframe_edges,
                                filepath_dataframe_labels)
            end
        end
    else
        println("ERROR : $collection_name - collection not avaiable. Retry with a new collection.")
    end
end

"""
    slide_cell_segmentation_with_download(collection_name::AbstractString)

Function for performing cell segmentation on histopathological slides present in
the `JHistint_DB` database associated with the collection name provided as an
argument. The function downloads the segmented slide, which is placed in the
same directory as the original slide. After generating the segmented slide,
the function proceeds with constructing and saving the corresponding graph
and adjacency matrix.

# Arguments
- `collection_name::AbstractString` = TCGA data collection for which to
perform cell segmentation.

# Notes
The function utilizes the `JHistint_DB` database for performing cell segmentation
on the histopathological slides associated with the provided collection name.
It generates a segmented slide and constructs a corresponding graph and adjacency matrix.
The output files are saved in a user-defined directory. The function may take a
considerable amount of time to complete, depending on the size of the slides
and the complexity of the segmentation algorithm. For each slide in the database,
cell segmentation is performed using the `apply_segmentation_with_download` function,
and the path where the result is saved is stored in the database using the
`load_seg_slide` function. The segmentation process is similar to that described
in the `slide_cell_segmentation_without_download` function, with the added step
of downloading the segmented image and placing it in the same directory as the
original slide. The segmentation process is defined in 6 steps:
- LOAD SLIDE ... (slide_id)
- APPLY SEGMENTATION ... (slide_id)
- BUILD SEGMENTED SLIDE ... (slide_id)
- BUILD GRAPH ... (slide_id)
- BUILD & SAVE ADJACENCY MATRIX ... (slide_id)
- SAVE SEGMENTED SLIDE ... (slide_id)
- J-SPACE features ... (slide_id)
The adjacency matrix is saved in text format in the same directory as both
the original and segmented images.
Finally, a confirmation message is printed for each segmented slide.
```julia
# Examples with valid input
julia> JHistint.slide_cell_segmentation_with_download("acc")
julia> JHistint.slide_cell_segmentation_with_download("bLca")
```
```julia
# Examples with invalid input
julia> JHistint.slide_cell_segmentation_with_download("ac")
julia> JHistint.slide_cell_segmentation_with_download("")
```
"""
function slide_cell_segmentation_with_download(collection_name::AbstractString)
    # Check the value of the parameter
    filepath_collection_list = joinpath(@__DIR__, "..", "collection", "collectionlist.json")
    download_collection_values(filepath_collection_list)
    collection_list = extract_collection_values(filepath_collection_list)

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
                println("")
                # J-SPACE Interface
                filepath_file_JSPACE = replace(filepath_matrix, ".txt" => "_Files_JSpace")
                if isdir(filepath_file_JSPACE)
                    # do nothing
                else
                    mkdir(filepath_file_JSPACE) # create directory for saving J_Space files
                end

                filepath_plot_JSPACE = replace(filepath_matrix, ".txt" => "_Plots_JSpace")
                if isdir(filepath_plot_JSPACE)
                    # do nothing
                else
                    mkdir(filepath_plot_JSPACE) # create directory for saving J_Space plots
                end

                filepath_reference_JSPACE = replace(filepath_matrix, ".txt" => "_reference.fasta")
                filepath_dataframe_labels = replace(filepath_matrix, r"....$" => "_dataframe_labels.csv")
                filepath_dataframe_edges = replace(filepath_matrix, r"....$" => "_dataframe_edges.csv")
                Start_J_Space(filepath_reference_JSPACE,
                                filepath_matrix,
                                filepath_file_JSPACE,
                                filepath_plot_JSPACE,
                                slide_id,
                                filepath_dataframe_edges,
                                filepath_dataframe_labels)
            end
        end
    else
        println("ERROR : $collection_name - collection not avaiable. Retry with a new collection.")
    end
end

"""
    start_segmentation_SOPHYSM_tessellation(filepath_input::AbstractString,
                                    filepath_output::AbstractString,
                                    thresholdGray::Float64,
                                    thresholdMarker::Float64,
                                    min_threshold::Float32,
                                    max_threshold::Float32)

Initiates the SOPHYSM segmentation process for histological image using
tessellation process.

# Arguments
- `filepath_input::AbstractString`: The file path to the input histological
image to be segmented.
- `filepath_output::AbstractString`: The file path where the segmented results
and related data will be saved.
- `thresholdGray::Float64`: The grayscale threshold used for initial
image processing.
- `thresholdMarker::Float64`: The marker threshold for identifying
cellular structures.
- `min_threshold`: Minimal threshold for considering segments area.
- `max_threshold`: Maximal threshold for considering segments area.
"""
function start_segmentation_SOPHYSM_tessellation(filepath_input::AbstractString,
                                    filepath_output::AbstractString,
                                    thresholdGray::Float64,
                                    thresholdMarker::Float64,
                                    min_threshold::Float32,
                                    max_threshold::Float32)
    apply_segmentation_SOPHYSM_tessellation(filepath_input,
                                filepath_output,
                                thresholdGray,
                                thresholdMarker,
                                min_threshold,
                                max_threshold)
end

"""
    start_segmentation_SOPHYSM_graph(filepath_input::AbstractString,
                                    filepath_output::AbstractString,
                                    thresholdGray::Float64,
                                    thresholdMarker::Float64,
                                    min_threshold::Float32,
                                    max_threshold::Float32)

Initiates the SOPHYSM segmentation process for histological image using
graph construction.

# Arguments
- `filepath_input::AbstractString`: The file path to the input histological
image to be segmented.
- `filepath_output::AbstractString`: The file path where the segmented results
and related data will be saved.
- `thresholdGray::Float64`: The grayscale threshold used for initial
image processing.
- `thresholdMarker::Float64`: The marker threshold for identifying
cellular structures.
- `min_threshold`: Minimal threshold for considering segments area.
- `max_threshold`: Maximal threshold for considering segments area.
"""
function start_segmentation_SOPHYSM_graph(filepath_input::AbstractString,
                                    filepath_output::AbstractString,
                                    thresholdGray::Float64,
                                    thresholdMarker::Float64,
                                    min_threshold::Float32,
                                    max_threshold::Float32)
    apply_segmentation_SOPHYSM_graph(filepath_input,
                                filepath_output,
                                thresholdGray,
                                thresholdMarker,
                                min_threshold,
                                max_threshold)
end
end

### end of file -- JHistint.jl
