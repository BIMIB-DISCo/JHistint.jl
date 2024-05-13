### -*- Mode: Julia -*-

### JHistint -- Julia Histopathology Interface.
### JHistint.jl

module JHistint

### Packages
using Logging
using HTTP
using JSON
using ZipFile
using TOML
using SQLite
using FileIO
using Dates
### Exported Functions
export test
export async_download_single_slide_from_collection
export download_single_slide_from_collection
export download_single_collection
export download_all_collection

### Included Files
include("DirectoryManager.jl")
include("JHistintLogger.jl")
include("apiManager.jl")
include("dbManager.jl")
include("zipManager.jl")

function test()
    jh_open_logger()
    jh_log_message("@info", "test")
    jh_close_logger()
end

### Main Functions
"""
    download_single_slide_from_collection(collection_name::AbstractString, path_to_save::AbstractString)

Function for downloading histological slides in SOPYHSM_app associated
with a collection available in TCGA.

# Arguments
- `collection_name::AbstractString` = Collection of TCGA data to download the
histological slides.
- `path_to_save::AbstractString` = Local folder path for saving
histological slides.
"""
function download_single_slide_from_collection(collection_name::AbstractString, 
                                               path_to_save::AbstractString)
    jh_open_logger()
    if Sys.iswindows() && path_to_save[1] == '/'
        path_to_save = path_to_save[2:end]
    end
    DirectoryManager.set_environment()
    # Check the value of the parameter
    filepath_collection_list = joinpath(DirectoryManager.CONFIG_DIR, 
                                        "collections", 
                                        "collectionlist.json")
    download_collection_values(filepath_collection_list)

    # List with all possible collection values
    collection_list = extract_collection_values(filepath_collection_list)
    collection_name = lowercase(collection_name)

    if collection_name in collection_list
        # Project Management (TCGA-OR-A5J1, TCGA-OR-A5J2, etc.)
        filepath_collection = joinpath(DirectoryManager.CONFIG_DIR, 
                                       "collections",
                                       "$(collection_name).json")
        download_project_infos(filepath_collection, collection_name)
        project_id = extract_project_id(filepath_collection)
        filepath_case = joinpath(DirectoryManager.CONFIG_DIR, 
                                 "cases", 
                                 "$(collection_name).json")
        casesID_values, casesNAME_values = getCasesForProject(filepath_case, project_id)

        # Slides Management
        if isdir(joinpath(path_to_save, "$collection_name"))
            jh_log_message("@info", "Updating data...")
        else
            mkdir(joinpath(path_to_save, "$collection_name"))
            jh_log_message("@info", "Create new collection folder")
        end

        i, j = casesID_values[1], casesNAME_values[1]
        single_casesID_values, single_casesNAME_values = getCasesForProject(filepath_case, i)
        x, y = single_casesID_values[1], single_casesNAME_values[1]
        if !isdir(joinpath(path_to_save, "$(collection_name)", "$j"))
            jh_log_message("@info", "Create new slide folder")
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
        jh_log_message("@info", "DOWNLOAD Slide complete: CASE NAME = $j - SLIDE ID = $y")
    else
        jh_log_message("@info", "Collection selected doesn't exist")
    end
    jh_close_logger()
end

"""
    async_download_single_slide_from_collection(collection_name::AbstractString, path_to_save::AbstractString)

Function for asyncronous download of histological slides in SOPYHSM_app associated
with a collection available in TCGA.

# Arguments
- `collection_name::AbstractString` = Collection of TCGA data to download the
histological slides.
- `path_to_save::AbstractString` = Local folder path for saving
histological slides.
"""
function async_download_single_slide_from_collection(collection_name::AbstractString, path_to_save::AbstractString)
    task = @task download_single_slide_from_collection(collection_name, path_to_save)
    schedule(task)
    return task
end

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
    download_all_collection(path_to_save::AbstractString)

Function for downloading histological slides associated with all
collections available in TCGA.

# Arguments
- `path_to_save::AbstractString` = Local folder path for saving
histological slides.

```julia
# Examples with valid input
julia> JHistint.download_all_collection("C:\\...")
```
"""
function download_all_collection(path_to_save::AbstractString)
    DirectoryManager.set_environment()
    # Collection Management (acc, blca, etc.)
    filepath_collection_list = joinpath(DirectoryManager.CONFIG_DIR, "collections", "collectionlist.json")
    download_collection_values(filepath_collection_list)
    collection_list = extract_collection_values(filepath_collection_list)

    for collection_name in collection_list
        # Project Management (TCGA-OR-A5J1, TCGA-OR-A5J2, etc.)
        filepath_collection_list = joinpath(DirectoryManager.CONFIG_DIR, "collections", "$(collection_name).json")
        download_project_infos(filepath_collection_list, collection_name)
        project_id = extract_project_id(filepath_collection_list)
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

end 
### end of module -- JHistint.jl
