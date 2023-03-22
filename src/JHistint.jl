module JHistint
include("functions.jl")

# Packages
using HTTP
using JSON
using ZipFile
using TOML

# Line of code to input the name of collection from command line. Not usable in Package
# collection_name = select_collection_name(collection_list)

function download_single_collection(collection_name::AbstractString)
    # Check the value of the parameter
    filepath_collection = joinpath(@__DIR__, "..", "collection", "collectionlist.jsn")
    download_collection_values(filepath_collection)
    collection_list = print_collection_values(filepath_collection)

    if lowercase(collection_name) in collection_list
        
        # Project Management (TCGA-OR-A5J1, TCGA-OR-A5J2, etc.)
        filepath_collection = joinpath(@__DIR__, "..", "collection", "$(collection_name).jsn")
        download_project_infos(filepath_collection, collection_name)
        project_id = extract_project_id(filepath_collection)
        filepath_case = joinpath(@__DIR__, "..", "case", "$(collection_name).jsn")
        download_project_infos(filepath_collection, collection_name)
        casesID_values, casesNAME_values = getCasesForProject(filepath_case, project_id)

        # Slides Management
        if isdir(joinpath(@__DIR__, "..", "slides", "svs", "$collection_name"))
            println("Update data ...")
        else
            mkdir(joinpath(@__DIR__, "..", "slides", "svs", "$collection_name"))
        end
        for (i, j) in zip(casesID_values, casesNAME_values)
            single_casesID_values, single_casesNAME_values = getCasesForProject(filepath_case, i)
            for (x, y) in zip(single_casesID_values, single_casesNAME_values)
                if isdir(joinpath(@__DIR__, "..", "slides", "svs", "$(collection_name)", "$j"))
                else
                    mkdir(joinpath(@__DIR__, "..", "slides", "svs", "$(collection_name)", "$j"))
                end
                filepath_slides = joinpath(@__DIR__, "..", "slides", "svs", "$(collection_name)", "$j", "$(y).zip")
                link_slides = "https://api.digitalslidearchive.org/api/v1/folder/$x/download"
                println("DOWNLOADING ... CASE NAME = $j - CASE ID = $i - SINGLE CASE NAME = $y")
                download_zip(link_slides, filepath_slides)
            end
        end
    else
        println("$collection_name collection not avaiable. Using collection_name in Config.toml")

        # Line of code for definition of "collection_name" from Config.toml
        filepath_config = joinpath(@__DIR__, "..", "Config.toml")
        config = TOML.parsefile(filepath_config)
        collection_name = config["collection_name"]
        if collection_name != "default" && lowercase(collection_name) in collection_list

            # Project Management (TCGA-OR-A5J1, TCGA-OR-A5J2, etc.)
            filepath_collection = joinpath(@__DIR__, "..", "collection", "$(collection_name).jsn")
            download_project_infos(filepath_collection, collection_name)
            project_id = extract_project_id(filepath_collection)
            filepath_case = joinpath(@__DIR__, "..", "case", "$(collection_name).jsn")
            download_project_infos(filepath_collection, collection_name)
            casesID_values, casesNAME_values = getCasesForProject(filepath_case, project_id)

            # Slides Management
            if isdir(joinpath(@__DIR__, "..", "slides", "svs", "$collection_name"))
                println("Update data ...")
            else
                mkdir(joinpath(@__DIR__, "..", "slides", "svs", "$collection_name"))
            end
            for (i, j) in zip(casesID_values, casesNAME_values)
                single_casesID_values, single_casesNAME_values = getCasesForProject(filepath_case, i)
                for (x, y) in zip(single_casesID_values, single_casesNAME_values)
                    if isdir(joinpath(@__DIR__, "..", "slides", "svs", "$(collection_name)", "$j"))
                    else
                        mkdir(joinpath(@__DIR__, "..", "slides", "svs", "$(collection_name)", "$j"))
                    end
                    filepath_slides = joinpath(@__DIR__, "..", "slides", "svs", "$(collection_name)", "$j", "$(y).zip")
                    link_slides = "https://api.digitalslidearchive.org/api/v1/folder/$x/download"
                    println("DOWNLOADING ... CASE NAME = $j - CASE ID = $i - SINGLE CASE NAME = $y")
                    download_zip(link_slides, filepath_slides)
                end
            end
        else 
            println("$collection_name collection not avaiable. Change collection_name field in Config.toml")
        end
    end
end

function download_all_collection()
    # Collection Management (acc, blca, etc.)
    filepath_collection = joinpath(@__DIR__, "..", "collection", "collectionlist.jsn")
    download_collection_values(filepath_collection)
    collection_list=print_collection_values(filepath_collection)
    
    for collection_name in collection_list
        # Project Management (TCGA-OR-A5J1, TCGA-OR-A5J2, etc.)
        #filepath_collection = "../collection/$collection_name.jsn"
        filepath_collection = joinpath(@__DIR__, "..", "collection", "$(collection_name).jsn")
        download_project_infos(filepath_collection, collection_name)
        project_id = extract_project_id(filepath_collection)
        filepath_case = joinpath(@__DIR__, "..", "case", "$(collection_name).jsn")
        download_project_infos(filepath_collection, collection_name)
        casesID_values, casesNAME_values = getCasesForProject(filepath_case, project_id)

        # Slides Management
        if isdir(joinpath(@__DIR__, "..", "slides", "svs", "$collection_name"))
        #if isdir("../slides/svs/$collection_name")
            println("Update data ...")
        else
            mkdir(joinpath(@__DIR__, "..", "slides", "svs", "$collection_name"))
            #mkdir("../slides/svs/$collection_name")
        end
        for (i, j) in zip(casesID_values, casesNAME_values)
            single_casesID_values, single_casesNAME_values = getCasesForProject(filepath_case, i)
            for (x, y) in zip(single_casesID_values, single_casesNAME_values)
                if isdir(joinpath(@__DIR__, "..", "slides", "svs", "$(collection_name)", "$j"))
                #if isdir("../slides/svs/$collection_name/$j")
                else
                    mkdir(joinpath(@__DIR__, "..", "slides", "svs", "$(collection_name)", "$j"))
                    #mkdir("../slides/svs/$collection_name/$j")
                end
                #filepath_slides = "../slides/svs/$collection_name/"*j*"/"*y*".zip"
                filepath_slides = joinpath(@__DIR__, "..", "slides", "svs", "$(collection_name)", "$j", "$(y).zip")
                link_slides = "https://api.digitalslidearchive.org/api/v1/folder/$x/download"
                println("DOWNLOADING ... CASE NAME = $j - CASE ID = $i - SINGLE CASE NAME = $y")
                download_zip(link_slides, filepath_slides)
            end
        end
    end
end
end
