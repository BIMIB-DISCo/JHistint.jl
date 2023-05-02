export download_collection_values
export extract_collection_values
export download_project_infos
export extract_project_id
export download_zip
export getCasesForProject

"""
    download_collection_values(filepath::AbstractString)

Function for downloading data from collections available in TCGA.

# Arguments
- `filepath::AbstractString` = Path where to save the obtained `.json` file from the API available in CDSA.

# Notes
The API requires the definition of `parentType` and `parentId`. `parentId` specifies the identifier of the collection. The collection of images associated with TCGA is identified by the code:
`5b9ef8e3e62914002e454c39`. The use of `limit=0` sets the absence of limits in the queried file, ensuring the complete download of the file. The API belongs to the category for managing the folders
stored in the repository. The downloaded file is `.json`.
"""
function download_collection_values(filepath::AbstractString)
    # ID for TCGA Collection
    idTCGA = "5b9ef8e3e62914002e454c39"
    # Download collection file as JSON file from the server
    url = "https://api.digitalslidearchive.org/api/v1/folder?parentType=collection&parentId=$idTCGA&limit=0&sort=lowerName&sortdir=1"
    response = HTTP.get(url)
    if response.status == 200
         open(filepath, "w") do file
        write(file, response.body)
    end
    else
        println("Error: HTTP request returned status code $(response.status)")
    end
end

"""
    extract_collection_values(filepath::AbstractString)

Function to extract the values of data collections from the `.json` file downloaded by the `download_collection_values` function.

# Arguments
- `filepath::AbstractString` = Path where the `collectionlist.json` file is stored.

# Return value
- `collection_values::Array{String}` = List of data collections available in TCGA.
"""
function extract_collection_values(filepath::AbstractString)
    # Read the collection file and insert into a list the data of the collection
    json_string = read(filepath, String)
    json_object = JSON.parse(json_string)

    collection_values = []
    # For every elements of the JSON object
    for item in json_object
        # Check the presence of "name" value
        if haskey(item, "name")
            # Add the value of "name" to the list
            push!(collection_values, item["name"])
        end
    end
    # The names of the collection are avaiable on GitHub repository or Documentation
    # count=0
    # println("Choose the imaging collection from the list above: ")
    # for i in collection_values
    #    count += 1
    #    println(count," - TCGA-",i)
    # end
    return collection_values
end

# Funzione per acquisizione in input della collection_name da linea di comando.
# function select_collection_name(collection_values::Vector{Any})
#    while true
#        print("Insert the number of the collection of interest: ")
#        num = parse(Int, readline())
#        if 1 ≤ num ≤ length(collection_values)
#            println("Collection selected: TCGA-",collection_values[num])
#            return collection_values[num]
#        end
#        println("Error: Number not in range")
#    end
# end

"""
    download_project_infos(filepath::AbstractString, collection_name::AbstractString)

Function to download metadata associated with the selected collection at startup.

# Arguments
- `filepath::AbstractString` = Path to save the `.json` file associated with the collection. The file is indicated
with the wording `collection_name.json`.
- `collection_name::AbstractString` = Name of the collection from which to download the slides.

# Notes
The API requires the definition of `parentType`, `parentId`, and `name`. The `name` attribute identifies the name of the collection from which you want to retrieve data (e.g., chol, esca, etc.).
The API belongs to the category for managing the folders stored in the repository. The downloaded file is `.json`.
"""
function download_project_infos(filepath::AbstractString, collection_name::AbstractString)
    # Download project file as JSON file from the server
    url = "https://api.digitalslidearchive.org/api/v1/folder?parentType=collection&parentId=5b9ef8e3e62914002e454c39&name=$collection_name&sort=lowerName&sortdir=1"
    response = HTTP.get(url)
    if response.status == 200
         open(filepath, "w") do file
        write(file, response.body)
    end
    else
        println("Error: HTTP request returned status code $(response.status)")
    end
end

"""
    extract_project_id(filepath::AbstractString)

Function to extract the `id` value from the metadata of the collection selected at startup.

# Arguments
- `filepath::AbstractString` = Path where the `collection_name.json` file is stored.

# Return value
- `project_id` = `id` of the collection.
"""
function extract_project_id(filepath::AbstractString)
    # Read the project info file and insert into a variable the data of the project ID
    json_string = read(filepath, String)
    json_object = JSON.parse(json_string)
    project_id=""
    # For every elements of the JSON object
    for item in json_object
        # Check the presence of "_id" value
        if haskey(item, "_id")
            # Add the value of "name" to the variable
            project_id = item["_id"]
        end
    end
    return project_id
end

"""
    getCasesForProject(filepath_case::AbstractString, project_id::AbstractString)

Function to download metadata associated with the cases of the selected collection at startup.

# Arguments
- `filepath::AbstractString` = Path where to save the `.json` file associated with the cases of the collection. The file is indicated with the term `collection_name.json`.
- `project_id::AbstractString` = `id` of the collection.

# Return values
- `casesID_values` = List of `id` of all the cases in the collection.
- `casesNAME_values` = List of `name` of all the cases in the collection.

# Notes
The API requires the definition of `parentType` and `parentId`. The `parentType` attribute is set to `folder` given the structure of the repository.
The `parentId` is set by defining the identifier of the chosen collection. The downloaded file is `.json`.
"""
function getCasesForProject(filepath_case::AbstractString, project_id::AbstractString)
    # Download case file as JSON file from the server
    url = "https://api.digitalslidearchive.org/api/v1/folder?parentType=folder&parentId=$project_id&limit=0&sort=lowerName&sortdir=1"
    response = HTTP.get(url)
    if response.status == 200
         open(filepath_case, "w") do file
        write(file, response.body)
    end
    else
        println("Error: HTTP request returned status code $(response.status)")
    end
    # Read the case file and insert into a list the data of the cases name
    json_string = read(filepath_case, String)
    json_object = JSON.parse(json_string)

    casesID_values = []
    casesNAME_values = []
    # For every elements of the JSON object
    for item in json_object
        # Check the presence of "_id" value
        if haskey(item, "_id")
            # Add the value of "_id" to the list
            push!(casesID_values, item["_id"])
            push!(casesNAME_values, item["name"])
        end
    end
    return casesID_values, casesNAME_values
end

"""
    download_zip(link::AbstractString, filepath::AbstractString)

Function for downloading histological slides in `.zip` format associated with the cases of the selected collection at startup.

# Arguments
- `link::AbstractString` = URL to access the API for slide download.
- `filepath::AbstractString` = Path to save the `.zip` file.
"""
function download_zip(link::AbstractString, filepath::AbstractString)
    response = HTTP.get(link)
    if response.status == 200
         open(filepath, "w") do file
        write(file, response.body)
    end
    else
        println("Error: HTTP request returned status code $(response.status)")
    end
end
