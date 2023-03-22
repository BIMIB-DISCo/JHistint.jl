export download_collection_values
export print_collection_values
export download_project_infos
export extract_project_id
export download_zip
export getCasesForProject
export select_collection_name

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

function print_collection_values(filepath::AbstractString)
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
    # The names of the collection are avaiable on GitHub repository
    # count=0
    # println("Choose the imaging collection from the list above: ")
    # for i in collection_values
    #    count += 1
    #    println(count," - TCGA-",i)
    # end    
    return collection_values
end

function select_collection_name(collection_values::Vector{Any})
    while true
        print("Insert the number of the collection of interest: ")
        num = parse(Int, readline())
        if 1 ≤ num ≤ length(collection_values)
            println("Collection selected: TCGA-",collection_values[num])
            return collection_values[num]
        end
        println("Error: Number not in range")
    end 
end

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

# Download a slides in .zip format given the API URL and the folder path
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
