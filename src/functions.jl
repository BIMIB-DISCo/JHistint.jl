export download_collection_values
export extract_collection_values
export download_project_infos
export extract_project_id
export download_zip
export getCasesForProject

"""
    download_collection_values(filepath::AbstractString)

Funzione per il download dei dati delle collezioni disponibili nel TCGA. 

# Argomenti
- `filepath::AbstractString` = Percorso in cui salvare il file .json ottenuto dalle API disponibili nel CDSA.

# Note
L'API richiede la definizione del `parentType` e del `parentId`. Il `parentId` specifica l'identificativo della collezione. La collezione di immagini associate al TCGA è identificata dal codice: `5b9ef8e3e62914002e454c39`. L'utilizzo del `limit=0` imposta l'assenza di limiti nel file interrogato, garantendo il download del file in modo completo. L'API appartiene alla categoria per gestire le folder memorizzate nel repository. Il file scaricato è `.json`.
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

Funzione per estrarre i valori delle collezioni di dati dal file `.json` scaricato dalla funzione `download_collection_values`.

# Argomenti
- `filepath::AbstractString` = Percorso in cui è memorizzato il file `collectionlist.json`.

# Valore di ritorno
- `collection_values::Array{String}` = Lista di collezioni di dati disponibili nel TCGA.
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

Funzione per il download dei metadati associati alla collezione selezionata all'avvio.

# Argomenti
- `filepath::AbstractString` = Percorso in cui salvare il file `.json` associato alla collezione. Il file è indicato
con la dicitura `collection_name.json`.
- `collection_name::AbstractString` = Nome della collezione di cui scaricare le slides.

# Note
L'API richiede la definizione del `parentType`, del `parentId` e del `name`. L'attributo `name` identifica il nome della collezione di cui si desidera prelevare i dati (esempio: chol, esca, etc.). L'API appartiene alla categoria per gestire le folder memorizzate nel repository. Il file scaricato è `.json`.
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

Funzione per estrarre il valore del `id` dai metadati della collezione selezionata all'avvio.

# Argomenti
- `filepath::AbstractString` =  Percorso in cui è memorizzato il file `collection_name.json`.

# Valori di ritorno
- `project_id` =  `id` della collezione.
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

Funzione per il download dei metadati associati ai casi della collezione selezionata all'avvio.

# Argomenti
- `filepath::AbstractString` =  Percorso in cui salvare il file `.json` associato ai casi della collezione. Il file è indicato
con la dicitura `collection_name.json`.
- `project_id::AbstractString` =  `id` della collezione.

# Valori di ritorno
- `casesID_values` =  Lista di `id` di tutti i casi della collezione.
- `casesNAME_values` = Lista di `name` di tutti i casi della collezione.

# Note
L'API richiede la definizione del `parentType` e del `parentId`. L'attributo `parentType` è impostato a `folder` data la struttura del repository. Il `parentId` è configurato definendo l'identificativo della collezione scelta. Il file scaricato è `.json`.
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

Funzione per il download delle slides istologiche in formato `.zip` associate ai casi della collezione selezionata all'avvio.

# Argomenti
- `link::AbstractString` = URL per l'accesso alla API per il download delle slides.
- `filepath::AbstractString` =  Percorso in cui salvare il file `.zip`.
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
