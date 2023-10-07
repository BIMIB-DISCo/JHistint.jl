### -*- Mode: Julia -*-

### DB Manager -- JHistint
### dbManager.jl

### Exported Functions
export insert_record_DB
export insert_record_DB_SOPHYSM
export query_extract_slide_svs
export load_seg_slide

"""
    insert_record_DB(col_name::AbstractString,
                          cas_name::AbstractString,
                          tcga_case_id::AbstractString,
                          sin_cas_name::AbstractString,
                          tcga_slide_id::AbstractString,
                          link_slide::AbstractString,
                          filepath_zip::AbstractString,
                          filepath_svs::AbstractString)

Function for storing in the `JHistint_DB` database the information associated
with each slide downloaded from the Cancer Digital Slide Archive (CDSA).

# Argomenti
- `col_name::AbstractString` = Collection name.
- `cas_name::AbstractString` = Case name, which corresponds to the `CASE-NAME`
displayed by the package.
- `tcga_case_id::AbstractString` = ID used by TCGA to identify the case.
- `sin_cas_name::AbstractString` = Name of the individual slide, which
corresponds to the `SLIDE-ID` displayed by the package.
- `tcga_slide_id::AbstractString` = ID used by TCGA to identify the slide.
- `link_slide::AbstractString` = Link to the APIs for downloading the slide.
- `filepath_zip::AbstractString` = Path where the `.zip` file is stored.
- `filepath_svs::AbstractString` = Path where the `.tif` file is stored.

# Notes
The `JHistint_DB` database is used for storing the information associated with
each slide downloaded from the CDSA.
The function takes a dictionary containing the information associated with the
slide and stores it in the database.
Data available in the `JHistint_DB` database for each slide:
- `collection_name TEXT` = Name of the collection.
- `case_name TEXT` = Name of the case.
- `TCGA_caseID TEXT` = ID used by TCGA to identify the case.
- `slide_ID TEXT` = Name of the individual slide case.
- `TCGA_slideID TEXT UNIQUE` = ID used by TCGA to identify the slide, `UNIQUE`
prevents duplicates from being generated.
- `slide_path_folder_zip TEXT` = Path where the `.zip` file is stored.
- `slide_path_folder_svs TEXT` = Path where the `.tif` file is stored.
- `slide_path_api TEXT` = Link to the API for downloading the slide.
- `slide_path_folder_seg TEXT` = Path where the segmented `.tif` file is stored.
- `slide_svs BLOB` = Histopathological slide (image).
- `slide_info_TSS TEXT` = Slide information - Tissue Source Site.
- `slide_info_participant_code TEXT` = Slide information - Participant Code,
alphanumeric string.
- `slide_info_sample_type TEXT` = Slide information - Sample Type. The values
associated with tumor samples are in the range 01-09. 10-19 indicates the range
for non-diseased normal samples. 20-29 indicates samples currently under control.
- `slide_info_vial TEXT` = Slide information - Vial. Related to the ordering of
the sample in the sequence of samples. The values range from A-Z.
- `slide_info_portion TEXT` = Slide information - Portion. Related to the
ordering of the analyzed portions associated with a sample. Takes values in the
range 01-99.
- `slide_info_type TEXT` = Slide information - Image Type. The possible
values are TS (Top Slide), BS (Bottom Slide), and MS (Middle Slide).
The alphanumeric value indicates the ordering of the slide.
- `slide_path_folder_matrix TEXT` = Path where the adjacency matrix `.txt`
file is stored.
- `matrix_data BLOB` = Adjacency matrix.
"""
function insert_record_DB(col_name::AbstractString,
                          cas_name::AbstractString,
                          tcga_case_id::AbstractString,
                          sin_cas_name::AbstractString,
                          tcga_slide_id::AbstractString,
                          link_slide::AbstractString,
                          filepath_zip::AbstractString,
                          filepath_svs::AbstractString)
    # Extract data form sin_cas_name
    res = split(sin_cas_name, "-")
    data = []
    for i in res
        push!(data, i)
    end
    tcga, TSS, participant_code, sample_type_vial, portion, type = data
    sample_type = sample_type_vial[1:2]
    vial = sample_type_vial[3:3]
    # Extract data from image slide
    # Setting for DEMO
    filepath_config = joinpath(@__DIR__, "..", "Config.toml")
    config = TOML.parsefile(filepath_config)
    demo = config["demo"]
    if (demo == 1)
        if (sin_cas_name == "TCGA-18-3406-01A-01-BS1")
            filepath_svs = joinpath(@__DIR__, "..", "input_example_demo",
                                    "slideExample1", "SlideExample_mini_1.tif")
        elseif (sin_cas_name == "TCGA-18-3406-01A-01-TS1")
            filepath_svs = joinpath(@__DIR__, "..", "input_example_demo",
                                    "slideExample2", "SlideExample_mini_2.tif")
        else
            filepath_svs = joinpath(@__DIR__, "..", "input_example_demo",
                                    "slideExample3", "SlideExample_mini_3.tif")
        end
    end
    svs_image = read(filepath_svs)

    # Connect to DB
    db = SQLite.DB(joinpath(@__DIR__, "..", "JHistint_DB"))
    # Create a Table(Slide)
    SQLite.execute(db, "CREATE TABLE IF NOT EXISTS Slide(collection_name TEXT,
                                        case_name TEXT,
                                        TCGA_caseID TEXT,
                                        slide_ID TEXT,
                                        TCGA_slideID TEXT UNIQUE,
                                        slide_path_folder_zip TEXT,
                                        slide_path_folder_svs TEXT,
                                        slide_path_api TEXT,
                                        slide_path_folder_seg TEXT,
                                        slide_svs BLOB,
                                        slide_info_TSS TEXT,
                                        slide_info_participant_code TEXT,
                                        slide_info_sample_type TEXT,
                                        slide_info_vial TEXT,
                                        slide_info_portion TEXT,
                                        slide_info_type TEXT,
                                        slide_path_folder_matrix TEXT,
                                        matrix_data BLOB)")
     stmt = SQLite.Stmt(db, "
        INSERT OR REPLACE INTO Slide (collection_name,
                           case_name,
                           TCGA_caseID,
                           slide_ID,
                           TCGA_slideID,
                           slide_path_folder_zip,
                           slide_path_folder_svs,
                           slide_path_api,
                           slide_svs,
                           slide_info_TSS,
                           slide_info_participant_code,
                           slide_info_sample_type,
                           slide_info_vial,
                           slide_info_portion,
                           slide_info_type) VALUES
                           (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)")

    DBInterface.execute(stmt, [col_name,
                               cas_name,
                               tcga_case_id,
                               sin_cas_name,
                               tcga_slide_id,
                               filepath_zip,
                               filepath_svs,
                               link_slide,
                               svs_image,
                               TSS,
                               participant_code,
                               sample_type,
                               vial,
                               portion,
                               type])
    SQLite.close(db)
end

"""
    insert_record_DB_SOPHYSM(col_name::AbstractString,
                          cas_name::AbstractString,
                          tcga_case_id::AbstractString,
                          sin_cas_name::AbstractString,
                          tcga_slide_id::AbstractString,
                          link_slide::AbstractString,
                          filepath_zip::AbstractString,
                          filepath_svs::AbstractString,
                          path_download_db::AbstractString)

Function for storing in the `JHistint_DB` database the information associated
with each slide downloaded from the Cancer Digital Slide Archive (CDSA).
The function is different from the standard model. It stores only informations
of the downloaded Slides.

# Argomenti
- `col_name::AbstractString` = Collection name.
- `cas_name::AbstractString` = Case name, which corresponds to the `CASE-NAME`
displayed by the package.
- `tcga_case_id::AbstractString` = ID used by TCGA to identify the case.
- `sin_cas_name::AbstractString` = Name of the individual slide, which
corresponds to the `SLIDE-ID` displayed by the package.
- `tcga_slide_id::AbstractString` = ID used by TCGA to identify the slide.
- `link_slide::AbstractString` = Link to the APIs for downloading the slide.
- `filepath_zip::AbstractString` = Path where the `.zip` file is stored.
- `filepath_svs::AbstractString` = Path where the `.tif` file is stored.
- `path_download_db::AbstractString` = Path where the DB file is stored.

# Notes
The `JHistint_DB` database is used for storing the information associated with
each slide downloaded from the CDSA.
The function takes a dictionary containing the information associated with
    the slide and stores it in the database.
Data available in the `JHistint_DB` database for each slide:
- `collection_name TEXT` = Name of the collection.
- `case_name TEXT` = Name of the case.
- `TCGA_caseID TEXT` = ID used by TCGA to identify the case.
- `slide_ID TEXT` = Name of the individual slide case.
- `TCGA_slideID TEXT UNIQUE` = ID used by TCGA to identify the slide,
`UNIQUE` prevents duplicates from being generated.
- `slide_path_folder_zip TEXT` = Path where the `.zip` file is stored.
- `slide_path_folder_svs TEXT` = Path where the `.tif` file is stored.
- `slide_path_api TEXT` = Link to the API for downloading the slide.
- `slide_path_folder_seg TEXT` = Path where the segmented `.tif` file is stored.
- `slide_svs BLOB` = Histopathological slide (image).
- `slide_info_TSS TEXT` = Slide information - Tissue Source Site.
- `slide_info_participant_code TEXT` = Slide information - Participant Code,
alphanumeric string.
- `slide_info_sample_type TEXT` = Slide information - Sample Type.
The values associated with tumor samples are in the range 01-09. 10-19
indicates the range for non-diseased normal samples. 20-29 indicates samples
    currently under control.
- `slide_info_vial TEXT` = Slide information - Vial. Related to the ordering
of the sample in the sequence of samples. The values range from A-Z.
- `slide_info_portion TEXT` = Slide information - Portion. Related to the
ordering of the analyzed portions associated with a sample. Takes values
in the range 01-99.
- `slide_info_type TEXT` = Slide information - Image Type. The possible
values are TS (Top Slide), BS (Bottom Slide), and MS (Middle Slide).
The alphanumeric value indicates the ordering of the slide.
- `slide_path_folder_matrix TEXT` = Path where the adjacency matrix
`.txt` file is stored.
- `matrix_data BLOB` = Adjacency matrix.
"""
function insert_record_DB_SOPHYSM(col_name::AbstractString,
                          cas_name::AbstractString,
                          tcga_case_id::AbstractString,
                          sin_cas_name::AbstractString,
                          tcga_slide_id::AbstractString,
                          link_slide::AbstractString,
                          filepath_zip::AbstractString,
                          filepath_svs::AbstractString,
                          path_download_db::AbstractString)
    # Extract data form sin_cas_name
    res = split(sin_cas_name, "-")
    data = []
    for i in res
        push!(data, i)
    end
    tcga, TSS, participant_code, sample_type_vial, portion, type = data
    sample_type = sample_type_vial[1:2]
    vial = sample_type_vial[3:3]
    # Extract data from image slide
    svs_image = read(filepath_svs)

    # Connect to DB
    db = SQLite.DB(joinpath(path_download_db, "JHistint_DB"))
    # Create a Table(Slide)
    SQLite.execute(db, "CREATE TABLE IF NOT EXISTS Slide(collection_name TEXT,
                                        case_name TEXT,
                                        TCGA_caseID TEXT,
                                        slide_ID TEXT,
                                        TCGA_slideID TEXT UNIQUE,
                                        slide_path_folder_zip TEXT,
                                        slide_path_folder_svs TEXT,
                                        slide_path_api TEXT,
                                        slide_svs BLOB,
                                        slide_info_TSS TEXT,
                                        slide_info_participant_code TEXT,
                                        slide_info_sample_type TEXT,
                                        slide_info_vial TEXT,
                                        slide_info_portion TEXT,
                                        slide_info_type TEXT,
                                        slide_path_folder_matrix TEXT,
                                        matrix_data BLOB)")
     stmt = SQLite.Stmt(db, "
        INSERT OR REPLACE INTO Slide (collection_name,
                           case_name,
                           TCGA_caseID,
                           slide_ID,
                           TCGA_slideID,
                           slide_path_folder_zip,
                           slide_path_folder_svs,
                           slide_path_api,
                           slide_svs,
                           slide_info_TSS,
                           slide_info_participant_code,
                           slide_info_sample_type,
                           slide_info_vial,
                           slide_info_portion,
                           slide_info_type) VALUES
                           (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)")

    DBInterface.execute(stmt, [col_name,
                               cas_name,
                               tcga_case_id,
                               sin_cas_name,
                               tcga_slide_id,
                               filepath_zip,
                               filepath_svs,
                               link_slide,
                               svs_image,
                               TSS,
                               participant_code,
                               sample_type,
                               vial,
                               portion,
                               type])
    SQLite.close(db)
end


"""
    query_extract_slide_svs(collection_name::AbstractString)

The function queries the `JHistint_DB` and extracts the list of slides
associated with the collection name provided as an argument.

# Arguments
- `collection_name::AbstractString`: Name of the slide collection to search
for in the `JHistint_DB`.

# Return value
- `slide_list`: List of tuples, each of which contains the ID of the slide,
the `.svs` file of the slide, and the path of the
folder containing the `.svs` file.
"""
function query_extract_slide_svs(collection_name::AbstractString)
    db = SQLite.DB(joinpath(@__DIR__, "..", "JHistint_DB"))
    stmt = SQLite.Stmt(db, "SELECT * FROM Slide WHERE collection_name LIKE '%' || ? || '%'")
    results = DataFrame(DBInterface.execute(stmt, [collection_name]))
    slide_list = []
    if !isempty(results)
        for row in eachrow(results)
            slide = (row.slide_ID, row.slide_svs, row.slide_path_folder_svs)
            push!(slide_list, slide)
        end
    end
    SQLite.close(db)
    return slide_list
end

"""
    load_seg_slide(filepath_seg::AbstractString, filepath_matrix::AbstractString, matrix::Matrix{Int64}, slide_id::AbstractString)

The function updates the `JHistint_DB` with the path of the segmented image file,
the path of the adjacency matrix file in text format, and the matrix itself.

# Arguments
- `filepath_seg::AbstractString`: Path of the segmented
image file to add to the DB.
- `filepath_matrix::AbstractString`: Path of the adjacency matrix file.
- `matrix::Matrix{Int64}`: Adjacency matrix.
- `slide_id::AbstractString`: ID of the slide to update with
the segmented image information.
"""
function load_seg_slide(filepath_seg::AbstractString, filepath_matrix::AbstractString, matrix::Matrix{Int64}, slide_id::AbstractString)
    db = SQLite.DB(joinpath(@__DIR__, "..", "JHistint_DB"))
    # seg_image = read(filepath_seg)
    stmt = SQLite.Stmt(db, "
       UPDATE Slide SET slide_path_folder_seg = ?, slide_path_folder_matrix = ?,
                        matrix_data = ? WHERE slide_ID = ?")
    DBInterface.execute(stmt, [filepath_seg, filepath_matrix, matrix, slide_id])
    SQLite.close(db)
end

### end of file -- dbManager.jl
