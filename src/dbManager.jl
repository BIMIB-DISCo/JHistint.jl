export insert_record_DB
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

Funzione per la memorizzazione nel database `JHistint_DB` delle informazioni associate ad ogni slide scaricata dal CDSA.

# Argomenti
- `col_name::AbstractString` = Nome della collezione.
- `cas_name::AbstractString` = Nome del caso, coincide con il `CASE-NAME` stampato a video dal package.
- `tcga_case_id::AbstractString` = ID utilizzato dal TCGA per identificare il caso.
- `sin_cas_name::AbstractString` = Nome del singolo caso, coincide con lo `SLIDE-ID` stampato a video dal package.
- `tcga_slide_id::AbstractString` = ID utilizzato dal TCGA per identificare la slide.
- `link_slide::AbstractString` = Link alle API per il download della slide.
- `filepath_zip::AbstractString` = Percorso in cui è memorizzato il file `.zip`.
- `filepath_svs::AbstractString` = Percorso in cui è memorizzato il file `.svs`.

# Note
Dati disponibili nel database `JHistint_DB` per ogni slide:
- `collection_name TEXT` = Nome della collezione.
- `case_name TEXT` = Nome del caso.
- `TCGA_caseID TEXT` = ID utilizzato dal TCGA per identificare il caso.
- `slide_ID TEXT` = Nome del singolo caso.
- `TCGA_slideID TEXT UNIQUE` = ID utilizzato dal TCGA per identificare la slide, `UNIQUE` evita la generazione di duplicati.
- `slide_path_folder_zip TEXT` = Percorso in cui è memorizzato il file `.zip`.
- `slide_path_folder_svs TEXT` = Percorso in cui è memorizzato il file `.tif`.
- `slide_path_api TEXT` = Link alle API per il download della slide.
- `slide_path_folder_seg TEXT` = Percorso in cui è memorizzato il file `.tif` segmentata.
- `slide_svs BLOB` = Slide istopatologica (immagine).
- `slide_info_TSS TEXT` = Informazioni sulla slide - Tissue Source Site.
- `slide_info_participant_code TEXT` = Informazioni sulla slide - Codice associato al Participant, stringa alfanumerica.
- `slide_info_sample_type TEXT` = Informazioni sulla slide - Sample Type. I valori associati ai campioni aventi tumori sono nell'intervallo 01 - 09. 10 - 19 indica l'intervallo dedicato a campioni normali non malati. 20 - 29 indica campioni attualmente sotto controllo.
- `slide_info_vial TEXT` = Informazioni sulla slide - Vial. Relativo all'ordinamento del campione nella sequenza di campioni. I valori variano tra A - Z.
- `slide_info_portion TEXT` = Informazioni sulla slide - Portion. Relativo all'ordinamento delle porzioni analizzate associate ad un campione. Assume valori nell'intervallo 01-99.
- `slide_info_type TEXT` = Informazioni sulla slide - Tipo di immagine. I valori assumbili sono TS (Top Slide), BS (Bottom Slide) e MS (Middle Slide). Il valore alfanumerico indica l'ordinamento della slide.
- `slide_path_folder_matrix TEXT` = Percorso in cui è memorizzato il file `.txt` della matrice di adiacenza.
- `matrix_data BLOB` = Matrice di adiacenza.
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
    filepath_svs = joinpath(@__DIR__, "..", "output_example", "example.tif")
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
    query_extract_slide_svs(collection_name::AbstractString)

La funzione interroga il `JHistint_DB` ed estrae la lista di slide associate al nome della collezione fornita come argomento.

# Argomenti
- `collection_name::AbstractString`: Nome della collezione di slide da ricercare nel `JHistint_DB`.

# Valore di ritorno
- `slide_list`: Lista di tuple, ognuna delle quali contiene l'ID della slide, il file `.svs` della slide e il percorso della cartella contenente il file `.svs`.
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

La funzione aggiorna il `JHistint_DB` con il percorso del file dell'immagine segmentata, il percorso
del file con la matrice di adiacenza in formato testo e la matrice stessa.

# Argomento
- `filepath_seg::AbstractString`: Percorso del file dell'immagine segmentata da aggiungere al DB.
- `filepath_matrix::AbstractString`: Percorso del file della matrice di adiacenza.
- `matrix::Matrix{Int64}`: Matrice di adiacenza.
- `slide_id::AbstractString`: ID della slide da aggiornare con le informazioni dell'immagine segmentata.
"""
function load_seg_slide(filepath_seg::AbstractString, filepath_matrix::AbstractString, matrix::Matrix{Int64}, slide_id::AbstractString)
    db = SQLite.DB(joinpath(@__DIR__, "..", "JHistint_DB"))
    # seg_image = read(filepath_seg)
    stmt = SQLite.Stmt(db, "
       UPDATE Slide SET slide_path_folder_seg = ?, slide_path_folder_matrix = ?, matrix_data = ?
                          WHERE slide_ID = ?")
    DBInterface.execute(stmt, [filepath_seg, filepath_matrix, matrix, slide_id])
    SQLite.close(db)
end
