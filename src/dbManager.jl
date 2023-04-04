export insert_record_DB

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
- `slide_path_folder_svs TEXT` = Percorso in cui è memorizzato il file `.svs`.
- `slide_path_api TEXT` = Link alle API per il download della slide.
- `slide_svs BLOB` = Slide istolpatologica (immagine).
- `slide_info_TSS TEXT` = Informazioni sulla slide - Tissue Source Site.
- `slide_info_participant_code TEXT` = Informazioni sulla slide - Codice associato al Participant, stringa alfanumerica.
- `slide_info_sample_type TEXT` = Informazioni sulla slide - Sample Type. I valori associati ai campioni aventi tumori sono nell'intervallo 01 - 09. 10 - 19 indica l'intervallo dedicato a campioni normali non malati. 20 - 29 indica campioni attualmente sotto controllo.
- `slide_info_vial TEXT` = Informazioni sulla slide - Vial. Relativo all'ordinamento del campione nella sequenza di campioni. I valori variano tra A - Z.
- `slide_info_portion TEXT` = Informazioni sulla slide - Portion. Relativo all'ordinamento delle porzioni analizzate associate ad un campione. Assume valori nell'intervallo 01-99.
- `slide_info_type TEXT` = Informazioni sulla slide - Tipo di immagine. I valori assumbili sono TS (Top Slide), BS (Bottom Slide) e MS (Middle Slide). Il valore alfanumerico indica l'ordinamento della slide.
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
    # filepath_svs = "C:/Users/nicom/Desktop/JHistint.jl/slides/svs/acc/TCGA-OR-A5J1/TCGA-OR-A5J1-01A-01-TS1/Cattura.PNG"
    # svs_image = read(filepath_svs)

    # Connect to DB
    db = SQLite.DB("JHistint_DB")
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
                                        slide_info_type TEXT)")

     stmt = SQLite.Stmt(db, "
        INSERT OR REPLACE INTO Slide (collection_name,
                           case_name,
                           TCGA_caseID,
                           slide_ID,
                           TCGA_slideID,
                           slide_path_folder_zip,
                           slide_path_folder_svs,
                           slide_path_api,
                           slide_info_TSS,
                           slide_info_participant_code,
                           slide_info_sample_type,
                           slide_info_vial,
                           slide_info_portion,
                           slide_info_type) VALUES
                           (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)")
    DBInterface.execute(stmt, [col_name,
                               cas_name,
                               tcga_case_id,
                               sin_cas_name,
                               tcga_slide_id,
                               filepath_zip,
                               filepath_svs,
                               link_slide,
                               TSS,
                               participant_code,
                               sample_type,
                               vial,
                               portion,
                               type])
    # Show tables in the database
    # SQLite.tables(db) # Empty Database with no tables
end
