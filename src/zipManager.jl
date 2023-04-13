export extract_slide

"""
    extract_slide(filepath_zip::AbstractString)

Funzione per estrarre il contenuto dei file .zip scaricati dal CDSA.

# Argomenti
- `filepath_zip::AbstractString` = Percorso in cui è salvato il file .zip relativo al singolo caso.
"""
function extract_slide(filepath_zip::AbstractString)
    dest_path = ""
    # Open .zip file
    archive = ZipFile.Reader(filepath_zip)
    # Extract .zip content
    for entry in archive.files
        # Build destination path - Change here for image type
        file_name = replace(entry.name, ".svs" => ".jpeg")
        dest_path = joinpath(dirname(filepath_zip), file_name)
        if !isdir(dirname(dest_path))
            mkdir(dirname(dest_path))
        end
        # Open destination file and insert the slide
        open(dest_path, "w") do dest_file
            write(dest_file, read(entry))
        end
    end
    # Close .zip file
    close(archive)
    return dest_path
end
