export extract_slide

"""
    extract_slide(filepath_zip::AbstractString)

Function to extract the contents of `.zip` files downloaded from CDSA.

# Arguments
- `filepath_zip::AbstractString` = Path where the `.zip` file for the individual case is saved.
"""
function extract_slide(filepath_zip::AbstractString)
    dest_path = ""
    # Open .zip file
    archive = ZipFile.Reader(filepath_zip)
    # Extract .zip content
    for entry in archive.files
        # Build destination path - Change here for image type
        file_name = replace(entry.name, ".svs" => ".tif")
        dest_path = joinpath(dirname(filepath_zip), file_name)
        # dest_path = joinpath(dirname(filepath_zip), entry.name)
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
