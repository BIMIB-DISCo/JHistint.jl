using HTTP

filepath_slides = joinpath(@__DIR__, "file.zip")
link = "https://api.digitalslidearchive.org/api/v1/folder/5b9f4d2ce62914002e95db26/download"

HTTP.open(:GET, link) do http
    open(filepath_slides, "w") do file
        write(file, http)
        close(file)
    end

end
