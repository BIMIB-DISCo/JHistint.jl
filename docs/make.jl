using Documenter
using JHistint
using HTTP
using JSON
using ZipFile
using TOML
using SQLite
using DataFrames
using Images
using ImageSegmentation
using ImageView
using FileIO
using Random

makedocs(
    sitename = "JHistint Documentation",
    format = Documenter.HTML(),
    modules = [JHistint],
    authors = "Niccolò Mandelli",
)
