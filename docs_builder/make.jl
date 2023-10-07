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
using ImageMagick
using ImageView
using FileIO
using Random
using IndirectArrays
using Graphs
using LightGraphs
using SimpleWeightedGraphs
using CSV
using J_Space
using Luxor
using Karnak
using MetaGraphs

makedocs(
    sitename = "JHistint Documentation",
    format = Documenter.HTML(),
    modules = [JHistint],
    authors = "Niccolò Mandelli",
    pages = [
        "Home" => "index.md"
        ]
)
