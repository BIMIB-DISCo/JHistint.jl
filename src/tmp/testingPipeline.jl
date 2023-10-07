# Packages
using HTTP
using JSON
using ZipFile
using TOML
using SQLite

using Images
using ImageSegmentation
using ImageMagick
using ImageView
using FileIO

using IndirectArrays
using SimpleWeightedGraphs
using CSV
using DelimitedFiles
using Distances
using Base

using MetaGraphs
using Graphs
using DataFrames
using Random
using NetworkLayout
using CairoMakie
using GraphMakie

# plotting
using GraphPlot
using Karnak
using Colors
using Luxor
using Compose
using LinearAlgebra

# deprecated
using Meshes
# new
using VoronoiCells
using GeometryBasics
using Plots

"""
`SegmentedImage` type contains the index-label mapping, assigned labels,
segment mean intensity and pixel count of each segment.
struct SegmentedImage{T<:AbstractArray,U<:Union{Colorant,Real}}
    image_indexmap::T
    segment_labels::Vector{Int}
    segment_means::Dict{Int,U}
    segment_pixel_count::Dict{Int,Int}
end
"""
@static if Base.VERSION <= v"1.0.5"
    # https://github.com/JuliaLang/julia/pull/29442
    _oneunit(::CartesianIndex{N}) where {N} = _oneunit(CartesianIndex{N})
    _oneunit(::Type{CartesianIndex{N}}) where {N} = CartesianIndex(ntuple(x -> 1, Val(N)))
else
    const _oneunit = Base.oneunit
end

# Once Base has colon defined here we can replace this
_colon(I::CartesianIndex{N}, J::CartesianIndex{N}) where N =
    CartesianIndices(map((i,j) -> i:j, Tuple(I), Tuple(J)))

# struct Point{T}
#    x::T
#    y::T
#     z::T
#end

"""
    G, vert_map = region_adjacency_graph(seg, weight_fn)

Constructs a region adjacency graph (RAG) from the `SegmentedImage`. It returns the RAG
along with a Dict(label=>vertex) map. `weight_fn` is used to assign weights to the edges.

    weight_fn(label1, label2)

Returns a real number corresponding to the weight of the edge between label1 and label2.

"""
function region_adjacency_graph(s::SegmentedImage, weight_fn::Function)

    function neighbor_regions!(df_cartesian_indices::AbstractArray, G::SimpleWeightedGraph, visited::AbstractArray, s::SegmentedImage, I::CartesianIndex)
        # n = Set{Int} - visited = Array - s = segmented image - p = CartesianIndex which define neighbors
        # R contains each possible index in s
        R = CartesianIndices(axes(s.image_indexmap))
        # I1 contains a Vector of only 1 with dimension equal to visited
        I1 = _oneunit(CartesianIndex{ndims(visited)})
        # Ibegin and Iend contains the first and last index of R
        Ibegin, Iend = first(R), last(R)
        # t is only a empty Vector with dimension equal to visited
        t = Vector{CartesianIndex{ndims(visited)}}()
        # add index I to Vector t
        push!(t, I)
        while !isempty(t)
            # Extract last element of t and save it in temp
            temp = pop!(t)
            # set index temp to true
            visited[temp] = true
            # _colon build an object CartesianIndices which include all the index from I to J (range) :
            # _colon(I::CartesianIndex{N}, J::CartesianIndex{N}) where N =
            #    CartesianIndices(map((i,j) -> i:j, Tuple(I), Tuple(J)))
            for J in _colon(max(Ibegin, temp-I1), min(Iend, temp+I1))
                if s.image_indexmap[temp] != s.image_indexmap[J]
                    # if s.image_indexmap[J] ∉ n
                        # If the values are different, it means they have two different colorings for the two points,
                        # therefore a neighbor has been identified, which is pushed into n.
                        # push!(n,s.image_indexmap[J])
                    if !Graphs.has_edge(G, vert_map[s.image_indexmap[I]], vert_map[s.image_indexmap[J]])
                        Graphs.add_edge!(G, vert_map[s.image_indexmap[I]], vert_map[s.image_indexmap[J]], weight_fn(s.image_indexmap[I], s.image_indexmap[J]))
                        # push!(added_indices, s.image_indexmap[J])
                        # push!(n,s.image_indexmap[J])
                    end
                elseif !visited[J]
                    # If they are equal, I place them in t, so that,
                    # as long as t is not empty, I can explore all the neighbors
                    # that have the same color.
                    push!(t,J)
                end
            end
        end
    end
    # Start
    visited  = fill(false, axes(s.image_indexmap))                           # Array to mark the pixels that are already visited
    # G        = SimpleWeightedGraph()                                         # The region_adjacency_graph
    # vert_map = Dict{Int,Int}()                                               # Map that stores (label, vertex) pairs
    # Build object for label (vertex) dataframe
    df_label = DataFrame()
    df_noisy_label = DataFrame()
    df_total_label = DataFrame()
    df_cartesian_indices = []
    df_color_indices = []
    # add vertices to graph
    # Graphs.add_vertices!(G,length(s.segment_labels))
    # setup `vert_map`
    # for (i,l) in enumerate(s.segment_labels)
    #    vert_map[l] = i
    # end
    # add edges to graph
    # For each CartesianIndices in s where the image_indexmap represent the image wich is a Multidimensional Array
    # The index of s.image_indexmap represent the pixel position in the segmented image
    # The value of s.image_indexmap represent the pixel color in the segmented image
    for p in CartesianIndices(axes(s.image_indexmap))
        # check if p of the segmented image s is not visited
        if !visited[p]
            push!(df_cartesian_indices, p)
        end
    end

    for i in s.segment_labels
        push!(df_color_indices, s.segment_means[i])
    end

    # Filter noisy segments from segmentation
    # Building object for DataFrame df_label
    df_label_filtered = Int[]
    df_cartesian_indices_filtered = CartesianIndex[]
    df_color_indices_filtered = []
    df_area_filtered = Int[]

    df_label_filtered_extra = Int[]
    df_cartesian_indices_filtered_extra = CartesianIndex[]
    df_color_indices_filtered_extra = []
    df_area_filtered_extra = Int[]

    df_label_total = Int[]
    df_cartesian_indices_total = CartesianIndex[]
    df_color_indices_total = []
    df_area_total = Int[]
    for i in 1:length(s.segment_labels)
        if(s.segment_pixel_count[i] > 3000)
            push!(df_label_filtered, s.segment_labels[i])
            push!(df_cartesian_indices_filtered, df_cartesian_indices[i])
            push!(df_color_indices_filtered, df_color_indices[i])
            push!(df_area_filtered, s.segment_pixel_count[i])
        end
        if(s.segment_pixel_count[i] <= 3000 &&
            s.segment_pixel_count[i] > 300)
            push!(df_label_filtered_extra, s.segment_labels[i])
            push!(df_cartesian_indices_filtered_extra, df_cartesian_indices[i])
            push!(df_color_indices_filtered_extra, df_color_indices[i])
            push!(df_area_filtered_extra, s.segment_pixel_count[i])
        end
        if(s.segment_pixel_count[i] > 300)
            push!(df_label_total, s.segment_labels[i])
            push!(df_cartesian_indices_total, df_cartesian_indices[i])
            push!(df_color_indices_total, df_color_indices[i])
            push!(df_area_total, s.segment_pixel_count[i])
        end
    end
    df_label.label = df_label_filtered
    df_label.position_label = df_cartesian_indices_filtered
    df_label.color_label = df_color_indices_filtered
    df_label.area = df_area_filtered

    df_noisy_label.label = df_label_filtered_extra
    df_noisy_label.position_label = df_cartesian_indices_filtered_extra
    df_noisy_label.color_label = df_color_indices_filtered_extra
    df_noisy_label.area = df_area_filtered_extra

    df_total_label.label = df_label_total
    df_total_label.position_label = df_cartesian_indices_total
    df_total_label.color_label = df_color_indices_total
    df_total_label.area = df_area_total

    df_label, df_noisy_label, df_total_label
end

"""
    compute_centroid_total_cells(s::SegmentedImage, df_label::DataFrame)

# Arguments:

# Return value:

# Notes:
"""
function compute_centroid_total_cells(s::SegmentedImage, df_label::DataFrame)
    # Start
    visited  = fill(false, axes(s.image_indexmap)) # Array to mark the pixels that are already visited
    df_centroids = CartesianIndex[]
    df_label_list = Int[]
    df_cells_indices = []
    for p in CartesianIndices(axes(s.image_indexmap))
        # 0.6854407198549426
        if s.segment_means[s.image_indexmap[p]] != 0
            if !visited[p]
                push!(df_cells_indices, p)
                visited[p] = true
                for j in CartesianIndices(axes(s.image_indexmap))
                    if s.image_indexmap[p] == s.image_indexmap[j] && !visited[j]
                        push!(df_cells_indices, j)
                        visited[j] = true
                    end
                end
                if(s.segment_pixel_count[s.image_indexmap[p]] > 300) # manual threshold for cells dimension in pixel (delete noise)
                    centroid = div(s.segment_pixel_count[s.image_indexmap[p]], 2)
                    if(!(df_cells_indices[centroid] in df_centroids))
                        push!(df_label_list, s.image_indexmap[p])
                        push!(df_centroids, df_cells_indices[centroid])
                    end
                end
                df_cells_indices = []
            end
        end
    end
    nuclei_list = df_label[!, :label]
    position = 99
    df_centroid_ordered = [CartesianIndex(0,0,0) for _ in 1:length(nuclei_list)]
    for index in 1:length(df_label_list)
        pointer = df_label_list[index]
        for j in 1:length(nuclei_list)
            if(nuclei_list[j] == pointer)
                position = j
            end
        end
        df_centroid_ordered[position] = df_centroids[index]
    end
    df_label.centroid = df_centroid_ordered
    df_label
end

"""
    compute_centroid_cells(s::SegmentedImage, df_label::DataFrame)

# Arguments:

# Return value:

# Notes:
"""
function compute_centroid_cells(s::SegmentedImage, df_label::DataFrame)
    # Start
    visited  = fill(false, axes(s.image_indexmap)) # Array to mark the pixels that are already visited
    df_centroids = CartesianIndex[]
    df_cells_indices = []
    df_label_list = Int[]
    for p in CartesianIndices(axes(s.image_indexmap))
        # 0.6854407198549426
        if s.segment_means[s.image_indexmap[p]] != 0
            if !visited[p]
                push!(df_cells_indices, p)
                visited[p] = true
                for j in CartesianIndices(axes(s.image_indexmap))
                    if s.image_indexmap[p] == s.image_indexmap[j] && !visited[j]
                        push!(df_cells_indices, j)
                        visited[j] = true
                    end
                end
                if(s.segment_pixel_count[s.image_indexmap[p]] > 3000) # manual threshold for cells dimension in pixel (delete noise)
                    centroid = div(s.segment_pixel_count[s.image_indexmap[p]], 2)
                    if(!(df_cells_indices[centroid] in df_centroids))
                        push!(df_label_list, s.image_indexmap[p])
                        push!(df_centroids, df_cells_indices[centroid])
                    end
                end
                df_cells_indices = []
            end
        end
    end
    nuclei_list = df_label[!, :label]
    position = 99
    df_centroid_ordered = [CartesianIndex(0,0,0) for _ in 1:length(nuclei_list)]
    for index in 1:length(df_label_list)
        pointer = df_label_list[index]
        for j in 1:length(nuclei_list)
            if(nuclei_list[j] == pointer)
                position = j
            end
        end
        df_centroid_ordered[position] = df_centroids[index]
    end
    df_label.centroid = df_centroid_ordered
    df_label
   df_label
end

"""
    compute_centroid_noisy_cells(s::SegmentedImage, df_label::DataFrame)

# Arguments:

# Return value:

# Notes:
"""
function compute_centroid_noisy_cells(s::SegmentedImage, df_label::DataFrame)
    # Start
    visited  = fill(false, axes(s.image_indexmap)) # Array to mark the pixels that are already visited
    df_centroids = CartesianIndex[]
    df_cells_indices = []
    df_label_list = Int[]
    for p in CartesianIndices(axes(s.image_indexmap))
        # 0.6854407198549426
        if s.segment_means[s.image_indexmap[p]] != 0
            if !visited[p]
                push!(df_cells_indices, p)
                visited[p] = true
                for j in CartesianIndices(axes(s.image_indexmap))
                    if s.image_indexmap[p] == s.image_indexmap[j] && !visited[j]
                        push!(df_cells_indices, j)
                        visited[j] = true
                    end
                end
                if(s.segment_pixel_count[s.image_indexmap[p]] <= 3000 &&
                    s.segment_pixel_count[s.image_indexmap[p]] > 300) # manual threshold for cells dimension in pixel (delete noise)
                    centroid = div(s.segment_pixel_count[s.image_indexmap[p]], 2)
                    if(!(df_cells_indices[centroid] in df_centroids))
                        push!(df_label_list, s.image_indexmap[p])
                        push!(df_centroids, df_cells_indices[centroid])
                    end
                end
                df_cells_indices = []
            end
        end
    end
    nuclei_list = df_label[!, :label]
    position = 99
    df_centroid_ordered = [CartesianIndex(0,0,0) for _ in 1:length(nuclei_list)]
    for index in 1:length(df_label_list)
        pointer = df_label_list[index]
        for j in 1:length(nuclei_list)
            if(nuclei_list[j] == pointer)
                position = j
            end
        end
        df_centroid_ordered[position] = df_centroids[index]
    end
    df_label.centroid = df_centroid_ordered
    df_label
   df_label
end

"""
    filter_dataframe_cells(df_label::DataFrame)

# Arguments:

# Return value:

# Notes:
"""
function filter_dataframe_cells(df_label::DataFrame)
    df_label_filtered = Int[]
    df_cartesian_indices_filtered = CartesianIndex[]
    df_color_indices_filtered = []
    df_area_filtered = Int[]
    df_centroid_filtered = CartesianIndex[]
    df_filtered = DataFrame()
    for i in eachrow(df_label)
        label = i.label
        position_label = i.position_label
        color_label = i.color_label
        area = i.area
        centroid = i.centroid
        if(area > 3000)
            push!(df_label_filtered, label)
            push!(df_cartesian_indices_filtered, position_label)
            push!(df_color_indices_filtered, color_label)
            push!(df_area_filtered, area)
            push!(df_centroid_filtered, centroid)
        end
    end
    df_filtered.label = df_label_filtered
    df_filtered.position_label = df_cartesian_indices_filtered
    df_filtered.color_label = df_color_indices_filtered
    df_filtered.area = df_area_filtered
    df_filtered.centroid = df_centroid_filtered
    df_filtered
end

"""
    filter_dataframe_extras(df_label::DataFrame)

# Arguments:

# Return value:

# Notes:
"""
function filter_dataframe_extras(df_label::DataFrame)
    df_label_filtered = Int[]
    df_cartesian_indices_filtered = CartesianIndex[]
    df_color_indices_filtered = []
    df_area_filtered = Int[]
    df_centroid_filtered = CartesianIndex[]
    df_filtered = DataFrame()
    for i in eachrow(df_label)
        label = i.label
        position_label = i.position_label
        color_label = i.color_label
        area = i.area
        centroid = i.centroid
        if(area <= 3000 && area > 300)
            push!(df_label_filtered, label)
            push!(df_cartesian_indices_filtered, position_label)
            push!(df_color_indices_filtered, color_label)
            push!(df_area_filtered, area)
            push!(df_centroid_filtered, centroid)
        end
    end
    df_filtered.label = df_label_filtered
    df_filtered.position_label = df_cartesian_indices_filtered
    df_filtered.color_label = df_color_indices_filtered
    df_filtered.area = df_area_filtered
    df_filtered.centroid = df_centroid_filtered
    df_filtered
end


"""
    add_column_is_cell(df_labels::DataFrame, df_noisy_labels::DataFrame, df_total_labels::DataFrame)

# Arguments:

# Return value:

# Notes:
"""
function add_column_is_cell(df_labels::DataFrame, df_noisy_labels::DataFrame, df_total_labels::DataFrame)
    nuclei_list = df_labels[!, :label]
    extra_list = df_noisy_labels[!, :label]
    total_list = df_total_labels[!, :label]
    is_cell_list = Bool[]
    for label_index in total_list
        if label_index in nuclei_list
            push!(is_cell_list, true)
        end
        if label_index in extra_list
            push!(is_cell_list, false)
        end
    end
    df_total_labels.is_cell = is_cell_list
    df_total_labels
end

"""
    build_dataframe_edges_from_grid(edge_list::Vector{Any}, df_total_labels::DataFrame)

# Arguments:

# Return value:

# Notes:
"""
function build_dataframe_edges_from_grid(edge_list::Vector{Any}, df_total_labels::DataFrame)
    # Build edge_list with label identifier
    cell_label_list = df_total_labels[!, :label]
    cell_area_list = df_total_labels[!, :area]
    edge_label = []
    edge_weight = []
    count = 1
    for pair_label in edge_list
        first_value = pair_label[1]
        second_value = pair_label[2]
        origin = cell_label_list[first_value]
        destination = cell_label_list[second_value]
        push!(edge_label, (origin, destination))
        if(cell_area_list[first_value] > cell_area_list[second_value])
            push!(edge_weight, (cell_area_list[first_value] - cell_area_list[second_value]))
        else
            push!(edge_weight, (cell_area_list[second_value] - cell_area_list[first_value]))
        end
        count = count + 1
    end
    # Build df_edges
    df = DataFrame()
    df_origin_column = []
    df_destination_column = []
    df_weight_column = []
    for pair_label in edge_label
        push!(df_origin_column, pair_label[1])
        push!(df_destination_column, pair_label[2])
    end
    for weight in edge_weight
        push!(df_weight_column, weight)
    end
    df.origin = df_origin_column
    df.destination = df_destination_column
    df.weight = df_weight_column
    df
end

"""
    build_graph_from_grid(df_labels::DataFrame,
                               df_noisy_labels::DataFrame,
                               df_total_labels::DataFrame,
                               w::Int64, h::Int64)

# Arguments:

# Return value:

# Notes:
"""
function build_graph_from_grid(df_labels::DataFrame,
                               df_noisy_labels::DataFrame,
                               df_total_labels::DataFrame,
                               w::Int64, h::Int64,
                               filepath_total_tess::AbstractString,
                               filepath_cell_tess::AbstractString)
    nuclei_list = df_labels[!, :centroid]
    extra_list = df_noisy_labels[!, :centroid]
    total_list = df_total_labels[!, :centroid]
    nuclei_label_list = df_labels[!, :label]
    cell_slot = length(nuclei_list)
    extra_slot = length(extra_list)
    total_slot = length(total_list)
    # VORONOICELLS IMPLEMENTATION
    position_array = Vector{GeometryBasics.Point2{Float64}}()
    cell_position_array = Vector{GeometryBasics.Point2{Float64}}()
    extra_position_array = Vector{GeometryBasics.Point2{Float64}}()
    for vertex in nuclei_list
        coordinates_str = match(r"\((.*)\)", string(vertex)).captures[1]
        coordinates = parse.(Int, split(coordinates_str, ", "))
        x, y ,z = coordinates
        point = GeometryBasics.Point2(y, -x)
        push!(cell_position_array, point)
    end
    for vertex in extra_list
        coordinates_str = match(r"\((.*)\)", string(vertex)).captures[1]
        coordinates = parse.(Int, split(coordinates_str, ", "))
        x, y ,z = coordinates
        point = GeometryBasics.Point2(y, -x)
        push!(extra_position_array, point)
    end
    for vertex in total_list
        coordinates_str = match(r"\((.*)\)", string(vertex)).captures[1]
        coordinates = parse.(Int, split(coordinates_str, ", "))
        x, y ,z = coordinates
        point = GeometryBasics.Point2(y, -x)
        push!(position_array, point)
    end
    # Check Position in Rectangle
    # s = ""
    # for e in position_array
    #    show(e)
    #    show(h+1)
    #    show(w+1)
    #    println()
    #    if(!(e[1] > 0 && e[1] < h))
    #        s = e[1]
    #        show("error: $s")
    #        throw("find : $s")
    #    end
    #    if(!(e[2] > 0 && e[2] < w))
    #        s = e[2]
    #        show("error: $s")
    #        throw("find : $s")
    #    end
    # send
    h = h + 1
    w = w + 1
    rect = VoronoiCells.Rectangle(GeometryBasics.Point2(0, 0),
                                  GeometryBasics.Point2(h, -w))
    edges = []
    tess = VoronoiCells.voronoicells(position_array, rect; edges)
    Plots.scatter(cell_position_array, markersize = 4, label = "Nuclei Centroid")
    # Plot centroid of extra cell
    # Plots.scatter(extra_position_array, markersize = 4, label = "Extra Centroid")
    Plots.annotate!([(cell_position_array[n][1],
                    cell_position_array[n][2],
                    Plots.text(nuclei_label_list[n])) for n in 1:cell_slot])
    # Plot number of extra cell
    # Plots.annotate!([(extra_position_array[n][1],
    #                extra_position_array[n][2],
    #                Plots.text(n)) for n in 1:extra_slot])
    plot_tessellation = Plots.plot!(tess, legend = :topleft)
    savefig(plot_tessellation, filepath_total_tess)
    df_edges = build_dataframe_edges_from_grid(edges, df_total_labels)
    # gui()
    rect = VoronoiCells.Rectangle(GeometryBasics.Point2(0, 0),
                                  GeometryBasics.Point2(h, -w))
    tess = VoronoiCells.voronoicells(cell_position_array, rect)
    Plots.scatter(cell_position_array, markersize = 4, label = "Nuclei Centroid")
    Plots.annotate!([(cell_position_array[n][1],
                    cell_position_array[n][2],
                    Plots.text(nuclei_label_list[n])) for n in 1:cell_slot])
    plot_tessellation = Plots.plot!(tess, legend = :topleft)
    savefig(plot_tessellation, filepath_cell_tess)
    return df_edges, edges

    # MATRIX IMPLEMENTATION
    # matrix = false(w, h)
    # show(matrix)
    # for p in CartesianIndex(matrix)
    #     if(p in nuclei_list)
    #         matrix[p] = true
    #     end
    # end

    # MESHES IMPLEMENTATION
    # grid = Meshes.CartesianGrid((h, w), (0.,0.), (1.,1.))
    # initialize
    # for p in grid
    #     show(p)
    #     println()
    # end
end

"""
    weighted_graph_to_adjacency_matrix(G::SimpleWeightedGraph{Int64, Float64}, n::Int64)

Converts a weighted graph represented as a `SimpleWeightedGraph` into an unweighted boolean adjacency matrix.

# Arguments:
- `G::SimpleWeightedGraph{Int64, Float64}`: Weighted graph represented as a `SimpleWeightedGraph` with integer vertex labels and floating-point edge weights.
- `n::Int64`: Number of nodes in the adjacency matrix.

# Return value:
- `adjacency_matrix`: `Matrix{Int64}` boolean adjacency matrix.

# Notes:
The function returns an `n` x `n` adjacency matrix representing the unweighted graph. If nodes `i` and `j` are adjacent,
the adjacency matrix will contain a value of `1` at position `(i,j)` and `(j,i)`. Otherwise, the adjacency matrix will contain a value of `0`.
"""
function weighted_graph_to_adjacency_matrix(G::SimpleWeightedGraph{Int64, Float64}, n::Int64)
    adjacency_matrix = zeros(Int, n, n)
    for i in 1:n
        for j in 1:n
            if Graphs.has_edge(G, i, j)
                adjacency_matrix[i, j] = 1
                adjacency_matrix[j, i] = 1
            end
        end
    end
    return adjacency_matrix
end

"""
    weighted_graph_to_adjacency_matrix_weight(G::SimpleWeightedGraph{Int64, Float64}, n::Int64)

# Arguments:

# Return value:

# Notes:
"""
function weighted_graph_to_adjacency_matrix_weight(G::SimpleWeightedGraph{Int64, Float64}, n::Int64)
    adjacency_matrix = zeros(Float32, n, n)
    for i in 1:n
        for j in 1:n
            adjacency_matrix[i, j] = -1
            adjacency_matrix[j, i] = -1
        end
    end

    for i in 1:n
        for j in 1:n
            if Graphs.has_edge(G, i, j)
                adjacency_matrix[i, j] = SimpleWeightedGraphs.get_weight(G, i, j)
                adjacency_matrix[j, i] = SimpleWeightedGraphs.get_weight(G, i, j)
            end
        end
    end
    return adjacency_matrix
end

"""
    tess_dataframe_to_adjacency_matrix_weight(df_total_labels::DataFrame, df_edges::DataFrame, edges::Vector{Any})

# Arguments:

# Return value:

# Notes:
"""
function tess_dataframe_to_adjacency_matrix_weight(df_total_labels::DataFrame, df_edges::DataFrame, edges::Vector{Any})
    cell_label_list = df_total_labels[!, :label]
    # origin_list = df_edges[!, :origin]
    # destination_list = df_edges[!, :destination]
    weight_list = df_edges[!, :weight]
    n = length(cell_label_list)
    adjacency_matrix = zeros(Int, n, n)
    for i in 1:n
        for j in 1:n
            adjacency_matrix[i, j] = -1
            adjacency_matrix[j, i] = -1
        end
    end
    count = 1
    for label_pair in edges
        adjacency_matrix[label_pair[1], label_pair[2]] =  weight_list[count]
        count = count + 1
    end
    return adjacency_matrix
end

"""
    build_dataframe_as_edgelist(mat::Matrix{Int64}, label_list::Vector{Int64})

# Arguments:

# Return value:

# Notes:
"""
function build_dataframe_as_edgelist(mat::Matrix{Int64}, label_list::Vector{Int64})
    df = DataFrame()
    df_origin_column = []
    df_destination_column = []
    df_weight_column = []
    n = size(mat, 1)
    for i in 1:n
        for j in 1:i
            if mat[i,j] != 0 && i in label_list && j in label_list
                push!(df_origin_column, i)
                push!(df_destination_column, j)
                push!(df_weight_column, mat[i,j])
            end
        end
    end
    df.origin = df_origin_column
    df.destination = df_destination_column
    df.weight = df_weight_column
    return df
end

"""
    build_dataframe_as_edgelist(mat::Matrix{Float32}, label_list::Vector{Int64})

# Arguments:

# Return value:

# Notes:
"""
function build_dataframe_as_edgelist(mat::Matrix{Float32}, label_list::Vector{Int64})
    df = DataFrame()
    df_origin_column = []
    df_destination_column = []
    df_weight_column = []
    n = size(mat, 1)
    for i in 1:n
        for j in 1:i
            if mat[i,j] != -1 && i in label_list && j in label_list
                push!(df_origin_column, i)
                push!(df_destination_column, j)
                push!(df_weight_column, mat[i,j])
            end
        end
    end
    df.origin = df_origin_column
    df.destination = df_destination_column
    df.weight = df_weight_column
    return df
end

"""
    save_adjacency_matrix(matrix::Matrix{Int64}, filepath_matrix::AbstractString)

Function to save an adjacency matrix represented as an integer matrix to a text file.

# Arguments:
- `matrix::Matrix{Int64}`: The integer matrix representing the adjacency matrix.
- `filepath_matrix::AbstractString`: The file path represented as a string indicating where to save the matrix.

# Notes
The function opens the file specified by the `filepath_matrix` path in write mode and writes the matrix in adjacency matrix format,
where each row represents the adjacent nodes of a node. The numbers in the matrix are separated by spaces.
"""
function save_adjacency_matrix(matrix::Matrix{Int64}, filepath_matrix::AbstractString)
    # Open txt file
    f = open(filepath_matrix, "w")
    # Write matrix dimensions
    n_rows, n_cols = size(matrix)
    matrix_string = ""
    for i in 1:n_rows
        write(f, " ")
        for j in 1:n_cols
            matrix_string = string("", matrix[i, j], "  ")
            write(f, matrix_string)
        end
        write(f, "\n")
    end
    close(f)
end

"""
    save_adjacency_matrix(matrix::Matrix{Float32}, filepath_matrix::AbstractString)

# Arguments:

# Return value:

# Notes:
"""
function save_adjacency_matrix(matrix::Matrix{Float32}, filepath_matrix::AbstractString)
    # Open txt file
    f = open(filepath_matrix, "w")
    # Write matrix dimensions
    n_rows, n_cols = size(matrix)
    matrix_string = ""
    for i in 1:n_rows
        write(f, " ")
        for j in 1:n_cols
            matrix_string = string("", matrix[i, j], "  ")
            write(f, matrix_string)
        end
        write(f, "\n")
    end
    close(f)
end

function metagraph_from_dataframe_JHistint(graph_type,
                                  df::DataFrame,
                                  origin::Symbol,
                                  destination::Symbol,
                                  weight::Symbol = Symbol(),
                                  edge_attributes::Union{Vector{Symbol}, Symbol}
                                  = Vector{Symbol}(),
                                  vertex_attributes::DataFrame = DataFrame(),
                                  vertex_id_col::Symbol = Symbol())
    # Map node names to vertex IDs
    nodes = sort!(unique!([df[:, origin]; df[:, destination]]))
    vertex_names = DataFrame(name=nodes, vertex_id = eachindex(nodes))
    # Merge in to original
    for c in [origin, destination]
        temp = rename(vertex_names, :vertex_id => Symbol(c, :_id))
        df = innerjoin(df, temp; on = c => :name)
    end
    # Merge additional attributes to names
    if vertex_attributes != DataFrame()
        idsym =
            vertex_id_col == Symbol() ?
            first(propertynames(vertex_attributes)) :
            vertex_id_col
        vertex_names = leftjoin(vertex_names,
                                vertex_attributes,
                                on = :name => idsym)
    end
    # Create Graph
    mg = graph_type(nrow(vertex_names))
    for r in eachrow(df)
        MetaGraphs.add_edge!(mg, r[Symbol(origin, :_id)], r[Symbol(destination, :_id)])
    end

    df[!, :Child] .= ""
    df[!, :Time] .= ""
    df[!, :Subpop_Child] .= ""
    # add times on graph
    times_df = df[:,[:Child, :Time]]
    vertex_names = leftjoin(vertex_names, times_df, on = :name => :Child)
    # add subpop_child on graph
    subpop_df = df[:,[:Child, :Subpop_Child]]
    vertex_names = leftjoin(vertex_names, subpop_df, on = :name => :Child)
    # Set vertex names and attributes
    attr_names = propertynames(vertex_names[!, Not(:vertex_id)])
    for r in eachrow(vertex_names)
        MetaGraphs.set_props!(mg, r[:vertex_id], Dict([a => r[a] for a in attr_names]))
    end
    # Set edge attributes
    if typeof(edge_attributes) == Symbol
        edge_attributes = Vector{Symbol}([edge_attributes])
    end
    origin_id = Symbol(origin, :_id)
    destination_id = Symbol(destination, :_id)
    for e in edge_attributes
        for r in eachrow(df)
            MetaGraphs.set_prop!(mg, MetaGraphs.Edge(r[origin_id], r[destination_id]), e, r[e])
        end
    end
    # Weight
    if weight != Symbol()
        for r in eachrow(df)
            set_prop!(mg,
                      MetaGraphs.Edge(r[origin_id], r[destination_id]),
                      :weight, r[weight])
        end
    end
    # Set name as index (Issue #9)
    set_indexing_prop!(mg, :name)
    return mg
end

function spatial_graph(path_dataframe_edges::String, path_dataframe_labels::String)
    df_edges_file = CSV.File(path_dataframe_edges)
    df_edges = DataFrame(df_edges_file)

    df_labels_file = CSV.File(path_dataframe_labels)
    df_labels = DataFrame(df_labels_file)
    G_meta = metagraph_from_dataframe_JHistint(MetaGraph, df_edges, :origin, :destination, :weight, :weight, df_labels, :label)
    return G_meta
end

function plot_lattice_JHistint(G::MetaGraph, Set_mut::Vector{Any}; dim::Int=3)
    driver_mut, labels, colors = get_drivermut_name_colors(G, Set_mut)
    mylayout = NetworkLayout.Spectral(dim=3)
    f, ax, p = graphplot(G,
                         layout = mylayout,
                         node_size = repeat([5], nv(G)),
                         edge_width=1.0,
                         node_color = colors)
    return f, ax, p, colors
end

function plot_lattice_metagraph(G::MetaGraph; dim::Int=3)
    # driver_mut, labels, colors = get_drivermut_name_colors(G, Set_mut)
    mylayout = NetworkLayout.Spectral(dim=2)
    f, ax, p = graphplot(G,
                         layout = mylayout,
                         node_size = repeat([5], nv(G)),
                         edge_width=1.0)
    return f, ax, p
end

function extract_vertex_position(G::MetaGraph)
    position_array = Luxor.Point[]
    for v in Graphs.vertices(g_meta)
        s = get_prop(g_meta, v, :centroid)
        coordinates_str = match(r"\((.*)\)", string(s)).captures[1]
        coordinates = parse.(Int, split(coordinates_str, ", "))
        x, y ,z = coordinates
        point = Luxor.Point(y, x)
        push!(position_array, point)
    end
    return position_array
end

function extract_vertex_color(G::MetaGraph)
    color_array = []
    for v in Graphs.vertices(g_meta)
        color_float = get_prop(g_meta, v, :color_label)
        #color_RGB = convert(RGB{N0f8}, HSL(color_float, color_float, color_float))
        push!(color_array, color_float)
    end
    return color_array
end


link = joinpath(@__DIR__, "..", "input_example_demo", "slideExample1", "SlideExample_mini_1.tif")
# load slide
println("LOAD SLIDE ...")
svs_image = read(link)
img = ImageMagick.load_(svs_image)

# define img dimensions
width = size(img, 1)
height = size(img, 2)

bw = Gray.(img) .> 0.15
dist = 1 .- distance_transform(feature_transform(bw))
markers = label_components(dist .< -0.3)
# watershed
println("APPLY SEGMENTATION ...")
segments = watershed(dist, markers)
# build segmented slide
println("BUILD SEGMENTED SLIDE ...")
labels = labels_map(segments)
colored_labels = IndirectArray(labels, distinguishable_colors(maximum(labels)))
masked_colored_labels = colored_labels .* (1 .- bw)
# build graph
println("BUILD GRAPH ...")
weight_fn(i,j) = euclidean(segment_pixel_count(segments,i), segment_pixel_count(segments,j))

df_labels = DataFrame()
df_noisy_labels = DataFrame()
df_total_labels = DataFrame()
df_edges = DataFrame()

df_labels, df_noisy_labels, df_total_labels = region_adjacency_graph(segments, weight_fn)
# define centroid
# df_labels = compute_centroid_cells(segments, df_labels)
# df_noisy_labels = compute_centroid_noisy_cells(segments, df_noisy_labels)
df_total_labels = compute_centroid_total_cells(segments, df_total_labels)
df_labels = filter_dataframe_cells(df_total_labels)
df_noisy_labels = filter_dataframe_extras(df_total_labels)
df_total_labels = add_column_is_cell(df_labels, df_noisy_labels, df_total_labels)
nvertices = length(vert_map)

# build graph - grid with cells
filepath_total_tess = replace(link, r"....$" => "_total_tessellation.png")
filepath_cell_tess = replace(link, r"....$" => "_cell_tessellation.png")
df_edges, edges = build_graph_from_grid(df_labels, df_noisy_labels, df_total_labels, width, height, filepath_total_tess, filepath_cell_tess)

# save dataframe label as .CSV
filepath_dataframe_labels = replace(link, r"....$" => "_dataframe_labels.csv")
CSV.write(filepath_dataframe_labels, df_total_labels)
# build and save adjacency matrix
println("BUILD ADJACENCY MATRIX ...")
# matrix = weighted_graph_to_adjacency_matrix_weight(G, nvertices)
matrix = tess_dataframe_to_adjacency_matrix_weight(df_total_labels, df_edges, edges)
filepath_matrix = replace(link, r"....$" => ".txt")
save_adjacency_matrix(matrix, filepath_matrix)
# build and save dataframe edgelist as .CSV
println("BUILD DATAFRAME EDGELIST ...")
# df_edges = build_dataframe_as_edgelist(matrix, df_labels.label)
filepath_dataframe_edges = replace(link, r"....$" => "_dataframe_edges.csv")
CSV.write(filepath_dataframe_edges, df_edges)

filepath_img_graph_vertex = replace(link, r"....$" => "_graph_vertex.png")
filepath_img_graph_edges = replace(link, r"....$" => "_graph_edges.png")
filepath_png = replace(link, r"....$" => "_seg-0.png")
filepath_seg = replace(link, r"....$" => "_seg.png")
save(filepath_seg, masked_colored_labels)
# svs_png = read(filepath_png)
# img = ImageMagick.load_(svs_png)
img_graph = Luxor.readpng(filepath_png)
w = img_graph.width
h = img_graph.height

# J-Space (temp metagraph for building img on graph)
g_meta = spatial_graph(filepath_dataframe_edges, filepath_dataframe_labels)
# Plot metagraph on slide
# Image with Vertices
@png begin
    Luxor.placeimage(img_graph, 0, 0, 0.8, centered=true)
    sethue("slateblue")
    Karnak.fontsize(7)
    drawgraph(g_meta,
        layout = extract_vertex_position(g_meta) .+ Karnak.Point(-w/2, -h/2),
        vertexlabels = [get_prop(g_meta, v, :name) for v in Graphs.vertices(g_meta)],
        vertexfillcolors = extract_vertex_color(g_meta),
        edgelines=:none
    )
end w h filepath_img_graph_vertex
# Image with Edges
@png begin
    Luxor.placeimage(img_graph, 0, 0, 0.8, centered=true)
    sethue("slateblue")
    Karnak.fontsize(7)
    drawgraph(g_meta,
        layout = extract_vertex_position(g_meta) .+ Karnak.Point(-w/2, -h/2),
        vertexlabels = [get_prop(g_meta, v, :name) for v in Graphs.vertices(g_meta)],
        vertexfillcolors = extract_vertex_color(g_meta),
    )
end w h filepath_img_graph_edges

#@svg begin
    # Luxor.Drawing(w, h)
#    Luxor.placeimage(image, 0, 0, 0.8)
#    sethue("slateblue")
#    move(0, 0)
#    Karnak.fontsize(8)
#    drawgraph(g_meta,
#        layout = extract_vertex_position(g_meta),
        # layout = [get_prop(g_meta, v, :position_label) for v in vertices(g_meta)],
#        vertexlabels = [get_prop(g_meta, v, :name) for v in vertices(g_meta)],
        # vertexlabels = [(degree(g_meta, n) == 1) ? find(n) : "" for n in vertices(g_meta)],
        # vertexfillcolors =
        #     [RGB(rand()/2, rand()/2, rand()/2)
        #       for i in 1:nv(g_meta)],
#        vertexfillcolors = extract_vertex_color(g_meta),
        # edgelines=:none
        # edgelabels = [G.weights[e.src, e.dst] for e in edges(G)]
#    )
#end 1000 1000 filepath_svg
