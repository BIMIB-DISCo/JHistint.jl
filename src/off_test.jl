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

struct Point{T}
    x::T
    y::T
    z::T
end

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
    G        = SimpleWeightedGraph()                                         # The region_adjacency_graph
    vert_map = Dict{Int,Int}()                                               # Map that stores (label, vertex) pairs
    # Build object for label (vertex) dataframe
    df_label = DataFrame()
    df_cartesian_indices = []
    df_color_indices = []
    # add vertices to graph
    Graphs.add_vertices!(G,length(s.segment_labels))
    # setup `vert_map`
    for (i,l) in enumerate(s.segment_labels)
        vert_map[l] = i
    end
    # add edges to graph
    # For each CartesianIndices in s where the image_indexmap represent the image wich is a Multidimensional Array
    # The index of s.image_indexmap represent the pixel position in the segmented image
    # The value of s.image_indexmap represent the pixel color in the segmented image
    for p in CartesianIndices(axes(s.image_indexmap))
        # check if p of the segmented image s is not visited
        if !visited[p]
            push!(df_cartesian_indices, p)
            # n = Set{Int16}()
            # Call neighbor_regions where :
            # n = Set{Int} - visited = Array - s = segmented image - p = CartesianIndex which define neighbors
            try
                neighbor_regions!(df_cartesian_indices, G, visited, s, p)
            catch oom
                if isa(oom, OutOfMemoryError)
                    # n = Set{Int}()
                    GC.gc()
                    println(">>> OOM")
                    exit()
                end
            end
            # for i in n
            #    Graphs.add_edge!(G, vert_map[s.image_indexmap[p]], vert_map[i], weight_fn(s.image_indexmap[p], i))
            # end
        end
    end

    for i in s.segment_labels
        push!(df_color_indices, s.segment_means[i])
    end
    df_label.label = s.segment_labels
    df_label.position_label = df_cartesian_indices
    df_label.color_label = df_color_indices
    G, vert_map, df_label
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
    build_dataframe_as_edgelist(mat::Matrix{Int64})

# Arguments:

# Return value:

# Notes:
"""
function build_dataframe_as_edgelist(mat::Matrix{Int64})
    df = DataFrame()
    df_origin_column = []
    df_destination_column = []
    df_weight_column = []
    n = size(mat, 1)
    for i in 1:n
        for j in 1:i
            if mat[i,j] != 0
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
    build_dataframe_as_edgelist(mat::Matrix{Float32})

# Arguments:

# Return value:

# Notes:
"""
function build_dataframe_as_edgelist(mat::Matrix{Float32})
    df = DataFrame()
    df_origin_column = []
    df_destination_column = []
    df_weight_column = []
    n = size(mat, 1)
    for i in 1:n
        for j in 1:i
            if mat[i,j] != -1
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
    println(mg)
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

function vertex_color(v)
    props = g.vprops[v]
    color = get_prop(props, :color, "grey")
    return color
end

function extract_position(G::MetaGraph)
    position_array = Point{}[]
    for v in vertices(g_meta)
        s = get_prop(g_meta, v, :position_label)
        coordinates_str = match(r"\((.*)\)", string(s)).captures[1]
        coordinates = parse.(Int, split(coordinates_str, ", "))
        x, y ,z = coordinates
        point = Point{Float64}(x, y, z)
        push!(position_array, point)
    end
    return position_array
end

link = joinpath(@__DIR__, "..", "input_example_demo", "slideExample1", "SlideExample_mini_1.tif")
# load slide
println("LOAD SLIDE ...")
svs_image = read(link)
img = ImageMagick.load_(svs_image)
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
G, vert_map, df_labels = region_adjacency_graph(segments, weight_fn)
nvertices = length(vert_map)
# save dataframe label as .CSV
filepath_dataframe_labels = replace(link, r"....$" => "_dataframe_labels.csv")
CSV.write(filepath_dataframe_labels, df_labels)
# build and save adjacency matrix
println("BUILD ADJACENCY MATRIX ...")
matrix = weighted_graph_to_adjacency_matrix_weight(G, nvertices)
filepath_matrix = replace(link, r"....$" => ".txt")
save_adjacency_matrix(matrix, filepath_matrix)
# build and save dataframe edgelist as .CSV
println("BUILD DATAFRAME EDGELIST ...")
df_edges = build_dataframe_as_edgelist(matrix)
filepath_dataframe_edges = replace(link, r"....$" => "_dataframe_edges.csv")
CSV.write(filepath_dataframe_edges, df_edges)
filepath_svg = replace(link, r"....$" => ".svg")
# filepath_png = replace(link, r"....$" => ".png")
# svs_png = read(filepath_png)
# img = ImageMagick.load_(svs_png)

# J-Space
g_meta = spatial_graph(filepath_dataframe_edges, filepath_dataframe_labels)

# Plot metagraph on slide
# plot_lattice_JHistint(g_meta)
point_array = extract_position(g_meta)
show(point_array)
@svg begin
    # Drawing(M)
    background("black")
    sethue("grey40")
    Karnak.fontsize(8)
    drawgraph(g_meta,
        layout = point_array,
        # layout = [get_prop(g_meta, v, :position_label) for v in vertices(g_meta)],
        vertexlabels = [get_prop(g_meta, v, :name) for v in vertices(g_meta)],
        # vertexlabels = [(degree(g_meta, n) == 1) ? find(n) : "" for n in vertices(g_meta)],
        # vertexfillcolors =
        #     [RGB(rand()/2, rand()/2, rand()/2)
        #       for i in 1:nv(g_meta)],
        # vertexfillcolors = [RGB(get_prop(g_meta, v, :color_label), get_prop(g_meta, v, :color_label), get_prop(g_meta, v, :color_label)) for v in vertices(g_meta)]
        vertexfillcolors = [normalize_hue(get_prop(g_meta, v, :color_label)) for v in vertices(g_meta)],
        # edgelabels = [G.weights[e.src, e.dst] for e in edges(G)]
    )
end 600 400 filepath_svg



# modify function in Start_J_Space
# API for DataFrames to include in J-Space, new file julia for DataFrames API
#label_1 = 12
#label_2 = 11
#df_edges = find_edges_from_vertex(filepath_dataframe_edges, label_1, label_2)

#df_neighbor = find_neighbors_from_vertex(filepath_dataframe_edges, filepath_dataframe_label, label_1)

#CI_example = CartesianIndex(95, 1, 1)
#df_vertex = find_vertex_from_CI(filepath_dataframe_label, CI_example)

#color_example = 0.6854407198549426
#df_vertex = find_vertex_from_color(filepath_dataframe_label, color_example)

#df_vertex = find_vertex_from_label(filepath_dataframe_label, label_1)

#degree = count_vertex_degree(filepath_dataframe_edges, label_1)

#edges = count_edges(filepath_dataframe_edges)

#vertices = count_vertices(filepath_dataframe_label)
