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


"""
    c = segment_pixel_count(seg, l)

Returns the count of pixels that are assigned label `l`. If no label is
supplied, it returns a Dict(label=>pixel_count) of all the labels.
"""
segment_pixel_count(seg::SegmentedImage, l::Int16) = seg.segment_pixel_count[l]
segment_pixel_count(seg::SegmentedImage) = seg.segment_pixel_count


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
    # Build object for dataframe
    df = DataFrame()
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
    df.label = s.segment_labels
    df.position_label = df_cartesian_indices
    df.color_label = df_color_indices
    G, vert_map, df
end

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

function build_dataframe_as_edgelist(mat::Matrix{Int64})
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

# J-SPACE FUNCTIONS
# in DataFrameGraphBridge.jl
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
    show(vertex_names)
    show(attr_names)
    for r in eachrow(vertex_names)
        MetaGraphs.set_props!(mg, r[:vertex_id], Dict([a => r[a] for a in attr_names]))
    end
    # Set edge attributes
    if typeof(edge_attributes) == Symbol
        edge_attributes = Vector{Symbol}([edge_attributes])
    end
    show(edge_attributes)
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


# change signatures (NOT IMPORT)
function spatial_graph(path_dataframe::String, path_matrix::String, seed::Int, n_cell::Int)
    mat = readdlm(path_matrix, Float32)
    df = CSV.File(path_dataframe)
    n = size(mat, 1) # Numero di nodi
    # n = 2 * count(x -> x != -1, mat)
    graph = SimpleWeightedGraph(n)
    for i in 1:n
        for j in 1:i
            if mat[i,j] != -1
                Graphs.add_edge!(graph, i, j, mat[i,j])
            end
        end
    end
    abstract_g = Graphs.Graph(graph)
    show(abstract_g)
    G_meta = metagraph_from_dataframe(SimpleWeightedGraph, df, "origin", "destination", "weight")
    return G_meta
end

# in J-Space.jl
function spatial_graph(path_dataframe_edges::String, path_dataframe_labels::String)
    df_edges_file = CSV.File(path_dataframe_edges)
    df_edges = DataFrame(df_edges_file)

    df_labels_file = CSV.File(path_dataframe_labels)
    df_labels = DataFrame(df_labels_file)
    G_meta = metagraph_from_dataframe_JHistint(MetaGraph, df_edges, :origin, :destination, :weight, :weight, df_labels, :label)
    return G_meta
end

function plot_lattice_JHistint(G::MetaGraph)
    mylayout = Spectral()
    f, ax, p = graphplot(G,
                         layout = mylayout,
                         node_size = 0.0,
                         edge_width=1.0)
    hidedecorations!(ax)
    hidespines!(ax)
    return f, ax, p
end

# in DataFrameAPI.jl
function find_edges_from_vertex(path_dataframe_edges::String, first_vertex::Int64, second_vertex::Int64)
    df_edges_file = CSV.File(path_dataframe_edges)
    df_edges = DataFrame(df_edges_file)
    # Apply filter
    filtered_edges = filter(row -> (row.origin == first_vertex && row.destination == second_vertex) ||
                                   (row.origin == second_vertex && row.destination == first_vertex), df_edges)
    new_df_edges = DataFrame(filtered_edges)
    return new_df_edges
end

function find_neighbors_from_vertex(path_dataframe_edges::String, path_dataframe_labels::String, vertex::Int64)
    df_edges_file = CSV.File(path_dataframe_edges)
    df_edges = DataFrame(df_edges_file)
    df_labels_file = CSV.File(path_dataframe_labels)
    df_labels = DataFrame(df_labels_file)
    # Extract neighbors vertex based on df_edges
    neighbors = unique(vcat(filter(row -> row.origin == vertex, df_edges).destination,
                            filter(row -> row.destination == vertex, df_edges).origin))
    # Apply filter row in df_labels based on label field
    neighbor_labels = filter(row -> row.label in neighbors, df_labels)
    return neighbor_labels
end

function find_vertex_from_CI(path_dataframe_labels::String, CI::CartesianIndex)
    df_labels_file = CSV.File(path_dataframe_labels)
    df_labels = DataFrame(df_labels_file)
    position = string(CI)
    filtered_labels = filter(row -> row.position_label == position, df_labels)
    df_label_result = DataFrame(filtered_labels)
    return df_label_result
end

function find_vertex_from_color(path_dataframe_labels::String, color::Float64)
    df_labels_file = CSV.File(path_dataframe_labels)
    df_labels = DataFrame(df_labels_file)
    filtered_labels = filter(row -> row.color_label == color, df_labels)
    df_label_result = DataFrame(filtered_labels)
    return df_label_result
end

function find_vertex_from_label(path_dataframe_labels::String, label::Int64)
    df_labels_file = CSV.File(path_dataframe_labels)
    df_labels = DataFrame(df_labels_file)
    filtered_labels = filter(row -> row.label == label, df_labels)
    df_label_result = DataFrame(filtered_labels)
    return df_label_result
end

function count_vertex_degree(path_dataframe_labels::String, vertex::Int64)
    df_labels_file = CSV.File(path_dataframe_labels)
    df_labels = DataFrame(df_labels_file)
    outgoing_degree = nrow(filter(row -> row.origin == vertex, df_labels))
    incoming_degree = nrow(filter(row -> row.destination == vertex, df_labels))
    return outgoing_degree + incoming_degree
end

function count_edges(path_dataframe_edges::String)
    df_edges_file = CSV.File(path_dataframe_edges)
    df_edges = DataFrame(df_edges_file)
    num_records = nrow(df_edges)
    return num_records
end

function count_vertices(path_dataframe_labels::String)
    df_labels_file = CSV.File(path_dataframe_labels)
    df_labels = DataFrame(df_labels_file)
    num_records = nrow(df_labels)
    return num_records
end

link = joinpath(@__DIR__, "..", "input_example_demo", "slideExample1", "SlideExample_mini_1.tif")
svs_image = read(link)
println("LOAD SLIDE ... ")
img = ImageMagick.load_(svs_image)

println("LOAD SLIDE OK ... ")
bw = Gray.(img) .> 0.15
dist = 1 .- distance_transform(feature_transform(bw))
markers = label_components(dist .< -0.3)

println("APPLY SEGMENTATION ... ")
segments = watershed(dist, markers)
labels = labels_map(segments)
colored_labels = IndirectArray(labels, distinguishable_colors(maximum(labels)))
masked_colored_labels = colored_labels .* (1 .- bw)

# second method
# segmented_slide = map(i->get_random_color(i), labels_map(segments)) .* (1 .-bw)

println("BUILD GRAPH ... ")
# Calcola distanza euclidea tra il numero di pixel nei segmenti i e j
weight_fn(i,j) = euclidean(segment_pixel_count(segments,i), segment_pixel_count(segments,j))
df_label = DataFrame()
G, vert_map, df_label = region_adjacency_graph(segments, weight_fn)
filepath_dataframe_label = replace(link, ".tif" => "_dataframe_label.csv")
CSV.write(filepath_dataframe_label, df_label)

println("BUILD & SAVE ADJACENCY MATRIX ... ")
nvertices = Graphs.nv(G)
matrix = weighted_graph_to_adjacency_matrix_weight(G, nvertices)
df_edges = build_dataframe_as_edgelist(matrix)
filepath_dataframe_edges = replace(link, ".tif" => "_dataframe_edges.csv")
CSV.write(filepath_dataframe_edges, df_edges)

filepath_matrix = replace(link, ".tif" => "_matrix.txt")
save_adjacency_matrix(matrix, filepath_matrix)

# J-Space
g_meta = spatial_graph(filepath_dataframe_edges, filepath_dataframe_label)
plot_lattice_JHistint(g_meta)

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
