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

using Distances


using LinearAlgebra
using Base



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

    function neighbor_regions!(G::SimpleWeightedGraph, visited::AbstractArray, s::SegmentedImage, I::CartesianIndex)
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
        added_indices = []
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
                    if s.image_indexmap[J] ∉ added_indices
                        # If the values are different, it means they have two different colorings for the two points,
                        # therefore a neighbor has been identified, which is pushed into n.
                        # push!(n,s.image_indexmap[J])
                        Graphs.add_edge!(G, vert_map[s.image_indexmap[I]], vert_map[s.image_indexmap[J]], weight_fn(s.image_indexmap[I], s.image_indexmap[J]))
                        push!(added_indices, s.image_indexmap[J])
                    end
                elseif !visited[J]
                    # If they are equal, I place them in t, so that,
                    # as long as t is not empty, I can explore all the neighbors
                    # that have the same color.
                    push!(t,J)
                end
                GC.gc()
            end
        end
        G
    end
    # Start
    visited  = fill(false, axes(s.image_indexmap))  # Array to mark the pixels that are already visited
    G        = SimpleWeightedGraph()                # The region_adjacency_graph
    vert_map = Dict{Int,Int}()                      # Map that stores (label, vertex) pairs

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
            # initialize n, fondamental to define the neighbor of p
            # n = Set{Int}()
            # Call neighbor_regions where :
            # n = Set{Int} - visited = Array - s = segmented image - p = CartesianIndex which define neighbors
            try
                G = neighbor_regions!(G, visited, s, p)
            catch oom
                if isa(oom, OutOfMemoryError)
                    # n = Set{Int}()
                    GC.gc()
                    println(">>> OOM")
                    exit()
                end
            # for i in n
            #     Graphs.add_edge!(G, vert_map[s.image_indexmap[p]], vert_map[i], weight_fn(s.image_indexmap[p], i))
            # end
            end
        end
    end
    G, vert_map
end

# DIRECT METHOD WITH CONNECTED COMPONENTS RETURNED BY WATERSHED ALGORITHM
function region_adjacency_graph_label(s::SegmentedImage, weight_fn::Function)

    function neighbor_regions!(n::Set{Int}, visited::AbstractArray, s::SegmentedImage, I::Int64)
        # n = Set{Int} - visited = Array - s = segmented image - p = CartesianIndex which define neighbors
        # R contains each possible index in s
        # R = CartesianIndices(axes(s.image_indexmap))
        # I1 contains a Vector of only 1 with dimension equal to visited
        # I1 = _oneunit(CartesianIndex{ndims(visited)})
        # Ibegin and Iend contains the first and last index of R
        # Ibegin, Iend = first(R), last(R)
        # t is only a empty Vector with dimension equal to visited
        t = Vector{SegmentedImage.segment_labels}()
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
            for J in s.segment_labels
                if s.segment_means[temp] != s.segment_means[J]
                    # Se sono diversi i valori allora vuol dire che hanno
                    # due colorazioni diverse per i due punti, di conseguenza
                    # si è individuato un neighbor, che viene pushato in n
                    push!(n,s.segments_labels)
                elseif !visited[J]
                    # Se sono uguali invece metto in t, in modo tale che,
                    # finché t non è vuoto posso esplorare tutti i vicini
                    # che hanno colore uguale
                    push!(t,J)
                end
            end
        end
        n
    end

    # Start
    G        = SimpleWeightedGraph() # The region_adjacency_graph
    # add vertices to graph
    Graphs.add_vertices!(G,maximum(s.segment_labels))

    for p in s.segment_labels
        show(p)
    end
    G, vert_map
end

# RECURSIVE METHOD (CAUSE STACK OVERFLOW)
function region_adjacency_graph_recursive(s::SegmentedImage, weight_fn::Function)

    function calculate_neighbor_recursive(n::Set{Int}, s::SegmentedImage, G::SimpleWeightedGraph, I1::CartesianIndex, Ibegin::CartesianIndex, Iend::CartesianIndex, temp::CartesianIndex, I::CartesianIndex, visited::AbstractArray)
        visited[temp] = true
        for J in _colon(max(Ibegin, temp-I1), min(Iend, temp+I1))
            if s.image_indexmap[temp] != s.image_indexmap[J]
                push!(n,s.image_indexmap[J])
                # Graphs.add_edge!(G, vert_map[s.image_indexmap[I]], vert_map[J], weight_fn(s.image_indexmap[I], J))
            elseif !visited[J]
                calculate_neighbor_recursive(n, s, G, I1, Ibegin, Iend, J, I, visited)
            end
        end
        G, n
    end

    function neighbor_regions!(n::Set{Int}, visited::AbstractArray, s::SegmentedImage, I::CartesianIndex, G::SimpleWeightedGraph, vert_map::Dict{Int,Int})
        R = CartesianIndices(axes(s.image_indexmap))
        I1 = _oneunit(CartesianIndex{ndims(visited)})
        Ibegin, Iend = first(R), last(R)

        # set index temp to true
        visited[I] = true
        for J in _colon(max(Ibegin, I-I1), min(Iend, I+I1))
            if s.image_indexmap[I] != s.image_indexmap[J]
                push!(n,s.image_indexmap[J])
                # Graphs.add_edge!(G, vert_map[s.image_indexmap[I]], vert_map[J], weight_fn(s.image_indexmap[I], J))
            elseif !visited[J]
                calculate_neighbor_recursive(n, s, G, I1, Ibegin, Iend, J, I, visited)
            end
        end
        G, n
    end

    # Start
    visited  = fill(false, axes(s.image_indexmap))  # Array to mark the pixels that are already visited
    G        = SimpleWeightedGraph()                # The region_adjacency_graph
    vert_map = Dict{Int,Int}()                      # Map that stores (label, vertex) pairs

    # add vertices to graph
    Graphs.add_vertices!(G,length(s.segment_labels))

    # setup `vert_map`
    for (i,l) in enumerate(s.segment_labels)
        vert_map[l] = i
    end

    # add edges to graph
    for p in CartesianIndices(axes(s.image_indexmap))
        if !visited[p]
            n = Set{Int}()
            G, n = neighbor_regions!(n, visited, s, p, G, vert_map)
            for i in n
                Graphs.add_edge!(G, vert_map[s.image_indexmap[p]], vert_map[i], weight_fn(s.image_indexmap[p], i))
            end
        end
    end
    G, vert_map
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

link = joinpath(@__DIR__, "..", "output_example", "example.tif")
svs_image = read(link)
println("LOAD SLIDE ... ")
img = ImageMagick.load_(svs_image)
println("LOAD SLIDE OK ... ")
bw = Gray.(img) .> 0.21
dist = 1 .- distance_transform(feature_transform(bw))
markers = label_components(dist .< -0.3)
println("APPLY SEGMENTATION ... ")
segments = watershed(dist, markers)
# ncolors = maximum(markers)
# randcolors = rand(RGB{N0f8}, ncolors)
# imgi = IndirectArray(markers, randcolors)
# imshow(imgi)


labels = labels_map(segments)
colored_labels = IndirectArray(labels, distinguishable_colors(maximum(labels)))
masked_colored_labels = colored_labels .* (1 .- bw)
# imshow(masked_colored_labels)
# second method
# segmented_slide = map(i->get_random_color(i), labels_map(segments)) .* (1 .-bw)
# generate_graph(segments)
# ::SegmentedImage{Array{Int64, 3}, Float64}

println("BUILD GRAPH ... ")
# generate graph for adjacency_matrix (J-Space)
# Calcola distanza euclidea tra il numero di pixel nei segmenti i e j
weight_fn(i,j) = euclidean(segment_pixel_count(segments,i), segment_pixel_count(segments,j))
G, vert_map = region_adjacency_graph(segments, weight_fn)

nvertices = Graphs.nv(G)
println("BUILD & SAVE ADJACENCY MATRIX ... ")
matrix = weighted_graph_to_adjacency_matrix(G, nvertices)
filepath_matrix = replace(link, ".tif" => "1.txt")
save_adjacency_matrix(matrix, filepath_matrix)
