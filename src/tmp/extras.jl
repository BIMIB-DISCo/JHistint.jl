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
