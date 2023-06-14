export get_random_color
export weighted_graph_to_adjacency_matrix
export weighted_graph_to_adjacency_matrix_weight
export save_adjacency_matrix
export region_adjacency_graph
export apply_segmentation_without_download
export apply_segmentation_with_download
export apply_segmentation_SOPHYSM

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

"""
    get_random_color(seed)

Function to return a random 8-bit RGB format color, using a specified seed.

# Arguments
- `seed`: An integer used to initialize the random number generator. If two calls to the function use the same seed, the same color will be generated.

# Return value
The function returns a random 8-bit RGB format color.
"""
function get_random_color(seed)
    Random.seed!(seed)
    rand(RGB{N0f8})
end

"""
    weighted_graph_to_adjacency_matrix(G::SimpleWeightedGraph{Int64, Float64}, n::Int64)

Converts a weighted graph represented as a `SimpleWeightedGraph` into an unweighted adjacency matrix.

# Arguments:
- `G::SimpleWeightedGraph{Int64, Float64}`: Weighted graph represented as a `SimpleWeightedGraph` with integer vertex labels and floating-point edge weights.
- `n::Int64`: Number of nodes in the adjacency matrix.

# Return value:
- `adjacency_matrix`: `Matrix{Int64}` adjacency matrix.

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

function weighted_graph_to_adjacency_matrix_weight(G::SimpleWeightedGraph{Int64, Float64}, n::Int64)
    adjacency_matrix = zeros(Int, n, n)
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
    apply_segmentation_without_download(slide_info::Tuple{String, Vector{UInt8}, String})

The function performs the segmentation of a histological image, generates its corresponding graph, and translates it into a symmetric adjacency matrix with only 0s and 1s.

# Arguments
- `slide_info::Tuple{String, Vector{UInt8}, String}`: A tuple containing the slide ID, the image obtained from the DB, and the path of the original image file.

# Return values
- `filepath_matrix`: The path where the adjacency matrix is stored in `.txt` format.
- `matrix`: The adjacency matrix constructed from the segmentation.

# Notes
The function uses the watershed segmentation algorithm to segment the image into different groups of pixels. Segmentation is performed using a feature transformation of the image (`feature_transform`) and labeling of connected components. The distance between the different regions is then calculated, and an adjacency graph of the regions is constructed using the `region_adjacency_graph` function. The resulting graph is then transformed into an adjacency matrix using the `weighted_graph_to_adjacency_matrix` function and saved to the path of the original image.
"""
function apply_segmentation_without_download(slide_info::Tuple{String, Vector{UInt8}, String})
    # initialization
    svs_image = slide_info[2]
    slide_id = slide_info[1]
    # start segmentation process
    println("LOAD SLIDE ... ($slide_id)")
    img = ImageMagick.load_(svs_image)
    bw = Gray.(img) .> 0.15
    img = nothing
    dist = 1 .- distance_transform(feature_transform(bw))
    bw = nothing
    markers = label_components(dist .< -0.3)

    println("APPLY SEGMENTATION ... ($slide_id)")
    segments = watershed(dist, markers)
    dist = nothing
    markers = nothing
    GC.gc()

    println("BUILD GRAPH ... ($slide_id)")
    weight_fn(i,j) = euclidean(segment_pixel_count(segments,i), segment_pixel_count(segments,j))
    G, vert_map = region_adjacency_graph(segments, weight_fn)
    nvertices = length(vert_map)
    vert_map = nothing
    segments = nothing
    GC.gc()

    println("BUILD & SAVE ADJACENCY MATRIX ... ($slide_id)")
    matrix = weighted_graph_to_adjacency_matrix(G, nvertices)
    filepath_matrix = replace(slide_info[3], ".tif" => ".txt")
    G = nothing
    nvertices = nothing
    save_adjacency_matrix(matrix, filepath_matrix)

    return filepath_matrix, matrix
end

"""
    apply_segmentation_with_download(slide_info::Tuple{String, Vector{UInt8}, String})

The function performs segmentation of a histological image, saves the segmented image in `.tif` format, generates the corresponding graph, and translates it into an adjacency matrix.

# Arguments
- `slide_info::Tuple{String, Vector{UInt8}, String}`: Tuple containing the slide ID, the image itself obtained from the DB, and the path of the original image file.

# Return values
- `filepath_seg`: The path where the segmented image is stored in `.tif` format.
- `filepath_matrix`: The path where the graph is stored in `.txt` format.
- `matrix`: The adjacency matrix constructed from the segmentation.

# Notes
The function uses the watershed segmentation algorithm to segment the image into different groups of pixels. Segmentation is performed using an image feature transformation (`feature_transform`) and connected component labeling. The distance between different regions is then calculated, and an adjacency graph of the regions is constructed using the `region_adjacency_graph` function.
The obtained graph is transformed into an adjacency matrix using the `weighted_graph_to_adjacency_matrix` function, which is saved in the path of the original image.
Finally, a segmented `.tif` image is saved, and the path of the segmented slide file is returned.
"""
function apply_segmentation_with_download(slide_info::Tuple{String, Vector{UInt8}, String})
    # initialization
    svs_image = slide_info[2]
    slide_id = slide_info[1]
    # start segmentation process
    println("LOAD SLIDE ... ($slide_id)")
    img = ImageMagick.load_(svs_image)
    bw = Gray.(img) .> 0.15
    dist = 1 .- distance_transform(feature_transform(bw))
    markers = label_components(dist .< -0.3)

    println("APPLY SEGMENTATION ... ($slide_id)")
    segments = watershed(dist, markers)

    println("BUILD SEGMENTED SLIDE ... ($slide_id)")
    labels = labels_map(segments)
    colored_labels = IndirectArray(labels, distinguishable_colors(maximum(labels)))
    masked_colored_labels = colored_labels .* (1 .- bw)
    imshow(masked_colored_labels)

    println("BUILD GRAPH ... ($slide_id)")
    weight_fn(i,j) = euclidean(segment_pixel_count(segments,i), segment_pixel_count(segments,j))
    G, vert_map = region_adjacency_graph(segments, weight_fn)
    nvertices = length(vert_map)

    println("BUILD & SAVE ADJACENCY MATRIX ... ($slide_id)")
    matrix = weighted_graph_to_adjacency_matrix(G, nvertices)
    filepath_matrix = replace(slide_info[3], ".tif" => ".txt")
    save_adjacency_matrix(matrix, filepath_matrix)

    println("SAVE SEGMENTED SLIDE ... ($slide_id)")
    filepath_seg = replace(slide_info[3], ".tif" => "_seg.tif")
    save(filepath_seg, masked_colored_labels)
    return filepath_seg, filepath_matrix, matrix
end

function apply_segmentation_SOPHYSM(filepath_input::AbstractString, filepath_output::AbstractString, thresholdGray::Float64, thresholdMarker::Float64)
    # load slide
    svs_image = read(filepath_input)
    img = ImageMagick.load_(svs_image)
    bw = Gray.(img) .> thresholdGray
    dist = 1 .- distance_transform(feature_transform(bw))
    markers = label_components(dist .< thresholdMarker)
    # watershed
    segments = watershed(dist, markers)
    # build segmented slide
    labels = labels_map(segments)
    colored_labels = IndirectArray(labels, distinguishable_colors(maximum(labels)))
    masked_colored_labels = colored_labels .* (1 .- bw)
    # build graph
    weight_fn(i,j) = euclidean(segment_pixel_count(segments,i), segment_pixel_count(segments,j))
    df = DataFrame()
    G, vert_map, df = region_adjacency_graph(segments, weight_fn)
    nvertices = length(vert_map)
    # build and save adjacency matrix
    matrix = weighted_graph_to_adjacency_matrix_weight(G, nvertices)
    filepath_matrix = replace(filepath_output, r"....$" => ".txt")
    save_adjacency_matrix(matrix, filepath_matrix)
    # save segmented slide
    filepath_seg = replace(filepath_output, r"....$" => "_seg.png")
    save(filepath_seg, masked_colored_labels)
end
