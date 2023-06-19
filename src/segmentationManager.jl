export get_random_color
export apply_segmentation_without_download
export apply_segmentation_with_download
export apply_segmentation_SOPHYSM

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
    svs_image = slide_info[2]
    slide_id = slide_info[1]

    println("LOAD SLIDE ... ($slide_id)")
    img = ImageMagick.load_(svs_image)
    bw = Gray.(img) .> 0.15
    dist = 1 .- distance_transform(feature_transform(bw))
    markers = label_components(dist .< -0.3)

    println("APPLY SEGMENTATION ... ($slide_id)")
    segments = watershed(dist, markers)

    println("BUILD GRAPH ... ($slide_id)")
    weight_fn(i,j) = euclidean(segment_pixel_count(segments,i), segment_pixel_count(segments,j))
    df = DataFrame()
    G, vert_map, df = region_adjacency_graph(segments, weight_fn)
    nvertices = length(vert_map)

    println("BUILD & SAVE ADJACENCY MATRIX ... ($slide_id)")
    matrix = weighted_graph_to_adjacency_matrix(G, nvertices)
    filepath_matrix = replace(slide_info[3], ".tif" => ".txt")
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
    svs_image = slide_info[2]
    slide_id = slide_info[1]

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
    df = DataFrame()
    G, vert_map, df = region_adjacency_graph(segments, weight_fn)
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
    df_labels = DataFrame()
    G, vert_map, df_labels = region_adjacency_graph(segments, weight_fn)
    nvertices = length(vert_map)
    # save dataframe label as .CSV
    filepath_dataframe_labels = replace(filepath_output, r"....$" => "_dataframe_labels.csv")
    CSV.write(filepath_dataframe_labels, df_labels)
    # build and save adjacency matrix
    matrix = weighted_graph_to_adjacency_matrix_weight(G, nvertices)
    filepath_matrix = replace(filepath_output, r"....$" => ".txt")
    save_adjacency_matrix(matrix, filepath_matrix)
    # build and save dataframe edgelist as .CSV
    df_edges = build_dataframe_as_edgelist(matrix)
    filepath_dataframe_edges = replace(filepath_output, r"....$" => "_dataframe_edges.csv")
    CSV.write(filepath_dataframe_edges, df_edges)
    # save segmented slide
    filepath_seg = replace(filepath_output, r"....$" => "_seg.png")
    save(filepath_seg, masked_colored_labels)
end
