export get_random_color
export weighted_graph_to_adjacency_matrix
export save_adjacency_matrix
export apply_segmentation_without_download
export apply_segmentation_with_download

"""
    get_random_color(seed)

Funzione per restituire un colore casuale in formato RGB a 8 bit, utilizzando un seme specificato.

# Argomenti
- `seed` = Un intero utilizzato per inizializzare il generatore di numeri casuali. Se due chiamate alla funzione utilizzano lo stesso seme, verrà generato lo stesso colore.

# Valore di ritorno
La funzione restituisce un colore casuale in formato RGB a 8 bit.
"""
function get_random_color(seed)
    Random.seed!(seed)
    rand(RGB{N0f8})
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
    apply_segmentation_without_download(slide_info::Tuple{String, Array{ColorTypes.RGB{FixedPointNumbers.N0f8}, 3}, String})

La funzione esegue la segmentazione di un'immagine istologica.

# Argomenti
- `slide_info::Tuple{String, Array{ColorTypes.RGB{FixedPointNumbers.N0f8}, 3}, String}`: Tupla contenente l'id della slide, l'immagine stessa ottenuto dal DB e il percorso del file dell'immagine originale.

# Valori di ritorno
- `filepath_seg`: Il percorso all'interno del quale è memorizzata la slide segmentata.
- `segmented_slide`: La slide istologica segmentata in formato `.tif`.

# Note
La funzione utilizza l'algoritmo di segmentazione watershed per segmentare l'immagine in diversi raggruppamenti di pixel. La segmentazione viene eseguita utilizzando una trasformazione delle caratteristiche dell'immagine (`feature_transform`) e l'etichettatura dei componenti connessi. Viene quindi calcolata la distanza tra le diverse regioni e viene costruito un grafo di adiacenza delle regioni, utilizzando la funzione `region_adjacency_graph`. Viene inoltre assegnato un colore casuale a ciascuna regione segmentata.
Infine, viene salvata un'immagine `.tif` segmentata e viene restituito il percorso del file della slide segmentata e la slide stessa.
"""
function apply_segmentation_without_download(slide_info::Tuple{String, Vector{UInt8}, String})
    # initialization
    svs_image = slide_info[2]
    slide_id = slide_info[1]
    # start segmentation process
    println("LOAD SLIDE ... ($slide_id)")
    img = ImageMagick.load_(svs_image)
    bw = Gray.(img) .> 0.21
    dist = 1 .- distance_transform(feature_transform(bw))
    markers = label_components(dist .< -0.00001)

    println("APPLY SEGMENTATION ... ($slide_id)")
    segments = watershed(dist, markers)

    println("BUILD GRAPH ... ($slide_id)")
    weight_fn(i,j) = euclidean(segment_pixel_count(segments,i), segment_pixel_count(segments,j))
    G, vert_map = region_adjacency_graph(segments, weight_fn)
    nvertices = length(vert_map)

    println("BUILD & SAVE ADJACENCY MATRIX ... ($slide_id)")
    matrix = weighted_graph_to_adjacency_matrix(G, nvertices)
    filepath_matrix = replace(slide_info[3], ".tif" => ".txt")
    save_adjacency_matrix(matrix, filepath_matrix)
    
    return filepath_matrix, matrix
end

function apply_segmentation_with_download(slide_info::Tuple{String, Vector{UInt8}, String})
    # initialization
    svs_image = slide_info[2]
    slide_id = slide_info[1]
    # start segmentation process
    println("LOAD SLIDE ... ($slide_id)")
    img = ImageMagick.load_(svs_image)
    bw = Gray.(img) .> 0.21
    dist = 1 .- distance_transform(feature_transform(bw))
    markers = label_components(dist .< -0.00001)

    println("APPLY SEGMENTATION ... ($slide_id)")
    segments = watershed(dist, markers)

    println("BUILD SEGMENTED SLIDE ... ($slide_id)")
    labels = labels_map(segments)
    colored_labels = IndirectArray(labels, distinguishable_colors(maximum(labels)))
    masked_colored_labels = colored_labels .* (1 .- bw)

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
