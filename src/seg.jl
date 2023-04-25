using Images
using ImageView
using ImageSegmentation
using FileIO
using Random
using ImageMagick
using ImageFiltering
using ImageMorphology
# lgg lusc
using LightGraphs
using SimpleWeightedGraphs
using IndirectArrays
using Graphs



function get_random_color(seed)
    Random.seed!(seed)
    rand(RGB{N0f8})
end

function save_adjacency_matrix(matrix::Matrix{Int64})
    # Apri il file di testo
    f = open("adjacency_matrix.txt", "w")
    # Scrivi le dimensioni della matrice (numero di righe e colonne)
    n_rows, n_cols = size(matrix)
    matrix_string = ""
    # Scrivi la matrice riga per riga
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

link = "C:/Users/nicom/Desktop/segmentation/TCGA-OR-A5J1-01A-01-TS1.CFE08710-54B8-45B0-86AE-500D6E36D8A5_001.tif"
img_byte = read(link)
img = ImageMagick.load_(img_byte)
bw = Gray.(img) .> 0.20
dist = 1 .- distance_transform(feature_transform(bw))
markers = label_components(dist .< -0.0001)
segments = watershed(dist, markers)
labels = labels_map(segments)

# first method
# colored_labels = IndirectArray(labels, distinguishable_colors(maximum(labels)))
# masked_labels = colored_labels .* (1 .- bw)

# second method
# segmented_slide = map(i->get_random_color(i), labels_map(segments)) .* (1 .-bw)
# generate_graph(segments)
# ::SegmentedImage{Array{Int64, 3}, Float64}

# generate graph for adjacency_matrix (J-Space)
weight_fn(i,j) = euclidean(segment_pixel_count(segments,i), segment_pixel_count(segments,j))
G, vert_map = region_adjacency_graph(segments, weight_fn)
nvertices = length(vert_map)
matrix = weighted_graph_to_adjacency_matrix(G, nvertices)
save_adjacency_matrix(matrix)

# imshow(segmented_slide)
save("C:/Users/nicom/Desktop/img.tif", masked_labels)
