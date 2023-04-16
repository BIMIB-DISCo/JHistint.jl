export get_random_color
export apply_segmentation

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

"""
    apply_segmentation(slide_info::Tuple{String, Array{ColorTypes.RGB{FixedPointNumbers.N0f8}, 3}, String})

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
function apply_segmentation(slide_info::Tuple{String, Vector{UInt8}, String})
    svs_image = slide_info[2]
    println("1")
    img = ImageMagick.load_(svs_image)
    # imshow(img)
    println("1")
    bw = Gray.(img) .> 0.21
    dist = 1 .- distance_transform(feature_transform(bw))
    markers = label_components(dist .< -0.3)
    segments = watershed(dist, markers)
    println("1")
    # segmented_slide = map(i->get_random_color(i), labels_map(segments)) .* (1 .-bw)
    labels = labels_map(segments)
    # colored_labels = IndirectArray(labels, distinguishable_colors(maximum(labels)))
    # masked_colored_labels = colored_labels .* (1 .- bw)
    # weight_fn(i,j) = euclidean(segment_pixel_count(segments,i), segment_pixel_count(segments,j))
    # G, vert_map = region_adjacency_graph(segments, weight_fn)
    filepath_seg = replace(slide_info[3], ".tif" => "_seg.tif")
    # filepath_seg = replace(slide_info[3], ".svs" => "_seg.tif")
    # println(G)
    # println(vert_map)
    println("1")
    save(filepath_seg, labels)
    return filepath_seg
end
