export get_random_color
export weighted_graph_to_adjacency_matrix
export save_adjacency_matrix
export apply_segmentation_without_download
export apply_segmentation_with_download


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

function region_adjacency_graph(s::SegmentedImage)
    function neighbor_regions!(n::Set{Int}, visited::AbstractArray, s::SegmentedImage, I::CartesianIndex)
        R = CartesianIndices(axes(s.image_indexmap))
        I1 = _oneunit(CartesianIndex{ndims(visited)})
        Ibegin, Iend = first(R), last(R)
        t = Vector{CartesianIndex{ndims(visited)}}()
        push!(t, I)
        while !isempty(t)
            temp = pop!(t)
            visited[temp] = true
            for J in _colon(max(Ibegin, temp-I1), min(Iend, temp+I1))
                if s.image_indexmap[temp] != s.image_indexmap[J]
                    push!(n,s.image_indexmap[J])
                elseif !visited[J]
                    push!(t,J)
                end
            end
        end
        n
    end

    visited  = fill(false, axes(s.image_indexmap))                           # Array to mark the pixels that are already visited
    G        = SimpleWeightedGraph()                                         # The region_adjacency_graph
    vert_map = Dict{Int,Int}()                                               # Map that stores (label, vertex) pairs
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
            neighbor_regions!(n, visited, s, p)
            for i in n
                Graphs.add_edge!(G, vert_map[s.image_indexmap[p]], vert_map[i])
            end
        end
    end
    G, vert_map
end

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
    weighted_graph_to_adjacency_matrix(G::SimpleWeightedGraph{Int64, Float64}, n::Int64)

Converte un grafo pesato rappresentato come un `SimpleWeightedGraph` in una matrice di adiacenza non pesata.

# Argomenti:
- `G::SimpleWeightedGraph{Int64, Float64}`: Grafo pesato rappresentato come un `SimpleWeightedGraph` con etichette di vertice intere e pesi dei bordi in virgola mobile.
- `n::Int64`: Numero di nodi nella matrice di adiacenza.

# Valore di ritorno
- `adjacency_matrix`: matrice di adiacenza `Matrix{Int64}`.

# Note
La funzione restituisce una matrice di adiacenza di dimensioni `n` x `n` rappresentante il grafo non pesato.
Se i nodi `i` e `j` sono adiacenti, la matrice di adiacenza conterrà un valore di `1` nella posizione `(i, j)` e `(j, i)`. Altrimenti, la matrice di adiacenza conterrà un valore di `0`.
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
    save_adjacency_matrix(matrix::Matrix{Int64}, filepath_matrix::AbstractString)

Salva una matrice di adiacenza rappresentata come una matrice di interi su un file di testo.

# Argomenti:
- `matrix::Matrix{Int64}`: La matrice di interi che rappresenta la matrice di adiacenza.
- `filepath_matrix::AbstractString`: Il percorso file rappresentato come una stringa che indica dove salvare la matrice.

# Note
La funzione apre il file specificato dal percorso `filepath_matrix` in modalità scrittura e scrive la matrice in formato di matrice di adiacenza,
dove ogni riga rappresenta i nodi adiacenti di un nodo. I numeri nella matrice sono separati da spazi.
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

La funzione esegue la segmentazione di un'immagine istologica, genera il rispettivo grafo e lo traduce in una
matrice di adiacenza simmetrica e con soli 0 e 1.

# Argomenti
- `slide_info::Tuple{String, Vector{UInt8}, String}`: Tupla contenente l'id della slide, l'immagine stessa ottenuta dal DB e il percorso del file dell'immagine originale.

# Valori di ritorno
- `filepath_matrix`: Il percorso all'interno del quale è memorizzata la matrice in formato `.txt`.
- `matrix`: La matrice di adiacenza costruita dalla segmentazione.

# Note
La funzione utilizza l'algoritmo di segmentazione watershed per segmentare l'immagine in diversi raggruppamenti di pixel. La segmentazione viene eseguita utilizzando una trasformazione delle caratteristiche dell'immagine (`feature_transform`) e l'etichettatura dei componenti connessi. Viene quindi calcolata la distanza tra le diverse regioni e viene costruito un grafo di adiacenza delle regioni, utilizzando la funzione `region_adjacency_graph`.
Il grafo ottenuto viene trasformato in una matrice di adiacenza con la funzione `weighted_graph_to_adjacency_matrix` che sarà salvata nel percorso della immagine originale.
"""
function apply_segmentation_without_download(slide_info::Tuple{String, Vector{UInt8}, String})
    # initialization
    svs_image = slide_info[2]
    slide_id = slide_info[1]
    # start segmentation process
    println("LOAD SLIDE ... ($slide_id)")
    img = ImageMagick.load_(svs_image)
    bw = Gray.(img) .> 0.21
    img = nothing
    dist = 1 .- distance_transform(feature_transform(bw))
    bw = nothing
    markers = label_components(dist .< -0.00001)

    println("APPLY SEGMENTATION ... ($slide_id)")
    segments = watershed(dist, markers)
    dist = nothing
    markers = nothing
    GC.gc()

    println("BUILD GRAPH ... ($slide_id)")
    # weight_fn(i,j) = euclidean(segment_pixel_count(segments,i), segment_pixel_count(segments,j))
    G, vert_map = region_adjacency_graph(segments)
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

La funzione esegue la segmentazione di un'immagine istologica, salva l'immagine segmentata in formato `.tif`,
genera il rispettivo grafo e lo traduce in una matrice di adiacenza.

# Argomenti
- `slide_info::Tuple{String, Vector{UInt8}, String}`: Tupla contenente l'id della slide, l'immagine stessa ottenuta dal DB e il percorso del file dell'immagine originale.

# Valori di ritorno
- `filepath_seg`: Il percorso all'interno del quale è memorizzata l'immagine segmentata in formato `.tif`.
- `filepath_matrix`: Il percorso all'interno del quale è memorizzata il grafo in formato `.txt`.
- `matrix`: La matrice di adiacenza costruita dalla segmentazione.

# Note
La funzione utilizza l'algoritmo di segmentazione watershed per segmentare l'immagine in diversi raggruppamenti di pixel. La segmentazione viene eseguita utilizzando una trasformazione delle caratteristiche dell'immagine (`feature_transform`) e l'etichettatura dei componenti connessi. Viene quindi calcolata la distanza tra le diverse regioni e viene costruito un grafo di adiacenza delle regioni, utilizzando la funzione `region_adjacency_graph`.
Il grafo ottenuto viene trasformato in una matrice di adiacenza con la funzione `weighted_graph_to_adjacency_matrix` che sarà salvata nel percorso della immagine originale.
Infine, viene salvata un'immagine `.tif` segmentata e viene restituito il percorso del file della slide segmentata.
"""
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
