#using Pkg
#Pkg.add("Gtk")
export build_GUI

global str = ""
function on_button_clicked(w)
  println("The button has been clicked")
  println("selezione: "*str)
end


function save_text(widget, data)
    text = data["text"].text
    # qui si può inserire il codice per salvare il testo
    println("Testo salvato: ", text)
end


function build_GUI(collection_values::Array)
    win = GtkWindow("JHistint GUI - Julia Histopathology Interface", 400, 400)
    box = GtkBox(:v)  # :h makes a horizontal layout, :v a vertical layout

    for item in collection_values
        println("TCGA - ",item)
        label = GtkLabel("TCGA - $item")
        GAccessor.selectable(label,true)
        push!(box, label)
    end

    button = GtkButton("Download Slides")
    label = GtkLabel("Select a collection (i.e. acc, blca, brca etc.):")
    entry = GtkEntry()

    push!(box, label)
    push!(box, entry)
    push!(box, button)
    push!(win, box)    
    showall(win)

    str = get_gtk_property(entry,:text,String)
    #signal_connect(save_text, button , "clicked", data)
    signal_connect(on_button_clicked, button, "clicked")
    return str
end