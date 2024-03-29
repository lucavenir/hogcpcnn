# Executable script meant to process raw inputs from the datasets
# The goal is to write standard input files, easy to interpret.
# The reason behind this pre-processing phase is that in this project
# multiple datasets are involved, maybe each with different "notations".
# Furthermore, we want to address problem such as non-occurrent genes in some files.

# We chose the following format for our files:
#       - G = (V,E) will be encoded with just a list of edges, like: "i j"
#       - Vertex labeling will have the "number name" format
#       - Looks like samples already have the format we want: "name heat"

import os

# The OUTPUT files of this script should be the following:
input_folder = '../in/'
adjacency_matrix = 'adjacency_matrix.txt'
vertex_labels = 'vertex_labels.txt'
samples = 'heat.txt'
rawinputs_folder = '../rawinput/'

def hinthi2012():
    # Processing HINT+HI2012 database.

    # Referencing to the global variables:
    global input_folder
    global adjacency_matrix
    global vertex_labels
    global rawinputs_folder
    global samples

    # File locations
    dataset_name = 'HINT+HI2012'
    edges_file = 'hint+hi2012_edge_file.txt'
    index_file = 'hint+hi2012_index_file.txt'
    heat_file = 'mutsigcv_expr_filtered.txt'
    input_folder += dataset_name+'/'

    try:
        os.mkdir(input_folder)
    except FileExistsError:
        pass

    filter = set()
    # With the edge list, we should be able to construct the graph later.
    vertices_dict = {}  # We want to re-map the vertices
    with open(rawinputs_folder+dataset_name+'/'+edges_file, 'r') as infp:
        lines = infp.readlines()
        outfp = open(input_folder+adjacency_matrix, 'w')

        i = 0
        for l in lines:
            s = l.split(" ")
            j = int(s[0])-1
            k = int(s[1])-1
            if j not in vertices_dict:
                vertices_dict[j] = i
                i += 1
            if k not in vertices_dict:
                vertices_dict[k] = i
                i += 1
            outfp.write(str(vertices_dict[j])+' '+str(vertices_dict[k])+'\n')
        outfp.close()

    # After this analysis, we store the actual number of nodes found
    # (here, n="nodes that are in a connected component")
    n = i

    # Now, processing the vertices' labels
    counter = 0
    vertices_names_dict = {}
    with open(rawinputs_folder+dataset_name+'/'+index_file, 'r') as infp:
        lines = infp.readlines()

        for l in lines:
            s = l.split(" ")
            vertex_number = int(s[0])-1
            vertex_name = str(s[1].split("\t")[0])  # Removing the '\t0\n'
            if vertex_number in vertices_dict:
                vertices_names_dict[vertex_name] = vertices_dict[vertex_number]
            else:
                counter += 1

    if counter > 0:
        print("Nel file index compaiono "+str(counter)+" nodi che non ci sono nel grafo")

    # Writing down the label files, ordering them by node.
    with open(input_folder+vertex_labels, 'w') as outfp:
        outfp.write(str(n)+'\n')
        outfp.write(str(n)+'\n')  # Writing it down twice IS necessary.
        for item in sorted(vertices_names_dict.items(), key=lambda kv:kv[1]):
            outfp.write(str(item[1])+' '+item[0]+'\n')

    # Finally, processing the heat file
    counter = 0  # DEBUG: perché ci sono 290 proteine escluse?
    heat_dict = {}
    with open(rawinputs_folder+heat_file, 'r') as infp:
        lines = infp.readlines()

        for l in lines:
            s = l.split(" ")
            vertex_name = str(s[0])
            heat_value = float(s[1])

            if vertex_name in vertices_names_dict:
                heat_dict[vertices_names_dict[vertex_name]] = heat_value
            else:
                counter += 1 if heat_value>0 else 0
        outfp.close()

    # TODO: da chiedere al professore ---> cosa fare di queste proteine mutate
    # in qualche sample, ma non presenti nel grafo di interazioni prot-prot?
    if counter > 0:
        print("WARNING: Nel file heat compaiono " +
        str(counter) +
        " proteine (con heat>0) i cui nodi sono isolati o non presenti nel grafo")

    # Writing down the label files, ordering them by node.
    with open(input_folder+samples, 'w') as outfp:
        sorted_vertices = sorted(heat_dict.items(), key=lambda kv:kv[0])
        for i in range(n):
            if i in heat_dict:
                outfp.write(str(i)+' '+str(heat_dict[i])+'\n')
            else:
                outfp.write(str(i)+' '+str(0)+'\n')

    # # TODO: implementare preventivamente alcuni filtri con la seguente idea.
    # in qualche modo una lista di nodi va esclusa dai giochi.
    # l'idea è di mettere tali nodi in un set ed evitare semplicemente di inserirli
    # dentro i heat.txt, adjacency_matrix.txt e vertex_labels.txt.
    # L'implementazione di quest'idea quindi andrebbe fatta SOPRA
    # filtering = set() è un modo per creare un set vuoto di nodi da escludere

def irefindex():
    # Processing iRedIndex database.

    # Referencing to the global variables:
    global input_folder
    global adjacency_matrix
    global vertex_labels
    global samples
    global rawinputs_folder

    # File locations
    dataset_name = 'iRefIndex'

    edges_file = 'irefindex_edge_file.txt'
    index_file = 'irefindex_index_file.txt'
    heat_file = 'mutsigcv_expr_filtered.txt'
    input_folder += dataset_name+'/'

    try:
        os.mkdir(input_folder)
    except FileExistsError:
        pass

    filter = set()  # TODO: later, we'll filter something out. Maybe.
    # With the edge list, we should be able to construct the graph later.
    vertices_dict = {}  # We want to re-map the vertices
    with open(rawinputs_folder+dataset_name+'/'+edges_file, 'r') as infp:
        lines = infp.readlines()
        outfp = open(input_folder+adjacency_matrix, 'w')

        i = 0
        for l in lines:
            s = l.split(" ")
            j = int(s[0])-1
            k = int(s[1])-1
            if j not in vertices_dict:
                vertices_dict[j] = i
                i += 1
            if k not in vertices_dict:
                vertices_dict[k] = i
                i += 1
            outfp.write(str(vertices_dict[j])+' '+str(vertices_dict[k])+'\n')
        outfp.close()

    # After this analysis, we store the actual number of nodes found
    # (here, n="nodes that are in a connected component")
    n = i

    # Now, processing the vertices' labels
    counter = 0
    vertices_names_dict = {}
    with open(rawinputs_folder+dataset_name+'/'+index_file, 'r') as infp:
        lines = infp.readlines()

        for l in lines:
            s = l.split(" ")
            vertex_number = int(s[0])-1
            vertex_name = str(s[1])[:-1]  # Removing the '\n'.
            if vertex_number in vertices_dict:
                vertices_names_dict[vertex_name] = vertices_dict[vertex_number]
            else:
                counter += 1

    if counter > 0:
        print("Nel file index compaiono "+str(counter)+" nodi che non ci sono nel grafo")

    # Writing down the label files, ordering them by node.
    with open(input_folder+vertex_labels, 'w') as outfp:
        outfp.write(str(n)+'\n')
        for item in sorted(vertices_names_dict.items(), key=lambda kv:kv[1]):
            outfp.write(str(item[1])+' '+item[0]+'\n')

    # Finally, processing the heat file
    counter = 0  # DEBUG: perché ci sono X proteine escluse?
    heat_dict = {}
    with open(rawinputs_folder+heat_file, 'r') as infp:
        lines = infp.readlines()

        for l in lines:
            s = l.split(" ")
            vertex_name = str(s[0])
            heat_value = float(s[1])

            if vertex_name in vertices_names_dict:
                heat_dict[vertices_names_dict[vertex_name]] = heat_value
            else:
                counter += 1 if heat_value>0 else 0
        outfp.close()

    # TODO: da chiedere al professore ---> cosa fare di queste proteine mutate
    # in qualche sample, ma non presenti nel grafo di interazioni prot-prot?
    if counter > 0:
        print("WARNING: Nel file heat compaiono " +
        str(counter) +
        " proteine (con heat>0) i cui nodi sono isolati o non presenti nel grafo")

    # Writing down the label files, ordering them by node.
    with open(input_folder+samples, 'w') as outfp:
        sorted_vertices = sorted(heat_dict.items(), key=lambda kv:kv[0])
        for i in range(n):
            if i in heat_dict:
                outfp.write(str(i)+' '+str(heat_dict[i])+'\n')
            else:
                outfp.write(str(i)+' '+str(0)+'\n')


def multinet():
    # Processing Multinet database.

    # Referencing to the global variables:
    global input_folder
    global adjacency_matrix
    global vertex_labels
    global samples
    global rawinputs_folder

    # File locations
    dataset_name = 'Multinet'

    edges_file = 'multinet_edge_file.txt'
    index_file = 'multinet_index_file.txt'
    heat_file = 'mutsigcv_expr_filtered.txt'
    input_folder += dataset_name+'/'

    try:
        os.mkdir(input_folder)
    except FileExistsError:
        pass

    filter = set()  # TODO: later, we'll filter something out. Maybe.
    # With the edge list, we should be able to construct the graph later.
    vertices_dict = {}  # We want to re-map the vertices
    with open(rawinputs_folder+dataset_name+'/'+edges_file, 'r') as infp:
        lines = infp.readlines()
        outfp = open(input_folder+adjacency_matrix, 'w')

        i = 0
        for l in lines:
            s = l.split(" ")
            j = int(s[0])-1
            k = int(s[1])-1
            if j not in vertices_dict:
                vertices_dict[j] = i
                i += 1
            if k not in vertices_dict:
                vertices_dict[k] = i
                i += 1
            outfp.write(str(vertices_dict[j])+' '+str(vertices_dict[k])+'\n')
        outfp.close()

    # After this analysis, we store the actual number of nodes found
    # (here, n="nodes that are in a connected component")
    n = i

    # Now, processing the vertices' labels
    counter = 0
    vertices_names_dict = {}
    with open(rawinputs_folder+dataset_name+'/'+index_file, 'r') as infp:
        lines = infp.readlines()

        for l in lines:
            s = l.split(" ")
            vertex_number = int(s[0])-1
            vertex_name = str(s[1])[:-1]  # Removing the '\n'.
            if vertex_number in vertices_dict:
                vertices_names_dict[vertex_name] = vertices_dict[vertex_number]
            else:
                counter += 1

    if counter > 0:
        print("Nel file index compaiono "+str(counter)+" nodi che non ci sono nel grafo")

    # Writing down the label files, ordering them by node.
    with open(input_folder+vertex_labels, 'w') as outfp:
        outfp.write(str(n)+'\n')
        for item in sorted(vertices_names_dict.items(), key=lambda kv:kv[1]):
            outfp.write(str(item[1])+' '+item[0]+'\n')

    # Finally, processing the heat file
    counter = 0  # DEBUG: perché ci sono X proteine escluse?
    heat_dict = {}
    with open(rawinputs_folder+heat_file, 'r') as infp:
        lines = infp.readlines()

        for l in lines:
            s = l.split(" ")
            vertex_name = str(s[0])
            heat_value = float(s[1])

            if vertex_name in vertices_names_dict:
                heat_dict[vertices_names_dict[vertex_name]] = heat_value
            else:
                counter += 1 if heat_value>0 else 0
        outfp.close()

    # TODO: da chiedere al professore ---> cosa fare di queste proteine mutate
    # in qualche sample, ma non presenti nel grafo di interazioni prot-prot?
    if counter > 0:
        print("WARNING: Nel file heat compaiono " +
        str(counter) +
        " proteine (con heat>0) i cui nodi sono isolati o non presenti nel grafo")

    # Writing down the label files, ordering them by node.
    with open(input_folder+samples, 'w') as outfp:
        sorted_vertices = sorted(heat_dict.items(), key=lambda kv:kv[0])
        for i in range(n):
            if i in heat_dict:
                outfp.write(str(i)+' '+str(heat_dict[i])+'\n')
            else:
                outfp.write(str(i)+' '+str(0)+'\n')

if __name__ == "__main__":
    # Here, uncomment the function corresponding to the database you want to decode.
    #hinthi2012()
    #irefindex()
    multinet()
