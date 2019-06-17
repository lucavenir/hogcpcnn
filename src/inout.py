import argparse
import os
import time
import datetime
# import networkx as nx
import numpy as np


def parse_command_line():
    """
    py:function:: parse_command_line()

    Command line parsing:
    algorithms here should be queried via command line, only.
    """
    # initializing 'parameters' with hard-coded defaults.
    current_time = time.time()
    timestamp = datetime.datetime.fromtimestamp(current_time).strftime('%Y%m%d_%H:%M:%S')
    default_folder = "run__" + timestamp

    # Initializing the dictionary containing the parameters with default values
    parameters = {
        'dataset': 'HINT+HI2012',
        'motif': 'u_triangle',
        'w': False,
        'timestamp': timestamp  # for now, this will be ignored
    }

    parser = argparse.ArgumentParser(description='Parser to receive user preferences and/or parameters.')

    # Input folder
    parser.add_argument('-d', type=str, help='Name of the dataset used')
    # Motif searched
    parser.add_argument('-m', type=str, help='Motif type')
    # Computing W only?
    parser.add_argument('-w', action="store_true", help='Type -w to compute the W matrix only')

    args = parser.parse_args()

    if args.d:
        parameters['dataset'] = args.d

    if args.m:
        parameters['motif'] = args.m

    if args.w:
        parameters['w'] = True

    return parameters

def load_input(dataset_name='HINT+HI2012', motif_name=None):
    """
    py:function:: load_input(dataset_name='HINT+HI2012')

    Parsing the input from the /in/ folder, assuming that the input has a
    'standard' format and that it has been already pre-processed with the script.

    :param dataset_name: String containing the name of the dataset to be loaded.

    :return: a dictionary, containing all the input we need.
    """

    inputs = {}
    # For compatibility issues, we're extracting the absolute path of the project.
    abs_path = os.path.dirname(os.path.abspath(__file__))
    project_path = abs_path[:-3]  # This will work just bc of the name of this dir
    input_path = project_path+'in/'+dataset_name+'/'
    if motif_name != None:
        input_path += 'temp/'+motif_name+'/'

    # First, we extract the vertex labeling, since it's easy and it gives us n;
    # (n = |V|)
    file_name = 'vertex_labels.txt'
    with open(input_path+file_name) as infp:
        lines = infp.readlines()
        k = int(lines[0])
        v_labeling = np.empty(k, dtype=object)
        for l in lines[1:]:
            s = l.split(" ")
            idx = int(s[0])
            # To extract the label is a bit tricky, but nothing to it honestly.
            name = str(s[1].split("\n")[0])  # We want to get rid of the '\n'
            v_labeling[idx] = name

    if motif_name == None:
        file_name = 'adjacency_matrix.txt'
        with open(input_path+file_name) as infp:
            lines = infp.readlines()
            adjacency_matrix = np.zeros((k,k), dtype=np.uint8)
            for l in lines:
                s = l.split(" ")
                i = int(s[0])
                j = int(s[1])
                adjacency_matrix[i,j] = 1
                adjacency_matrix[j,i] = 1  # assuming undirected graphs
    else:
        file_name = 'w.txt'
        with open(input_path+file_name) as infp:
            lines = infp.readlines()
            w = np.zeros((k,k), dtype=np.uint16)
            for l in lines:
                s = l.split(",")
                edge = s[0]
                weight = int(s[1])
                s = edge.split(" ")
                i = int(s[0])
                j = int(s[1])
                w[i,j] = weight
                w[j,i] = weight  # assuming undirected graphs

    file_name = 'heat.txt'
    with open(input_path+file_name) as infp:
        lines = infp.readlines()
        heat = np.zeros(k, dtype=np.float64)
        for l in lines:
            s = l.split(" ")
            v = int(s[0])
            heat_value = float(s[1])
            heat[v] = heat_value


    """
    # IDEA: questo codice metteva in conto la possibilità di caricare in memoria
    più matrici con diversi motif, in futuro. Non credo sia utile ma intanto lascio qui.

    temp_path = input_path+'temp/'
    if os.path.exists(temp_path):
        files = [
            f
            for f in os.listdir(temp_path)
            if os.path.isfile(os.path.join(temp_path, f))
        ]
        w = [
            np.zeros((n,n), dtype=np.uint16)
            for i in range(len(files))
        ]

        for i, file_name in enumerate(files):
            temp_path = input_path+'temp/'
            with open(temp_path+file_name) as infp:
                lines = infp.readlines()
                for l in lines:
                    s = l.split(",")
                    weight = int(s[1])
                    s = s[0].split(" ")
                    v = int(s[0])
                    u = int(s[1])
                    w[i][v,u] = weight
                    w[i][u,v] = weight

            motif_name = str(file_name[:-4])
            inputs[motif_name] = w[i]
    """
    inputs['v_labels'] = v_labeling
    inputs['heat'] = heat
    if motif_name != None:
        inputs['w'] = w
    else:
        inputs['adj'] = adjacency_matrix

    return inputs

def write_to_temp(inputs, w, dataset_name='HINT+HI2012', motif_name='u_triangle'):
    # For compatibility issues, we're extracting the absolute path of the project.
    abs_path = os.path.dirname(os.path.abspath(__file__))
    project_path = abs_path[:-3]  # This will work just bc of the name of this dir

    temp_path = project_path+'/in/'+dataset_name+'/temp/'
    try:
        os.mkdir(temp_path)
    except FileExistsError:
        pass

    temp_path += motif_name+'/'
    try:
        os.mkdir(temp_path)
    except FileExistsError:
        pass

    filename = 'unfiltered.txt'
    with open(temp_path+filename, 'w') as outfp:
        for i, row in enumerate(w):
            for j, el in enumerate(row[i+1:]):
                if el!=0:
                    outfp.write(str(i)+' '+str(i+j+1)+','+str(el)+'\n')

    vertices_dict = {}  # We want to re-map the vertices
    reverse_vertices_dict = {}
    filename = 'w.txt'
    with open(temp_path+filename, 'w') as outfp:
        k = 0
        for i, row in enumerate(w):
            for j, el in enumerate(row[i+1:]):
                if el!=0:
                    if i not in vertices_dict:
                        vertices_dict[i] = k
                        reverse_vertices_dict[k] = i
                        k += 1
                    if i+j+1 not in vertices_dict:
                        vertices_dict[i+j+1] = k
                        reverse_vertices_dict[k] = i+j+1
                        k += 1
                    outfp.write(
                        str(vertices_dict[i])+' '+str(vertices_dict[i+j+1])+','
                    )
                    outfp.write(str(el)+'\n')

    # Storing the actual number of nodes found
    # (here, n="nodes that are in a connected component")
    n = k
    print("Filtered out "+str(len(w) - n)+" disconnected nodes")

    # Now, writing down the vertices' labels
    filename = 'vertex_labels.txt'
    with open(temp_path+filename, 'w') as outfp:
        outfp.write(str(n)+'\n')
        for i in range(n):
            vertex_number = i
            vertex_name = inputs['v_labels'][reverse_vertices_dict[i]]
            outfp.write(
                str(vertex_number)+' '+str(vertex_name)+'\n'
            )

    # Finally, writing down the new heat file
    filename = 'heat.txt'
    with open(temp_path+filename, 'w') as outfp:
        for i in range(n):
            vertex_number = i
            heat = inputs['heat'][reverse_vertices_dict[i]]
            outfp.write(
                str(vertex_number)+' '+str(heat)+'\n'
            )


def read_from_temp(dataset_name='HINT+HI2012', motif_name='u_triangle'):
    """
    py:function:: read_from_temp(dataset_name='HINT+HI2012', motif='u_triangle')

    Parsing the input from the /in/dataset_name/temp/ folder, assuming that we
    can find the w matrix stored, there. The format should be linked to whatever
    the 'write_to_temp' wrote. This function will extract the W matrix corresponding
    to the motif given.

    :param dataset_name: String containing the name of the dataset to be loaded.
    :param motif: String containing the name of the motif wanted.

    :return: the numpy matrix containing W.
    """

    # For compatibility issues, we're extracting the absolute path of the project.
    abs_path = os.path.dirname(os.path.abspath(__file__))
    project_path = abs_path[:-3]  # WARNING This will work bc of the name of this dir
    input_path = project_path+'/in/'+dataset_name+'/'

    file_name = 'vertex_labels.txt'
    with open(input_path+file_name) as infp:
        lines = infp.readlines()
        n = int(lines[0])

    w = np.zeros((n,n), dtype=np.uint16)
    file_name = '/temp/'+'w_'+motif_name+'.txt'
    with open(input_path+file_name) as infp:
        lines = infp.readlines()
        m = len(lines)
        for l in lines:
            s = l.split(",")
            edge = s[0]
            weight = int(s[1])
            s = edge.split(" ")
            i = int(s[0])
            j = int(s[1])
            w[i,j] = weight
            w[j,i] = weight  # assuming undirected graphs

    return w

def write_outputs(subgraphs=None):
    """
        py::funcion write_outputs()
        We expect the outputs of our procedures to be... subgraphs!
        Since we already have the graph saved in our /in/ folders,
        we just want to save every subgraph as a list of vertices.

        For efficiency purposes, we can store those vertices inside tuples.

        :param subgraphs: List of tuples: each tuple are the vertices that induce a subgraph
    """

    output_path = '../out/'+dataset_name+'/'
    try:
        os.mkdir(output_path)
    except FileExistsError:
        pass

    filename = 'out.txt'
    with open(output_path+filename) as outfp:
        # TODO: vorremmo scrivere I NOMI delle proteine, non i loro numeri.
        for s in subgraphs:
            outfp.write(s)
            outfp.write('\n')
