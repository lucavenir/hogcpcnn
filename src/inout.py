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

    current_time = time.time()
    timestamp = datetime.datetime.fromtimestamp(current_time).strftime('%Y%m%d_%H:%M:%S')
    default_folder = "run__" + timestamp

    # Initializing the dictionary containing the parameters with default values
    parameters = {
        'dataset': 'HINT+HI2012',
        'motif': 'triangle',
        'k': 4,
        'method': 'HOTNET2',
        'soft': False,
        'delta': 0.000496,  # this default value comes from HotNet2: it doesn't necessary makes sense
        'timestamp': timestamp  # for now, this will be ignored
    }

    parser = argparse.ArgumentParser(description='Parser to receive user preferences and/or parameters.')

    # Input folder
    parser.add_argument('-db', type=str, help='Name of the dataset used')
    # Motif searched
    parser.add_argument('-m', type=str, help='Motif type')
    # Diffusion process algorithm selection
    parser.add_argument('-k', type=int, help='Clique size, when selected')
    parser.add_argument('-f', type=str, help='Method to compute f')
    # Delta parameter selection
    parser.add_argument('-d', type=float, help='Cut-off value for edges in the graph H')
    # Soft version of the problem?
    parser.add_argument('-s', action='store_true', help='Type -s to enable the soft version')
    args = parser.parse_args()

    if args.d:
        parameters['delta'] = float(args.d)

    if args.db:
        parameters['dataset'] = args.db

    if args.m:
        parameters['motif'] = args.m

    if args.k:
        parameters['k'] = args.k

    if args.f:
        parameters['method'] = args.f

    if args.s:
        parameters['soft'] = True

    return parameters

def load_input(dataset_name='HINT+HI2012', motif_name=None, soft=False, k=4):
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
        if motif_name == 'clique':
            motif_path = input_path + 'temp/'+str(k)+motif_name
        else:
            motif_path = input_path + 'temp/'+motif_name

        if soft:
            motif_path += '_s'
        motif_path += '/'


    # First, we extract the vertex labeling, since it's easy and it gives us n;
    # (n = |V|)
    file_name = 'vertex_labels.txt'
    file_path = input_path if motif_name==None else motif_path
    with open(file_path+file_name) as infp:
        lines = infp.readlines()
        n = int(lines[0])
        k = int(lines[1])
        v_labeling = np.empty(k, dtype=object)
        rev_labeling = {}
        for l in lines[2:]:
            s = l.split(" ")
            idx = int(s[0])
            # To extract the label is a bit tricky, but nothing to it honestly.
            name = str(s[1].split("\n")[0])  # We want to get rid of the '\n'
            v_labeling[idx] = name
            rev_labeling[name] = idx

    file_name = 'adjacency_matrix.txt'
    with open(input_path+file_name) as infp:
        lines = infp.readlines()
        adjacency_matrix = np.zeros((n,n), dtype=np.uint8)
        for l in lines:
            s = l.split(" ")
            i = int(s[0])
            j = int(s[1])
            adjacency_matrix[i,j] = 1
            adjacency_matrix[j,i] = 1  # assuming undirected graphs

    if motif_name!=None:
        file_name = 'w.txt'
        with open(motif_path+file_name) as infp:
            lines = infp.readlines()
            w = np.zeros((k,k), dtype=np.float64)
            for l in lines:
                s = l.split(",")
                edge = s[0]
                weight = float(s[1])
                s = edge.split(" ")
                i = int(s[0])
                j = int(s[1])
                w[i,j] = weight
                w[j,i] = weight  # assuming undirected graphs

    file_name = 'heat.txt'
    file_path = input_path if motif_name==None else motif_path
    with open(file_path+file_name) as infp:
        lines = infp.readlines()
        heat = np.zeros(k, dtype=np.float64)
        for l in lines:
            s = l.split(" ")
            v = int(s[0])
            heat_value = float(s[1])
            heat[v] = heat_value

    file_name = 'classics.txt'
    file_path = project_path+'in/'
    with open(file_path+file_name) as infp:
        classics_list = np.zeros(k, dtype=np.bool)
        lines = infp.readlines()
        for l in lines:
            gene = l.split('\n')[0]
            try:
                classics_list[rev_labeling[gene]] = True
            except KeyError:
                pass

    inputs['classics'] = classics_list
    inputs['rev_labels'] = rev_labeling
    inputs['v_labels'] = v_labeling
    inputs['heat'] = heat
    inputs['adj'] = adjacency_matrix
    if motif_name != None:
        inputs['w'] = w

    return inputs

def write_transition_matrix(inputs, w, dataset_name='HINT+HI2012', motif_name='triangle', soft=False, k=4):
    # For compatibility issues, we're extracting the absolute path of the project.
    abs_path = os.path.dirname(os.path.abspath(__file__))
    project_path = abs_path[:-3]  # This will work just bc of the name of this dir

    temp_path = project_path+'/in/'+dataset_name+'/temp/'
    try:
        os.mkdir(temp_path)
    except FileExistsError:
        pass

    if motif_name=='clique':
        temp_path += str(k)+motif_name
    else:
        temp_path += motif_name

    if soft:
        temp_path += '_s'
    temp_path += '/'

    try:
        os.mkdir(temp_path)
    except FileExistsError:
        pass

    k = len(w)
    n = len(inputs['adj'])
    filename = 'w.txt'
    with open(temp_path+filename, 'w') as outfp:
        for i, row in enumerate(w):
            for j, el in enumerate(row):
                if el!=0:
                    outfp.write(
                        str(i)+' '+str(j)+','
                    )
                    outfp.write(str(el)+'\n')

    # Now, writing down the vertices' labels
    filename = 'vertex_labels.txt'
    with open(temp_path+filename, 'w') as outfp:
        outfp.write(str(n)+'\n')
        outfp.write(str(k)+'\n')
        for i in range(k):
            vertex_number = i
            vertex_name = inputs['v_labels'][i]
            outfp.write(
                str(vertex_number)+' '+str(vertex_name)+'\n'
            )

    # Finally, writing down the new heat file
    filename = 'heat.txt'
    with open(temp_path+filename, 'w') as outfp:
        for i in range(k):
            vertex_number = i
            heat = inputs['heat'][i]
            outfp.write(
                str(vertex_number)+' '+str(heat)+'\n'
            )

def write_strong_ccs(s_cc_list, v_labels, parameters):
    # For compatibility issues, we're extracting the absolute path of the project.
    abs_path = os.path.dirname(os.path.abspath(__file__))
    project_path = abs_path[:-3]  # WARNING: This will work just bc of the name of this dir ('src')

    out_path = project_path+'/out/'+parameters['dataset']+'/'
    try:
        os.mkdir(out_path)
    except FileExistsError:
        pass

    if parameters['soft']:
        motif_name += '_s'

    out_path += 'delta='+str(parameters['delta'])+'/'
    try:
        os.mkdir(out_path)
    except FileExistsError:
        pass

    if parameters['motif'] == 'clique':
        filename = str(parameters['k'])+parameters['motif']+'_delta='+str(parameters['delta'])+'.txt'
    else:
        filename = parameters['motif']+'_delta='+str(parameters['delta'])+'.txt'

    with open(out_path+filename, 'w') as outfp:
        m = len(s_cc_list)
        outfp.write("Found "+str(m)+" strongly connected components.\n")
        for i,set in enumerate(s_cc_list):
            outfp.write('#'+str(i+1)+' (len: '+ str(len(set)) +'): {')
            set_list = sorted(
                list(set),
                key=lambda el:v_labels[el]
            )
            for el in set_list:
                outfp.write(' '+str(v_labels[el]))
            outfp.write(' }\n')
