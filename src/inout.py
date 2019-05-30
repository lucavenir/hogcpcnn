import argparse
# import os
# import time
# import datetime
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
        'database': 'HINT+HI2012',
        'delta': 0.8,  # TBD
        #'timeout': 10e12,
        'timestamp': timestamp,  # for now, this will be ignored
        # TODO: Set a warning to the user if no output folder is given;
        # it should say something like...
        # "WARNING: the default output folder is used and will be overwritten."
        'output_folder': 'current_run',
        'verbosity': 0  # by default, no debug messages
    }

    parser = argparse.ArgumentParser(description='Parser to receive user preferences and/or parameters.')

    # Input folder
    parser.add_argument('--i', type=str, help='Path of the folder containing the input')
    # Delta parameter
    parser.add_argument('--d', type=float, help='TBD: not used yet')
    # Output folder
    parser.add_argument('--o', type=str, help='Path of the desired output folder')
    # Verbosity of the debug # TODO:
    parser.add_argument('--v', type=int, help='How much ya wanna know?')

    args = parser.parse_args()

    if args.i:
        parameters['input_folder'] = args.i

    if args.d:
        parameters['delta'] = float(args.d)

    if args.o:
        parameters['output_folder'] = args.o

    return parameters

def load_input(dataset_name='HINT+HI2012'):
    """
    py:function:: load_input()

    Parsing the input from the /in/ folder, assuming that the input has a
    'standard' format and that it has been already pre-processed with the script.

    :param dataset: String containing the name of the dataset to be loaded.

    :return: a dictionary, containing all the input we need.
    """

    inputs = {}
    input_path = '../in/'+dataset_name+'/'

    # First, we extract the vertex labeling, since it's easy and it gives us n;
    # (n = |V|)
    file_name = 'vertex_labels.txt'
    with open(input_path+file_name) as infp:
        lines = infp.readlines()
        n = int(lines[0])
        v_labeling = np.empty(n, dtype=object)
        for l in lines[1:]:
            s = l.split(" ")
            idx = int(s[0])
            # To extract the label is a bit tricky, but nothing to it honestly.
            name = str(s[1].split("\n")[0])  # We want to get rid of the '\n'
            v_labeling[idx] = (idx, name)

    file_name = 'adjacency_matrix.txt'
    with open(input_path+file_name) as infp:
        lines = infp.readlines()
        m = len(lines)
        adjacency_matrix = np.zeros((n,n), dtype=np.int8)
        for l in lines:
            s = l.split(" ")
            i = int(s[0])
            j = int(s[1])
            adjacency_matrix[i,j] = 1
            adjacency_matrix[j,i] = 1  # assuming undirected graphs

    file_name = 'heat.txt'
    with open(input_path+file_name) as infp:
        lines = infp.readlines()
        heat = np.zeros(n, dtype=np.float)
        for l in lines:
            s = l.split(" ")
            v = int(s[0])
            heat_value = float(s[1])
            heat[v] = heat_value

    inputs['v_dict'] = v_labeling
    inputs['adj'] = adjacency_matrix
    inputs['heat'] = heat

    return inputs

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
