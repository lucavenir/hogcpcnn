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
import string
import random

# The OUTPUT files of this script should be the following:
input_folder = '../in/'
adjacency_matrix = 'adjacency_matrix.txt'
vertex_labels = 'vertex_labels.txt'
samples = 'heat.txt'

def random_string(str_len=5):
    """
        Generates a random string of fixed length
    """
    letters = string.ascii_uppercase
    return ''.join(random.choice(letters) for i in range(str_len))

def random_graph(size='small', edge_density='low'):
    """
        Generates a random graph of a given vertex and edge size.
    """
    # Referencing to the global variables:
    global input_folder
    global adjacency_matrix
    global vertex_labels
    global samples

    # Graph size
    n = 300 if size=='small' else 5000 if size=='medium' else 10000
    p = 0.05 if edge_density=='low' else 0.10 if edge_density=='medium' else 0.15

    # File locations
    dataset_name = 'random_'+size+'size_'+edge_density+'edgedensity'
    input_folder += dataset_name+'/'

    # Generating the vertices first
    vertices = [i for i in range(n)]

    # Then, the edges
    edges = []
    for i in range(n):
        for j in range(n):
            if random.random() < p:
                edges.append((i,j))

    # Finally, the heat values (at random, with an arbitrary bound)
    heat = [
        (v, random.random()*n/1000 if random.random()*n/1000>p/2 else 0)
        for v in vertices
    ]

    # Creating the input folder we need
    try:
        os.mkdir(input_folder)
    except FileExistsError:
        pass

    with open(input_folder+vertex_labels, 'w') as outfp:
        outfp.write(str(n)+'\n')
        for v in vertices:
            outfp.write(str(v)+' '+random_string()+'\n')

    with open(input_folder+samples, 'w') as outfp:
        for t in heat:
            outfp.write(str(t[0])+' '+str(t[1])+'\n')

    with open(input_folder+adjacency_matrix, 'w') as outfp:
        for e in edges:
            outfp.write(str(e[0])+' '+str(e[1])+'\n')

if __name__ == "__main__":
    random_graph('small', 'low')
