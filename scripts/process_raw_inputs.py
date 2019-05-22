# Executable script meant to process raw inputs from the datasets
# The goal is to write standard input files, easy to interpret.
# The reason behind this pre-processing phase is that in this project
# multiple datasets are involved, each with different "notations".

# We want a "good" input such as the following:
#       - VERTEXES: they should be just.. numbers!
#           - the i-th gene should be the (i-1)-th entry of the matrix
#           (the i-1 happens because of indexing)
#       - ADJACENCY MATRIX: each row should be "i j"
#           (indicating that the edge between i and j exists)
#       - VERTEX LABELS: "i name_i"
#           (vertex i's name is "name_i")
#       - SAMPLES: It depends on the version we're using.
#               DETERMINISTIC: sample_id [mutated_gene_1, ..., mutated_gene_m]
#               PROBABILISTIC: sample_id [(gene, p_muted), ...]
#           (for the probabilistic version, the list consist of proabilities > 0)
#           (when not present, p=0)

# The OUTPUT files of this script should be the following:
input_folder = '../in/'
adjacency_matrix = 'adjacency_matrix.txt'
vertex_labels = 'vertex_labels.txt'
samples = 'samples.txt'

def hinthi2012():
    # Processing HINT+HI2012 database.
    # hint+hi2012_edge_file.txt has (almost) the same structure we desire.
    # All we need to do is to get rid of the 1s and 0s, useless for our purposes.

    # File location
    dataset_name = 'HINT+HI2012'
    rawinputs_folder = '../rawinput/'
    edges_file = 'hint+hi2012_edge_file.txt'
    index_file = 'hint+hi2012_index_file.txt'
    samples_file = 'mutation_frequency_expr_filtered.txt'
    input_folder += dataset_name+'/'

    # Processing the edge matrix first
    with open(rawinputs_folder+dataset_name+'/'+edges_file, 'r') as infp:
        lines = fp.readlines()
        outfp = open(input_folder+adjacency_matrix, 'w')
        for l in lines:
            s=l.split(" ")
            outfp.write(s[0]+' '+s[1]+'\n')

        outfp.close()

    # Then, processing the vertexes' labels
    with open(rawinputs_folder+dataset_name+'/'+index_file, 'r') as infp:
        lines = fp.readlines()
        outfp = open(input_folder+vertex_labels, 'w')
        for l in lines:
            s=l.split(" ")
            outfp.write(s[0]+' '+s[1].split("\t")[0]+'\n')

        outfp.close()

    # Then, again, processing the samples
    with open(rawinputs_folder+samples_file, 'r') as infp:
        lines = fp.readlines()
        outfp = open(input_folder+samples, 'w')
        for l in lines:
            # TODO


if __name__ == "__process_raw_inputs__":
    # Here, uncomment the function corresponding to the database you want to decode.
    hinthi2012()
    #irefindex()
    #multinet()
