# TBD
# import pandas as pd
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
        'input_folder': '../in/random_input_example',
        'delta': 0.8,  # TBD
        #'timeout': 10e12,
        'timestamp': timestamp,  # for now, this will be ignored
        # TODO: Set a warning to the user if no output folder is given;
        # it should say something like...
        # "WARNING: the default output folder is used and will be overwritten."
        'output_folder': '../out/current_run'
    }

    parser = argparse.ArgumentParser(description='Parser to receive user preferences and/or parameters.')

    # Input folder
    parser.add_argument('--i', type=str, help='Path of the folder containing the input')
    # Delta parameter
    parser.add_argument('--d', type=float, help='TBD: not used yet')
    # Output folder
    parser.add_argument('--o', type=str, help='Path of the desired output folder')

    args = parser.parse_args()

    if args.i:
        parameters['input_folder'] = args.i

    if args.d:
        parameters['delta'] = float(args.d)

    if args.o:
        parameters['output_folder'] = args.o

    return parameters

def load_input():
    # TODO: dopo molti ripensamenti, la struttura dati migliore è una matrice
    # NxN, con "indici" da 0 a N-1. Chiaramente, abbiamo N molto alti e con
    # "buchi" in mezzo (es. ci sono nodi che compaiono in altri file ma
    # non compaiono nella rete). Soluzione: creare una mappa su un file temporaneo
    # che permette di indicizzare i nodi da 1 a TOT, senza avere buchi in mezzo.
    # QUINDI: cartella "temp" da creare in ogni istanza (e.g. "HINT+HI2012").
    # DA IMPLEMENTARE IN "process_raw_inputs.py"
    pass
