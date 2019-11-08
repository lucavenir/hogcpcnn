#!/usr/bin/env python3
import os
import time
import datetime
import paramiko
import getpass

'''
    --- WARNING ---
    FIRST TIME INFO, MANDATORY TODOs

    (1)     Se non hai mai avviato questo script, come prima cosa esegui il login:
            >>> ssh username@login.dei.unipd.it
            >>> git clone <url>

            (sì, la repository nel nostro spazio dei è necessaria)
            tl; dr: clonare la repository dentro la 'home' del nostro spazio sul DEI.

    (2)     Per evitare inutili disagi, eseguire questo script sempre da dentro la cartella src/.
'''

# Parametri per l'esecuzione

# input files
databases = ['HINT+HI2012']
motifs = [
    'clique'
]  # Which motif should we search inside the graph?
soft = [False]  # Which version should we run?
ks = [4,5]  # Cliques only: size of the searched clique
deltas = [
    0.0005,
    0.0055,
    0.0267,
    0.0675,
    0.1345,
    0.1000,
    0.22612
]
alphas = [
    0.2,
    0.4,
    0.5,
    0.6,
    0.8
]

# Scripts to be executed
print("Excecuting...")
for d in deltas:
    for db in databases:
        for m in motifs:
            for s in soft:
                for k in ks:
                    for a in alphas:
                        # Create a local file that will be sent to the server (the infamous '.job' file)
                        # Main output folder:

                        # Formatting/constructing the instruction to be given:
                        instruction = "python3 ../src/__init__.py"

                        # Options to be added:
                        instruction += " -db " + str(db)
                        instruction += " -m " + str(m)
                        instruction += " -a " + str(a)
                        instruction += ' -d' + str(d)
                        if s:
                            instruction += " -s"
                        if m=='clique':
                            instruction += ' -k ' + str(k)

                        print("Executing the following instruction: ")
                        print(instruction)
                        os.system(instruction)
