# Executable script meant to plot the TPR/FPR analysis.
# The goal is to write down easy-to-understand plots of the TPR and FPR.
# One single plot will have every possible value of delta considered.

# In our case, such delta values considered are:
#   - 0.000496 (paper value)
#   - 0.0055 (approx. 10 times more)
#   - 0.067 (approx. 10 times more)
#   - 0.019196354 (geometric mean between the last two)
#   - TBD: delta value to be computed for each motif like they did on the paper

import matplotlib.pyplot as plt
import os

# The INPUT files of this script should be the following:
input_folder = '../out/'
database = 'HINT+HI2012'
motifs = [
    'triangle',
    'tailed_triangle',
    'quadrangle',
    'tailed_quadrangle',
    '4clique',
    '5clique',
    'no'
]
colours = [
    'b',
    'g',
    'r',
    'c',
    'm',
    'y',
    'k',
]
delta = 0.000496

# The OUTPUT files of this script should be the following:
output_folder = '../out/plots/'

def plotting():
    # Referencing to the global variables:
    global input_folder
    global database
    global motifs
    global delta
    global output_folder
    global colours

    # File locations
    input_folder += database+'/delta='+delta+'/'

    try:
        os.mkdir(output_folder)
    except FileExistsError:
        pass

    tpr_values = dict()
    fpr_values = dict()
    tprfpr_ks = dict()
    for motif in motifs:
        tprfpr_file = 'tprfpr_'+motif+'_delta='+delta+'.txt'
        tprfpr_ks[motif] = list()
        tpr_values[motif] = list()
        fpr_values[motif] = list()

        with open(input_folder+tprfpr_file, 'r') as infp:
            lines = infp.readlines()
            for l in lines:
                s = l.split(' ')
                k = int(s[0].split('=')[1].split(':')[0])
                s1 = s.split(',')
                tpr = float(s1[0])
                fpr = float(s1[1])

                tprfpr_ks[motif].append(k)
                tpr_values[motif].append(tpr)
                fpr_values[motif].append(fpr)

    plt.title('Results with delta='+delta)
    plt.ylabel('TPR')
    plt.xlabel('FPR')
    plt.axis([0,1,0,1])

    for i, motif in enumerate(motifs):
        tpr = tpr_values[motif]
        fpr = fpr_values[motif]
        k = tprfpr_ks[motif]
        plt.plot(fpr, tpr, fmt='-'+colours[i], label=motif)

    plt.legend()
    plt.show()

if __name__ == "__main__":
    plotting()
