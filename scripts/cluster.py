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
    'quadrangle'
]  # Which motif should we search inside the graph?
soft = [False]  # Which version should we run?
# ks = [4,5]  # Cliques only: size of the searched clique
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

server = "login.dei.unipd.it"

# Crea un'istanza del client SSH
ssh = paramiko.SSHClient()
# <<Required, since "login.dei.unipd.it" is not a "well known ssh host">> (???)
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

# Per evitare di scrivere il proprio username, se la macchina dalla quale lo sia avvia è "nota"
if getpass.getuser() == 'venir':
    pwd = str(getpass.getpass(prompt="Password: "))
    username = "venirluca"
elif getpass.getuser() == 'iltuonomesultuoPC':  # dumbdumb
    pwd = str(getpass.getpass(prompt="Password: "))
    username = "iltuonomeDEI"
else:
    username = str(input("Please, enter your username (@" + server + "): "))
    pwd = str(getpass.getpass(prompt=username+"@" + server + ": "))

# Connesione al server
ssh.connect(server, username=username, password=pwd)

# Apertura di un'istanza per il file transfer (con secure file transfer protocol)
sftp = ssh.open_sftp()  # per trasferire i file

# Cartella remota del nostro progetto
remote_path = "/home/" + username + "/Thesis/hogcpcnn/"

# Cartella locale del nostro progetto
local_path = os.path.dirname(os.getcwd()) + "/"

print("Acquiring remote path..." + remote_path)  # /home/venirluca/Thesis/hogcpcnn/
print("Acquiring local path..." + local_path)  # /home/venir/Documents/Thesis/hogcpcnn/

# Files to be uploaded
files = ['src/__init__.py', 'src/lib.py', 'src/inout.py']
print("Excecuting...")
with open("commands.job", "w") as fp:
    for d in deltas:
        for db in databases:
            for m in motifs:
                for s in soft:
                    #for k in ks:
                    for a in alphas:
                        # Create a local file that will be sent to the server (the infamous '.job' file)
                        # Main output folder:

                        fp.write("#!/bin/bash \n")

                        # Formatting/constructing the instruction to be given:
                        instruction = "time python3 " + remote_path + "src/__init__.py"

                        # Options to be added:
                        instruction += " -db " + str(db)
                        instruction += " -m " + str(m)
                        instruction += " -a " + str(a)
                        instruction += ' -d' + str(d)
                        if s:
                            instruction += " -s"
                        if m=='clique':
                            instruction += ' -k ' + str(k)

                        # Saving the output to a log file:
                        if m=='clique':
                            output_logfilename = 'db='+str(d) + '_' + 'm='+str(k)+str(m)
                        else:
                            output_logfilename = 'db='+str(d) + '_' + 'm='+str(m)
                        if s:
                            output_logfilename += '_s'
                        output_logfilename += '_results.log'
                        # instruction += ' > '
                        # instruction += remote_path + "out/" + current_folder +'/'+ output_logfilename
                        instruction += '\n'
                        fp.write(instruction)

print("Copying files")
current_time = time.time()
timestamp = datetime.datetime.fromtimestamp(current_time).strftime('%Y%m%d_%H:%M:%S')
current_folder = "run__" + timestamp

# Dato che la cartella corrente è un timestamp, siamo sicuri di poterla creare sempre (in remoto)
sftp.mkdir(remote_path + "out/" + current_folder)

for file in files:
    file_remote = remote_path + file
    file_local = local_path + file

    print(file_local + ' >>> ' + file_remote)
    try:
        sftp.remove(file_remote)
    except IOError:
        print(file + " was not on the cluster")
        pass

    sftp.put(file_local, file_remote)

# Put the file on the current folder on the cluster and delete the local one
print(local_path + 'scripts/commands.job' ' >>> ' + remote_path + 'out/' + current_folder + '/commands.job')
sftp.put(local_path + 'scripts/commands.job', remote_path + 'out/' + current_folder + '/commands.job')
os.remove("commands.job")

# Give this job to the cluster
ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command("export SGE_ROOT=/usr/share/gridengine \n" +  #  necessary
                                                     "cd {0}out/{1} \n".format(remote_path, current_folder) +  # also necessary
                                                     "qsub -cwd commands.job")
# dev'essere "tutto assieme"
# o si dimentica dell'export.
# una singola chiamata di exec_command fa dimenticare tutto quello che è stato fatto prima
# qsub -cwd == "current working directory". DEVE ESSERE MESSO PRIMA!!!

# Print output and errors
print(ssh_stdout.read().decode('utf-8'))
print(ssh_stderr.read().decode('utf-8'))

time.sleep(5-.250)

sftp.close()
ssh.close()
