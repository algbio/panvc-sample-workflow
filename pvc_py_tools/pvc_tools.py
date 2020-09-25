# Copyright (c) Daniel Valenzuela, Tuukka Norri 2019–2020
# Licenced under the MIT licence.

import os.path
import sys
import gzip
from pathlib import Path
from datetime import datetime

    
def PVC_read_len_from_reads(reads_file):
    start = datetime.now()
    print (f"{start} Read len from: {reads_file}")
    file_path, file_extension = os.path.splitext(reads_file)
    read_len = 0
    if file_extension == ".gz":
        with gzip.open(reads_file, 'r') as myfile:
            a_line = myfile.readline()
            a_read = myfile.readline().decode().rstrip('\n')
    else:
        with open(reads_file, 'r') as myfile:
            a_line = myfile.readline()
            a_read = myfile.readline().rstrip('\n')
    read_len =len(a_read)
    end = datetime.now()
    print(f"{end} Time taken (PVC_read_len_from_reads): {end - start}")
    return read_len

def PVC_save_var(var_value, var_name, target_dir):
    file_name = target_dir + "/" + var_name + ".txt" 
    with open(file_name, "w") as g:
        g.write(str(var_value))

def PVC_load_var(var_name, target_dir):
    file_name = target_dir + "/" + var_name + ".txt" 
    assert os.path.isfile(file_name), f"Expected {file_name} to be a file"
    with open(file_name, 'r') as myfile:
        ans = myfile.read()
    return int(ans)

def PVC_get_chr_list(pgindex_dir):
    chr_list_file = pgindex_dir + "/chr_list.txt"
    chr_list = [line.rstrip('\n') for line in open(chr_list_file)]
    return chr_list

def PVC_get_first_chr(chr_list_file):
    with open(chr_list_file) as f:
        line = f.readline()
        first_chr = str(int(line))
    return first_chr

def call_or_die(command):
    import subprocess
    start = datetime.now()
    print(f"{start} About to call:")
    print(command)
    result = subprocess.run(command, shell=True, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, executable="/bin/bash")
    end = datetime.now()
    print(f"{end} Time taken (call_or_die): {end - start}")
    if (result.returncode == 0):
        print ("Stdout:\n" + result.stdout.decode())  
        print ("Stderr:\n" + result.stderr.decode())  
    elif (command.startswith("grep") and result.returncode == 1):
        assert(result.stdout.decode() == "")
        assert(result.stderr.decode() == "")
    else:
        print ("Return code not zero:" + str(result.returncode))
        print ("Stdout:\n" + result.stdout.decode())  
        print ("Stderr:\n" + result.stderr.decode())  
        exit(result.returncode)
    print("--------------")
