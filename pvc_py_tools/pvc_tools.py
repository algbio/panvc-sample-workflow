import os.path
import sys
import gzip
from config import *
from pathlib import Path
from datetime import datetime

def PVC_make_std_ref(args):
    start = datetime.now()
    print(f"{start} Making std_ref.fa")
    output_file = args.output_dir + "/std_ref.fa"
    with open(args.chr_list_file) as f, open(output_file, "w") as out_file:
        for line in f:
            chr_id = str(int(line))
            print ("processing chr: " + chr_id)
            out_file.write(">" + chr_id + "\n")
            src_file = args.output_dir + "/" + chr_id + "/recombinant.n1.gapped" 
            assert(Path(src_file).is_file())
            with open(src_file) as g:
                for seq in g:
                    out_file.write(seq.replace("-","").rstrip() + "\n")
    # Index here in order to allow multiple instances of panvc_call_variants operate on the same index.
    end = datetime.now()
    print(f"{end} Time taken (PVC_make_std_ref): {end - start}")
    PVC_index_ref(output_file)
    
def PVC_index_ref(output_file):
    start = datetime.now()
    print(f"{start} Indexing reference")
    if not os.path.isfile(output_file + ".bwt"):
        bwa_index_command = BWA_BIN + " index " + output_file
        call_or_die(bwa_index_command)

    if not os.path.isfile(output_file + ".fai"):
        faidx_command = SAMTOOLS_BIN + " faidx " + output_file
        call_or_die(faidx_command)

    output_ext_removed = os.path.splitext(output_file)[0]
    output_dict = output_ext_removed + ".dict"
    if not os.path.isfile(output_dict):
        gatk_dict_command = GATK_BIN + " CreateSequenceDictionary --REFERENCE=" + output_file + " --OUTPUT=" + output_dict
        call_or_die(gatk_dict_command)
    end = datetime.now()
    print(f"{end} Time taken (PVC_index_ref): {end - start}")

    
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

## This validates that the samples names are the same for all the chromosomes. 
## If this assumption is not true, the pipeline breaks.
def PVC_validate_names(pgindex_dir):
    import filecmp
    start = datetime.now()
    print (f"{start} Validating read names")
    chr_list = PVC_get_chr_list(pgindex_dir)
    first_chr = chr_list[0]
    for curr_chr in chr_list:
        f1 = pgindex_dir + "/" + first_chr + "/n_refs.txt"
        f2 = pgindex_dir + "/" + curr_chr+ "/n_refs.txt"
        if (not filecmp.cmp(f1, f1)):
            print ("Inconsistent input: different chromosmes have different number of genomes in the pan-genome")
            exit(33)
        samples_1 = pgindex_dir + "/" + first_chr + "/names.plain"
        samples_2 = pgindex_dir + "/" + curr_chr+ "/names.plain"
        if (not PVC_compare_samples_names(samples_1, samples_2)):
            print ("Inconsistent input: samples have different names in different chromosmes in the pan-genome")
            print ("File 1: " + samples_1)
            print ("File 2: " + samples_2)
            exit(33)
    end = datetime.now()
    print(f"{end} Time taken (PVC_validate_names): {end - start}")

def PVC_compare_samples_names(file_1, file_2):
    start = datetime.now()
    print (f"{start} Comparing sample names")
    f1 = open(file_1,"r")
    names_1 = f1.readlines()
    f2 = open(file_2,"r")
    names_2 = f2.readlines()
    retval = True
    if len(names_1) != len(names_2):
        retval = False
    for i in range(1, len(names_1)):
        if names_1[i] != names_2[i]:
            retval = False
            break
    end = datetime.now()
    print(f"{end} Time taken (PVC_compare_samples_names): {end - start}")
    return retval

def PVC_get_first_chr(chr_list_file):
    with open(chr_list_file) as f:
        line = f.readline()
        first_chr = str(int(line))
    return first_chr

def PVC_get_system_ulimit():
    nfiles = call_and_get_result("ulimit -n")
    if (nfiles == "unlimited"):
        nfiles = "10000"
    return int(nfiles)

def PVC_merge_reads(reads_1, reads_2, output_folder):
    start = datetime.now()
    print (f"{start} Merging reads")
    #TODO: replace bash_command with a native python function. A system call with pipes can be dangerous.
    assert(Path(reads_1).is_file())
    assert(Path(reads_2).is_file())
    output = output_folder + "/reads_ALL_RENAMED.fq.gz"
    if os.path.isfile(output):
        return output
    file_path, file_extension = os.path.splitext(reads_1)
    bash_command = ""
    if file_extension == ".gz":
        bash_command = "awk '{if ((NR-1)%4==0) print \"@\"1+(NR-1)/4; else print }' <(zcat "+ reads_1 +") <(zcat " + reads_2 + ") | gzip > " + output
    else:
        bash_command = "awk '{if ((NR-1)%4==0) print \"@\"1+(NR-1)/4; else print }' " + reads_1 + " " + reads_2 + " | gzip > " + output
    call_or_die(bash_command)
    print ("Done")
    end = datetime.now()
    print(f"{end} Time taken (PVC_merge_reads): {end - start}")
    return output

def PVC_splitA2M(A2M_FILENAME, target_folder):
    start = datetime.now()
    print (f"{start} splitting A2M")
    n_refs = 0
    f = open(A2M_FILENAME,"r")
    out_num = 0
    valid_a2m_input = False
    for line in f:
        if ";" in line :
            continue
        if ">" in line :
            out_num = out_num + 1
            #close?
            out_file = open(target_folder + '/recombinant.n'+str(out_num)+'.gapped', 'w')
            valid_a2m_input = True
            continue
        if (not valid_a2m_input):
            print("Input file is not a proper A2M file. Headers ('>name') were missing")
            exit(33)
        out_file.write(line)
        n_refs = n_refs + 1
    end = datetime.now()
    print(f"{end} Time taken (PVC_splitA2M): {end - start}")
    return n_refs

def PVC_gapped_len(filename):
    f = open(filename, "r")
    line = f.readline().rstrip()
    return len(line)

# what about newlines? shall I remove them?

def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def ensure_clean_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        import shutil
        shutil.rmtree(directory)
        os.makedirs(directory)

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

def call_and_get_result(command):    
    import subprocess
    start = datetime.now()
    print(f"{start} About to call:")
    print (command)
    result = subprocess.run(command, shell=True, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, executable="/bin/bash")
    end = datetime.now()
    print(f"{end} Time taken (call_and_get_result): {end - start}")
    if (result.returncode != 0):
        print ("Return code not zero:" + str(result.returncode))
        print ("Stdout:\n" + result.stdout.decode())  
        print ("Stderr:\n" + result.stderr.decode())
        exit(result.returncode)
    print ("--------------")
    return result.stdout.decode().rstrip()

def file_exists_or_die(filename):
    if not os.path.isfile(filename):
        print("******")
        print("Cannot find " + filename)
        print("Aborting execution")
        print("******")
        sys.exit(1);


