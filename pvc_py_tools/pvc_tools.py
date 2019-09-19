import os.path
import sys
import gzip
from config import *
from pathlib import Path

def PVC_make_std_ref(args):
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
    bwa_index_command = BWA_BIN + " index " + output_file
    call_or_die(bwa_index_command)

    faidx_command = SAMTOOLS_BIN + " faidx " + output_file
    call_or_die(faidx_command)

    output_dict = args.output_dir + "/std_ref.dict"
    gatk_dict_command = GATK_BIN + " CreateSequenceDictionary --REFERENCE=" + output_file + " --OUTPUT=" + output_dict
    call_or_die(gatk_dict_command)

    
def PVC_read_len_from_reads(reads_file):
    print ("Read len from: " + reads_file)
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
    return read_len

def PVC_save_var(var_value, var_name, target_dir):
    file_name = target_dir + "/" + var_name + ".txt" 
    with open(file_name, "w") as g:
        g.write(str(var_value))

def PVC_load_var(var_name, target_dir):
    file_name = target_dir + "/" + var_name + ".txt" 
    assert(os.path.isfile(file_name))
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

def PVC_compare_samples_names(file_1, file_2):
    f1 = open(file_1,"r")
    names_1 = f1.readlines()
    f2 = open(file_2,"r")
    names_2 = f2.readlines()
    if len(names_1) != len(names_2):
        return False
    for i in range(1, len(names_1)):
        if names_1[i] != names_2[i]:
            return False
    return True

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
    return output

def PVC_sequence_num_to_name(SAMPLES_FILENAME, n_chroms , PLOIDY, CHR_ID, SEQ_ID):
    assert (SEQ_ID >= 1)
    hap_suffix = ""
    if (SEQ_ID != 1) :
        offset = SEQ_ID - 2
        target_line = 1 + (offset//PLOIDY)
        hap_suffix = "_" + str(offset % PLOIDY)
    else:
        target_line = 0
        hap_suffix = ""

    with open(SAMPLES_FILENAME, 'r') as f:
        curr_name = f.readline()

    # Open and read file into buffer
    f = open(SAMPLES_FILENAME,"r")
    sample_names = f.readlines()

    chr_suffix = ""
    if (n_chroms != 1) :
        chr_suffix="_CHR"+str(CHR_ID)

    #we are currently deleting whitespaces from the samples names;
    #if a replacement is needed just change the next line.
    SPACEREPLACEMENT=""
    #answer="PG_REF_"+str(SEQ_ID)+"_CHR"+str(CHR_ID)
    answer = sample_names[target_line].strip().replace(" ", SPACEREPLACEMENT)+hap_suffix+chr_suffix
    return answer

def PVC_splitA2M(A2M_FILENAME, target_folder):
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
    print ("About to call:")
    print (command)
    print ("--------------")
    result = subprocess.run(command, shell=True, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, executable="/bin/bash")
    if (result.returncode == 0):
        return
    if (command.startswith("grep") and result.returncode == 1):
        assert(result.stdout.decode() == "")
        assert(result.stderr.decode() == "")
        return
    else:
        print ("Return code not zero:" + str(result.returncode))
        print ("Stdout:\n" + result.stdout.decode())  
        print ("Stderr:\n" + result.stderr.decode())  
        exit(result.returncode)
def call_and_get_result(command):    
    import subprocess
    print ("About to call:")
    print (command)
    print ("--------------")
    result = subprocess.run(command, shell=True, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, executable="/bin/bash")
    if (result.returncode == 0):
        return result.stdout.decode().rstrip()
    else:
        print ("Return code not zero:" + str(result.returncode))
        print ("Stdout:\n" + result.stdout.decode())  
        print ("Stderr:\n" + result.stderr.decode())
        exit(result.returncode)

def file_exists_or_die(filename):
    if not os.path.isfile(filename):
        print("******")
        print("Cannot find " + filename)
        print("Aborting execution")
        print("******")
        sys.exit(1);


