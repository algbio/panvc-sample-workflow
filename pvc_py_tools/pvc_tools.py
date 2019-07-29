import os.path
import sys
import gzip
from config import *
from pathlib import Path
#from config import *

def PVC_make_std_ref(args):
    output_file = args.output_dir + "/std_ref.fa"
    out_file = open(output_file, "w")
    with open(args.chr_list_file) as f:
        for line in f:
            chr_id = str(int(line))
            print ("processing chr: " + chr_id)
            out_file.write(">" + chr_id + "\n")
            src_file = args.output_dir + "/" + chr_id + "/recombinant.n1.gapped" 
            assert(Path(src_file).is_file())
            with open(src_file) as g:
                for seq in g:
                    out_file.write(seq.replace("-","").rstrip() + "\n")

    
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

def PVC_validate_names(pgindex_dir):
    #TODO: implement me
    """
    F_CHR=$( head -n1 ${CHR_LIST_FILE})
    if [ ${DEBUG_MODE} -eq 1 ]; then
    cat ${CHR_LIST_FILE} | while read CHR_ID
    do
      utils_assert_files_are_equal ${PG_INDEX_OUTPUT_DIR}/${F_CHR}/n_refs.txt ${PG_INDEX_OUTPUT_DIR}/${CHR_ID}/n_refs.txt
      utils_assert_namefiles_are_almost_equal ${PG_INDEX_OUTPUT_DIR}/${F_CHR}/names.plain ${PG_INDEX_OUTPUT_DIR}/${CHR_ID}/names.plain
    done
    fi
    """
    pass

#TODO: receive output dir as well, we do not want this file to lie in the input folder, which is the current behavior
def PVC_merge_reads(reads_1, reads_2):
    ## assert reads_1, 2, are valid files.
    basedir = os.path.dirname(reads_1)
    output = basedir + "/reads_ALL_RENAMED.fq.gz"
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

def PVC_sequence_num_to_name(SAMPLES_FILENAME, CHRS_FILENAME, PLOIDITY, CHR_ID, SEQ_ID):
    assert (SEQ_ID >= 1)
    hap_suffix = ""
    if (SEQ_ID != 1) :
        offset = SEQ_ID - 2
        target_line = 1 + (offset//PLOIDITY)
        hap_suffix = "_" + str(offset % PLOIDITY)
    else:
        target_line = 0
        hap_suffix = ""

    with open(SAMPLES_FILENAME, 'r') as f:
        curr_name = f.readline()

    # Open and read file into buffer
    f = open(SAMPLES_FILENAME,"r")
    sample_names = f.readlines()

    total_chroms = sum(1 for line in open(CHRS_FILENAME))
    
    chr_suffix = ""
    if (total_chroms != 1) :
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
    for line in f:
        if ";" in line :
            continue
        if ">" in line :
            out_num = out_num + 1
            #close?
            out_file = open(target_folder + '/recombinant.n'+str(out_num)+'.gapped', 'w')
            continue
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
    ret_code = subprocess.call([command], shell=True, executable="/bin/bash")
    if ret_code != 0:
        print ("Return code not zero:" + str(ret_code))
        sys.exit(ret_code)

def file_exists_or_die(filename):
    if not os.path.isfile(filename):
        print("******")
        print("Cannot find " + filename)
        print("Aborting execution")
        print("******")
        sys.exit(1);


