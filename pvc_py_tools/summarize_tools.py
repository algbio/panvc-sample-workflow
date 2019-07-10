import os.path
import sys
import subprocess
from config import *

def file_exists_or_die(filename):
    if not os.path.isfile(filename):
        print("******")
        print("Cannot find " + filename)
        print("Aborting execution")
        print("******")
        sys.exit(1);


## To understand the next two functions see: 
## https://www.biostars.org/p/138116/#138118
## and the SAM/BAM file specification.

def get_NMAPS_FROM_SINGLE_END_READS(bam_filename):
    command_args = [SAMTOOLS_PATH, "view", "-F", "0x904", "-c", bam_filename]
    result = subprocess.run(command_args, stdout=subprocess.PIPE)
    output = result.stdout.decode('utf-8').rstrip()
    return output

def get_NMAPS_FROM_PAIRED_END_READS(bam_filename):
    command_str = " ".join([SAMTOOLS_PATH, "view", "-F", "0x4", bam_filename, " | cut -f 1 | sort | uniq | wc -l "])
    result = subprocess.run(command_str, shell=True, stdout=subprocess.PIPE)
    output = result.stdout.decode('utf-8').rstrip()
    return output
    

def count_reads(reads_filename):
    command_args = ["cat", reads_filename, "| awk 'NR % 4 == 1' | tr ' ' '_' | cut -c 2- | sort | uniq | wc -l"]
    if reads_filename.endswith(".gz"):
        command_args[0] = "zcat"

    command_str = " ".join(command_args)
    #print ("About to run:", command_str)
    result = subprocess.run(command_str, shell=True, stdout=subprocess.PIPE)
    output = result.stdout.decode('utf-8').rstrip()
    return output

def get_NREADS(reads_1, reads_2, input_reads_are_paired):
    answer = count_reads(reads_1)
    if input_reads_are_paired :
        answer = answer + count_reads(reads_2)
    return answer

def SummarizeAdhocRef(adhocref_filename, ref_filename):
    input_file = open(adhocref_filename)
    adhoc_seq = input_file.readline().rstrip('\n')
    input_file.close()
    
    input_file = open(ref_filename)
    ref_seq = input_file.readline().rstrip('\n')
    input_file.close()

    assert(len(ref_seq) == len(adhoc_seq))
    diff = 0
    for (x,y) in zip (ref_seq, adhoc_seq):
        if x != y:
            diff = diff + 1
    
    #answer = ref_seq + "\n" + adhoc_seq + "\n"
    answer = ""
    answer += "Number of positions changed to generate ad-hoc ref: " + str(diff) + "\n"
    return answer




def write_summary(working_dir, pgindex_dir, input_reads_are_paired, reads_1, reads_2):
    if (not input_reads_are_paired): assert reads_2 == ""
    
    BAMFILE_to_adhoc = working_dir + "/ext_vc_files/aligned_deduplicated.bam"
    
    #TODO: verify the following is the file I want. 
    #E.g. we may want to  deduplicate, and see if it will make any difference?
    BAMFILE_to_pg = working_dir + "/sam_files/all_sorted.bam"
    ## Assert file exists
    file_exists_or_die(BAMFILE_to_pg)
    file_exists_or_die(BAMFILE_to_adhoc)
    
    N_INPUT_READS = get_NREADS(reads_1, reads_2, input_reads_are_paired)

    ## Even if the reads were paired-end, we did treat them (and renamed) as if they were
    ## singl end:
    N_READS_TO_PG = get_NMAPS_FROM_SINGLE_END_READS(BAMFILE_to_pg)
    
    # The ext. vc. tool however, respected the original reads:
    if input_reads_are_paired :
        N_READS_TO_ADHOC = get_NMAPS_FROM_PAIRED_END_READS(BAMFILE_to_adhoc)
    else :
        N_READS_TO_ADHOC = get_NMAPS_FROM_SINGLE_END_READS(BAMFILE_to_adhoc)

    
    summary_filename = working_dir + "/summary.txt"
    summary_file = open(summary_filename, "w")
    summary_file.write("INPUT READS : " + str(N_INPUT_READS) + "\n")
    summary_file.write("\n")
    summary_file.write("BAM file with reads aligned to pan genome : " + BAMFILE_to_pg + "\n")
    summary_file.write("Number of reads aligned to pan genome     : " + str(N_READS_TO_PG) + "\n") 
    summary_file.write("Percentage of reads aligned to pan genome : " + str(100*float(N_READS_TO_PG)/float(N_INPUT_READS)) + "\n") 
    summary_file.write("\n")
    summary_file.write("BAM file with reads aligned to ad-hoc ref : " + BAMFILE_to_adhoc + "\n")
    summary_file.write("Number of reads aligned to ad-hoc ref     : " + str(N_READS_TO_ADHOC) + "\n") 
    summary_file.write("Percentage of reads aligned to ad-hoc ref : " + str(100*float(N_READS_TO_ADHOC)/float(N_INPUT_READS)) + "\n") 
    summary_file.write("\n")
    
    ## for each chromosome:
    chr_file_list = pgindex_dir + "/chr_list.txt"
    for line in open(chr_file_list):
        chr_id = line.rstrip('\n')
        adhoc_ref_aligned_to_ref = working_dir + "/adhoc_ref_files/" + chr_id + "/adhoc_reference.aligned_to_ref"
        ref_file = pgindex_dir + "/" + chr_id + "/recombinant.n1.gapped"
        adhocref_summary = SummarizeAdhocRef(adhoc_ref_aligned_to_ref, ref_file)
        summary_file.write("Ad-hoc ref summary for chromosome "+ chr_id + "\n")
        summary_file.write(adhocref_summary)
        summary_file.write("Original reference file: " +  ref_file)
        summary_file.write("Ad-hoc   reference file: " +  adhocref_summary)
    
    summary_file.write("\n")
    score_pgm = working_dir + "/logs/scores_matrix.pgm"
    summary_file.write("A picture of the scores matrix is available at: " + score_pgm + "\n")

    summary_file.close()




