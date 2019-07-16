#!/usr/bin/python3
import argparse
from pvc_tools import *
from pathlib import Path

def align_main():
    parser = argparse.ArgumentParser(description='Align reads to a PanVC-indexed pan-genome.')
    parser.add_argument("--pgindex_dir",          required=True, help="PanVC index directory")
    parser.add_argument("-r1", "--reads_file_1",  required=True, help="reads file 1 (aka)    ")
    parser.add_argument("-r2", "--reads_file_2",  default="",    help="reads file 2 (optional)")
    parser.add_argument("-o",  "--sam_output_dir", required=True, help="Output folder to store variant calling results")

    parser.add_argument("--debug", action='store_true', help="Will run extra checks, degrading performance.")
    parser.add_argument("-t", "--n_threads",  type=int, default=1,    help="Number of threads")
    ## TODO: check if we really want this to be a parameter HERE. (it will be for the index though.)
    parser.add_argument("-d", "--max_edit_distance",  type=int, default=10,    help="Max edit distance")
    parser.add_argument("-p", "--ploidity",  type=int, default=2,    help="Ploidity of the genomes")
    args = parser.parse_args()
    PVC_align(args)

def PVC_align(args):
    pgindex_dir = args.pgindex_dir;
    reads_1 = args.reads_file_1
    reads_2 = args.reads_file_2
    output_folder = args.sam_output_dir

    paired_flag = (args.reads_file_2 != "")

    debug_mode = args.debug
    n_threads = args.n_threads
    max_edit_distance = args.max_edit_distance
    ploidity = args.ploidity

    ensure_dir(output_folder)

    reads_all = reads_1
    if (paired_flag):
        reads_all = PVC_merge_reads(reads_1, reads_2)

    n_refs = PVC_nrefs_from_index(pgindex_dir)
    read_len = PVC_read_len_from_reads(reads_all)
    read_len_file = output_folder + "/read_len.txt"
    with open(read_len_file, 'w') as f:
        f.write(str(read_len))

    chr_list_file = pgindex_dir + "/chr_list.txt"
    assert(os.path.isfile(chr_list_file))


    print ("Reference contains : " + str(n_refs) + " references")
    print ("Read len           : " + str(read_len))
    ##TODO:  do I want the old "reuse sam option?"

    sequence_all_file = pgindex_dir + "/recombinant.all.fa"

    chic_align_flags ="--secondary=LZ --threads=" + str(n_threads) + " --kernel-options=--n-ceil=C," + str(max_edit_distance) + ",0"
    alignment_command = ""
    if debug_mode:
        sam_all_plain = output_folder + "/mapped_reads_all.sam"
        alignment_command = chic_align_bin + " " + chic_align_flags + " " +sequence_all_file + " " + reads_all + " --output="+ sam_all_plain 
        samtools_view_command = samtools_bin + " view -Sb " + sam_all_plain + " > " + output_folder + "/all_mapped.bam"
        call_or_die(alignment_command)
        call_or_die(samtools_view_command)
    else:
        alignment_command = chic_align_bin + " " + chic_align_flags + " " +sequence_all_file + " " + reads_all +" | " + samtools_bin + " view -Sb -  > " + output_folder + "/all_mapped.bam"
        call_or_die(alignment_command)

    sort_command = samtools_bin + " sort -@ " + str(n_threads) + " " + output_folder+"/all_mapped.bam" + " " + output_folder + "/all_sorted"
    call_or_die(sort_command)

    split_command = bamtools_bin + " split -reference -in "+output_folder+"/all_sorted.bam"
    call_or_die(split_command)

    with open(chr_list_file) as f:
        for line in f:
            chr_id = int(line)
            print ("processing chr: " + str(chr_id))
            samples_name_file = pgindex_dir + "/" + str(chr_id) + "/names.plain"
            Path(output_folder + "/" + str(chr_id)).mkdir()

            for curr_ref in range(1,n_refs+1):
                print ("Ref id: " + str(curr_ref))
                curr_fasta_name = PVC_sequence_num_to_name(samples_name_file, chr_list_file, ploidity, chr_id, curr_ref)
                bam_file = output_folder + "/all_sorted.REF_" + curr_fasta_name  + ".bam"
                sam_file = output_folder + "/" + str(chr_id)  + "/mapped_reads_to" + str(curr_ref) + ".sam.gz"
                if not os.path.isfile(bam_file):
                    assert not os.path.exists(sam_file)
                    Path(sam_file).touch()
                else:
                    command = samtools_bin + " view -@ " + str(n_threads) + " " + bam_file + " | gzip > " + sam_file
                    call_or_die(command)

if __name__ == "__main__":
    align_main()

