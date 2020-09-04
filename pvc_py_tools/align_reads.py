#!/usr/bin/python3
import argparse
import glob
import re
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
    parser.add_argument("-p", "--ploidy",  type=int, default=2,    help="Ploidy of the genomes")
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
    max_edit_distance = PVC_load_var("max_edit_distance", pgindex_dir)
    ploidy = args.ploidy

    ensure_dir(output_folder)

    reads_all = reads_1
    if (paired_flag):
        reads_all = PVC_merge_reads(reads_1, reads_2, output_folder)

    n_refs = PVC_load_var("n_refs", f"{pgindex_dir}/1") # FIXME: Allow varying number of references per chromosome.
    max_read_len = PVC_load_var("max_read_len", pgindex_dir)
    read_len = PVC_read_len_from_reads(reads_all)
    assert(read_len <= max_read_len)
    
    PVC_save_var(read_len, "read_len", output_folder)
    chr_list = PVC_get_chr_list(pgindex_dir)

    print ("Reference contains : " + str(n_refs) + " references")
    print ("Read len           : " + str(read_len))
    ##TODO(Future):  We used to have a "reuse sam" flag to avoid reruning the whole pipeline. 
    ## Now this option is gone. The real fix is not to bring it back, but to migrate to snake-make so that we can re-run the pipeline from any intermediate 
    ## stage if some changes were done.

    align_cmd = f"{chic_align_bin} --secondary-report=NONE --threads={n_threads} --kernel-options=--n-ceil=C,{max_edit_distance},0 --max-ed={max_edit_distance} --split-output-by-reference --output-prefix={output_folder} --samtools-path={SAMTOOLS_BIN} {pgindex_dir}/recombinant.all.fa {reads_all}"
    call_or_die(align_cmd)

    # FIXME: add memory limit to sort.
    # FIXME: glob is not the best approach. Use PVC_sequence_num_to_name instead?
    unsorted_bams = glob.glob(f"{output_folder}/all_mapped.REF_*.bam")
    print("Got unsorted bams:", unsorted_bams)
    for unsorted_bam in unsorted_bams:
        rname = unsorted_bam[len(f"{output_folder}/all_mapped.REF_"):]
        rname = rname[:-(len(".bam"))]
        sort_command = f"{SAMTOOLS_BIN} sort -@ {n_threads} --output-fmt BAM -o {output_folder}/all_sorted.REF_{rname}.bam {unsorted_bam}"
        call_or_die(sort_command)
        os.remove(unsorted_bam)

    for chr_id in chr_list:
        samples_name_file = pgindex_dir + "/" + chr_id + "/names.plain"
        Path(output_folder + "/" + chr_id).mkdir()

        for curr_ref in range(1,n_refs+1):
            curr_fasta_name = PVC_sequence_num_to_name(samples_name_file, len(chr_list), ploidy, chr_id, curr_ref)
            bam_file = output_folder + "/all_sorted.REF_" + curr_fasta_name  + ".bam"
            sam_file = output_folder + "/" + chr_id  + "/mapped_reads_to" + str(curr_ref) + ".sam.gz"
            if not os.path.isfile(bam_file):
                assert not os.path.exists(sam_file)
                Path(sam_file).touch()
            else:
                command = SAMTOOLS_BIN + " view -@ " + str(n_threads) + " " + bam_file + " | gzip > " + sam_file
                call_or_die(command)

if __name__ == "__main__":
    align_main()

