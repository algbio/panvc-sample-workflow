#!/usr/bin/python3

# Copyright (c) Daniel Valenzuela, Tuukka Norri 2019–2020
# Licenced under the MIT licence.

import argparse
import distutils.spawn as ds
import glob
import os
import re
from pvc_tools import PVC_save_var, PVC_read_len_from_reads, call_or_die
from pathlib import Path


def run_panvc_aligner(reads_all, pgindex_dir, chr_list, ploidy, n_refs, max_read_len, max_edit_distance, n_threads, output_folder):
    samtools_path = ds.find_executable("samtools")

    read_len = PVC_read_len_from_reads(reads_all)
    assert(read_len <= max_read_len)
    
    PVC_save_var(read_len, "read_len", output_folder)

    print ("Reference contains : " + str(n_refs) + " references")
    print ("Read len           : " + str(read_len))
    ##TODO(Future):  We used to have a "reuse sam" flag to avoid reruning the whole pipeline. 
    ## Now this option is gone. The real fix is not to bring it back, but to migrate to snake-make so that we can re-run the pipeline from any intermediate 
    ## stage if some changes were done.

    align_cmd = f"chic-align --secondary-report=NONE --threads={n_threads} --kernel-options=--n-ceil=C,{max_edit_distance},0 --max-ed={max_edit_distance} --split-output-by-reference --output-prefix={output_folder} --samtools-path={samtools_path} --verbose=3 {pgindex_dir}/recombinant.all.fa {reads_all}"
    call_or_die(align_cmd)

    # FIXME: add memory limit to sort.
    # FIXME: glob is not the best approach. Use PVC_sequence_num_to_name instead? (or replace with sample indices)
    unsorted_bams = glob.glob(f"{output_folder}/all_mapped.REF_*.bam")
    print("Got unsorted bams:", unsorted_bams)
    for unsorted_bam in unsorted_bams:
        rname = unsorted_bam[len(f"{output_folder}/all_mapped.REF_"):]
        rname = rname[:-(len(".bam"))]
        sort_command = f"samtools sort -@ {n_threads} --output-fmt BAM -o {output_folder}/all_sorted.REF_{rname}.bam {unsorted_bam}"
        call_or_die(sort_command)
        os.remove(unsorted_bam)

    for chr_id in chr_list:
        samples_name_file = pgindex_dir + "/" + chr_id + "/names.plain"
        #Path(output_folder + "/" + chr_id).mkdir()

        for curr_ref in range(1,n_refs+1):
            curr_fasta_name = str(curr_ref) #PVC_sequence_num_to_name(samples_name_file, len(chr_list), ploidy, chr_id, curr_ref)
            bam_file = output_folder + "/all_sorted.REF_" + curr_fasta_name  + ".bam"
            sam_file = output_folder + "/" + chr_id  + "/mapped_reads_to" + str(curr_ref) + ".sam.gz"
            if not os.path.isfile(bam_file):
                assert not os.path.exists(sam_file)
                Path(sam_file).touch()
            else:
                command = "samtools view -@ " + str(n_threads) + " " + bam_file + " | gzip > " + sam_file
                call_or_die(command)
