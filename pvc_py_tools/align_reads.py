#!/usr/bin/python3

# Copyright (c) Daniel Valenzuela, Tuukka Norri 2019–2020
# Licenced under the MIT licence.

import argparse
import distutils.spawn as ds
import glob
import os
import re
import sys
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

    align_cmd = f"chic-align --secondary-report=NONE --threads={n_threads} --kernel-options=--n-ceil=C,{max_edit_distance},0 --max-ed={max_edit_distance} --split-output-by-reference --output-prefix={output_folder} --samtools-path={samtools_path} --no-samtools-threads --verbose=3 {pgindex_dir}/recombinant.all.fa {reads_all}"
    call_or_die(align_cmd)

    # FIXME: create a new workflow step for this.
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

    # FIXME: CHIC outputs a BAM file for only each reference to which reads were aligned, but Snakemake expects a predetermined list of outputs. The following file is a marker but creating the remaining files would be an option.
    Path(f"{output_folder}/chic_aligner_did_finish").touch()


def convert_panvc_output(chr_list, n_refs, n_threads, output_folder):
    for chr_id in chr_list:
        Path(f"{output_folder}/{chr_id}").mkdir(parents = True, exist_ok = True)
        for curr_ref in range(1, n_refs + 1):
            input_bam = f"{output_folder}/all_sorted.REF_pg_ref_{chr_id}_{curr_ref}.bam"
            output_sam = f"{output_folder}/{chr_id}/mapped_reads_to{curr_ref}.sam.gz"
            if not os.path.isfile(input_bam):
                sys.stderr.write(f"{input_bam} does not exist; creating an empty file for {chr_id}:{curr_ref}.\n")
                assert not os.path.exists(output_sam)
                # FIXME: this does not generate a valid gzip archive.
                Path(output_sam).touch()
            else:
                command = f"samtools view -@ {n_threads} {input_bam} | gzip > {output_sam}"
                call_or_die(command)
