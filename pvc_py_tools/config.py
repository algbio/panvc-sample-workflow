#!/usr/bin/python3
#import sys
import os
#PANVC_DIR = sys.path[0]
#Alternative approach:
CURR_PATH = os.path.dirname(os.path.realpath(__file__))
PANVC_DIR = CURR_PATH + "/.."
chic_index_bin = PANVC_DIR + "/components/pan_genome_index_real/CHIC/src/chic_index"
vcf2names_bin = PANVC_DIR + "/components/pan_genome_index_real/ext/vcflib/bin/vcfsamplenames"
vcf2fastas_bin = PANVC_DIR + "/components/pan_genome_index_real/ext/vcf2multialign"
chic_align_bin = PANVC_DIR + "/components/pan_genome_index_real/CHIC/src/chic_align"
samtools_bin = PANVC_DIR + "/ext_var_call_pipelines/ext/samtools-0.1.19/samtools"
bamtools_bin = PANVC_DIR + "/components/pan_genome_index_real/ext/bamtools/build/src/toolkit/bamtools"
VCFCHECK_BIN = PANVC_DIR + "/components/pan_genome_index_real/ext/vcflib/bin/vcfcheck"
HP_BIN = PANVC_DIR +  "/components/lightweight_heaviest_paths/src/heaviest_path"
MATRIX_PRINT_BIN = PANVC_DIR +  "/components/lightweight_heaviest_paths/src/print_matrix"
SPLIT_AND_SORT_BIN = PANVC_DIR + "/components/lightweight_heaviest_paths/src/split_and_sort"
PILEUP_BIN = PANVC_DIR + "/components/lightweight_heaviest_paths/src/pileup"
## VC    
BWA_BIN = PANVC_DIR + "/ext_var_call_pipelines/ext/bwa-0.7.17/bwa"
SAMTOOLS_BIN = PANVC_DIR + "/ext_var_call_pipelines/ext/samtools-0.1.19/samtools"
BCFTOOLS_BIN = PANVC_DIR + "/ext_var_call_pipelines/ext/samtools-0.1.19/bcftools/bcftools"
VCFUTILS_BIN = PANVC_DIR + "/ext_var_call_pipelines/ext/samtools-0.1.19/bcftools/vcfutils.pl"
GATK_BIN = PANVC_DIR + "/ext_var_call_pipelines/ext/GATK/gatk-4.1.0.0/gatk"
    
