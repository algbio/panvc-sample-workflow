#!/usr/bin/python3
import sys
PANVC_DIR = sys.path[0]
#Alternative approach:
#CURR_PATH = os.path.dirname(os.path.realpath(__file__))
#PANVC_PATH = CURR_PATH + "/.."
chic_index_bin = PANVC_DIR + "/components/pan_genome_index_real/CHIC/src/chic_index"
vcf2names_bin = PANVC_DIR + "/components/pan_genome_index_real/ext/vcflib/bin/vcfsamplenames"
vcf2fastas_bin = PANVC_DIR + "/components/pan_genome_index_real/ext/vcf2multialign"
chic_align_bin = PANVC_DIR + "/components/pan_genome_index_real/CHIC/src/chic_align"
samtools_bin = PANVC_DIR + "/ext_var_call_pipelines/ext/samtools-0.1.19/samtools"
bamtools_bin = PANVC_DIR + "/components/pan_genome_index_real/ext/bamtools/build/src/toolkit/bamtools"
