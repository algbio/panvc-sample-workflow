import distutils.spawn as ds
import os


CURR_PATH = os.path.dirname(os.path.realpath(__file__))
PANVC_DIR = CURR_PATH + "/.."
chic_index_bin		= PANVC_DIR + "/components/pan_genome_index_real/CHIC/src/chic_index"
vcf2names_bin		= PANVC_DIR + "/components/pan_genome_index_real/ext/vcflib/bin/vcfsamplenames"
vcfbreakmulti_bin	= PANVC_DIR + "/components/pan_genome_index_real/ext/vcflib/bin/vcfbreakmulti"
vcfcombine_bin		= PANVC_DIR + "/components/pan_genome_index_real/ext/vcflib/bin/vcfcombine"
vcfcreatemulti_bin	= PANVC_DIR + "/components/pan_genome_index_real/ext/vcflib/bin/vcfcreatemulti"
vcf2fastas_bin		= PANVC_DIR + "/components/pan_genome_index_real/ext/vcf2multialign"
CHIC_ALIGN_BIN		= PANVC_DIR + "/components/pan_genome_index_real/CHIC/src/chic_align"
VCFCHECK_BIN		= PANVC_DIR + "/components/pan_genome_index_real/ext/vcflib/bin/vcfcheck"
HP_BIN				= PANVC_DIR +  "/components/lightweight_heaviest_paths/src/heaviest_path"
MATRIX_PRINT_BIN	= PANVC_DIR +  "/components/lightweight_heaviest_paths/src/print_matrix"
SPLIT_AND_SORT_BIN	= PANVC_DIR + "/components/lightweight_heaviest_paths/src/split_and_sort"
PILEUP_BIN			= PANVC_DIR + "/components/lightweight_heaviest_paths/src/pileup"

BWA_BIN				= ds.find_executable("bwa")
SAMTOOLS_BIN		= ds.find_executable("samtools")
BCFTOOLS_BIN		= ds.find_executable("bcftools")
GATK_BIN			= ds.find_executable("gatk")
