# Copyright (c) Daniel Valenzuela, Tuukka Norri 2019–2020
# Licenced under the MIT licence.

import distutils.spawn as ds
import os


CURR_PATH = os.path.dirname(os.path.realpath(__file__))
PANVC_DIR = CURR_PATH + "/.."
CHIC_INDEX_BIN		= PANVC_DIR + "/components/CHIC/src/chic_index"
CHIC_ALIGN_BIN		= PANVC_DIR + "/components/CHIC/src/chic_align"
HP_BIN				= PANVC_DIR + "/components/lightweight_heaviest_paths/src/heaviest_path"
MATRIX_PRINT_BIN	= PANVC_DIR + "/components/lightweight_heaviest_paths/src/print_matrix"
SPLIT_AND_SORT_BIN	= PANVC_DIR + "/components/lightweight_heaviest_paths/src/split_and_sort"
PILEUP_BIN			= PANVC_DIR + "/components/lightweight_heaviest_paths/src/pileup"

BWA_BIN				= ds.find_executable("bwa")
SAMTOOLS_BIN		= ds.find_executable("samtools")
BCFTOOLS_BIN		= ds.find_executable("bcftools")
GATK_BIN			= ds.find_executable("gatk")
