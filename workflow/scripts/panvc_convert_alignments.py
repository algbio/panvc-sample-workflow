# Copyright (c) Tuukka Norri 2020
# Licenced under the MIT licence.


import sys

sys.path.append(f"{snakemake.scriptdir}/../../pvc_py_tools")
from align_reads import convert_panvc_output


sam_output_dir = f"{snakemake.params.output_root}/sam_files"
convert_panvc_output(
	snakemake.params.chromosome_list,
	snakemake.threads,
	sam_output_dir
)
