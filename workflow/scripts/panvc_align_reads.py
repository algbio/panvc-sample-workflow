# Copyright (c) Tuukka Norri 2020-2021
# Licenced under the MIT licence.


import sys
import distutils.spawn as ds

sys.path.append(f"{snakemake.params.basedir}/pvc_py_tools")
from pvc_tools import PVC_read_len_from_reads

read_len = PVC_read_len_from_reads(snakemake.input.reads_all)
assert(read_len <= snakemake.params.max_read_len)

samtools_path = ds.find_executable("samtools")

import snakemake.shell as shell # If imported earlier, the module masks the global variable “snakemake”.
shell((
	"mkdir -p {snakemake.output} &&"
	" chic-align"
	" --secondary-report=NONE"
	" --threads={snakemake.threads}"
	" --kernel-options=--n-ceil=C,{snakemake.params.max_edit_distance},0"
	" --max-ed={snakemake.params.max_edit_distance}"
	" --split-output-by-reference"
	" --output-prefix={snakemake.output}"
	" --samtools-path={samtools_path}"
	" --no-samtools-threads"
	" --verbose=3"
	" {snakemake.input.recombinant_all} {snakemake.input.reads_all}"
))
