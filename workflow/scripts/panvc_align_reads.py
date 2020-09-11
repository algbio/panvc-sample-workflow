import sys

sys.path.append(f"{snakemake.scriptdir}/../../pvc_py_tools")
from align_reads import run_panvc_aligner

sam_output_dir = f"{snakemake.params.output_root}/sam_files"
run_panvc_aligner(
	snakemake.input.reads_file_1,
	snakemake.input.reads_file_2,
	snakemake.params.index_root,
	snakemake.params.chromosome_list,
	snakemake.params.ploidy,
	snakemake.params.n_refs,
	snakemake.params.max_read_len,
	snakemake.params.max_edit_distance,
	snakemake.threads,
	sam_output_dir
)
