import sys
sys.path.append(f"{snakemake.scriptdir}/../../pvc_py_tools")
from sam_to_positions import SamToPos

SamToPos(
	snakemake.params.sam_output_dir[0],
	snakemake.input.recombinant_all[0],
	snakemake.params.chromosome_list,
	snakemake.params.sensibility,
	snakemake.params.n_refs,
	snakemake.log[0]
)
