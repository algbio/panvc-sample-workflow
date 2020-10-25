# Copyright (c) Tuukka Norri 2020
# Licenced under the MIT licence.

import argparse
import os
import sys
import tempfile
import yaml

CURRENT_PATH = os.path.dirname(__file__) or "."
sys.path.append(f"{CURRENT_PATH}/pvc_py_tools")
from pvc_tools import PVC_load_var


def write_config(
	output_config,
	index_root,
	output_root,
	reads_file_1,
	reads_file_2,
	vc_methods,
	ploidy,
	ploidy_file,
	max_memory_MB,
	should_store_merged_reads_with_originals,
	should_run_baseline_vc,
	benchmark_dir = "benchmark-call",
	tempdir = tempfile.gettempdir()
):
	with open(output_config, "x") as f:
		config = {}
		
		config["index_root"] = index_root
		config["output_root"] = output_root
		
		config["reads_file_1"] = reads_file_1
		config["reads_file_2"] = reads_file_2
		
		workflow = ["pg"]
		if should_run_baseline_vc:
			workflow.append("baseline")
		config["workflow"] = workflow
		
		variant_caller = []
		if "GATK" in vc_methods:
			variant_caller.append("gatk")
		if "SAMTOOLS" in vc_methods:
			variant_caller.append("samtools")
		config["variant_caller"] = variant_caller
		
		chromosome_list = []
		for line in open(f"{index_root}/chr_list.txt"):
			chrom = line.rstrip("\n")
			chromosome_list.append(chrom)
		config["chromosome_list"] = chromosome_list
		
		config["ploidy"] = ploidy
		config["ploidy_file"] = ploidy_file
		
		config["sensibility"] = 5
		
		config["n_refs"] = int(PVC_load_var("n_refs", f"{index_root}/{chromosome_list[0]}")) # FIXME: Allow varying the number of reference sequences per chromosome.
		config["max_edit_distance"] = int(PVC_load_var("max_edit_distance", index_root))
		config["max_read_len"] = int(PVC_load_var("max_read_len", index_root))
		
		config["max_memory_MB"] = max_memory_MB
		
		reads_1_dir  = os.path.dirname(reads_file_1)
		reads_2_dir  = os.path.dirname(reads_file_2)
		reads_1_base = os.path.basename(reads_file_1)
		reads_2_base = os.path.basename(reads_file_2)
		reads_common_prefix = os.path.commonprefix([reads_1_base, reads_2_base])
		if should_store_merged_reads_with_originals:
			config["reads_all_path"] = f"{reads_1_dir}/{reads_common_prefix}_ALL_RENAMED.fq.gz"
		else:
			config["reads_all_path"] = f"{output_root}/reads_ALL_RENAMED.fq.gz"

		if os.path.exists(config["reads_all_path"]):
			print(f"NOTE: {config['reads_all_path']} already exists.", file = sys.stderr)
		
		config["benchmark_dir"] = benchmark_dir
		config["tempdir"] = tempdir
		
		yaml.dump(config, f)


def main():
	parser = argparse.ArgumentParser(description = "Align reads to a PanVC-indexed pan-genome.")
	parser.add_argument("--pgindex-dir", required = True, help = "PanVC index directory")
	parser.add_argument("--benchmark-dir", default = "benchmark-call", help = "Variant calling benchmark directory")
	parser.add_argument("-r1", "--reads-file-1", required = True, help = "reads file 1")
	parser.add_argument("-r2", "--reads-file-2", required = True, help = "reads file 2")
	parser.add_argument("-o", "--output-dir", required = True, help = "Output directory for variant calling results")
	parser.add_argument("--vc-method", nargs = "+", choices = ["GATK", "SAMTOOLS"], default = ["GATK"], help = "Methods to use for variant calling.")
	parser.add_argument("--baseline-vc", action = "store_true", help = "Run the variant caller on the standard reference as a baseline.")
	parser.add_argument("-p", "--ploidy", type = int, default = 2, help = "Ploidy of the genomes")
	parser.add_argument("-f", "--ploidy-file", type = str, default = "GRCh38", help = "Ploidy for bcftools when using Samtools")
	parser.add_argument("--store-merged-reads-with-originals", action = "store_true", default = False, help = "Store merged reads file in the same directory as the original files.")
	parser.add_argument("-m", "--max-memory-MB", type = int, default = 8192, help = "Amount of RAM to use with most of the tools (MB)")
	parser.add_argument("-c", "--output-config", type = str, default = "panvc-config-align.yaml", help = "Configuration file output path")
	
	args = parser.parse_args()
	write_config(
		args.output_config,
		args.pgindex_dir,
		args.output_dir,
		args.reads_file_1,
		args.reads_file_2,
		args.vc_method,
		args.ploidy,
		args.ploidy_file,
		args.max_memory_MB,
		args.store_merged_reads_with_originals,
		args.baseline_vc,
		benchmark_dir = args.benchmark_dir
	)
	
if __name__ == "__main__":
	main()
