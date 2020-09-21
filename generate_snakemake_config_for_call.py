# Copyright (c) Tuukka Norri 2020
# Licenced under the MIT licence.

import argparse
import sys
import tempfile
import yaml

CURRENT_PATH = dirname(__file__)
sys.path.append(f"{CURRENT_PATH}/pvc_py_tools")
from pvc_tools import PVC_load_var


def main():
	parser = argparse.ArgumentParser(description = "Align reads to a PanVC-indexed pan-genome.")
	parser.add_argument("--pgindex-dir", required = True, help = "PanVC index directory")
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
	
	with open(args.output_config, "x") as f:
		config = {}
		
		pgindex_dir = args.pgindex_dir
		config["index_root"] = pgindex_dir
		config["output_root"] = args.output_dir
		
		config["reads_file_1"] = args.reads_file_1
		config["reads_file_2"] = args.reads_file_2
		
		workflow = ["pg"]
		if args.baseline_vc:
			workflow.append("baseline")
		config["workflow"] = workflow
		
		variant_caller = []
		if "GATK" in args.vc_method:
			variant_caller.append("gatk")
		if "SAMTOOLS" in args.vc_method:
			variant_caller.append("samtools")
		config["variant_caller"] = variant_caller
		
		chromosome_list = []
		for line in open(f"{pgindex_dir}/chr_list.txt"):
			chrom = line.rstrip("\n")
			chromosome_list.append(chrom)
		config["chromosome_list"] = chromosome_list
		
		config["ploidy"] = args.ploidy
		config["ploidy_file"] = args.ploidy_file
		
		config["sensibility"] = 5
		
		config["n_refs"] = int(PVC_load_var("n_refs", f"{pgindex_dir}/1")) # FIXME: Allow varying the number of reference sequences per chromosome.
		config["max_edit_distance"] = int(PVC_load_var("max_edit_distance", pgindex_dir))
		config["max_read_len"] = int(PVC_load_var("max_read_len", pgindex_dir))
		
		config["max_memory_MB"] = args.max_memory_MB
		
		if args.store_merged_reads_with_originals:
			base = os.path.dirname(args.reads_file_1)
			config["reads_all_path"] = f"{base}/reads_ALL_RENAMED.fq.gz"
		else:
			config["reads_all_path"] = f"{pgindex_dir}/reads_ALL_RENAMED.fq.gz"
		
		tempdir = tempfile.gettempdir()
		config["tempdir"] = tempdir
		
		yaml.dump(config, f)

if __name__ == "__main__":
	main()
