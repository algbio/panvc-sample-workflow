# Copyright (c) Tuukka Norri 2020
# Licenced under the MIT licence.

import argparse
import re
import sys
import tempfile
import yaml

parse_fname_re = re.compile(r"([^/]+)[.]a2m$")


def count_sequences(path):
	count = 0
	with open(path, "r") as f:
		for line in f:
			if line.startswith(">"):
				count += 1
	return count


def chr_names(paths):
	global parse_fname_re
	for path in paths:
		m = parse_fname_re.search(path)
		if not m:
			print(f"ERROR: Unable to extract chromosome name from “{path}”.", file = sys.stderr)
			sys.exit(1)
		
		chr_name = m.group(1)
		yield (chr_name, path)


def write_config(config_output_path, pgindex_dir, input_a2m_paths, max_edit_distance, max_read_length, max_memory_MB, a2m_path_prefix = "", benchmark_dir = "benchmark-index", tempdir = tempfile.gettempdir()):
	with open(config_output_path, "x") as f:
		config = {}
		
		config["index_root"] = pgindex_dir
		config["input_a2m"] = [(chr_name, f"{a2m_path_prefix}/{path}" if 0 < len(a2m_path_prefix) else path) for chr_name, path in chr_names(input_a2m_paths)]
		config["n_refs"] = count_sequences(input_a2m_paths[0])
		
		config["max_edit_distance"] = max_edit_distance
		config["max_read_length"] = max_read_length
		
		config["max_memory_MB"] = max_memory_MB
		config["benchmark_dir"] = benchmark_dir
		config["tempdir"] = tempdir
		
		yaml.dump(config, f)


def main():
	parser = argparse.ArgumentParser(description = "Align reads to a PanVC-indexed pan-genome.")
	parser.add_argument("--input-a2m", nargs = "+", help = "Input A2M files. File name without suffix will be used as the chromosome name.")
	parser.add_argument("--pgindex-dir", required = True, help = "PanVC index directory")
	parser.add_argument("-c", "--output-config", type = str, default = "panvc-config-index.yaml", help = "Configuration file output path")
	parser.add_argument("-d", "--max-edit-distance", type = int, default = 10, help = "Maximum edit distance allowed by the index")
	parser.add_argument("-l", "--max-read-length", type = int, default = 120, help = "Maximum read length that will be possible to align using the index")
	parser.add_argument("-m", "--max-memory-MB", type = int, default = 8192, help = "Amount of RAM to use with most of the tools (MB)")

	args = parser.parse_args()
	write_config(args.output_config, args.pgindex_dir, args.input_a2m, args.max_edit_distance, args.max_read_length, args.max_memory_MB)

if __name__ == "__main__":
	main()
