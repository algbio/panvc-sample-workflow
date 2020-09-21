# Copyright (c) Tuukka Norri 2020
# Licenced under the MIT licence.

from Bio import SeqIO
import pathlib


def split_a2m(input_path, output_prefix):
	global snakemake
	
	record_ids = []
	
	gapped_len = 0
	seq_count = 0
	
	for idx, record in enumerate(SeqIO.parse(input_path, "fasta"), 1):
		record_ids.append(record.id)
		
		with open(f"{output_prefix}/recombinant.n{idx}.gapped", "x") as f:
			for c in record.seq: # FIXME: this is likely inefficient.
				f.write(c)
			f.write("\n")
			
		gapped_len = len(record.seq)
		seq_count = idx
		
	assert 1 == len(snakemake.output.msa_len)
	with open(snakemake.output.msa_len[0], "x") as f:
		f.write(str(gapped_len))
	
	assert snakemake.config["n_refs"] == seq_count, "All chromosomes need to have the same number of reference sequences"
	assert 1 == len(snakemake.output.n_refs)
	with open(snakemake.output.n_refs[0], "x") as f:
		f.write(str(seq_count))
	
	return record_ids


dst_dir = f"{snakemake.config['index_root']}/{snakemake.wildcards.chr_id}"

# Create the directory.
pathlib.Path(dst_dir).mkdir(parents = True, exist_ok = True)

# Find the processed file.
input_path = None
for chr_id, path in snakemake.config["input_a2m"]:
	if chr_id == snakemake.wildcards.chr_id:
		input_path = path
		break
assert input_path is not None

# Gapped recombinants
record_ids = split_a2m(input_path, dst_dir)

# names.plain
assert 1 == len(snakemake.output.names_plain)
with open(snakemake.output.names_plain[0], "x") as f:
	for rec_id in record_ids:
		f.write(rec_id)
		f.write("\n")
