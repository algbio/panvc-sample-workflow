from Bio import SeqIO

assert 1 == len(snakemake.output)
with open(snakemake.output[0], "x") as output_file:
	for chr_id, path in snakemake.config["input_a2m"]:
		output_file.write(">")
		output_file.write(chr_id)
		output_file.write("\n")
		
		for record in SeqIO.parse(path, "fasta"):
			for c in record.seq:	# FIXME: likely inefficient.
				if '-' != c:
					output_file.write(c)
			break
		output_file.write("\n")
