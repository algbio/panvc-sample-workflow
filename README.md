PanVC Sample Workflow
=====================

This is a workflow for variant calling using short reads that utilizes [PanVC](https://gitlab.com/dvalenzu/PanVC/-/tree/PanVC-2.0-rc-tsnorri). The workflow uses an index generated from a multiple sequence alignment. The sequences may be e.g. founder sequences generated from the variants of a group of donors. Using the predicted DNA sequences of the donors is also possible.


Requirements
------------
 * Linux on x86-64 for using pre-built binaries
 * Python 3
 * [Snakemake](https://snakemake.readthedocs.io/)
 * [Conda](https://conda.io/)


Installing
----------
Snakemake is used to run the workflow. To simplify the installation of the tools that are used as part of the workflow, Conda environments (in `workflow/envs`) are used. To set up the environment, Snakemake can run Conda automatically.

One tool used by PanVC, `combine_msa_vcf` is currently not available via Conda. To install it, Snakemake will download a distribution as part of variant calling and extract it to a subdirectory of this repository.


### Installing on platforms other than Linux on x86-64

Currently, pre-built binaries of PanVC and `combine_msa_vcf` are only available for Linux on x86-64. To run the workflow on other architectures, tools from the following repositories need to be installed:

 * [PanVC](https://gitlab.com/dvalenzu/PanVC/-/tree/PanVC-2.0-rc-tsnorri)
 * `combine_msa_vcf` from [vcf2multialign](https://github.com/tsnorri/vcf2multialign)


Running
-------

Running the workflow involves the following steps:

 1. Install the required tools with Conda
 2. Generate indexing input for PanVC
 3. Prepare an index
 4. Run the short read alignment and variant calling workflow

### Installing the required tools with Conda

The required tools will be installed automatically when running Snakemake for either indexing or variant calling. However, to e.g. run multiple instances of Snakemake in parallel, the tools may be installed to e.g. `/path/to/conda/environment` as follows:

 1. Run `snakemake --snakefile Snakefile.install --cores 1 --printshellcmds --use-conda --conda-prefix /path/to/conda/environment conda_environment`
 2. Run `snakemake --snakefile Snakefile.install --cores 1 --printshellcmds --use-conda --conda-prefix /path/to/conda/environment conda_environment_gatk`
 3. Run `snakemake --snakefile Snakefile.install --cores 1 --printshellcmds get_combine_msa_vcf`

### Generating indexing input

Indexing input for PanVC consists of multiple alignments of sequences for one or more chromosomes. For each chromosome, the same number of input sequences must be used. Either predicted sequences of a cohort of donors or founder sequences may be used.

PanVC accepts (a subset of) A2M formatted multiple alignments as its inputs. In particular, the first sequence needs to be the reference sequences and only uppercase characters and dashes may be used in the sequences.

Suitable multiple alignments of either predicted or founder sequences may be generated from VCF files with [vcf2multialign](https://github.com/tsnorri/vcf2multialign).

### Preparing an index

An index may be generated from a set of A2M files as follows. The amount of memory required depends on the length and complexity of the inputs.

 1. Use `generate_snakemake_config_for_index.py` to generate a configuration file for Snakemake. Here, the maximum read length and maximum edit distance for each read may be specified. Please use `python3 generate_snakemake_config_for_index.py --help` for listing all available options. For instance, to prepare an index for two chromosomes, `chr1` and `chr2`, with maximum read length of 105 and maximum edit distance of 10, the following command could be used: `python3 generate_snakemake_config_for_index.py --input-a2m chr1.a2m chr2.a2m --pgindex-dir index-output --output-config config-index.yaml --max-edit-distance 10 --max-read-length 105 --max-memory-MB 20000`. (The script may also be loaded as a module to be used from Python.)
 2. Run Snakemake with the configuration file with e.g. `snakemake --configfile config-index.yaml --snakefile Snakefile.index --cores 16 --printshellcmds --use-conda --conda-prefix /path/to/conda/environment --resources mem_mb=20000`

### Running the short read alignment and variant calling workflow

The generated index may be used to align a set of short reads and call variants as follows. As in the previous step, the amount of memory required depends on the length and complexity of the inputs.

 1. Use `generate_snakemake_config_for_call.py` to generate a configuration file for Snakemake. Please use `python3 generate_snakemake_config_for_call.py --help` for listing all available options. For instance, to call variants from a human sample using GRCh38 reference with both GATK and Samtools (for e.g. comparison), the following command could be used: `python3 generate_snakemake_config_for_call.py --pgindex-dir index-output --output-config config-call.yaml -r1 reads1.fq.gz -r2 reads2.fq.gz --output-dir call-output --vc-method GATK SAMTOOLS --ploidy 2 --ploidy-file GRCh38 --max-memory-MB 200000`.
 2. (Optional) Generate merged reads if the same reads are aligned multiple times. This is done as part of the variant calling process but can also be done beforehand.
    1. Add `--store-merged-reads-with-originals` to the parameters of `generate_snakemake_config_for_call.py`. This causes the merged reads to be stored in the same directory as the original reads instead of using the output directory.
    2. Run Snakemake with the generated configuration and use `panvc_generate_reads_all` as the target name, e.g. `snakemake --configfile config-call.yaml --cores 1 --printshellcmds --use-conda --conda-prefix /path/to/conda/environment -- panvc_generate_reads_all`. Using more than one core may speed up the process in case the reads are compressed.
 3. Run Snakemake with the configuration file with e.g. `snakemake --configfile config-call.yaml --cores 16 --printshellcmds
 --use-conda --conda-prefix /path/to/conda/environment --resources mem_mb=200000`
