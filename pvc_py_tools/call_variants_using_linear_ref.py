#!/usr/bin/python3
from config import *
from pvc_tools import *
import argparse
import shutil
from pathlib import Path
import os
import tempfile

def VC_main():
    parser = argparse.ArgumentParser(description='Align reads to a PanVC-indexed pan-genome.')
    parser.add_argument("reference"             , help="Reference genome.")
    parser.add_argument("output_file"           , help="Name of the desired output file (vcf).")
    parser.add_argument("working_dir"           , help="Working directroy, where intermediate files live.")
    parser.add_argument("-r1", "--reads_file_1",  required=True, help="reads file 1           ")
    parser.add_argument("-r2", "--reads_file_2",  default="",    help="reads file 2 (optional)")
    parser.add_argument("-p", "--ploidy",  type=int, default=2,    help="Ploidy of the genomes")
    parser.add_argument("-f", "--ploidy-file", type = str, default = "GRCh38", help = "Ploidy for bcftools when using SAMTOOLS for alignment")

    parser.add_argument("--debug", action='store_true', help="Will run extra checks, degrading performance.")
    parser.add_argument("-t", "--n_threads",      type=int, default=1,    help="Number of threads")
    parser.add_argument("-m", "--max_memory_MB",  type=int, default=1000,    help="Available ram (MB)")
    parser.add_argument("--vc_base_method",  choices=["GATK", "SAMTOOLS"], default = "SAMTOOLS",    help="Method to call variants with the ad-hoc reference.")
    args = parser.parse_args()
    call_variants_with_standard_methods(args)

def call_variants_with_standard_methods(args):
    if (args.vc_base_method == "GATK"):
        BwaGATKVC(args)
    else:
        BwaSamtoolsVC(args)


def BwaSamtoolsVC(args):
    working_dir = args.working_dir
    reference = args.reference
    reads_file_1 = args.reads_file_1
    reads_file_2 = args.reads_file_2
    n_threads = args.n_threads
    output_file = args.output_file
    ploidy = args.ploidy
    ploidy_file = args.ploidy_file
    tempdir = tempfile.gettempdir()
    
    debug_mode = args.debug
    paired_flag = False if reads_file_2 == "" else True
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    assert(Path(reference).is_file())
    assert(Path(reads_file_1).is_file())
    if paired_flag:
        assert(Path(reads_file_2).is_file())

    print ("Bwa path: " + BWA_BIN)
    assert(Path(BWA_BIN).is_file())

    PVC_index_ref(reference)
    bwa_align_comand = f"{BWA_BIN} mem -t {n_threads} {reference} {reads_file_1} {reads_file_2} | {SAMTOOLS_BIN} view -b -@ {n_threads} -o {working_dir}/aligned_reads.bam"
    call_or_die(bwa_align_command)
    
    markdup_comm = GATK_BIN + f" --java-options '-Djava.io.tmpdir={tempdir} -Dsnappy.disable=true -Dsamjdk.snappy.disable=true'" \
                            + f" MarkDuplicates" \
                            + f" --INPUT {working_dir}/aligned_reads.bam"\
                            + f" --OUTPUT {working_dir}/aligned_deduplicated.bam"\
                            + f" --TMP_DIR={tempdir}"\
                            + f" --METRICS_FILE {working_dir}/metrics.txt"\
                            + f" --VALIDATION_STRINGENCY SILENT"\
                            + f" --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500"\
                            + f" --ASSUME_SORT_ORDER queryname"\
                            + f" --CREATE_MD5_FILE true" \
                            + f" --USE_JDK_DEFLATER true" \
                            + f" --USE_JDK_INFLATER true"
    call_or_die(markdup_comm)
    

    memory_per_thread = "%.0f" % (1024.0 * 1024.0 * args.max_memory_MB / n_threads)
    samtools_sort_command = f"{SAMTOOLS_BIN} sort -@ {n_threads} -m {memory_per_thread} --output-fmt BAM -o {working_dir}/sorted-alns2.bam {working_dir}/aligned_deduplicated.bam"
    call_or_die(samtools_sort_command)

    #TODO: rewrite without pipes. A system call with pipes can be dangerous.
    # FIXME Apparentlly ploidy should be specified for bcftools mpileup. However, there does not seem to be an option to do that and according to the manual, --samples-file’s second column is only handled by bcftools call.
    pileup_command = f"{BCFTOOLS_BIN} mpileup --output-type u --fasta-ref {reference} {working_dir}/sorted-alns2.bam | {BCFTOOLS_BIN} call --ploidy {ploidy_file} --output-type u --variants-only --multiallelic-caller > {working_dir}/var.raw.bcf"
    call_or_die(pileup_command)
    assert(Path(working_dir + "/var.raw.bcf").is_file())

    filter_command = BCFTOOLS_BIN + " view " + working_dir + "/var.raw.bcf | " + VCFUTILS_BIN + " varFilter -D100 > " + working_dir + "/var.flt.vcf"
    call_or_die(filter_command)

    shutil.copy(working_dir + "/var.flt.vcf", output_file)

def BwaGATKVC(args):
    max_memory_MB = str(args.max_memory_MB)
    working_dir = args.working_dir
    reference = args.reference
    reads_file_1 = args.reads_file_1
    reads_file_2 = args.reads_file_2
    n_threads = args.n_threads
    output_file = args.output_file
    ploidy = args.ploidy
    tempdir = tempfile.gettempdir()
    
    debug_mode = args.debug
    paired_flag = False if reads_file_2 == "" else True
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    assert(Path(reference).is_file())
    assert(Path(reads_file_1).is_file())
    if paired_flag:
        assert(Path(reads_file_2).is_file())

    assert(Path(BWA_BIN).is_file())
    
    print ("############################################")
    print (" Read alignment")
    print ("############################################")
    PVC_index_ref(reference)

    bwa_align_command = f"{BWA_BIN} mem -K 100000000 -v 3 -t {n_threads} -Y {reference} {reads_file_1} {reads_file_2} | {GATK_BIN} --java-options '-Xmx{max_memory_MB}M -Djava.io.tmpdir={tempdir} -Dsnappy.disable=true -Dsamjdk.snappy.disable=true' SortSam --TMP_DIR={tempdir} --USE_JDK_DEFLATER true --USE_JDK_INFLATER true --SORT_ORDER queryname --INPUT /dev/stdin --OUTPUT {working_dir}/aligned_reads.bam"
    call_or_die(bwa_align_command)
    
    # Running GATK can cause a “pure virtual method called” error related to a native library.
    # According to this thread, passing --USE_JDK_DEFLATER true --USE_JDK_INFLATER true can help: https://gatkforums.broadinstitute.org/gatk/discussion/12438/gatk-haplotypecaller-crashes-randomly

    ######
    ## The "uBAM" issue of GATK
    ## -> convert original fastq reads to uBAM
    
    fq2sam_command = f"{GATK_BIN} --java-options '-Xmx{max_memory_MB}M -Djava.io.tmpdir={tempdir} -Dsnappy.disable=true -Dsamjdk.snappy.disable=true' FastqToSam -F1 {reads_file_1}"
    if (paired_flag):
        fq2sam_command += f" -F2 {reads_file_2}"
    fq2sam_command += f" --SAMPLE_NAME sample1 -O {working_dir}/input_reads_unaligned.ubam --TMP_DIR={tempdir} --USE_JDK_DEFLATER true --USE_JDK_INFLATER true"
    call_or_die(fq2sam_command)

    ## merge above uBAM file with the output of BwaMem
    merge_command = GATK_BIN \
            + f" --java-options '-Djava.io.tmpdir={tempdir} -Dsnappy.disable=true -Dsamjdk.snappy.disable=true'"\
            + f" MergeBamAlignment"\
            + f" --TMP_DIR={tempdir}"\
            + f" --VALIDATION_STRINGENCY SILENT"\
            + f" --EXPECTED_ORIENTATIONS FR"\
            + f" --ATTRIBUTES_TO_RETAIN X0"\
            + f" --ALIGNED_BAM {working_dir}/aligned_reads.bam"\
            + f" --UNMAPPED_BAM {working_dir}/input_reads_unaligned.ubam"\
            + f" --OUTPUT {working_dir}/merged_reads.bam"\
            + f" --REFERENCE_SEQUENCE {reference}"\
            + f" --PAIRED_RUN true"\
            + f" --SORT_ORDER unsorted"\
            + f" --IS_BISULFITE_SEQUENCE false"\
            + f" --ALIGNED_READS_ONLY false"\
            + f" --CLIP_ADAPTERS false"\
            + f" --MAX_RECORDS_IN_RAM 2000000"\
            + f" --ADD_MATE_CIGAR true"\
            + f" --MAX_INSERTIONS_OR_DELETIONS -1"\
            + f" --PRIMARY_ALIGNMENT_STRATEGY MostDistant"\
            + f" --UNMAPPED_READ_STRATEGY COPY_TO_TAG"\
            + f" --ALIGNER_PROPER_PAIR_FLAGS true"\
            + f" --UNMAP_CONTAMINANT_READS true"\
            + f" --USE_JDK_DEFLATER true" \
            + f" --USE_JDK_INFLATER true"

    #Probably too much, relying on defaults
    #merge_comm+=' --java-options "-Dsamjdk.compression_level=${compression_level}"'
    #merge_comm+=' --PROGRAM_GROUP_VERSION "${bwa_version}"'
    #merge_comm+=' --PROGRAM_RECORD_ID bwamem'
    #merge_comm+=' --PROGRAM_GROUP_NAME bwamem'
    #merge_comm+=' --PROGRAM_GROUP_COMMAND_LINE '${BWA_COMMAND}
    call_or_die(merge_command)
    
    print ("############################################")
    print (" Mark duplicates:")
    print ("############################################")
    mark_duplicates_command = GATK_BIN\
            + f" --java-options '-Djava.io.tmpdir={tempdir} -Dsnappy.disable=true -Dsamjdk.snappy.disable=true'"\
            + f" MarkDuplicates"\
            + f" --TMP_DIR={tempdir}"\
            + f" --INPUT {working_dir}/merged_reads.bam"\
            + f" --OUTPUT {working_dir}/aligned_deduplicated.bam "\
            + f" --METRICS_FILE {working_dir}/metrics.txt "\
            + f" --VALIDATION_STRINGENCY SILENT "\
            + f" --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 "\
            + f" --ASSUME_SORT_ORDER queryname"\
            + f" --CREATE_MD5_FILE true" \
            + f" --USE_JDK_DEFLATER true" \
            + f" --USE_JDK_INFLATER true"
    call_or_die(mark_duplicates_command)

    sortfix_command = GATK_BIN \
            + f" --java-options '-Xms4000m -Djava.io.tmpdir={tempdir} -Dsnappy.disable=true -Dsamjdk.snappy.disable=true'"\
            + f" SortSam"\
            + f" --TMP_DIR={tempdir}"\
            + f" --INPUT {working_dir}/aligned_deduplicated.bam"\
            + f" --OUTPUT {working_dir}/aligned_deduplicated_sorted.bam"\
            + f" --SORT_ORDER coordinate"\
            + f" --CREATE_INDEX false"\
            + f" --CREATE_MD5_FILE false"\
            + f" --USE_JDK_DEFLATER true" \
            + f" --USE_JDK_INFLATER true"
    call_or_die(sortfix_command)

    sortfix_command_2 = GATK_BIN \
            + f" --java-options '-Xms500m -Djava.io.tmpdir={tempdir} -Dsnappy.disable=true -Dsamjdk.snappy.disable=true'"\
            + f" SetNmAndUqTags "\
            + f" --TMP_DIR={tempdir}"\
            + f" --INPUT {working_dir}/aligned_deduplicated_sorted.bam "\
            + f" --OUTPUT {working_dir}/sortedfixed.bam "\
            + f" --CREATE_INDEX true "\
            + f" --CREATE_MD5_FILE true "\
            + f" --REFERENCE_SEQUENCE {reference}" \
            + f" --USE_JDK_DEFLATER true" \
            + f" --USE_JDK_INFLATER true"
    call_or_die(sortfix_command_2)
    
    print ("############################################")
    print (" Base recalibration:")
    print ("############################################")
    # We Skip that for now, as we dont have a knownSites. 
    # For that case, the documentation suggested call SNP variants once, and use those as the known sites.
    # And repeat until convergence ... :S
    #java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R reference.fa -I realigned_reads.bam -knownSites dbsnp.vcf -knownSites gold_indels.vcf -o recal_data.table 
    
    print ("############################################")
    print (" Call variants:")
    print ("############################################")
    vc_command = GATK_BIN \
        + f" --java-options '-Xmx{max_memory_MB}M -Djava.io.tmpdir={tempdir} -Dsnappy.disable=true -Dsamjdk.snappy.disable=true'" \
        + f" HaplotypeCaller" \
        + f" -ploidy {ploidy}" \
        + f" -R {reference}" \
        + f" -I {working_dir}/sortedfixed.bam" \
        + f" -O {working_dir}/raw_variants.vcf" \
        + f" -pairHMM LOGLESS_CACHING" \
        + f" --use-jdk-deflater" \
        + f" --use-jdk-inflater"
    call_or_die(vc_command)
    ## -L ${interval_list} \ 
    ## -contamination ${default=0 contamination} ${true="-ERC GVCF" false="" make_gvcf}              
    
    shutil.copy(f"{working_dir}/raw_variants.vcf", output_file)
    
    #TODO: replace with a native python function. A system call with pipes can be dangerous.
    if (paired_flag):
        histogram_command = SAMTOOLS_BIN + " view -F 0x4 " + working_dir + "/sortedfixed.bam  | cut -f 5 | sort | uniq -c | awk '{print $2 \" \" $1}' > " + working_dir + "/histogram.txt"
        call_or_die(histogram_command)

if __name__ == "__main__":
    VC_main()

