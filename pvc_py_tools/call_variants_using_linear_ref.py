#!/usr/bin/python3
from config import *
from pvc_tools import *
import argparse
import shutil
from pathlib import Path
import os

def VC_main():
    parser = argparse.ArgumentParser(description='Align reads to a PanVC-indexed pan-genome.')
    parser.add_argument("reference"             , help="Reference genome.")
    parser.add_argument("output_file"           , help="Name of the desired output file (vcf).")
    parser.add_argument("working_dir"           , help="Working directroy, where intermediate files live.")
    parser.add_argument("-r1", "--reads_file_1",  required=True, help="reads file 1           ")
    parser.add_argument("-r2", "--reads_file_2",  default="",    help="reads file 2 (optional)")

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
    max_mem_bytes = str(args.max_memory_MB) + "000000"
    working_dir = args.working_dir
    reference = args.reference
    reads_file_1 = args.reads_file_1
    reads_file_2 = args.reads_file_2
    n_threads = args.n_threads
    output_file = args.output_file
    ploidy = args.ploidy
    
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
    bwa_index_command = BWA_BIN + " index " + reference
    call_or_die(bwa_index_command)

    bwa_align_command = BWA_BIN + " mem -t" + str(n_threads) + " " + reference + " "  + reads_file_1 + " " + reads_file_2 + " > " + working_dir + "/aligned.sam"
    call_or_die(bwa_align_command)
    
    samtools_view_command = SAMTOOLS_BIN + " view -Sb " + working_dir + "/aligned.sam > " + working_dir + "/aligned.bam"
    call_or_die(samtools_view_command)

    markdup_comm = GATK_BIN + " MarkDuplicates" \
                            + " --INPUT "       + working_dir + "/aligned.bam"\
                            + " --OUTPUT "       + working_dir + "/aligned_deduplicated.bam"\
                            + " --METRICS_FILE " + working_dir + "/metrics.txt"\
                            + " --VALIDATION_STRINGENCY SILENT"\
                            + " --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500"\
                            + " --ASSUME_SORT_ORDER queryname"\
                            + " --CREATE_MD5_FILE true"
    call_or_die(markdup_comm)
    

    samtools_sort_command = SAMTOOLS_BIN + " sort -m " + str(max_mem_bytes) + " -f " + working_dir + "/aligned_deduplicated.bam " + working_dir + "/sorted-alns2.bam"
    call_or_die(samtools_sort_command)

    #TODO: rewrite without pipes. A system call with pipes can be dangerous.
    samtools_pileup_command = SAMTOOLS_BIN + " mpileup -uf " + reference + " " + working_dir + "/sorted-alns2.bam | " + BCFTOOLS_BIN + " view --ploidy " + ploidy + " -bvcg - > " + working_dir + "/var.raw.bcf"
    call_or_die(samtools_pileup_command)
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
    
    debug_mode = args.debug
    paired_flag = False if reads_file_2 == "" else True
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    assert(Path(reference).is_file())
    assert(Path(reads_file_1).is_file())
    if paired_flag:
        assert(Path(reads_file_2).is_file())

    assert(Path(BWA_BIN).is_file())
    
    reference_basename = os.path.splitext(reference)[0]
    reference_dict = reference_basename + ".dict"

    print ("############################################")
    print (" Index reference") 
    print ("############################################")
    bwa_index_command = BWA_BIN + " index " + reference
    call_or_die(bwa_index_command)
    
    reference_is_indexed = False  ##TODO(optimization): actual check. This would be automatic using snake-make
    if (not reference_is_indexed):
        bwa_align_command = BWA_BIN + " mem -t" + str(n_threads) + " " + reference + " "  + reads_file_1 + " " + reads_file_2 + " > " + working_dir + "/aligned.sam"
        #TODO(scalability): -a allows whole genome. should it be a parameter or should it be infered from the context
        call_or_die(bwa_align_command)

    faidx_command = SAMTOOLS_BIN + " faidx " + reference
    call_or_die(faidx_command)
    
    if (Path(reference_dict).exists()):
        assert(Path(reference_dict).is_file)
        Path(reference_dict).unlink()
    
    gatk_dict_command = GATK_BIN + " CreateSequenceDictionary --REFERENCE=" + reference + " --OUTPUT=" + reference_dict
    call_or_die(gatk_dict_command)

    print ("############################################")
    print (" Read alignment")
    print ("############################################")
    
    bwa_align_command = BWA_BIN + " mem -K 100000000 -v 3 -t 16 -Y " + reference + " " + reads_file_1 + " " + reads_file_2 + " > " + working_dir + "/aligned_reads.sam"
    call_or_die(bwa_align_command)

######
## The "uBAM" issue of GATK
## -> convert original fastq reads to uBAM

    fq2sam_command =  GATK_BIN + " --java-options -Xmx" + max_memory_MB + "M  FastqToSam" 
    fq2sam_command = fq2sam_command + " -F1 " + reads_file_1
    if (paired_flag):
        fq2sam_command = fq2sam_command + " -F2 " + reads_file_2
    fq2sam_command = fq2sam_command + " --SAMPLE_NAME sample1 -O " + working_dir + "/input_reads_unaligned.ubam"

    call_or_die(fq2sam_command)

## merge above uBAM file with the output of BwaMem
    merge_command = GATK_BIN \
            + " MergeBamAlignment"\
            + " --VALIDATION_STRINGENCY SILENT"\
            + " --EXPECTED_ORIENTATIONS FR"\
            + " --ATTRIBUTES_TO_RETAIN X0"\
            + " --ALIGNED_BAM " + working_dir + "/aligned_reads.sam"\
            + " --UNMAPPED_BAM " + working_dir + "/input_reads_unaligned.ubam"\
            + " --OUTPUT " + working_dir + "/merged_reads.bam"\
            + " --REFERENCE_SEQUENCE " + reference\
            + " --PAIRED_RUN true"\
            + " --SORT_ORDER unsorted"\
            + " --IS_BISULFITE_SEQUENCE false"\
            + " --ALIGNED_READS_ONLY false"\
            + " --CLIP_ADAPTERS false"\
            + " --MAX_RECORDS_IN_RAM 2000000"\
            + " --ADD_MATE_CIGAR true"\
            + " --MAX_INSERTIONS_OR_DELETIONS -1"\
            + " --PRIMARY_ALIGNMENT_STRATEGY MostDistant"\
            + " --UNMAPPED_READ_STRATEGY COPY_TO_TAG"\
            + " --ALIGNER_PROPER_PAIR_FLAGS true"\
            + " --UNMAP_CONTAMINANT_READS true"

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
            + " MarkDuplicates"\
            + " --INPUT " + working_dir + "/merged_reads.bam"\
            + " --OUTPUT " + working_dir + "/aligned_deduplicated.bam "\
            + " --METRICS_FILE " + working_dir + "/metrics.txt "\
            + " --VALIDATION_STRINGENCY SILENT "\
            + " --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 "\
            + " --ASSUME_SORT_ORDER queryname"\
            + " --CREATE_MD5_FILE true"
    call_or_die(mark_duplicates_command)

    sortfix_command = GATK_BIN\
            + " --java-options  -Xms4000m SortSam"\
            + " --INPUT " + working_dir + "/aligned_deduplicated.bam"\
            + " --OUTPUT /dev/stdout"\
            + " --SORT_ORDER coordinate"\
            + " --CREATE_INDEX false"\
            + " --CREATE_MD5_FILE false"\
            + " | "\
            + GATK_BIN\
            + " --java-options -Xms500m"\
            + " SetNmAndUqTags "\
            + " --INPUT /dev/stdin "\
            + " --OUTPUT " + working_dir + "/sortedfixed.bam "\
            + " --CREATE_INDEX true "\
            + " --CREATE_MD5_FILE true "\
            + " --REFERENCE_SEQUENCE " + reference
    call_or_die(sortfix_command)
    
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
    java_opt = ""
    vc_command = GATK_BIN + " --java-options -Xmx" + max_memory_MB + "M " + java_opt + " HaplotypeCaller -ploidy " + ploidy  + " -R " + reference + " -I " + working_dir + "/sortedfixed.bam -O " + working_dir + "/raw_variants.vcf"
    call_or_die(vc_command)
## -L ${interval_list} \ 
## -contamination ${default=0 contamination} ${true="-ERC GVCF" false="" make_gvcf}              
    
    shutil.copy(working_dir + "/raw_variants.vcf", output_file)
    
    #TODO: replace with a native python function. A system call with pipes can be dangerous.
    if (paired_flag):
        histogram_command = SAMTOOLS_BIN + " view -F 0x4 " + working_dir + "/sortedfixed.bam  | cut -f 5 | sort | uniq -c | awk '{print $2 \" \" $1}' > " + working_dir + "/histogram.txt"
        call_or_die(histogram_command)

if __name__ == "__main__":
    VC_main()

