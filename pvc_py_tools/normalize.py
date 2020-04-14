#!/usr/bin/python3
import sys
from pvc_tools import *
from pathlib import Path
from Bio import SeqIO
from random import randint

PANVC_DIR = sys.path[0]
BIOTOOLS_DIR = "%s/components/normalize_vcf/ext/biotools" % PANVC_DIR
sys.path.append(BIOTOOLS_DIR)
import a2m_alignment_to_vcf as aatv
import combine_vcf

# return 1/0
def underlyng_seqs_equals(seq_1, seq_2):
  clean_1=""
  clean_2=""
  for c in seq_1:
    if(c != '-'):
      clean_1+=c
  for c in seq_2:
    if(c != '-'):
      clean_2+=c
  
  if (clean_1 == clean_2):
    return 1
  else:
    i = 0
    for c in clean_1:
      if(c != clean_2[i]):
        print('error in pos ' + str(i) + ' quitting,')
        print('because: ' + clean_1[i] + ' not eq ' + clean_2[i] + '. bye.')
        return 0
      else:
        i = i+1
    return 0

def apply_vcf(input_fasta,  input_vcf, output_alignment, debug_mode, secondary_vcf, iteration):
    assert(iteration == "First" or iteration == "Second")
    if (secondary_vcf is not "NULL"):
        overlapping_vars_file = open(secondary_vcf, "w")
    else:
        overlapping_vars_file = open("/dev/null", "w")

    #### IF FASTA:::
    print ("Input fasta:" + input_fasta)
    print ("Input VCF: " + input_vcf)
    records=list(SeqIO.parse(input_fasta, "fasta"));
    assert(len(records) == 1)
    my_seq=records[0].seq;
    #############
    #with open(input_fasta, 'r') as f:
    #    my_seq= f.readline()


    #print(seq_ref.id)
    #print(repr(seq_ref.seq))
    #print(len(seq_ref))

    ## OPEN VCF and modify A and B accordingly

    vcf_file = open(input_vcf);
    ok = 0;
    problem = 0;
    many_vars = 0;

    marker = 0;
    new_a = "";
    new_b = "";

    overlapping_vars = 0
    low_qual = 0
    undefined_gt = 0
    deletions = 0

    for line in vcf_file:
        if line.startswith("#"):
            continue;
        values = line.split('\t');
        if (len(values) != 10):
            print("We assumed genotyping one sample at the type.")
            print("If this changed, method to apply vcf must be revisited.")
            print("Aborting")
            exit(33)
      
        pos = int(values[1]) - 1;
        ref = values[3];
        #var = values[4];
        alts = values[4].split(",");
        qual = values[5];
        status = values[6];
        sample_data = values[9]
      
        if (qual == "."):
            qual = "0"
        if (float(qual) < 2.0 and status != "PASS"):
            low_qual = low_qual + 1;
            continue
        curr_len = len(ref)
        ref_original = my_seq[pos:pos+curr_len];
        if (pos < marker):
            sys.stderr.write('\nSkip a line (overlapping var)\n')
            #sys.stderr.write(line)
            overlapping_vars = overlapping_vars + 1
            line = line.replace("0/1", "1/1")
            overlapping_vars_file.write(line)
            continue
        else:
            overlapping_vars_file.write(line)

        
        #sample_data looks like "0/1:blabla"
        if iteration == "First":
            hap_i = 2
        else:
            hap_i = 0
        var_id = sample_data[hap_i]
        if var_id == "0":
            continue
        if var_id == ".":
            #var_id = "1"  # If we wanted to "rescue" a variant?
            continue
        #now var_id is typically 1, but sometimes could be 2 (or three if the vcf were a multisample)
        if (int(var_id) > len(alts)):
            if (debug_mode):
                print ("Current line appears to be malformed, or there is a problem in this parser.")
                exit(33)
            else:
                continue ## just ignoring the offending line
        var = alts[int(var_id)-1]
      
        if (str(ref_original) != ref):
            sys.stderr.write('In pos: '+str(pos)+' ref is '+ref+' and var is: '+ var+' extract ref is: '+ref_original);
            problem = problem + 1;
            sys.stderr.write('Problem with file:'+ input_vcf+' and with file:'+input_fasta+'\n');
            assert(False);
        # We are now in the ok case: 
        ok = ok + 1;
        new_a += str(my_seq[marker:pos]);
        new_b += str(my_seq[marker:pos]);
        
        token_a = ref;
        token_b = var;
      
        if (len(ref) > len(var)):
            token_b += "-" * (len(ref) - len(var));   # we do want an alignment again!
            deletions = deletions + 1;
        elif (len(ref) < len(var)):
            token_a += "-" * (len(var) - len(ref));
        else:
            assert(len(var) == len(ref));
      
        new_a += token_a;
        new_b += token_b;
      
        marker = pos + len(ref);

    # Tail of the file
    new_a += str(my_seq[marker:]);
    new_b += str(my_seq[marker:]);
    sys.stderr.write('Low Qual lines skipped: '+str(low_qual)+'\n')
    sys.stderr.write('Undefined GT (ie "./.") lines skipped: '+str(undefined_gt)+'\n')
    sys.stderr.write('Overlapping lines skipped: '+str(overlapping_vars)+'\n')
    sys.stderr.write('Many Variatios lines skipped: '+str(many_vars)+'\n')
    sys.stderr.write('Vars applied : '+str(ok)+'\n')
    #print 'Ok:', ok , 'problem:', problem, 'many variations', many_vars;


    ## The following is for debug_mode only, as it will slow the code.
    ## Also this verification is done afterwards in bash and in the projector code.
    if (debug_mode):
        if (underlyng_seqs_equals(my_seq, new_a) == 1):
            sys.stderr.write('apply_vcf: Underlyng seqs consistent, everithing ok.\n');
        else:
            sys.stderr.write('apply_vcf: FAILED, changed original seq.\n');
            sys.exit(21); 

    #print '>A';
    out_file = open(output_alignment, "w")
    out_file.write(new_a + "\n");
    #print '>B\'';
    out_file.write(new_b + "\n");
    out_file.close()
    

    vcf_file.close();
    overlapping_vars_file.close()

def break_multiallelic_vars(vcf_filename):
    import shutil
    tmp_filename = vcf_filename +".tmp_vars_" + str(randint(1,1000000))
    command = vcfbreakmulti_bin + " " + vcf_filename + " > " + tmp_filename
    call_or_die(command)
    shutil.copy(tmp_filename, vcf_filename)
    ## TODO: Move instead? or delete tmp_file?

def split_vcf_by_chrom(all_vcf_files, adhoc_ref_output_folder, chr_list):
    multi_output = {chr_id : open(adhoc_ref_output_folder + "/" + chr_id + "/vars.vcf", "w") for chr_id in chr_list}
    with open(all_vcf_files, 'r') as f:
        for line in f:
            if line[0] == "#":
                for a,b in multi_output.items():
                    b.write(line)
                continue
            for chr_id in chr_list:
                if"adhoc_ref_chr" + chr_id in line:
                    multi_output[chr_id].write(line)
    for a,b in multi_output.items():
        b.close()

def validate_equal_underlying_seqs(file_1, file_2):
    with open(file_1) as f:
        seq_1 = f.readline().replace("-","")
    with open(file_2) as f:
        seq_2 = f.readline().replace("-","")
    assert(seq_1 == seq_2)

def merge_homozygous_variants(vcf_1, vcf_2):
    output_1 = []
    output_2 = []
    input_1 = open(vcf_1, "r")
    input_2 = open(vcf_2, "r")
    
    line_2 = input_2.readline()
    line_1 = input_1.readline()
    
    while line_1 and line_2:
        if line_1[0] == "#":
            output_1.append(line_1)
            line_1 = input_1.readline()
            continue
        if line_2[0] == "#":
            output_2.append(line_2)
            line_2 = input_2.readline()
            continue
        pos_1 = int(line_1.split('\t')[1]) - 1;
        pos_2 = int(line_2.split('\t')[1]) - 1;
        if (pos_1 == pos_2):
            if (line_1 == line_2):
                line_1 = line_1.replace("0/1", "1/1")
                output_1.append(line_1)
                pass
            else:
                output_1.append(line_1)
                output_2.append(line_2)
            line_1 = input_1.readline()
            line_2 = input_2.readline()
            continue
        if (pos_1 < pos_2):
            output_1.append(line_1)
            line_1 = input_1.readline()
            continue
        else:
            output_2.append(line_2)
            line_2 = input_2.readline()
            continue
    
    ##Read remaining of file 
    while line_1:
        output_1.append(line_1)
        line_1 = input_1.readline()
    while line_2:
        output_2.append(line_2)
        line_2 = input_2.readline()
    #Write output to same files: 
    input_1.close()
    input_2.close()
    with open(vcf_1,"w") as output_file:
        for line in output_1:
            output_file.write(line)
    with open(vcf_2,"w") as output_file:
        for line in output_2:
            output_file.write(line)

#TODO: rename method, as now the postprocessing does more.
def remove_last_sample_from_vcf(input_vcf, debug_mode):
    output = []
    with open(input_vcf,"r") as input_file:
        for line in input_file:
            if line[0] == "#" and not line.startswith("#CHROM"):
                output.append(line)
                continue
            pos = line.rfind("\t")
            if (pos == -1):
                if (debug_mode):
                    print("Wrong vcf, aborting.\n")
                    print ("Offending line:\n")
                    print (line)
                    exit(33)
                else:
                    continue
            new_line = line[0:pos] + "\n"
            
            if not line.startswith("#CHROM"):
                new_line = new_line.replace("1/1", "0/1")
            
            output.append(new_line)
    with open(input_vcf,"w") as output_file:
        for line in output:
            output_file.write(line)


def concat_vcfs(adhoc_ref_output_folder, chr_list):
    #TODO: vcfconcat is using only the headers from the first vcf file, hence losing the  ##contig=<ID=X,length=400> for all X that are not in the first file.
    #We do not use these header info anyway and the field is optional so this not urgent.
    vcftools_path = PANVC_DIR + "/ext_var_call_pipelines/ext/vcftools_0.1.12b/perl/"
    vcfconcat = PANVC_DIR + "/ext_var_call_pipelines/ext/vcftools_0.1.12b/perl/vcf-concat"
    assert(Path(vcfconcat).is_file())
    input_vcfs = ""
    for chr_id in chr_list:
        input_vcfs = input_vcfs + " " + adhoc_ref_output_folder + "/" +chr_id + "/vars.vcf.normalized.vcf"
    
    target_vcf =  adhoc_ref_output_folder + "/all_vars_normalized.vcf"
    assert(not Path(target_vcf).exists())
    concat_command = "perl -I " + vcftools_path + " " + vcfconcat + " " + input_vcfs + " > " + target_vcf
    call_or_die(concat_command)


def normalize_vcf_haploid_snps_only(pgindex_dir, all_vcf_files, adhoc_ref_output_folder, debug_mode, ploidy):

    def read_file(f):
        return f.read().replace('\n', '')

    def read_path(path):
        with open(path, 'r') as f:
            return read_file(f)

    print("Warning: only haploid SNPs are currently handled.", file = sys.stderr)
    if 1 != ploidy:
        print("Warning: normalization not done.", file = sys.stderr)
        return False
    
    alignment_to_vcf_bin = PANVC_DIR + "/components/normalize_vcf/ext/biotools/a2m_alignment_to_vcf.py"
    
    # Split by chromosome.
    chr_list = PVC_get_chr_list(pgindex_dir)
    split_vcf_by_chrom(all_vcf_files, adhoc_ref_output_folder, chr_list)

    for chr_id in chr_list:
        print("Normalizing variants in chromosome %s" % chr_id, file = sys.stderr)

        # Get the original multialigned reference.
        original_ref_aln = "%s/%s/recombinant.n1.gapped" % (pgindex_dir, chr_id)
        assert(Path(original_ref_aln).is_file())
        with open(original_ref_aln, 'r') as original_ref_aln_file:
            original_ref_aln_seq = read_file(original_ref_aln_file)

            # Get the VCF file.
            curr_vcf_file = "%s/%s/vars.vcf" % (adhoc_ref_output_folder, chr_id)
            assert(Path(curr_vcf_file).is_file())

            # Get the current ad-hoc reference.
            curr_adhoc_ref_prefix = "%s/%s/adhoc_reference" % (adhoc_ref_output_folder, chr_id)
            adhoc_ref_aln = curr_adhoc_ref_prefix + ".aligned_to_ref"
            assert(Path(adhoc_ref_aln).is_file())
            adhoc_ref_aln_seq = read_path(adhoc_ref_aln)

            # Generate a variant file from the ad-hoc reference w.r.t. the original reference.
            adhoc_variants_path = "%s/%s/adhoc_reference_variants.vcf" % (adhoc_ref_output_folder, chr_id)
            assert(not Path(adhoc_variants_path).exists())
            with open(adhoc_variants_path, "w") as dst_file:
                aatv.process_sequences([original_ref_aln_seq, adhoc_ref_aln_seq], ["original", "chr%s" % chr_id], dst_file, "adhoc_ref_chr%s" % chr_id)

            # Combine the variant files, replace the REF column with the reference.
            # Since only SNPs are allowed, the gapped reference should not actually have gaps.
            output_vcf = curr_vcf_file + ".normalized.vcf"
            assert(not Path(output_vcf).exists())
            print("Combining %s and %s to %s" % (adhoc_variants_path, curr_vcf_file, output_vcf), file = sys.stderr)
            with open(adhoc_variants_path, "r") as adhoc_vars_file, open(curr_vcf_file, "r") as vars_file, open(output_vcf, "w") as dst_file:
                combine_vcf.combine_vcf([vars_file, adhoc_vars_file], dst_file, chr_id = "adhoc_ref_chr%s" % chr_id, reference_seq = original_ref_aln_seq)
    concat_vcfs(adhoc_ref_output_folder, chr_list)
    return True

    
def normalize_vcf(pgindex_dir, all_vcf_files, adhoc_ref_output_folder, debug_mode, ploidy):
    alignment_to_vcf_bin = PANVC_DIR + "/components/normalize_vcf/combine_msa_vcf" # FIXME: move to correct place
    
    # Split by chromosome.
    chr_list = PVC_get_chr_list(pgindex_dir)
    split_vcf_by_chrom(all_vcf_files, adhoc_ref_output_folder, chr_list)
    
    for chr_id in chr_list:
        print("Normalizing variants in chromosome %s" % chr_id, file = sys.stderr)

        # Get the original multialigned reference.
        original_ref_aln = "%s/%s/recombinant.n1.gapped" % (pgindex_dir, chr_id)
        assert(Path(original_ref_aln).is_file())
        
        # Get the VCF file.
        curr_vcf_file = "%s/%s/vars.vcf" % (adhoc_ref_output_folder, chr_id)
        assert(Path(curr_vcf_file).is_file())
        
        # Get the current ad-hoc reference.
        curr_adhoc_ref_prefix = "%s/%s/adhoc_reference" % (adhoc_ref_output_folder, chr_id)
        adhoc_ref_aln = curr_adhoc_ref_prefix + ".aligned_to_ref"
        assert(Path(adhoc_ref_aln).is_file())
        
        # Combine and normalize.
        output_vcf = curr_vcf_file + ".normalized.vcf"
        combine_command = "%s --ref='%s' --alt='%s' --variants='%s' --output-chr='%s' --ploidy='%d' > %s" % \
            (alignment_to_vcf_bin, original_ref_aln, adhoc_ref_aln, curr_vcf_file, chr_id, ploidy, output_vcf)
        call_or_die(combine_command)
    
    concat_vcfs(adhoc_ref_output_folder, chr_list)
    return True
