#!/usr/bin/python3
import sys
from pvc_tools import *
from pathlib import Path
from Bio import SeqIO
from random import randint

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

def apply_vcf(input_fasta,  input_vcf, output_alignment, debug_mode, secondary_vcf):
    if (secondary_vcf is not "NULL"):
        overlapping_vars_file = open(secondary_vcf, "w")
    else:
        #f = open(os.devnull, 'w')
        overlapping_vars_file = open("/dev/null", "w")

    #### IF FASTA:::
    print ("Input:" + input_fasta)
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

    overlapping_vars = 0;
    low_qual = 0;
    deletions = 0;

    for line in vcf_file:
      if line.startswith("#"):
        continue;
      values = line.split('\t');
      pos=int(values[1]) - 1;
      ref=values[3];
      var=values[4];
      qual=values[5];
      status=values[6];
      
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
        overlapping_vars_file.write(line)
        #assert(False);
        continue

      if "," in var: 
        many_vars = many_vars + 1;
        continue
      
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
    tmp_filename = "./tmp_vars_" + str(randint(1,1000000))
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


def projected_alignment_to_vcf(curr_aligned_vars, chr_id, pgindex_dir, curr_adhoc_ref_prefix, vars_output_filename):
    pass

def normalize_vcf(pgindex_dir, all_vcf_files, adhoc_ref_output_folder, debug_mode):
    
    assert(Path(all_vcf_files).is_file())
    assert(Path(adhoc_ref_output_folder).is_dir())
    assert(Path(pgindex_dir).is_dir())
    chr_list = PVC_get_chr_list(pgindex_dir)
    
    split_vcf_by_chrom(all_vcf_files, adhoc_ref_output_folder, chr_list)
    
    PANVC_DIR = sys.path[0]
    for chr_id in chr_list:
        print ("Normalizing vars in chr: " + chr_id)
        curr_vcf_file = adhoc_ref_output_folder + "/" + chr_id + "/vars.vcf"
        assert(Path(curr_vcf_file).is_file())
        
        break_multiallelic_vars(curr_vcf_file)
        
        curr_adhoc_ref_prefix = adhoc_ref_output_folder + "/" + chr_id + "/adhoc_reference" 
        
        curr_adhoc_ref_file_fasta = curr_adhoc_ref_prefix + ".fasta"
        assert(Path(curr_adhoc_ref_file_fasta).is_file())
        
        curr_aligned_vars =  curr_vcf_file + ".applied"
        assert(not Path(curr_aligned_vars).exists())
        secondary_vcf = "./tmp_secondary_vcf_" + str(randint(1,10000)) 
        assert(not Path(secondary_vcf).exists())
        apply_vcf(curr_adhoc_ref_file_fasta, curr_vcf_file, curr_aligned_vars, debug_mode, secondary_vcf)
        assert(Path(curr_aligned_vars).is_file())
        assert(Path(secondary_vcf).is_file())
         
        curr_adhoc_ref_aligned_to_ref = curr_adhoc_ref_prefix + ".aligned_to_ref"
        assert(Path(curr_adhoc_ref_aligned_to_ref).is_file())
            
        

        a1 = pgindex_dir + "/" + chr_id + "/recombinant.n1.gapped"
        a2 = curr_adhoc_ref_aligned_to_ref
        x1 = curr_vcf_file + ".applied_vars.seq1"
        x2 = curr_vcf_file + ".applied_vars.seq2"
        #TODO: replace with a native python function. A system call with pipes can be dangerous.
        command_head = "head -n1 " + curr_aligned_vars + " | tr -d '\n' > " + x1
        command_tail = "tail -n1 " + curr_aligned_vars + " | tr -d '\n' > " + x2
        call_or_die(command_head)
        call_or_die(command_tail)
        
        if debug_mode:
            validate_equal_underlying_seqs(a2, x1)
        output_prefix = curr_vcf_file
        project_command =  PANVC_DIR + "/components/normalize_vcf/projector/projector " + a1 + " " + a2 + " " + x1 + " " + x2 + " " + output_prefix
        call_or_die(project_command)

        new_ref = output_prefix + ".newrefgapped"
        gap_positions_bin = PANVC_DIR + "/components/pan_genome_index_real/store_gap_pos/src/store_gaps"
        assert(Path(new_ref).is_file())  # created also by projector
        assert(Path(gap_positions_bin).is_file())
        gap_pos_prefix = new_ref + ".gaps"
        command_gaps = gap_positions_bin + " " + new_ref + " " + gap_pos_prefix
        call_or_die(command_gaps)

        msa = output_prefix + ".msa"
        msa2vcf= PANVC_DIR + "/components/normalize_vcf/ext/jvarkit/dist/msa2vcf.jar"
        tmp_vcf = curr_vcf_file + ".normalized.tmp.vcf"
        command_msa2vcf = "java -jar " + msa2vcf + " -c Reference " + msa + " -R " + chr_id + " > " + tmp_vcf
        call_or_die(command_msa2vcf)
            
        normalized_vcf = curr_vcf_file + ".normalized.vcf"
        renormalizer = PANVC_DIR  + "/components/normalize_vcf/renormalizer/renormalizer"
        renormalizer_command = renormalizer + " " + tmp_vcf + " " + new_ref + ".gaps > " + normalized_vcf 
        call_or_die(renormalizer_command)
        if (debug_mode):
            print ("Validating vcf")
            validate_vcf_command = VCFCHECK_BIN + " -f " + pgindex_dir + "/std_ref.fa " + normalized_vcf 
            validation_result = call_and_get_result(validate_vcf_command)
            if (validation_result != ""):
                #Note: the original version of this code built the fai index, but
                #apparently vcfcheck builds it whent it is not found.
                print ("The normalized VCF does not pass the validation. This should not occur.")
                print ("Validation command and output follows:")
                print ("command: " + validate_vcf_command)
                print ("output : " + validation_result + "\n")
                exit(33)

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
