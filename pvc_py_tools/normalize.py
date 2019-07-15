#!/usr/bin/python3
import sys
#curr_dir = sys.path[0]
#sys.path.append(curr_dir + '../../pvc_py_tools')
sys.path.append('/home/local/dvalenzu/Repos/Code/PanVC/pvc_py_tools')
from pvc_tools import *
from pathlib import Path
from Bio import SeqIO

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

def apply_vcf(INPUT_FIE,  VCF_FILE_NAME, OUTPUT_FILENAME):
    #### IF FASTA:::
    print ("Input:" + INPUT_FIE)
    records=list(SeqIO.parse(INPUT_FIE, "fasta"));
    assert(len(records) == 1)
    my_seq=records[0].seq;
    #############
    #with open(INPUT_FIE, 'r') as f:
    #    my_seq= f.readline()


    #print(seq_ref.id)
    #print(repr(seq_ref.seq))
    #print(len(seq_ref))

    ## OPEN VCF and modify A and B accordingly

    vcf_file = open(VCF_FILE_NAME);
    ok = 0;
    problem = 0;
    many_vars = 0;

    marker = 0;
    new_a = "";
    new_b = "";

    redundant_vars = 0;
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
        sys.stderr.write('\nSkip a line (redundant var)\n')
        #sys.stderr.write(line)
        redundant_vars = redundant_vars + 1
        #assert(False);
        continue

      if "," in var: 
        many_vars = many_vars + 1;
        continue
      
      if (str(ref_original) != ref):
        sys.stderr.write('In pos: '+str(pos)+' ref is '+ref+' and var is: '+ var+' extract ref is: '+ref_original);
        problem = problem + 1;
        sys.stderr.write('Problem with file:'+ VCF_FILE_NAME+' and with file:'+INPUT_FIE+'\n');
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
    sys.stderr.write('Redundant lines skipped: '+str(redundant_vars)+'\n')
    sys.stderr.write('Many Variatios lines skipped: '+str(many_vars)+'\n')
    sys.stderr.write('Vars applied : '+str(ok)+'\n')
    #print 'Ok:', ok , 'problem:', problem, 'many variations', many_vars;


    ## The following is for debug only, as it will slow the code.
    ## Also this verification is done afterwards in bash and in the projector code.

    if (underlyng_seqs_equals(my_seq, new_a) == 1):
      sys.stderr.write('[apply_vcf.py] Underlyng seqs consistent, Im happy \n');
    else:
      sys.stderr.write('[apply_vcf.py] I FAILED, and changed original seq. CHECK apply_vcf.py \n');
      sys.exit(21); 

    #print '>A';
    out_file = open(OUTPUT_FILENAME, "w")
    out_file.write(new_a + "\n");
    #print '>B\'';
    out_file.write(new_b + "\n");
    out_file.close()
    

    vcf_file.close();

def normalize_vcf(pgindex_dir, all_vcf_files, adhoc_ref_output_folder):
    
    assert(Path(all_vcf_files).is_file())
    assert(Path(adhoc_ref_output_folder).is_dir())
    assert(Path(pgindex_dir).is_dir())
    chr_list_file = pgindex_dir + "/chr_list.txt"
    assert(Path(chr_list_file).is_file())
    
    # TODO: this grep could be done in one pass over all_vcf_files. 
    with open(chr_list_file) as f:
        for line in f:
            chr_id = str(int(line))
            print ("Splitting " + chr_id)
            curr_vcf_dir = adhoc_ref_output_folder + "/" + chr_id 
            curr_vcf_file = curr_vcf_dir + "/vars.vcf"
            assert(Path(curr_vcf_dir).is_dir())
            assert(not Path(curr_vcf_file).exists())
            command  = 'grep "adhoc_ref_chr' + chr_id + '" ' +  all_vcf_files + " > " + curr_vcf_file
            call_or_die(command)
    
    
    with open(chr_list_file) as f:
        for line in f:
            chr_id = str(int(line))
            print ("Normalizing vars in chr: " + chr_id)
            curr_vcf_file = adhoc_ref_output_folder + "/" + chr_id + "/vars.vcf"
            assert(Path(curr_vcf_file).is_file())
            
            curr_adhoc_ref_prefix = adhoc_ref_output_folder + "/" + chr_id + "/adhoc_reference" 
            curr_adhoc_ref_file_fasta = curr_adhoc_ref_prefix + ".fasta"
            curr_adhoc_ref_aligned_to_ref = curr_adhoc_ref_prefix + ".aligned_to_ref"
            assert(Path(curr_adhoc_ref_file_fasta).is_file())
            assert(Path(curr_adhoc_ref_aligned_to_ref).is_file())
            
            curr_aligned_vars =  curr_vcf_file + ".applied"
            assert(not Path(curr_aligned_vars).exists())
            apply_vcf(curr_adhoc_ref_file_fasta, curr_vcf_file, curr_aligned_vars)
            assert(Path(curr_aligned_vars).is_file())

            a1 = pgindex_dir + "/" + chr_id + "/recombinant.n1.gapped"
            a2 = curr_adhoc_ref_aligned_to_ref
            x1 = curr_vcf_file + ".applied.seq1"
            x2 = curr_vcf_file + ".applied.seq2"
            command_head = "head -n1 " + curr_aligned_vars + " | tr -d '\n' > " + x1
            command_tail = "tail -n1 " + curr_aligned_vars + " | tr -d '\n' > " + x2
            call_or_die(command_head)
            call_or_die(command_tail)
            #if [ ${DEBUG_MODE} -eq 1 ]; then
            #   ${DIR}/validate_equal_sequences.sh ${A2} ${X1}
            #fi
            output_prefix = curr_vcf_file
            project_command =  "/home/local/dvalenzu/Repos/Code/PanVC/components/normalize_vcf/projector/projector " + a1 + " " + a2 + " " + x1 + " " + x2 + " " + output_prefix
            call_or_die(project_command)

            new_ref = output_prefix + ".newrefgapped"
            gap_positions_bin = "/home/local/dvalenzu/Repos/Code/PanVC/components/pan_genome_index_real/store_gap_pos/src/store_gaps"
            assert(Path(new_ref).is_file())  # created also by projector
            assert(Path(gap_positions_bin).is_file())
            gap_pos_prefix = new_ref + ".gaps"
            command_gaps = gap_positions_bin + " " + new_ref + " " + gap_pos_prefix
            call_or_die(command_gaps)

            msa = output_prefix + ".msa"
            msa2vcf= "/home/local/dvalenzu/Repos/Code/PanVC/components/normalize_vcf/ext/jvarkit/dist/msa2vcf.jar"
            tmp_vcf = curr_vcf_file + ".normalized.tmp.vcf"
            command_msa2vcf = "java -jar " + msa2vcf + " -c Reference " + msa + " -R " + chr_id + " > " + tmp_vcf
            call_or_die(command_msa2vcf)
            
            normalized_vcf = curr_vcf_file + ".normalized.vcf"
            renormalizer = "/home/local/dvalenzu/Repos/Code/PanVC/components/normalize_vcf/renormalizer/renormalizer"
            renormalizer_command = renormalizer + " " + tmp_vcf + " " + new_ref + ".gaps > " + normalized_vcf 
            call_or_die(renormalizer_command)
            #TODO: the following validation:
            """
            VCFCHECK_BIN=${DIR}/../pan_genome_index_real/ext/vcflib/bin/vcfcheck
            utils_assert_file_exists ${VCFCHECK_BIN}
            if [ ${DEBUG_MODE} -eq 1 ]; then
                >&2 echo "Validating VCFs..."
                TMP_ORIGINAL_REF_FA="${DIR}/tmp_chr${CHR_ID}.fa"
                echo ">${CHR_ID}" > ${TMP_ORIGINAL_REF_FA}
                cat ${A1} | tr -d '-' | tr -d '\n' >> ${TMP_ORIGINAL_REF_FA}
                echo "tmp ref in: ${TMP_ORIGINAL_REF_FA}"
                SAMTOLS_BIN=${DIR}/../../ext_var_call_pipelines/ext/samtools-0.1.19/samtools
                utils_assert_file_exists ${SAMTOOLS_BIN}
                ${SAMTOOLS_BIN} faidx ${TMP_ORIGINAL_REF_FA} 
                utils_assert_file_exists ${TMP_ORIGINAL_REF_FA}.fai
                utils_assert_file_exists ${TMP_ORIGINAL_REF_FA}
                utils_assert_file_exists ${NORMALIZED_VCF}
                CHECK_OUTPUT=$(${VCFCHECK_BIN}  -f ${TMP_ORIGINAL_REF_FA} ${NORMALIZED_VCF})
                if [[ ${CHECK_OUTPUT} ]]; then
                  >&2 echo "Validation failed:"
                  >&2 echo "${CHECK_OUTPUT}"
                  utils_die "stopping here..."
                else
                  >&2 echo "Validation passed."
                fi

                rm -f ${TMP_ORIGINAL_REF_FA}*
            fi
            """

    vcftools_path = "/home/local/dvalenzu/Repos/Code/PanVC/ext_var_call_pipelines/ext/vcftools_0.1.12b/perl/"
    vcfconcat = "/home/local/dvalenzu/Repos/Code/PanVC/ext_var_call_pipelines/ext/vcftools_0.1.12b/perl/vcf-concat"
    assert(Path(vcfconcat).is_file())
    input_vcfs = ""
    with open(chr_list_file) as f:
        for line in f:
            chr_id = str(int(line))
            input_vcfs = input_vcfs + " " + adhoc_ref_output_folder + "/" +chr_id + "/vars.vcf.normalized.vcf"
    
    target_vcf =  adhoc_ref_output_folder + "/all_vars_normalized.vcf"
    assert(not Path(target_vcf).exists())
    concat_command = "perl -I " + vcftools_path + " " + vcfconcat + " " + input_vcfs + " > " + target_vcf
    call_or_die(concat_command)

def normalize_main():
    n_args = len(sys.argv);
    if(n_args != 4):
        print('Number of arguments: ' +str(n_args) +  ' arguments, is incorrect:')
        sys.exit();
    all_vcf_files = sys.argv[1];
    adhoc_ref_output_folder = sys.argv[2];
    pgindex_dir = sys.argv[3];
    normalize_vcf(pgindex_dir, all_vcf_files, adhoc_ref_output_folder)

if __name__ == "__main__":
    normalize_main()


