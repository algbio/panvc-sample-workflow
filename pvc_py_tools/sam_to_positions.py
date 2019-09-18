#!/usr/bin/python3
import sys
import re
#import StringIO
import gzip
import os 
from Bio import SeqIO

### Assumptuions:
### The sam record contains the pattern, we are not covered for the * case, so it will trigger an
### index out of bouds while verifying pattern.


op_finder = re.compile('\d\d*\D')

# SENSIBILITY = NONE == 0 
# Means that no error will be accepted.
# SENSIBILITY = ALL means max_error = infinity.
# Any integer is accepted and is the max error allowed.


# return 1 / 0 depending on if the string was successfully understood or no
def cigar_to_inervals(output, cigar_string, pattern, start_pos, reference, max_error):
  
  #print >> sys.stderr, 'Reference: ' + str(reference)
  #print >> sys.stderr, 'Reference length: ' + str(len(reference))
  count_errors = 0;
  original_cigar = cigar_string
  #print original_cigar
  #print pattern
  #print start_pos
  while (len(cigar_string) != 0):
    #process operations in CIGAR string:
      #print "sp:"+str(start_pos);
      op = op_finder.match(cigar_string);
      op_str = op.group();
      #cigar_string_range = op.span();
      cigar_string=cigar_string[op.end():]
      #print cigar_string;
      op_length = int((op_str[:-1]));
      op_type = (op_str[-1]);
      if(op_type == "M"):
        #validate is match:
        next_pos = start_pos + op_length;
        match_count = 0;
        report_pos = start_pos;
        for i in range(op_length):
          #print >> sys.stderr, 'Inspecting: ' + str(start_pos) + "+"+str(i) +":"+ str(start_pos+i)
          if (reference[start_pos + i] == pattern[i]):
            if (match_count == 0):
              report_pos = start_pos + i;
            match_count += 1;
          else:
            count_errors += 1; 
            if (count_errors > max_error):
              return True
            if (match_count != 0):
              line=(str(report_pos)+" "+str(match_count));
              output.write(line+"\n");
            ##start_pos = start_pos + i;
            match_count = 0;
        if (match_count!= 0):
          line=(str(report_pos)+" "+str(match_count));
          output.write(line+"\n");
        #else:
        #  print >> sys.stderr, original_cigar+'[XXXXXXX]CIGAR is a happy match'; 
        start_pos = start_pos + op_length;
        pattern = pattern[op_length:];
      elif (op_type == "D"):
        # Deletion from the reference:
        # We only advance the start position. dont reduce pattern, dont emit nothing.
        count_errors += op_length 
        if (count_errors > max_error):
          return True
        start_pos = start_pos + op_length;
      elif (op_type == "S"):
        # we dont move start position: 
        # we are processing an insertion in the pattern, we procced to remove this area from the pattern, and that is
        count_errors += op_length 
        if (count_errors > max_error):
          return True
        pattern = pattern[op_length:];
      elif (op_type == "I"):
        # we dont move start position: 
        # we are processing an insertion into  the referecne, we procced to remove this area from the pattern, and that is
        count_errors += op_length 
        if (count_errors > max_error):
          return True
        pattern = pattern[op_length:];
      elif (op_type == "H"):
        # we dont do nothing. The string in pattern was already "HARD CLIPPED".
        count_errors += op_length 
        if (count_errors > max_error):
          return True
      else: 
        print >> sys.stderr, original_cigar+'[XXXXXXX]CIGAR string unkkown.';
        return True;
  if (len(pattern) > 0):
    print("Error!. Inconsitency in sam format!");
    sys.exit(1);
    return True;
  return False;

def read_fasta_special(file_name, seq_name):
    fasta_sequences = SeqIO.parse(open(file_name),'fasta')
    for curr_fasta in fasta_sequences:
        name, sequence = curr_fasta.id, str(curr_fasta.seq)
        if (name == seq_name):
            # print >> sys.stderr, 'bingo'
            return sequence
        #else:
        #    print >> sys.stderr, 'nope: '+name+" != "+seq_name


def sam_process(input_file_name, output_obj, reference, max_error):
  f = gzip.open(input_file_name);
  not_maped= 0;
  cigar_unknown = 0;
  maped = 0;
  for line in f:
    values = line.decode().split('\t');  #### herea
    if (len(values) < 9):
      continue;
    if (values[1] == '4'):
      not_maped = not_maped + 1
      continue;
    cigar_string = values[5];
    pattern = values[9];
    start_pos = int(values[3]) - 1;
    #print >> sys.stderr, 'Processing:'
    #print >> sys.stderr, line
    error = cigar_to_inervals(output_obj, cigar_string, pattern, start_pos, reference, max_error);
    if (error):
      cigar_unknown = cigar_unknown + 1;
    else:
      maped = maped + 1;

  f.close();
  return (maped, cigar_unknown, not_maped);

def main():
    n_args = len(sys.argv);
    if(n_args != 7):
      print >> sys.stderr, 'Number of arguments:', n_args, 'arguments, is incorrect:'
      sys.exit();

    #print >> sys.stderr, 'Running:'
    #print >> sys.stderr, sys.argv[0] + ' ' + sys.argv[1]+ ' ' + sys.argv[2] + ' ' + sys.argv[3]

    SAM_FOLDER=sys.argv[1];
    REFERENCE_FILE_NAME=sys.argv[2];
    CHR_LIST_FILE_NAME=sys.argv[3];
    SENSIBILITY=sys.argv[4];
    N_REFS=sys.argv[5];
    LOG_FILE_NAME=sys.argv[6];

    SamToPos(SAM_FOLDER, REFERENCE_FILE_NAME, CHR_LIST_FILE_NAME, SENSIBILITY, N_REFS, LOG_FILE_NAME)

def SamToPos(SAM_FOLDER, REFERENCE_FILE_NAME, CHR_LIST_FILE_NAME, SENSIBILITY, N_REFS, LOG_FILE_NAME):
    log_file = open(LOG_FILE_NAME, 'w')
    if (SENSIBILITY == 'NONE'):
      #max_error = True;
      max_error = 0;
      log_file.write('[sam_to_positions]: ****** Only exact matches accepted. *******\n')
    elif (SENSIBILITY == 'ALL'):
      #max_error = False;
      max_error = float("inf")
      log_file.write('[sam_to_positions]: ****** All matches accepted. *******\n')
    else:
      #max_error = False;
      max_error = int(SENSIBILITY);
      log_file.write('[sam_to_positions]: ****** Max size of mismatch/indels is:'+ SENSIBILITY +'*******\n')

    fasta_sequences = SeqIO.parse(open(REFERENCE_FILE_NAME),'fasta')

    with open(CHR_LIST_FILE_NAME) as f:
        chr_list = f.readlines()
    chr_list = [x.strip() for x in chr_list] 
    n_chrs=len(chr_list)

    i = 0
    for curr_fasta in fasta_sequences:
        name, sequence = curr_fasta.id, str(curr_fasta.seq)
        ref_id= 1 +(i//n_chrs)
        chr_id=chr_list[i%n_chrs]
        ## This validation is no longer used, as it would require new parameter to this script.
        #sample_name = sequence_num_to_name(SAMPLES_FILENAME, PLOIDY, ref_id)
        #seq_name=sample_name+"_CHR"+str(chr_id)
        #print >> log_file, 'My name:'+seq_name;
        #print >> log_file, 'FA name:'+name;
        #assert(name == seq_name)
        CURR_SAM_FOLDER=SAM_FOLDER+"/"+str(chr_id)
        INPUT_SAM_FILENAME=CURR_SAM_FOLDER+"/mapped_reads_to"+str(ref_id)+".sam.gz"
        OUTPUT_FILENAME=CURR_SAM_FOLDER+"/mapped_reads_to"+str(ref_id)+".pos"
        output_object = open(OUTPUT_FILENAME, "w")
        (maped, cigar_unknown, not_maped) =  sam_process(INPUT_SAM_FILENAME, output_object, sequence, max_error);
        output_object.close()
        total = maped + cigar_unknown + not_maped;
        log_file.write('[sam_to_positions]: *************************\n')
        log_file.write('[sam_to_positions]: Reports of mapes reads REPORT(MAPED|CIGAR_UNKNOWN|NOT_MAPED|TOTAL)\n')
        log_file.write('REPORT('+str(maped)+'|'+str(cigar_unknown)+'|'+str(not_maped)+'|'+str(total)+')\n')
        log_file.write('[sam_to_positions]: *************************\n')
        i = i + 1
    log_file.write('NEW SAM TO POS DONE\n')
