#!/usr/bin/env python
import sys
import StringIO
from panvc_utils import sequence_num_to_name

n_args = len(sys.argv);
if(n_args != 6):
    print 'Number of arguments:', n_args, 'arguments, is incorrect:'
    print 'Usage:'
    print sys.argv[0], 'samples_filename PLOIDITY SEQ_ID'
    sys.exit(1);

SAMPLES_FILENAME=sys.argv[1];
CHRS_FILENAME=sys.argv[2];
PLOIDITY=int(sys.argv[3]);
CHR_ID=int(sys.argv[4]);
SEQ_ID=int(sys.argv[5]);
assert (SEQ_ID >= 1)

answer = sequence_num_to_name(SAMPLES_FILENAME, CHRS_FILENAME, PLOIDITY, CHR_ID, SEQ_ID)

print(answer)

