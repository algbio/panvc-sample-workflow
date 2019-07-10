#!/usr/bin/python3
# TODO: delete this file, as summrize_tools will be directly imported and used.
import sys
from summarize_tools import *

n_args = len(sys.argv);
if(n_args != 4 and n_args != 5):
    print('Number of arguments:', n_args, 'arguments, is incorrect:')
    print('Usage:')
    print(sys.argv[0], 'VC_WORKING_DIR PG_INDEX_DIR READS_1 (READS_2)')
    sys.exit(1);

working_dir = sys.argv[1];
pgindex_dir = sys.argv[2];
reads_1 = sys.argv[3]
reads_2 = ""
paired_flag = False
if (n_args == 5):
    reads_2 = sys.argv[4]
    paired_flag = True
write_summary(working_dir, pgindex_dir, paired_flag, reads_1, reads_2)
