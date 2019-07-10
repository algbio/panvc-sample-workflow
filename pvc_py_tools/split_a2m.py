#!/usr/bin/env python
import sys
import StringIO
from panvc_utils import splitA2M

n_args = len(sys.argv);
if(n_args != 2):
    print 'Number of arguments:', n_args, 'arguments, is incorrect:'
    print 'Usage:'
    print sys.argv[0], 'A2M_FILENAME'
    sys.exit(1);

A2M_FILENAME=sys.argv[1];

splitA2M(A2M_FILENAME)
