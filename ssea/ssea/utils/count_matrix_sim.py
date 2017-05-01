'''
SSEA: Sample Set Enrichment Analysis
'''
import sys
import os
import argparse
import logging
import numpy as np

from ssea.lib.countdata import CountMatrix

__author__ = "Matthew Iyer, Yashar Niknafs"
__copyright__ = "Copyright 2012-2017"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"

DEFAULT_LENGTH = 1000
DEFAULT_VALUE = 100

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--nrows', type=int, default=100)
    parser.add_argument('--ncols', type=int, default=500)
    args = parser.parse_args()

    header_fields = ['rowname', 'length']
    header_fields.extend(['S%d' % i for i in xrange(args.ncols)])
    print '\t'.join(header_fields)
    for i in xrange(args.nrows):
        fields = ['G%d' % i, DEFAULT_LENGTH]
        fields.extend([DEFAULT_VALUE for i in xrange(args.ncols)])
        print '\t'.join(map(str, fields))
    return 0

if __name__ == '__main__':
    sys.exit(main())
