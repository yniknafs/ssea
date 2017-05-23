'''
SSEA: Sample Set Enrichment Analysis
'''
import sys
import os
import argparse
import logging
import numpy as np

from ssea.lib.countdata import BigCountMatrix

__author__ = "Matthew Iyer, Yashar Niknafs"
__copyright__ = "Copyright 2012-2017"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"

DEFAULT_NA_VALUE = 'NA'

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('--na-value', dest='na_values',
                        default=[DEFAULT_NA_VALUE], action='append',
                        help='Value to interpret as missing/invalid '
                        'in weight matrix [default=%(default)s]')
    parser.add_argument('--tsv', dest='tsv_file')
    parser.add_argument('-o', '--output-dir', dest='output_dir')
    # parse args
    args = parser.parse_args()
    tsv_file = args.tsv_file
    output_dir = args.output_dir
    na_values = args.na_values
    # check args
    if not os.path.exists(tsv_file):
        parser.error('Input file "%s" not found' % (tsv_file))
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # convert matrix
    logging.info('Importing')
    bm = BigCountMatrix.from_tsv(tsv_file, output_dir,
                              na_values=na_values)
    logging.info("Estimating size factors")
    bm.estimate_size_factors('deseq')
    bm.close()

if __name__ == '__main__':
    sys.exit(main())
