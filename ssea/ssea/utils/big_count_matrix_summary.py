'''
Created on Jan 30, 2014

@author: mkiyer
'''
import logging
import argparse
import sys
import numpy as np

# local imports
from ssea.lib.countdata import BigCountMatrix

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('matrix_dir')
    args = parser.parse_args()
    # get args
    matrix_dir = args.matrix_dir
    bm = BigCountMatrix.open(matrix_dir)
    # get total counts per library
    logging.debug('Getting total counts per library')
    lib_sizes = np.empty(bm.shape[1], dtype=np.float)
    for j in xrange(bm.shape[1]):
        a = bm.counts_t[j,:]
        a = a[np.isfinite(a)]
        lib_sizes[j] = a.sum()
    lib_sizes /= 1.0e6
    # normalize
    logging.debug('Normalizing and summarizing counts per transcript')
    print '\t'.join(['transcript_id', 'exprtot', 'exprmean', 'exprmedian', 'exprmax',
                     'expr9999', 'expr999', 'expr99', 'expr95', 'expr90'])
    for i in xrange(bm.shape[0]):
        t_id = bm.rownames[i]
        length = bm.lengths[i]
        a = bm.counts[i,:]
        valid = np.isfinite(a)
        anorm = (a[valid] / lib_sizes[valid]) / length
        # get stats
        fields = [t_id, np.sum(anorm), np.mean(anorm),
                  np.median(anorm), np.max(anorm),
                  np.percentile(anorm, 99.99),
                  np.percentile(anorm, 99.9),
                  np.percentile(anorm, 99),
                  np.percentile(anorm, 95),
                  np.percentile(anorm, 90)]
        print '\t'.join(map(str, fields))
    bm.close()

if __name__ == '__main__':
    sys.exit(main())
