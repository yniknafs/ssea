'''
SSEA: Sample Set Enrichment Analysis
'''
import os
import logging
import collections
import numpy as np

from dtypes import FLOAT_DTYPE

__author__ = "Matthew Iyer"
__copyright__ = "Copyright 2012-2017"
__credits__ = ["Matthew Iyer"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


class CountMatrix(object):
    VERSION = '0.1.0'
    MEMMAP_FILE = 'counts.memmap'
    MEMMAP_T_FILE = 'counts.transpose.memmap'
    ROWNAMES_FILE = 'rownames.txt'
    COLNAMES_FILE = 'colnames.txt'
    SIZE_FACTORS_FILE = 'size_factors.txt'
    LENGTHS_FILE = 'lengths.txt'

    def __init__(self):
        self.rownames = []
        self.colnames = []
        self.lengths = []
        self.counts = None
        self.counts_t = None
        self.size_factors = None

    @property
    def shape(self):
        return (len(self.rownames), len(self.colnames))

    @staticmethod
    def open(input_dir):
        counts_file = os.path.join(input_dir, CountMatrix.MEMMAP_FILE)
        counts_t_file = os.path.join(input_dir, CountMatrix.MEMMAP_T_FILE)
        rownames_file = os.path.join(input_dir, CountMatrix.ROWNAMES_FILE)
        colnames_file = os.path.join(input_dir, CountMatrix.COLNAMES_FILE)
        size_factors_file = os.path.join(input_dir, CountMatrix.SIZE_FACTORS_FILE)
        lengths_file = os.path.join(input_dir, CountMatrix.LENGTHS_FILE)
        self = CountMatrix()
        self.matrix_dir = input_dir
        with open(rownames_file, 'r') as fileh:
            self.rownames = [line.strip() for line in fileh]
        with open(colnames_file, 'r') as fileh:
            self.colnames = [line.strip() for line in fileh]
        with open(lengths_file, 'r') as fileh:
            self.lengths = np.array([int(line.strip()) for line in fileh])
        with open(size_factors_file, 'r') as fileh:
            self.size_factors = np.array([float(line.strip()) for line in fileh])
        self.counts = np.memmap(counts_file, dtype=FLOAT_DTYPE, mode='r',
                                shape=(len(self.rownames), len(self.colnames)))
        self.counts_t = np.memmap(counts_t_file, dtype=FLOAT_DTYPE, mode='r',
                                  shape=(len(self.colnames), len(self.rownames)))
        return self

    def close(self):
        del self.counts
        del self.counts_t

    def copy(self, output_dir, rowsubset=None, colsubset=None):
        if (rowsubset is None) or (len(rowsubset) == 0):
            rowsubset = set(self.rownames)
        else:
            rowsubset = set(rowsubset)
        if (colsubset is None) or (len(colsubset) == 0):
            colsubset = set(self.colnames)
        else:
            colsubset = set(colsubset)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        counts_file = os.path.join(output_dir, CountMatrix.MEMMAP_FILE)
        counts_t_file = os.path.join(output_dir, CountMatrix.MEMMAP_T_FILE)
        rownames_file = os.path.join(output_dir, CountMatrix.ROWNAMES_FILE)
        colnames_file = os.path.join(output_dir, CountMatrix.COLNAMES_FILE)
        size_factors_file = os.path.join(output_dir, CountMatrix.SIZE_FACTORS_FILE)
        lengths_file = os.path.join(output_dir, CountMatrix.LENGTHS_FILE)
        # write rownames and colnames
        with open(rownames_file, 'w') as fileh:
            row_inds = []
            ind = 0
            for rowname in self.rownames:
                if rowname in rowsubset:
                    print >>fileh, rowname
                    row_inds.append(ind)
                ind += 1
        with open(colnames_file, 'w') as fileh:
            col_inds = []
            ind = 0
            for colname in self.colnames:
                if colname in colsubset:
                    print >>fileh, colname
                    col_inds.append(ind)
                ind += 1
        # write rows
        fp = np.memmap(counts_file, dtype=FLOAT_DTYPE, mode='w+',
                       shape=(len(row_inds),len(col_inds)))
        for i,ind in enumerate(row_inds):
            fp[i,:] = np.array(self.counts[ind,col_inds], dtype=np.float)
        del fp
        # write cols
        fp = np.memmap(counts_t_file, dtype=FLOAT_DTYPE, mode='w+',
                       shape=(len(col_inds),len(row_inds)))
        for j,ind in enumerate(col_inds):
            fp[j,:] = np.array(self.counts_t[ind,row_inds], dtype=np.float)
        del fp
        # size factors
        with open(size_factors_file, 'w') as fileh:
            for ind in col_inds:
                print >>fileh, self.size_factors[ind]
        # lengths
        with open(lengths_file, 'w') as fileh:
            for ind in row_inds:
                print >>fileh, self.lengths[ind]

    @staticmethod
    def from_tsv(input_file, output_dir, na_values=None):
        '''
        convert/copy a tab-delimited file containing numeric weight data
        to a CountMatrix (memmap files)

        input_file: string path to tab-delimited matrix file
        output_dir: output directory
        na_val: set of values corresponding to missing data
        '''
        if na_values is None:
            na_values = set(['NA'])
        else:
            na_values = set(na_values)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        counts_file = os.path.join(output_dir, CountMatrix.MEMMAP_FILE)
        counts_t_file = os.path.join(output_dir, CountMatrix.MEMMAP_T_FILE)
        rownames_file = os.path.join(output_dir, CountMatrix.ROWNAMES_FILE)
        colnames_file = os.path.join(output_dir, CountMatrix.COLNAMES_FILE)
        lengths_file = os.path.join(output_dir, CountMatrix.LENGTHS_FILE)

        self = CountMatrix()
        self.matrix_dir = output_dir
        # get rownames, colnames, and lengths
        with open(input_file, 'r') as fileh:
            header_fields = fileh.next().strip().split('\t')
            self.colnames = list(header_fields[2:])
            for line in fileh:
                rowname, length = line.strip().split('\t', 2)[:2]
                self.rownames.append(rowname)
                self.lengths.append(int(length))

        # write rownames, colnames, and lengths
        with open(rownames_file, 'w') as fileh:
            for rowname in self.rownames:
                print >>fileh, rowname
        with open(colnames_file, 'w') as fileh:
            for colname in self.colnames:
                print >>fileh, colname
        with open(lengths_file, 'w') as fileh:
            for length in self.lengths:
                print >>fileh, str(length)

        # create memmap files
        self.counts = np.memmap(counts_file, dtype=FLOAT_DTYPE, mode='w+',
                                shape=(len(self.rownames), len(self.colnames)))
        with open(input_file, 'r') as fileh:
            fileh.next() # skip header
            for i,line in enumerate(fileh):
                fields = line.strip().split('\t')
                # convert to floats and store
                counts = [(np.nan if x in na_values else float(x))
                          for x in fields[2:]]
                self.counts[i,:] = counts
        # write transpose
        self.counts_t = np.memmap(counts_t_file, dtype=FLOAT_DTYPE, mode='w+',
                                  shape=(len(self.colnames), len(self.rownames)))
        self.counts_t[:] = self.counts.T[:]
        return self

    def _estimate_size_factors_deseq(self):
        '''
        Implements the procedure as described in DESeq:

        Each column is divided by the geometric means of the rows.
        the median of these ratios (skipping genes with zero values)
        is used as the size factor for each column.
        '''
        nrows = len(self.rownames)
        ncols = len(self.colnames)
        # compute row geometric means
        geomeans = np.empty(nrows, dtype=np.float)
        for i in xrange(nrows):
            a = np.around(self.counts[i,:])
            valid = np.logical_and((a > 0), np.isfinite(a))
            if np.all(valid):
                geomeans[i] = np.exp(np.log(a[valid]).mean())
            else:
                geomeans[i] = np.nan
            if i % 10000 == 0:
                logging.debug("%d %d" % (i, np.isfinite(geomeans[:i]).sum()))
        # ignore rows of zero or nan
        valid_geomeans = np.logical_and((geomeans > 0), np.isfinite(geomeans))
        size_factors = np.empty(ncols, dtype=np.float)
        for j in xrange(ncols):
            a = np.around(self.counts_t[j,:])
            lib_valid = np.logical_and((a > 0), np.isfinite(a))
            valid = np.logical_and(lib_valid, valid_geomeans)
            if np.any(valid):
                size_factors[j] = np.median(a[valid] / geomeans[valid])
            else:
                size_factors[j] = np.nan
            #print j, self.colnames[j], size_factors[j], valid.sum(), a[valid].sum()
        return size_factors

    def estimate_size_factors(self, method='total'):
        '''
        estimate size factors (normalization factors) for libraries
        '''
        ncols = len(self.colnames)
        if method == 'deseq':
            size_factors = self._estimate_size_factors_deseq()
        elif method == 'median':
            size_factors = np.empty(ncols, dtype=np.float)
            for j in xrange(ncols):
                a = self.counts_t[j,:]
                valid = np.logical_and((a > 0),np.isfinite(a))
                size_factors[j] = np.median(a[valid])
            size_factors = size_factors / np.median(size_factors)
        elif method == 'total':
            size_factors = np.empty(ncols, dtype=np.float)
            for j in xrange(ncols):
                a = self.counts_t[j,:]
                a = a[np.isfinite(a)]
                size_factors[j] = a.sum()
            size_factors = size_factors / np.median(size_factors)
        self.size_factors = size_factors
        size_factors_file = os.path.join(self.matrix_dir, CountMatrix.SIZE_FACTORS_FILE)
        with open(size_factors_file, 'w') as fileh:
            for x in self.size_factors:
                print >>fileh, x
        return size_factors
