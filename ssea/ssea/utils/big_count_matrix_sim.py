'''
SSEA: Sample Set Enrichment Analysis
'''
import sys
import os
import argparse
import logging
from collections import namedtuple
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


ABUNDANCE_MEMMAP_FILE = 'abundances.memmap'
COUNTS_MEMMAP_FILE = 'counts.memmap'


class SimGene:
    def __init__(self, name, length, expr, fold_change, frac_de):
        self.name = name
        self.length = int(length)
        self.expr = float(expr)
        self.fold_change = float(fold_change)
        self.frac_de = float(frac_de)

    def __str__(self):
        return '{}_{}bp_{}expr_{}fc_{}de'.\
            format(self.name, self.length, self.expr, self.fold_change,
                   self.frac_de)

    @staticmethod
    def header():
        return '\t'.join(['name', 'length', 'expr', 'fold_change', 'frac_de'])


class SimSample:
    def __init__(self, name, total_counts, membership):
        self.name = name
        self.total_counts = int(total_counts)
        self.membership = int(membership)

    def __str__(self):
        return '{}_{}_{}'.format(self.name, self.total_counts, self.membership)

    @staticmethod
    def header():
        return '\t'.join(['name', 'total_counts', 'membership'])


def gen_samples(num_samples, num_cancers, total_count_range):
    samples = []
    num_normals = num_samples - num_cancers
    for n in xrange(num_cancers):
        total_counts = np.random.randint(total_count_range[0], total_count_range[1])
        s = SimSample('C{}'.format(n), total_counts, 1)
        samples.append(s)
    for n in xrange(num_normals):
        total_counts = np.random.randint(total_count_range[0], total_count_range[1])
        s = SimSample('N{}'.format(n), total_counts, 0)
        samples.append(s)
    return samples


def read_genes(genes_file, null_expr_frac):
    genes = []
    total_expr = 0.0
    with open(genes_file, 'r') as f:
        f.next()
        for line in f:
            fields = line.strip().split('\t')
            g = SimGene(*fields)
            genes.append(g)
            total_expr += g.expr
    # recalculate gene relative abundances
    for g in genes:
        g.expr = (g.expr / total_expr) * (1.0 - null_expr_frac)
    return genes


def gen_null_genes(num_null_genes, null_expr_frac, length_range):
    null_exprs = np.random.random(num_null_genes)
    null_exprs = (null_exprs / null_exprs.sum()) * null_expr_frac
    genes = []
    for i in xrange(num_null_genes):
        length = np.random.randint(length_range[0], length_range[1])
        g = SimGene('null{}'.format(i), length, null_exprs[i], 1.0, 0.0)
        genes.append(g)
    return genes


def sim_rel_abundances(g, membership):
    # choose de at random
    hit_inds = (membership > 0).nonzero()[0]
    np.random.shuffle(hit_inds)
    num_de = int(round(g.frac_de * hit_inds.shape[0]))
    de_inds = hit_inds[:num_de]
    de = np.zeros(membership.shape[0], dtype=np.int)
    de[de_inds] = 1
    # pick from base or de counts
    expr_base = (g.length * 1.0e-3) * g.expr
    expr_de = expr_base * g.fold_change
    exprs = np.select([de == 0, de == 1], [expr_base, expr_de])
    return exprs


def create_relative_abundance_matrix(samples, genes, output_dir):
    abundance_memmap_file = os.path.join(output_dir, ABUNDANCE_MEMMAP_FILE)
    m = np.memmap(abundance_memmap_file, dtype=np.float, mode='w+',
                  shape=(len(genes), len(samples)))
    membership = np.array([s.membership for s in samples], dtype=np.int)
    # apply fold change and differential expression frac to genes
    for i, g in enumerate(genes):
        exprs = sim_rel_abundances(g, membership)
        m[i,:] = exprs
    # renormalize each sample such that sample exprs sum to 1
    for j, s in enumerate(samples):
        total_exprs = m[:,j].sum()
        m[:,j] /= total_exprs
    return m


def sim_sample_counts(relexprs, total_counts, step=100000):
    counts = np.zeros(relexprs.shape[0], dtype=np.float)
    probs = relexprs.cumsum()
    c = 0
    while (c + step) <= total_counts:
        counts[:] += np.bincount(probs.searchsorted(np.random.random(step)), minlength=counts.shape[0])
        c += step
    if (total_counts - c) > 0:
        step = total_counts - c
        counts[:] += np.bincount(probs.searchsorted(np.random.random(step)), minlength=counts.shape[0])
    return counts


def sim_counts(abundance_mat, samples, output_dir):
    total_counts = np.array([s.total_counts for s in samples], dtype=np.int)
    counts_memmap_file = os.path.join(output_dir, COUNTS_MEMMAP_FILE)
    count_mat = np.memmap(counts_memmap_file, dtype=np.float, mode='w+',
                          shape=abundance_mat.shape)
    for j, s in enumerate(samples):
        logging.debug('Simulating reads for sample %d (%s)' % (j, str(s)))
        counts = sim_sample_counts(abundance_mat[:, j], s.total_counts)
        count_mat[:, j] = counts
    return count_mat


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', dest='verbose',
                        action="store_true",
                        default=False,
                        help='Enabled detailed logging '
                        '(for debugging)')
    parser.add_argument('-o', dest='output_dir')
    parser.add_argument('-n', '--num-samples', dest='num_samples',
                        type=int, default=0)
    parser.add_argument('-ncancer', '--num-cancers', dest='num_cancers',
                        type=int, default=0)
    parser.add_argument('--total-count-range',
                        type=str,
                        default='50000000,100000000',
                        help='X,Y (e.g. 50000000,100000000)')
    parser.add_argument('--samples', dest='samples_file')
    parser.add_argument('--num-null-genes', dest='num_null_genes',
                        type=int, default=1000)
    parser.add_argument('--null-expr-frac', dest='null_expr_frac',
                        type=float, default=0.9)
    parser.add_argument('--null-length-range', dest='null_length_range',
                        type=str,
                        default='250,25000',
                        help='X,Y (e.g. 250,25000)')
    parser.add_argument('genes_file')
    args = parser.parse_args()

    # setup logging
    if args.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(level=level,
                        format="%(asctime)s pid=%(process)d "
                               "%(levelname)s - %(message)s")

    # input args
    num_samples = args.num_samples
    num_cancers = args.num_cancers
    assert num_cancers <= num_samples
    num_normals = num_samples - num_cancers
    total_count_range = map(int, args.total_count_range.split(','))
    num_null_genes = args.num_null_genes
    null_expr_frac = args.null_expr_frac
    null_length_range = map(int, args.null_length_range.split(','))

    # samples
    samples = []
    samples.extend(gen_samples(num_samples, num_cancers, total_count_range))
    # samples specified in file
    if args.samples_file:
        with open(args.samples_file, 'r') as f:
            f.next()
            for line in f:
                fields = line.strip().split('\t')
                samples.append(SimSample(*fields))

    # read genes
    genes = read_genes(args.genes_file, null_expr_frac)
    # generate null genes
    genes.extend(gen_null_genes(num_null_genes, null_expr_frac,
                                null_length_range))

    abundance_mat = create_relative_abundance_matrix(samples, genes, args.output_dir)
    count_mat = sim_counts(abundance_mat, samples, args.output_dir)

    # output dir
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # print count matrix
    counts_tsv_file = os.path.join(args.output_dir, 'counts.tsv')
    with open(counts_tsv_file, 'w') as f:
        fields = ['rowname', 'length']
        for s in samples:
            fields.append(str(s))
        print >>f, '\t'.join(fields)
        for i in xrange(count_mat.shape[0]):
            g = genes[i]
            fields = [str(g), str(g.length)]
            fields.extend(map(str, count_mat[i,:]))
            print >>f, '\t'.join(fields)

    # print sample set file
    sample_set_file = os.path.join(args.output_dir, 'sample_set.smt')
    with open(sample_set_file, 'w') as f:
        fields = ['name', 'desc']
        fields.extend([str(s) for s in samples])
        print >>f, '\t'.join(fields)
        fields = [args.output_dir, 'simulated counts']
        fields.extend([str(s.membership) for s in samples])
        print >>f, '\t'.join(fields)



if __name__ == '__main__':
    sys.exit(main())
