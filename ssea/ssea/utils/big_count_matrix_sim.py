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


class SimGene:
    def __init__(self, name, length, fpkm, fold_change, frac_de):
        self.name = name
        self.length = int(length)
        self.fpkm = float(fpkm)
        self.fold_change = float(fold_change)
        self.frac_de = float(frac_de)

    def __str__(self):
        return '{}_{}bp_{}fpkm_{}fc_{}de'.\
            format(self.name, self.length, self.fpkm, self.fold_change,
                   self.frac_de)

    @staticmethod
    def header():
        return '\t'.join(['name', 'length', 'fpkm', 'fold_change', 'frac_de'])


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


def sim_gene_counts(g, membership, total_counts):
    # choose de at random
    hit_inds = (membership > 0).nonzero()[0]
    np.random.shuffle(hit_inds)
    num_de = int(round(g.frac_de * hit_inds.shape[0]))
    de_inds = hit_inds[:num_de]
    de = np.zeros(membership.shape[0], dtype=np.int)
    de[de_inds] = 1
    # pick from base or de counts
    counts_base = (total_counts * 1.0e-6) * (g.length * 1.0e-3) * g.fpkm
    counts_de = counts_base * g.fold_change
    counts = np.select([de == 0, de == 1], [counts_base, counts_de])
    # sample from poisson distribution
    counts = np.random.poisson(counts)
    return counts


def sim_counts(genes, samples):
    membership = np.array([s.membership for s in samples], dtype=np.int)
    total_counts = np.array([s.total_counts for s in samples], dtype=np.int)
    for g in genes:
        counts = sim_gene_counts(g, membership, total_counts)
        yield g, counts


def gen_samples(num_samples, num_hits, total_count_range):
    assert(num_samples >= num_hits)
    num_misses = num_samples - num_hits
    samples = []
    for i in xrange(num_samples):
        name = 's{}'.format(i)
        total_counts = np.random.randint(*total_count_range)
        membership = 0 if (i < num_misses) else 1
        samples.append(SimSample(name, total_counts, membership))
    return samples


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--num-samples', type=int, default=None)
    parser.add_argument('--num-hits', type=int, default=None)
    parser.add_argument('--total-count-range', type=str, default=None,
                        help='X,Y (e.g. 50000000,100000000)')
    parser.add_argument('--samples', dest='samples_file')
    parser.add_argument('-o', dest='output_dir')
    parser.add_argument('genes_file')
    args = parser.parse_args()

    # read samples
    samples = []
    if args.samples_file:
        with open(args.samples_file, 'r') as f:
            f.next()
            for line in f:
                fields = line.strip().split('\t')
                samples.append(SimSample(*fields))
    

    # read genes
    genes = []
    with open(args.genes_file, 'r') as f:
        f.next()
        for line in f:
            fields = line.strip().split('\t')
            genes.append(SimGene(*fields))


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
        for g, counts in sim_counts(genes, samples):
            fields = [str(g), str(g.length)]
            fields.extend(map(str, counts))
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
