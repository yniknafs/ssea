'''
SSEA: Sample Set Enrichment Analysis
'''
import os
import argparse
import logging
import pickle

from base import WEIGHT_METHODS, ParserError, SampleSet
from countdata import CountMatrix


__author__ = "Matthew Iyer, Yashar Niknafs"
__copyright__ = "Copyright 2012-2017"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


class Args:
    # default arguments
    PROG = 'ssea'
    DESCRIPTION = 'Sample Set Enrichment Analysis (RNA-Seq)'
    VERBOSE = False
    NUM_PROCESSES = 1
    OUTPUT_DIR = 'ssea_output'
    NUM_PERMUTATIONS = 1000
    WEIGHT_MISS_DEFAULT = 'log'
    WEIGHT_HIT_DEFAULT = 'log'
    WEIGHT_PARAM = 1.0
    NOISE_LOC = 1.0
    NOISE_SCALE = 1.0
    CONFIG_JSON_FILE = 'config.json'
    SAMPLE_SET_JSON_FILE = 'sample_set.json'
    RESULTS_JSON_FILE = 'results.json'
    OUTPUT_HISTS_FILE = 'hists.npz'
    SAMPLE_SET_INFO_FILE = 'sample_sets.tsv'
    LOG_DIR = 'log'
    TMP_DIR = 'tmp'

    @staticmethod
    def create():
        parser = argparse.ArgumentParser(description=Args.DESCRIPTION)
        parser.add_argument('-v', '--verbose', dest='verbose',
                            action="store_true",
                            default=Args.VERBOSE,
                            help='Enabled detailed logging '
                            '(for debugging)')
        parser.add_argument('-p', '--num-processes', type=int,
                            metavar='N',
                            dest='num_processes',
                            default=Args.NUM_PROCESSES,
                            help='Run in parallel with N '
                            'processes [default=%(default)s]')
        parser.add_argument("-o", "--output-dir", dest="output_dir",
                            metavar='DIR',
                            default=Args.OUTPUT_DIR,
                            help='directory where output files will be '
                            'stored (if already exists then --resume must '
                            'be specified) [default=%(default)s]')
        parser.add_argument('--perms', dest='perms',
                            type=int, metavar='N',
                            default=Args.NUM_PERMUTATIONS,
                            help='Number of permutations '
                            '[default=%(default)s]')
        parser.add_argument('--weight-miss', dest='weight_miss',
                         	choices=WEIGHT_METHODS.keys(),
                         	default=Args.WEIGHT_MISS_DEFAULT,
                         	help='Weighting method for elements not in set '
                         	'[default=%(default)s]')
        parser.add_argument('--weight-hit', dest='weight_hit',
                         	choices=WEIGHT_METHODS.keys(),
                         	default=Args.WEIGHT_HIT_DEFAULT,
                         	help='Weighting method for elements in set '
                         	'[default=%(default)s]')
        parser.add_argument('--weight-param', dest='weight_param', type=float,
                         	default=Args.WEIGHT_PARAM,
                         	help='Either log2(n + X) for log transform or '
                         	'pow(n,X) for exponential (root) transform '
                         	'[default=%(default)s]')
        parser.add_argument('--noise-loc', dest='noise_loc', type=float,
                         	default=Args.NOISE_LOC,
                         	help='noise parameter (advanced)'
                         	'[default=%(default)s]')
        parser.add_argument('--noise-scale', dest='noise_scale', type=float,
                         	default=Args.NOISE_SCALE,
                         	help='noise parameter (advanced)'
                         	'[default=%(default)s]')
        parser.add_argument('-s', '--sample-set', required=True,
                            dest='sample_set_file',
                            help='File containing sample set')
        parser.add_argument('-c', '--count-matrix', required=True,
                            dest='count_matrix_path',
                            help='Path to count matrix')
        return parser

    @staticmethod
    def load(filename):
        return pickle.load(open(filename))

    @staticmethod
    def dump(args, filename):
        pickle.dump(args, open(filename, 'wb'))

    @staticmethod
    def log(args, func=logging.info):
        logging.info('%s version %s' % (Args.PROG, __version__))
        spacer = '-' * 78
        fmt = '{:<35}{:<35}'
        func(spacer)
        func(fmt.format('verbose logging:', str(args.verbose)))
        func(fmt.format('num processes:', str(args.num_processes)))
        func(fmt.format('output directory:', str(args.output_dir)))
        func(fmt.format('permutations:', str(args.perms)))
        func(fmt.format('weight_miss:', str(args.weight_miss)))
        func(fmt.format('weight_hit:', str(args.weight_hit)))
        func(fmt.format('weight_param:', str(args.weight_param)))
        func(fmt.format('noise_loc:', str(args.noise_loc)))
        func(fmt.format('noise_scale:', str(args.noise_scale)))
        func(fmt.format('sample_set_file:', str(args.sample_set_file)))
        func(fmt.format('count_matrix_path:', str(args.count_matrix_path)))

    @staticmethod
    def parse():
        # parse command line arguments
        parser = Args.create()
        args = parser.parse_args()
        # process and check arguments
        if os.path.exists(args.output_dir):
        	parser.error("Output directory '%s' already exists" %
                         args.output_dir)
        args.num_processes = max(1, args.num_processes)
        args.perms = max(1, args.perms)
        # check weight methods
        if isinstance(args.weight_miss, basestring):
            args.weight_miss = WEIGHT_METHODS[args.weight_miss]
        if isinstance(args.weight_hit, basestring):
            args.weight_hit = WEIGHT_METHODS[args.weight_hit]
        if args.weight_param < 0.0:
            parser.error('weight param < 0.0 invalid')
        elif ((args.weight_miss == 'log' or args.weight_hit == 'log')):
            if args.weight_param < 1.0:
                parser.error('weight param %f < 1.0 not allowed with '
                             'log methods' % (args.weight_param))
        # sample set file
        if not os.path.exists(args.sample_set_file):
            parser.error('sample set file "%s" not found' %
                        (args.sample_set_file))
        # count matrix path
        if not os.path.exists(args.count_matrix_path):
            parser.error('count matrix path "%s" not found' %
                         (args.count_matrix_path))
        args.count_matrix_path = os.path.abspath(args.count_matrix_path)
        return args


class Results(object):
    TMP_DIR = 'tmp'
    LOG_DIR = 'log'
    STATUS_FILE = 'status.json'
    ARGS_FILE = 'args.pickle'
    SAMPLE_SET_FILE = 'sample_set.json'

    def __init__(self, output_dir):
        self.output_dir = output_dir
        self.tmp_dir = os.path.join(output_dir, Results.TMP_DIR)
        self.log_dir = os.path.join(output_dir, Results.LOG_DIR)
        self.args_file = os.path.join(output_dir, Results.ARGS_FILE)
        self.sample_set_file = os.path.join(output_dir, Results.SAMPLE_SET_FILE)


class Run(object):

    def log_args(self):
        Args.log(self.args)

    @staticmethod
    def create():
        self = Run()
        # parse command line args
        args = Args.parse()
        self.args = args
        self.results = Results(args.output_dir)

        # setup logging
        if args.verbose:
            level = logging.DEBUG
        else:
            level = logging.INFO
        logging.basicConfig(level=level,
                            format="%(asctime)s pid=%(process)d "
                                   "%(levelname)s - %(message)s")
        # open count matrix
        cm = CountMatrix.open(args.count_matrix_path)
        shape = cm.shape

        # open sample sets
        logging.info("Reading sample set")
        sample_set = SampleSet.parse_smt(args.sample_set_file)[0]
        # ensure every sample found in count matrix
        if not set(sample_set.value_dict).issubset(cm.colnames):
            raise ParserError('samples in sample set do not match sample '
                              'names in count matrix')
        # sample sets must have at least one hit and one miss
        if len(sample_set) < 2:
            raise ParserError("sample set size <2")
        nhits = sum(x == 1 for x in sample_set.value_dict.itervalues())
        nmisses = sum(x == 0 for x in sample_set.value_dict.itervalues())
        if nhits < 1 or nmisses < 1:
            raise ParserError("sample set invalid")
        logging.debug("Sample set %s size %d" %
                      (sample_set.name, len(sample_set)))

        # create output directories
        results = self.results
        if not os.path.exists(results.output_dir):
            logging.debug("Creating output directory '%s'" %
                          (results.output_dir))
            os.makedirs(results.output_dir)
        if not os.path.exists(results.tmp_dir):
            logging.debug("Creating tmp directory '%s'" % (results.tmp_dir))
            os.makedirs(results.tmp_dir)
        if not os.path.exists(results.log_dir):
            logging.debug("Creating log directory '%s'" % (results.log_dir))
            os.makedirs(results.log_dir)

        # write command line args
        Args.dump(args, results.args_file)
        # write sample set file
        with open(results.sample_set_file, 'w') as fp:
            print >>fp, sample_set.to_json()

        return self
