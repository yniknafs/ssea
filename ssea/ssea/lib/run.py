'''
SSEA: Sample Set Enrichment Analysis
'''
import os
import argparse
import logging

from base import WEIGHT_METHODS
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
        parser.add_argument('-s', '--sample-sets', required=True,
                            dest='sample_set_file',
                            help='File containing sample sets')
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
        func(fmt.format('sample_sets:', str(args.sample_set_file)))
        func(fmt.format('count_matrix:', str(args.count_matrix_path)))

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
    STATUS_FILE = 'status.json'
    ARGS_FILE = 'args.pickle'
    SAMPLE_SET_FILE = 'sample_sets.tsv'

    def __init__(self, output_dir):
        self.output_dir = output_dir
        self.tmp_dir = os.path.join(output_dir, Results.TMP_DIR)
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
        # setup logging
        if args.verbose:
            level = logging.DEBUG
        else:
            level = logging.INFO
        logging.basicConfig(level=level,
                            format="%(asctime)s pid=%(process)d "
                                   "%(levelname)s - %(message)s")

        self.args = args
        self.results = Results(args.output_dir)

        # create output directories
        results = self.results
        if not os.path.exists(results.output_dir):
            logging.debug("Creating output directory '%s'" %
                          (results.output_dir))
            os.makedirs(results.output_dir)
        if not os.path.exists(results.tmp_dir):
            logging.debug("Creating tmp directory '%s'" % (results.tmp_dir))
            os.makedirs(results.tmp_dir)

        # open count matrix
        cm = CountMatrix.open(args.count_matrix_path)
        if cm.size_factors is None:
            parser.error('No size factors found in count matrix')
        shape = cm.shape
        cm.close()

        return self

        # open sample sets
        sample_sets = read_sample_sets(config)
        # read current file if exists
        sample_set_info_file = os.path.join(config.output_dir,
                                            Config.SAMPLE_SET_INFO_FILE)
        ss_infos = []
        ss_names = set()
        max_ss_id = 0
        if os.path.exists(sample_set_info_file):
            for ss_info in read_sample_set_info(sample_set_info_file):
                ss_infos.append(ss_info)
                ss_names.add(ss_info.name)
                max_ss_id = max(max_ss_id, ss_info.id)
        # read new sample sets
        for sample_set in sample_sets:
            logging.info("Sample Set: %s" % (sample_set.name))
            if sample_set.name in ss_names:
                logging.warning('Sample Set "%s" already exists in output '
                                'directory' % (sample_set.name))
                continue
            # create sample set directory
            ss_id = (max_ss_id + 1)
            sample_set_dirname = 'ss_%d' % (ss_id)
            sample_set_path = os.path.join(config.output_dir, sample_set_dirname)
            assert not os.path.exists(sample_set_path)
            max_ss_id += 1
            if not os.path.exists(sample_set_path):
                logging.debug("Creating sample set directory '%s'" % (sample_set_path))
                os.makedirs(sample_set_path)
            # create temp directory
            tmp_dir = os.path.join(sample_set_path, Config.TMP_DIR)
            if not os.path.exists(tmp_dir):
                logging.debug("Creating tmp directory '%s'" % (tmp_dir))
                os.makedirs(tmp_dir)
            # create log directory
            log_dir = os.path.join(sample_set_path, Config.LOG_DIR)
            if not os.path.exists(log_dir):
                logging.debug("Creating log directory '%s'" % (log_dir))
                os.makedirs(log_dir)
            # write configuration file
            logging.debug("Writing configuration file")
            config_file = os.path.join(sample_set_path, Config.CONFIG_JSON_FILE)
            with open(config_file, 'w') as fp:
                print >>fp, config.to_json()
            # write sample set file
            logging.debug("Writing sample set file")
            sample_set_file = os.path.join(sample_set_path, Config.SAMPLE_SET_JSON_FILE)
            with open(sample_set_file, 'w') as fp:
                print >>fp, sample_set.to_json()
            # map work into chunks
            worker_file = os.path.join(sample_set_path, WORKER_FILE)
            write_worker_info(worker_file,
                              num_jobs=shape[0],
                              num_processes=config.num_processes)
            # add to sample set info list
            ss_info = SSInfo(id=ss_id,
                             dirname=sample_set_dirname,
                             name=sample_set.name,
                             desc=sample_set.desc,
                             size=len(sample_set))
            ss_infos.append(ss_info)
            # mark job status as ready
            JobStatus.set(sample_set_path, JobStatus.READY)
        # write sample set info list
        write_sample_set_info(ss_infos, sample_set_info_file)
        return 0
