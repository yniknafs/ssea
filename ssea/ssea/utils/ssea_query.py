'''
Created on Jan 30, 2014

@author: mkiyer
'''
import logging
import argparse
import os
import sys
import glob
from multiprocessing import Pool

# local imports
from ssea.lib.config import Config 
from ssea.lib.countdata import BigCountMatrix
from ssea.lib.base import Result, JobStatus

# # functions to test operator
# OP_TEST_FUNCS = {'<=': lambda a,b: a<=b,
#                  '>=': lambda a,b: a>=b,
#                  '<': lambda a,b: a<b,
#                  '>': lambda a,b: a>b}
#         # parse threshold arguments of the form 'attribute,value'
#         # for example: nominal_p_value,0.05
#         if args.thresholds is not None:
#             for arg in args.thresholds:
#                 m = re.match(r'(.+)([<>]=?)(.+)', arg)
#                 if m is None:
#                     parser.error('error parsing threshold argument "%s"' % (arg))
#                 attr, op, value = m.groups()
#                 if attr not in Result.FIELDS:
#                     parser.error('threshold attribute "%s" unknown' % (attr))
#                 if op not in OP_TEST_FUNCS:
#                     parser.error('unrecognized operator "%s"' % (op))
#                 value = float(value)
#                 self.thresholds.append((attr,op,value))       
#             
# def parse_and_filter_results(filename, thresholds):
#     passed = 0
#     failed = 0    
#     with open(filename, 'r') as fin:
#         for line in fin:
#             # load json document (one per line)
#             result = Result.from_json(line.strip())      
#             # apply thresholds
#             skip = False
#             for attr,op,threshold in thresholds:
#                 value = getattr(result, attr)
#                 if not OP_TEST_FUNCS[op](value, threshold):
#                     skip = True
#                     break
#             if skip:
#                 failed += 1
#                 continue
#             passed += 1
#             yield result
#     logging.info("Parsed %d results, %d passed and %d failed" %
#                  (passed+failed, passed, failed))
    
    
def check_path(path):
    if not os.path.exists(path):
        logging.error('Directory not found at path "%s"' % (path))
        return False
    elif not os.path.isdir(path):
        logging.error('Directory not valid at path "%s"' % (path))
        return False
    elif JobStatus.get(path) != JobStatus.DONE:
        logging.info('Result status not set to "DONE" at path "%s"' % (path))
        return False
    return True

def find_paths(input_dirs):
    paths = set()
    for input_dir in input_dirs:
        if not os.path.exists(input_dir):
            logging.error("input dir '%s' not found" % (input_dir))
        if not os.path.isdir(input_dir):
            logging.error("input dir '%s' not a valid directory" % (input_dir))
        input_dir = os.path.abspath(input_dir)
        for path in glob.iglob(os.path.join(input_dir, "*")):
            if os.path.exists(path) and os.path.isdir(path):
                if check_path(path):
                    paths.add(path)
    return paths

def parse_results(filename):
    with open(filename, 'r') as fp:
        for line in fp:
            result = Result.from_json(line.strip())
            yield result

def query_worker(args):    
    (input_path, matrix_dir, meta, fdr_threshold, 
     frac_threshold, fpr_threshold, prec_threshold) = args
    bm = BigCountMatrix.open(matrix_dir)
    ss_compname = os.path.basename(input_path)
    results_file = os.path.join(input_path, Config.RESULTS_JSON_FILE)
    i = 0
    sig = 0
    lines = []
    # extract data
    for res in parse_results(results_file):
        i += 1
        #if (i % 10000) == 0:
        #    logging.debug('%s %d results' % (ss_compname, i))
        transcript_id = bm.rownames[res.t_id]
        if (meta is not None) and (transcript_id not in meta):
            continue
        core_size = res.core_hits + res.core_misses
        if core_size == 0:
            prec = 0.0
        else:
            prec = res.core_hits / float(core_size)
        num_misses = res.core_misses + res.null_misses
        if num_misses == 0:
            fpr = 0.0
        else:
            fpr = res.core_misses / float(num_misses)
        if ((res.ss_fdr_q_value <= fdr_threshold) and 
            (abs(res.ss_frac) >= frac_threshold) and
            (fpr <= fpr_threshold) and
            (prec >= prec_threshold)):
            if meta is None:
                fields = [transcript_id]
            else:
                fields = list(meta[transcript_id])
            fields.extend([ss_compname, res.es, res.nes, 
                           res.ss_fdr_q_value, res.ss_frac, fpr, prec])
            lines.append('\t'.join(map(str,fields)))
            sig += 1
    bm.close()
    logging.debug('Found %d results for path %s' % (sig, input_path))
    return lines

def query_parallel(input_paths, matrix_dir, meta, meta_header_fields, 
                   fdr_threshold, frac_threshold, fpr_threshold, 
                   prec_threshold, num_processes):
    tasklist = []
    for input_path in input_paths:
        tasklist.append((input_path, matrix_dir, meta, fdr_threshold, 
                         frac_threshold, fpr_threshold, prec_threshold))
    header_fields = meta_header_fields + ['ss_compname', 'es', 'nes', 'fdr', 'frac', 'fpr', 'prec']
    print '\t'.join(header_fields)
    # create pool
    pool = Pool(processes=num_processes)
    result_iter = pool.imap_unordered(query_worker, tasklist)
    for lines in result_iter:
        print '\n'.join(lines)
    pool.close()
    pool.join()

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('--fdr', dest='fdr', type=float, default=1.0)
    parser.add_argument('--frac', dest='frac', type=float, default=0.0)
    parser.add_argument('--fpr', dest='fpr', type=float, default=1.0)
    parser.add_argument('--prec', dest='prec', type=float, default=0.0)
    parser.add_argument('--meta', dest='metadata_file', default=None)
    parser.add_argument('-p', dest='num_processes', type=int, default=1)
    parser.add_argument('-i', dest='input_paths_file', default=None)
    parser.add_argument('matrix_dir')
    args = parser.parse_args()
    # get args
    matrix_dir = args.matrix_dir
    fdr_threshold = args.fdr
    frac_threshold = args.frac
    fpr_threshold = args.fpr
    prec_threshold = args.prec
    num_processes = args.num_processes
    input_paths_file = args.input_paths_file
    metadata_file = args.metadata_file
    # check args
    fdr_threshold = max(0.0, min(fdr_threshold, 1.0))
    frac_threshold = max(0.0, min(frac_threshold, 1.0))
    input_paths = []
    if input_paths_file is None:
        logging.error('No input directories specified (use -i).. Exiting.')
        return 1
    if not os.path.exists(input_paths_file):
        logging.error('Input paths file "%s" not found' % (input_paths_file))
    else:
        with open(input_paths_file) as fileh:
            for line in fileh:
                path = line.strip()
                if path in input_paths:
                    continue
                if check_path(path):
                    input_paths.append(path)
    if len(input_paths) == 0:
        logging.error('No valid SSEA results directories found.. Exiting.')
        return 1
    meta = None
    if metadata_file is not None:
        logging.debug('Parsing transcript metadata')
        meta = {}
        with open(metadata_file) as f:
            meta_header_fields = f.next().strip().split()
            for line in f:
                fields = line.strip().split('\t')
                meta[fields[0]] = fields
        logging.debug('Found metadata for %d transcripts' % (len(meta)))
    else:
        meta = None
        meta_header_fields = ['transcript_id']
    # parse results
    logging.debug('SSEA results: %d' % (len(input_paths)))
    logging.debug('FDR threshold: %f' % (fdr_threshold))
    logging.debug('Frac threshold: %f' % (frac_threshold))
    logging.debug('Prec threshold: %f' % (prec_threshold))
    logging.debug('FPR threshold: %f' % (fpr_threshold))
    logging.debug('Num processes: %d' % (num_processes))
    query_parallel(input_paths, matrix_dir, meta, meta_header_fields, 
                   fdr_threshold, frac_threshold, fpr_threshold, 
                   prec_threshold, num_processes)

if __name__ == '__main__':
    sys.exit(main())