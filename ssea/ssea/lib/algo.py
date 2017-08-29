'''
Created on Oct 9, 2013

@author: mkiyer
'''
import os
import logging
import shutil
from collections import namedtuple
from multiprocessing import Process
import numpy as np

# local imports
import ssea.lib.cfisher as fisher
from ssea.lib.batch_sort import batch_sort, batch_merge
from ssea.lib.kernel import ssea_kernel, RandomState
from ssea.lib.base import Result, interp
from ssea.lib.countdata import BigCountMatrix

# results directories
TMP_DIR = 'tmp'

# temporary files
JSON_UNSORTED_SUFFIX = '.unsorted.json'
JSON_SORTED_SUFFIX = '.sorted.json'
NPY_HISTS_SUFFIX = '.hists.npz'

# improves readability of code
KernelResult = namedtuple('KernelResult', ('ranks',
                                           'norm_counts',
                                           'norm_counts_miss',
                                           'norm_counts_hit',
                                           'es_val',
                                           'es_rank',
                                           'es_run',
                                           'perms'))

# batch sort configuration
SORT_BUFFER_SIZE = 32000

# floating point precision to output in report
FLOAT_PRECISION = 8
SCIENTIFIC_NOTATION_PRECISION = 200

# to compute FDR q values we need to store normalized enrichment
# scores (NES). we accomplish this storing a histogram of NES values during
# runtime and then merging histograms from parallel processes prior to
# computing global statistics. the NES can range from [-inf, +inf] making
# it difficult for a histogram to store enough granularity in few bins.
# we accomplish this here by using log-space histogram bins that vary from
# 0.1 to 100. values outside this range are clipped and stored so that no
# values are lost
NUM_NES_BINS = 10001
NES_BINS = np.logspace(-1,2,num=NUM_NES_BINS,base=10)
LOG_NES_BINS = np.log10(NES_BINS)
NES_MIN = NES_BINS[0]
NES_MAX = NES_BINS[-1]


def _init_hists():
    '''returns a dictionary with histogram arrays initialized to zero'''
    return {'null_nes_pos': np.zeros(NUM_NES_BINS-1, dtype=np.float),
            'null_nes_neg': np.zeros(NUM_NES_BINS-1, dtype=np.float),
            'obs_nes_pos': np.zeros(NUM_NES_BINS-1, dtype=np.float),
            'obs_nes_neg': np.zeros(NUM_NES_BINS-1, dtype=np.float)}


def _cmp_json_nes(line):
    '''comparison function for batch_sort'''
    res = Result.from_json(line.strip())
    return abs(res.nes)

def ssea_kernel_debug(counts, size_factors, membership, rng, config):
    '''
    counts: numpy array of float values
    size_factors: normalization factors for counts
    membership: int array (0 or 1) with set membership
    rng: RandomState object
    config: Config object
    '''
    # run kernel to generate a range of observed enrichment scores
    resample_rand_seeds = np.empty(config.resampling_iterations, dtype=np.int)
    resample_count_ranks = np.empty((config.resampling_iterations, counts.shape[0]), dtype=np.int)
    resample_es_vals = np.zeros(config.resampling_iterations, dtype=np.float)
    resample_es_ranks = np.zeros(config.resampling_iterations, dtype=np.int)
    # save random number generator seed before running kernel
    resample_rand_seeds[0] = rng.seed
    k = ssea_kernel(counts, size_factors, membership, rng,
                    resample_counts=True,
                    permute_samples=False,
                    add_noise=True,
                    noise_loc=config.noise_loc,
                    noise_scale=config.noise_scale,
                    method_miss=config.weight_miss,
                    method_hit=config.weight_hit,
                    method_param=config.weight_param)
    k = KernelResult._make(k)
    resample_count_ranks[0] = k.ranks
    resample_es_vals[0] = k.es_val
    resample_es_ranks[0] = k.es_rank

def ssea_kernel_debug_run(config, sample_set, output_basename,
                startrow=None, endrow=None):
    '''
    main SSEA loop (single processor)

    config: Config object
    sample_set: SampleSet object
    output_basename: prefix for writing result files
    '''
    # initialize random number generator
    rng = RandomState()
    # open data matrix
    bm = BigCountMatrix.open(config.matrix_dir)
    # determine range of matrix to process
    if startrow is None:
        startrow = 0
    if endrow is None:
        endrow = bm.shape[0]
    assert startrow < endrow
    # get membership array for sample set
    membership = sample_set.get_array(bm.colnames)
    valid_samples = (membership >= 0)
    # setup histograms
    hists = _init_hists()
    # setup report file
    unsorted_json_file = output_basename + JSON_UNSORTED_SUFFIX
    outfileh = open(unsorted_json_file, 'wb')
    for i in xrange(startrow, endrow):
        logging.debug("\tRow: %d (%d-%d)" % (i, startrow, endrow))
        # read from memmap
        counts = np.array(bm.counts[i,:], dtype=np.float)
        # remove 'nan' values
        valid_inds = np.logical_and(valid_samples, np.isfinite(counts))
        # subset counts, size_factors, and membership array
        counts = counts[valid_inds]
        size_factors = bm.size_factors[valid_inds]
        valid_membership = membership[valid_inds]
        # write dummy results for invalid rows
        if (valid_inds.sum() == 0) or (np.all(counts == 0)):
            res = Result()
        else:
            # run ssea
            ssea_kernel_debug(counts, size_factors, valid_membership, rng, config)
        break

def ssea_run(counts, size_factors, membership, rng, config):
    '''
    counts: numpy array of float values
    size_factors: normalization factors for counts
    membership: int array (0 or 1) with set membership
    rng: RandomState object
    config: Config object
    '''
    # run kernel to generate a range of observed enrichment scores
    resample_rand_seeds = np.empty(config.resampling_iterations, dtype=np.int)
    resample_count_ranks = np.empty((config.resampling_iterations, counts.shape[0]), dtype=np.int)
    resample_es_vals = np.zeros(config.resampling_iterations, dtype=np.float)
    resample_es_ranks = np.zeros(config.resampling_iterations, dtype=np.int)
    for i in xrange(config.resampling_iterations):
        # save random number generator seed before running kernel
        resample_rand_seeds[i] = rng.seed
        k = ssea_kernel(counts, size_factors, membership, rng,
                        resample_counts=True,
                        permute_samples=False,
                        add_noise=True,
                        noise_loc=config.noise_loc,
                        noise_scale=config.noise_scale,
                        method_miss=config.weight_miss,
                        method_hit=config.weight_hit,
                        method_param=config.weight_param,
                        median_center=config.median_center)
        k = KernelResult._make(k)
        resample_count_ranks[i] = k.ranks
        resample_es_vals[i] = k.es_val
        resample_es_ranks[i] = k.es_rank
    # find the median ES value
    median_index = int(config.resampling_iterations / 2)
    median_index = resample_es_vals.argsort()[median_index]
    median_es_val = resample_es_vals[median_index]
    # choose whether to use the positive or negative side of the
    # distribution based on the median ES value
    if median_es_val == 0:
        return Result(), np.zeros((0,), dtype=np.float)
    elif median_es_val < 0:
        signfunc = np.less
    else:
        signfunc = np.greater
    # subset to include only the corresponding side of the distribution
    resample_sign_inds = signfunc(resample_es_vals, 0)
    resample_rand_seeds = resample_rand_seeds[resample_sign_inds]
    resample_count_ranks = resample_count_ranks[resample_sign_inds]
    resample_es_vals = resample_es_vals[resample_sign_inds]
    resample_es_ranks = resample_es_ranks[resample_sign_inds]
    # determine the median value
    median_index = int(resample_es_vals.shape[0] / 2)
    median_index = resample_es_vals.argsort()[median_index]
    # select the kernel run representing the median ES value
    rand_seed = resample_rand_seeds[median_index]
    es_val = resample_es_vals[median_index]
    es_rank = resample_es_ranks[median_index]
    ranks = resample_count_ranks[median_index]
    es_sign = cmp(es_val, 0)
    # permute samples and determine ES null distribution
    null_es_vals = np.zeros(config.perms, dtype=np.float)
    null_es_ranks = np.zeros(config.perms, dtype=np.float)
    i = 0
    j = 0
    while i < config.perms:
        j += 1
        k = ssea_kernel(counts, size_factors, membership, rng,
                        resample_counts=True,
                        permute_samples=True,
                        add_noise=True,
                        noise_loc=config.noise_loc,
                        noise_scale=config.noise_scale,
                        method_miss=config.weight_miss,
                        method_hit=config.weight_hit,
                        method_param=config.weight_param,
                        median_center=config.median_center)
        k = KernelResult._make(k)
        # logging.debug('es_val: %f, es_sign: %d' % (k.es_val, es_sign))
        # print 'ES SCORES', resample_es_vals
        # print 'ESVAL', k.es_val, 'PVAL', k, k.es_rank
        # print 'hi', k, k.es_rank
        # if j > 1000:
        #     null_es_vals[i] = 0.01 * es_sign
        #     null_es_ranks[i] = 0.01 * es_sign
        #     i += 1
        # else:
        if cmp(k.es_val, 0) == es_sign:
            null_es_vals[i] = k.es_val
            null_es_ranks[i] = k.es_rank
            i += 1
    # Subset the null ES scores to only positive or negative
    null_sign_inds = signfunc(null_es_vals, 0)
    null_es_vals = null_es_vals[null_sign_inds]
    null_es_ranks = null_es_ranks[null_sign_inds]
    # Adjust for variation in gene set size. Normalize ES(S,null)
    # and the observed ES(S), separately rescaling the positive and
    # negative scores by dividing by the mean of the ES(S,null) to
    # yield normalized scores NES(S,null)
    null_es_mean = np.fabs(null_es_vals.mean())
    null_nes_vals = null_es_vals / null_es_mean
    # Normalize the observed ES(S) by rescaling by the mean of
    # the ES(S,null) separately for positive and negative ES(S)
    nes_val = es_val / null_es_mean
    # estimate nominal p value for S from ES(S,null) by using the
    # positive or negative portion of the distribution corresponding
    # to the sign of the observed ES(S)
    p_value = (np.fabs(null_es_vals) >= np.fabs(es_val)).sum().astype(np.float)
    p_value /= null_es_vals.shape[0]

    # create result object for this SSEA test
    r = Result()
    r.rand_seed = int(rand_seed)
    r.es = round(es_val, FLOAT_PRECISION)
    r.es_rank = int(es_rank)
    r.nominal_p_value = round(p_value, SCIENTIFIC_NOTATION_PRECISION)
    r.nes = round(nes_val, FLOAT_PRECISION)

    # save resampled es range
    i = resample_es_vals.argmin()
    r.es_min = resample_es_vals[i]
    r.es_rank_min = resample_es_ranks[i]
    i = resample_es_vals.argmax()
    r.es_max = resample_es_vals[i]
    r.es_rank_max = resample_es_ranks[i]
    # save null distribution points
    r.null_es_mean = null_es_mean
    i = null_es_vals.argmin()
    r.null_es_min = null_es_vals[i]
    r.null_es_rank_min = null_es_ranks[i]
    i = null_es_vals.argmax()
    r.null_es_max = null_es_vals[i]
    r.null_es_rank_max = null_es_ranks[i]

    # membership arrays
    rank_inds = ranks.argsort()
    if es_val < 0:
        m_core = rank_inds >= es_rank
    else:
        m_core = rank_inds <= es_rank
    m_null = np.logical_not(m_core)
    m_hits = (membership > 0)
    m_misses = np.logical_not(m_hits)
    m_core_hits = np.logical_and(m_core, m_hits)

    # leading edge stats
    core_hits = m_core_hits.sum()
    core_misses = np.logical_and(m_core, m_misses).sum()
    null_hits = np.logical_and(m_null, m_hits).sum()
    null_misses = np.logical_and(m_null, m_misses).sum()
    num_hits = core_hits + null_hits
    num_misses = core_misses + null_misses

    # fold change
    norm_counts = counts / size_factors
    core_hits_mean = norm_counts[m_core_hits].mean() if core_hits > 0 else 0.0
    hits_mean = norm_counts[m_hits].mean() if num_hits > 0 else 0.0
    misses_mean = norm_counts[m_misses].mean() if num_misses > 0 else 1.0
    core_fold_change = np.log2((core_hits_mean + 1e-5 )/ (misses_mean + 1e-5))
    fold_change = np.log2((hits_mean + 1e-5) /( misses_mean + 1e-5))


    # fisher exact test (one-sided hypothesis that LE is enriched)
    fisher_p_value = fisher.pvalue(core_hits, core_misses, null_hits, null_misses).right_tail
    # odds ratio
    n = np.inf if null_hits == 0 else float(core_hits) / null_hits
    d = np.inf if null_misses == 0 else float(core_misses) / null_misses
    if np.isfinite(n) and np.isfinite(d):
        if n == 0 and d == 0:
            odds_ratio = np.nan
        else:
            odds_ratio = np.inf if d == 0 else (n / d)
    elif np.isfinite(d):
        odds_ratio = np.inf
    else:
        odds_ratio = np.nan if n == 0 else 0.0

    # create dictionary result
    r.core_hits = int(core_hits)
    r.core_misses = int(core_misses)
    r.null_hits = int(null_hits)
    r.null_misses = int(null_misses)
    r.fisher_p_value = np.round(fisher_p_value, SCIENTIFIC_NOTATION_PRECISION)
    r.odds_ratio = odds_ratio
    r.fold_change = fold_change
    r.core_fold_change = core_fold_change
    # return result and null distribution for subsequent fdr calculations
    return r, null_nes_vals


def ssea_serial(config, sample_set, output_basename,
                startrow=None, endrow=None):
    '''
    main SSEA loop (single processor)

    config: Config object
    sample_set: SampleSet object
    output_basename: prefix for writing result files
    '''
    # initialize random number generator
    rng = RandomState()
    # open data matrix
    bm = BigCountMatrix.open(config.matrix_dir)
    # determine range of matrix to process
    if startrow is None:
        startrow = 0
    if endrow is None:
        endrow = bm.shape[0]
    assert startrow < endrow
    # get membership array for sample set
    membership = sample_set.get_array(bm.colnames)
    valid_samples = (membership >= 0)
    # setup histograms
    hists = _init_hists()
    # setup report file
    unsorted_json_file = output_basename + JSON_UNSORTED_SUFFIX
    outfileh = open(unsorted_json_file, 'wb')
    for i in xrange(startrow, endrow):
        logging.debug("\tRow: %d (%d-%d)" % (i, startrow, endrow))
        # read from memmap
        counts = np.array(bm.counts[i,:], dtype=np.float)
        # remove 'nan' values
        valid_inds = np.logical_and(valid_samples, np.isfinite(counts))
        # subset counts, size_factors, and membership array
        counts = counts[valid_inds]
        size_factors = bm.size_factors[valid_inds]
        valid_membership = membership[valid_inds]
        # write dummy results for invalid rows
        if (valid_inds.sum() == 0) or (np.all(counts == 0)):
            res = Result()
        else:
            # run ssea
            res, null_nes_vals = ssea_run(counts, size_factors,
                                          valid_membership, rng, config)
            # update histograms
            null_keys = []
            obs_keys = []
            if res.es < 0:
                null_keys.append('null_nes_neg')
                obs_keys.append('obs_nes_neg')
            elif res.es > 0:
                null_keys.append('null_nes_pos')
                obs_keys.append('obs_nes_pos')
            for k in xrange(len(null_keys)):
                null_nes = np.clip(np.fabs(null_nes_vals), NES_MIN, NES_MAX)
                obs_nes = np.clip(np.fabs(res.nes), NES_MIN, NES_MAX)
                hists[null_keys[k]] += np.histogram(null_nes, NES_BINS)[0]
                hists[obs_keys[k]] += np.histogram(obs_nes, NES_BINS)[0]
        # save t_id
        res.t_id = i
        res.name = bm.rownames[i]
        # convert to json and write
        print >>outfileh, res.to_json()
    # close report file
    outfileh.close()
    # save histograms to a file
    output_hist_file = output_basename + NPY_HISTS_SUFFIX
    np.savez(output_hist_file, **hists)
    # cleanup
    bm.close()
    # sort output json file by abs(NES)
    logging.debug("Worker %s: sorting results" % (output_basename))
    # make tmp dir for sorting
    if os.path.exists(output_basename):
        shutil.rmtree(output_basename)
    os.makedirs(output_basename)
    # call batch sort python function
    sorted_json_file = output_basename + JSON_SORTED_SUFFIX
    batch_sort(input=unsorted_json_file,
               output=sorted_json_file,
               key=_cmp_json_nes,
               buffer_size=SORT_BUFFER_SIZE,
               tempdirs=[output_basename])
    # remove tmp dir
    shutil.rmtree(output_basename)
    # remove unsorted json file
    os.remove(unsorted_json_file)
    logging.debug("Worker %s: done" % (output_basename))
    return 0

def ssea_map(config, sample_set, worker_basenames, worker_chunks):
    '''
    parallel map step of SSEA run

    processes chunks of input matrix in parallel using multiprocessing
    '''
    # start worker processes
    procs = []
    for i in xrange(len(worker_basenames)):
        basename = worker_basenames[i]
        startrow, endrow = worker_chunks[i]
        args = (config, sample_set, basename, startrow, endrow)
        p = Process(target=ssea_serial, args=args)
        p.start()
        procs.append(p)
    # join worker processes (wait for processes to finish)
    for p in procs:
        p.join()
    return 0


def compute_qvalues(json_iterator, hists_file):
    '''
    computes fdr q values from json Result objects sorted
    by abs(NES) (low to high)

    json_iterator: iterator that yields json objects in sorted order
    hists_file: contains histogram data from null distribution
    '''
    # load histogram data
    hists = np.load(hists_file)
    # compute cumulative sums for fdr interpolation
    cdfs = {}
    for k in ('null_nes_neg', 'null_nes_pos', 'obs_nes_neg', 'obs_nes_pos'):
        h = hists[k]
        cdf = np.zeros(h.shape[0]+1, dtype=np.float)
        cdf[1:] = h.cumsum()
        cdfs[k] = cdf
    # keep track of minimum FDR and rank for positive
    # and negative NES separately
    NEG = 0
    POS = 1
    null_keys = ['null_nes_neg', 'null_nes_pos']
    obs_keys = ['obs_nes_neg', 'obs_nes_pos']
    tot_obs = [cdfs['obs_nes_neg'][-1], cdfs['obs_nes_pos'][-1]]
    cur_ranks = [tot_obs[0], tot_obs[1]]
    min_fdrs = [1.0, 1.0]
    # perform merge of sorted json files
    for line in json_iterator:
        # load json document (one per line)
        res = Result.from_json(line.strip())
        es = res.es
        log_nes_clip = np.log10(np.clip(abs(res.nes), NES_MIN, NES_MAX))
        if es != 0:
            if es < 0:
                sign_ind = NEG
                sign = -1.0
            else:
                sign_ind = POS
                sign = 1.0
            # For a given NES(S) = NES* >= 0, the FDR is the ratio of the
            # percentage of all permutations NES(S,null) >= 0, whose
            # NES(S,null) >= NES*, divided by the percentage of observed S with
            # NES(S) >= 0, whose NES(S) >= NES*, and similarly for
            # NES(S) = NES* <= 0.
            # to compute a sample set specific FDR q value we look at the
            # aggregated enrichment scores for all tests of that sample set
            # compute the cumulative sums of NES histograms use interpolation
            # to find fraction NES(null) >= NES* and account for the observed
            # permutation in the null set interpolate NES in log space
            null_nes_cumsum = cdfs[null_keys[sign_ind]]
            null_n = interp(log_nes_clip, LOG_NES_BINS, null_nes_cumsum)
            obs_nes_cumsum = cdfs[obs_keys[sign_ind]]
            obs_n = interp(log_nes_clip, LOG_NES_BINS, obs_nes_cumsum)
            n = 1.0 - (null_n / null_nes_cumsum[-1])
            d = 1.0 - (obs_n / obs_nes_cumsum[-1])
            #print 'SS_ID=%d ES=%f NES=%f n=%f (%f / %f) d=%f (%f / %f)' % (i, res.es, res.nes, n, null_n, null_nes_cumsum[-1], d, obs_n, obs_nes_cumsum[-1])
            # update json dict
            if (n <= 0.0) or (d <= 0.0):
                res.ss_fdr_q_value = 0.0
            else:
                res.ss_fdr_q_value = n / d
            #print 'SS_ID=%d ES=%f NES=%f fdr=%f minfdr=%f' % (i, res.es, res.nes, res.ss_fdr_q_value, min_fdrs[i])
            # compare with minimum FDR and adjust minimum FDR if necessary
            if res.ss_fdr_q_value < min_fdrs[sign_ind]:
                min_fdrs[sign_ind] = res.ss_fdr_q_value
            else:
                res.ss_fdr_q_value = min_fdrs[sign_ind]
            res.ss_rank = cur_ranks[sign_ind]
            res.ss_frac = sign * (1.0 - ((res.ss_rank - 1) / float(tot_obs[sign_ind])))
            cur_ranks[sign_ind] -= 1
        # convert back to json
        yield res.to_json()
        yield os.linesep
    # cleanup
    hists.close()


def ssea_reduce(input_basenames, output_json_file, output_hist_file):
    '''
    reduce step of SSEA run

    merges the null distribution histograms from individual worker
    processes from the map step, and then merges the sorted results
    json files from individual worker processes and computes global
    fdr q-value statistics
    '''
    # merge NES histograms
    logging.debug("Merging %d worker histograms" % (len(input_basenames)))
    hists = _init_hists()
    json_iterables = []
    for i in xrange(len(input_basenames)):
        # create sorted json file streams
        json_file = input_basenames[i] + JSON_SORTED_SUFFIX
        json_fileh = open(json_file, 'rb', 64*1024)
        json_iterables.append(json_fileh)
        # aggregate numpy arrays
        hist_file = input_basenames[i] + NPY_HISTS_SUFFIX
        npzfile = np.load(hist_file)
        for k in hists.iterkeys():
            hists[k] += npzfile[k]
        npzfile.close()
    np.savez(output_hist_file, **hists)
    # perform merge of sorted json files
    try:
        with open(output_json_file, 'wb', 64*1024) as output:
            iterator = batch_merge(_cmp_json_nes, *json_iterables)
            output.writelines(compute_qvalues(iterator, output_hist_file))
    finally:
        for iterable in json_iterables:
            try:
                iterable.close()
            except Exception:
                pass
    # remove worker files
    logging.debug("Removing temporary files")
    for i in xrange(len(input_basenames)):
        hist_file = input_basenames[i] + NPY_HISTS_SUFFIX
        os.remove(hist_file)
        json_file = input_basenames[i] + JSON_SORTED_SUFFIX
        os.remove(json_file)
    return 0

def parse_results(filename):
    with open(filename, 'r') as fp:
        for line in fp:
            result = Result.from_json(line.strip())
            yield result

def print_tsv(input_json_file, output_tsv_file):
    headero = ['id', 'log2FC', 'log2FC_le', 'es', 'nes', 'pval', 'fdr', 'percentile']
    with open(output_tsv_file, 'wb', 64*1024) as f:
        print >>f, '\t'.join(headero)
        for r in parse_results(input_json_file):
            lineo = [r.name, r.fold_change, r.core_fold_change, r.es,
                    r.nes, r.nominal_p_value, r.ss_fdr_q_value, r.ss_frac]
            print >>f, '\t'.join(map(str,lineo))
