#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>

#include "rng.h"

long* base_norm_counts;

void swap(double *p,double *q) {
   double t;

   t=*p;
   *p=*q;
   *q=t;
}

void sort(double a[],int n) {
   int i,j,temp;

   for(i=0;i<n-1;i++) {
      for(j=0;j<n-i-1;j++) {
         if(a[j]>a[j+1])
            swap(&a[j],&a[j+1]);
      }
   }
}

static PyObject* py_median_center_counts(PyObject* self, PyObject* args) {

    PyObject* counts_numpy_array;
    PyObject* size_factors_numpy_array;
    int counts_arraylen;
    int random_seed_int;
    int resample_flag;
    int add_noise_flag;
    double noise_loc;
    double noise_scale;

    if (!PyArg_ParseTuple(args, "OIOIIIdd", &counts_numpy_array, &counts_arraylen, &size_factors_numpy_array, &random_seed_int, &resample_flag, &add_noise_flag, &noise_loc, &noise_scale)) {
        return NULL;
    }

    double* counts = (double*)PyArray_DATA(counts_numpy_array);


    double counts_clone[counts_arraylen];
    for(int i = 0; i < counts_arraylen; i++) {
        counts_clone[i] = counts[i];
    }

    sort(counts_clone,counts_arraylen);
    int n = (counts_arraylen+1) / 2 - 1;
    double med = counts_clone[n];

    double* size_factors = (double*)PyArray_DATA(size_factors_numpy_array);

    for (int i = 0; i < counts_arraylen; i++) {
          counts[i] = fabs(log2(counts[i]/med));
        }


    return Py_BuildValue("i", random_seed_int);
}

static PyObject* py_normalize_counts(PyObject* self, PyObject* args) {

    PyObject* counts_numpy_array;
    PyObject* size_factors_numpy_array;
    int counts_arraylen;
    int random_seed_int;
    int resample_flag;
    int add_noise_flag;
    double noise_loc;
    double noise_scale;

    if (!PyArg_ParseTuple(args, "OIOIIIdd", &counts_numpy_array, &counts_arraylen, &size_factors_numpy_array, &random_seed_int, &resample_flag, &add_noise_flag, &noise_loc, &noise_scale)) {
        return NULL;
    }

    double* counts = (double*)PyArray_DATA(counts_numpy_array);


    double* size_factors = (double*)PyArray_DATA(size_factors_numpy_array);

    for (int i = 0; i < counts_arraylen; i++) {
        if (resample_flag == 1) {
            counts[i] = lcg_poisson(&random_seed_int, counts[i]);
        }

        counts[i] = counts[i] / size_factors[i];

        if (add_noise_flag == 1) {
            counts[i] += lcg_poisson(&random_seed_int, noise_loc);
            counts[i] += lcg_double(&random_seed_int) * noise_scale;
        }
    }


    return Py_BuildValue("i", random_seed_int);
}

static PyObject* py_power_transform(PyObject* self, PyObject* args) {

    PyObject* norm_counts_numpy_array;
    int norm_counts_arraylen;

    PyObject* output_numpy_array;

    int method;
    double param;
    int UNWEIGHTED;
    int WEIGHTED;
    int EXP;
    int LOG;

    if (!PyArg_ParseTuple(args, "OIOIdIIII", &norm_counts_numpy_array, &norm_counts_arraylen, &output_numpy_array, &method, &param, &UNWEIGHTED, &WEIGHTED, &EXP, &LOG)) {
        return NULL;
    }

    double* norm_counts = (double*)PyArray_DATA(norm_counts_numpy_array);
    double* output_array = (double*)PyArray_DATA(output_numpy_array);

    if (method == UNWEIGHTED) {
        for (int i = 0; i < norm_counts_arraylen; i++) {
            output_array[i] = 1.0;
        }
    } else if (method == WEIGHTED) {
        for (int i = 0; i < norm_counts_arraylen; i++) {
            output_array[i] = norm_counts[i];
        }
    } else if (method == EXP) {
        for (int i = 0; i < norm_counts_arraylen; i++) {
            output_array[i] = pow(norm_counts[i], param);
        }
    } else if (method == LOG) {
        for (int i = 0; i < norm_counts_arraylen; i++) {
            if (norm_counts[i] < 0) {
                output_array[i] = -log2(fabs(norm_counts[i]) + param);
            } else {
                output_array[i] = log2(fabs(norm_counts[i]) + param);
            }
        }
    } else {
        return NULL; // Error
    }

    return Py_BuildValue("i", 0);
}

static PyObject* py_shufflei(PyObject* self, PyObject* args) {

    PyObject* input_numpy_array;
    int arraylen;
    int random_seed_int;

    if (!PyArg_ParseTuple(args, "OII", &input_numpy_array, &arraylen, &random_seed_int)) {
        return NULL;
    }

    long* array = (long*)PyArray_DATA(input_numpy_array);

    for (int i = arraylen - 1; i >=0; i--) {
        const int j = lcg_range(&random_seed_int, 0, i);

        // Intentionally not using XOR swap to allow compiler to optimize as appropriate
        long temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }

    return Py_BuildValue("i", random_seed_int);
}

static PyObject* py_random_walk(PyObject* self, PyObject* args) {

    PyObject* weights_miss_numpy_array;
    PyObject* weights_hit_numpy_array;
    PyObject* es_run_numpy_array;
    PyObject* membership_numpy_array;
    PyObject* ranks_numpy_array;
    PyObject* perm_numpy_array;
    int nsamples;

    if (!PyArg_ParseTuple(args, "OOIOOOO", &weights_miss_numpy_array, &weights_hit_numpy_array, &nsamples, &es_run_numpy_array, &membership_numpy_array, &ranks_numpy_array, &perm_numpy_array)) {
        return NULL;
    }

    double* weights_miss = (double*)PyArray_DATA(weights_miss_numpy_array);
    double* weights_hit = (double*)PyArray_DATA(weights_hit_numpy_array);
    double* es_run = (double*)PyArray_DATA(es_run_numpy_array);

    long* membership = (long*)PyArray_DATA(membership_numpy_array);
    long* ranks = (long*)PyArray_DATA(ranks_numpy_array);
    long* perm = (long*)PyArray_DATA(perm_numpy_array);

    if (nsamples == 0) {
        return Py_BuildValue("(f,I,O)", 0.0, 0, es_run_numpy_array);
    }

    double* phit = (double*)malloc(sizeof(double) * nsamples);
    double* pmiss= (double*)malloc(sizeof(double) * nsamples);

    for (int i = 0; i < nsamples; i++) {
        const long r = ranks[i];
        const double wt_miss = weights_miss[r];
        const double wt_hit = weights_hit[r];

        const long p = perm[r];

        if (i == 0) {
            phit[i] = 0.0;
            pmiss[i] = 0.0;
        } else {
            phit[i] = phit[i-1];
            pmiss[i] = pmiss[i-1];
        }

        if (membership[p] == 1) {
            phit[i] += wt_hit;
        } else {
            pmiss[i] += wt_miss;
        }
    }

    const long last = nsamples - 1;

    if ((phit[last] > 0) || (pmiss[last] > 0)) {
        if (phit[last] == 0) {

            free(phit);
            free(pmiss);
            for (int i = 0; i < nsamples; i++) {
                es_run[i] = -1.0;
            }

            return Py_BuildValue("(f,I,O)", -1.0, last, es_run_numpy_array);
        } else if (pmiss[last] == 0) {

            free(phit);
            free(pmiss);
            for (int i = 0; i < nsamples; i++) {
                es_run[i] = 1.0;
            }

            return Py_BuildValue("(f,I,O)", 1.0, 0, es_run_numpy_array);
        } else {

            float es_val = 0.0;
            int es_rank = 0;

            for (int i = 0; i < nsamples; i++) {
                es_run[i] = (phit[i] / phit[last]) - (pmiss[i] / pmiss[last]);
                if (fabs(es_run[i]) >= fabs(es_val)) {
                    es_val = es_run[i];
                    es_rank = i;
                }
            }

            free(phit);
            free(pmiss);
            return Py_BuildValue("(f,I,O)", es_val, es_rank, es_run_numpy_array);
        }
    }

    free(phit);
    free(pmiss);
    return Py_BuildValue("(f,I,O)", 0.0, 0, es_run_numpy_array);
}

static int reverse_argsort_compare(const void* a, const void* b) {
    const long aa = *((long *)a);
    const long bb = *((long *)b);

    if (base_norm_counts[aa] < base_norm_counts[bb]) {
        return 1;
    } else if (base_norm_counts[aa] == base_norm_counts[bb]) {
        return 0;
    } else {
        return -1;
    }
}

static PyObject* py_argsort(PyObject* self, PyObject* args) {
    PyObject* norm_counts_numpy_array;
    PyObject* ranks_numpy_array;
    int array_len;

    if (!PyArg_ParseTuple(args, "OIO", &norm_counts_numpy_array, &array_len, &ranks_numpy_array)) {
        return NULL;
    }

    long* norm_counts = PyArray_DATA(norm_counts_numpy_array);
    long* ranks = PyArray_DATA(ranks_numpy_array);

    base_norm_counts = norm_counts;
    qsort(ranks, array_len, sizeof(long), reverse_argsort_compare);

    return Py_BuildValue("O", ranks_numpy_array);
}

static PyMethodDef CCkernelMethods[] = {
    {"c_normalize_counts", py_normalize_counts, METH_VARARGS},
    {"c_power_transform", py_power_transform, METH_VARARGS},
    {"c_shufflei", py_shufflei, METH_VARARGS},
    {"c_random_walk", py_random_walk, METH_VARARGS},
    {"c_argsort", py_argsort, METH_VARARGS},
    {"c_median_center_counts", py_median_center_counts, METH_VARARGS},
    {NULL, NULL}
};

PyMODINIT_FUNC
initckernel(void)
{
    (void) Py_InitModule("ckernel", CCkernelMethods);
    import_array();
}
