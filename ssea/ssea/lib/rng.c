/* file: rng.c */
/* linear congruential random number generator */
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "rng.h"

/* a.k.a. RAND_MAX */
#define MODULUS ((1U << 31) - 1)

/* seed generator using time */
int lcg_init_state()
{
    int x;
    x = ((int) time((time_t *) NULL)) % MODULUS;
    return(x);
}

/* integer pseudo-random number based on seed */
inline int lcg_rand(int seed)
{
    return (seed * 1103515245 + 12345) & MODULUS;
}

double lcg_double(int *seedp)
{
    int seed = lcg_rand(*seedp);
    *seedp = seed;
    return(((double) seed) / MODULUS);
}

int lcg_range(int *seedp, const int a, const int b)
{
    int seed = lcg_rand(*seedp);
    *seedp = seed;
    double x = ((double) seed) / ((unsigned int) MODULUS + 1);
    return(a + (int)((b - a + 1) * x));
}

/*
* copied from numpy.random source code
* log-gamma function to support some of these distributions. The
* algorithm comes from SPECFUN by Shanjie Zhang and Jianming Jin and their
* book "Computation of Special Functions", 1996, John Wiley & Sons, Inc.
*/
static double loggam(double x)
{
    double x0, x2, xp, gl, gl0;
    long k, n;

    static double a[10] = {8.333333333333333e-02,-2.777777777777778e-03,
         7.936507936507937e-04,-5.952380952380952e-04,
         8.417508417508418e-04,-1.917526917526918e-03,
         6.410256410256410e-03,-2.955065359477124e-02,
         1.796443723688307e-01,-1.39243221690590e+00};
    x0 = x;
    n = 0;
    if ((x == 1.0) || (x == 2.0))
    {
        return 0.0;
    }
    else if (x <= 7.0)
    {
        n = (long)(7 - x);
        x0 = x + n;
    }
    x2 = 1.0/(x0*x0);
    xp = 2*M_PI;
    gl0 = a[9];
    for (k=8; k>=0; k--)
    {
        gl0 *= x2;
        gl0 += a[k];
    }
    gl = gl0/x0 + 0.5*log(xp) + (x0-0.5)*log(x0) - x0;
    if (x <= 7.0)
    {
        for (k=1; k<=n; k++)
        {
            gl -= log(x0-1.0);
            x0 -= 1.0;
        }
    }
    return gl;
}

/* adapted from numpy.random source code */
long lcg_poisson_mult(int *seedp, const double lam)
{
    long X = 0;
    double enlam = exp(-lam);
    double prod = 1.0;
    double U;

    while (1)
    {
        U = lcg_double(seedp);
        prod *= U;
        if (prod > enlam)
        {
            X += 1;
        }
        else
        {
            /* update seed before returning */
            return X;
        }
    }
}

/* adapted from numpy.random source code */
#define LS2PI 0.91893853320467267
#define TWELFTH 0.083333333333333333333333
long lcg_poisson_ptrs(int *seedp, const double lam)
{
    const double slam = sqrt(lam);
    const double loglam = log(lam);
    const double b = 0.931 + 2.53*slam;
    const double a = -0.059 + 0.02483*b;
    const double invalpha = 1.1239 + 1.1328/(b-3.4);
    const double vr = 0.9277 - 3.6224/(b-2);

    while (1)
    {
        const double U = lcg_double(seedp) - 0.5;
        const double V = lcg_double(seedp);
        const double us = 0.5 - fabs(U);
        const long k = (long)floor((2*a/us + b)*U + lam + 0.43);
        if ((us >= 0.07) && (V <= vr))
        {
            return k;
        }
        if ((k < 0) ||
            ((us < 0.013) && (V > us)))
        {
            continue;
        }
        if ((log(V) + log(invalpha) - log(a/(us*us)+b)) <=
            (-lam + k*loglam - loggam(k+1)))
        {
            return k;
        }
    }
}

/* adapted from numpy.random source code */
long lcg_poisson(int *seedp, const double lam)
{
    if (lam >= 10)
    {
        return lcg_poisson_ptrs(seedp, lam);
    }
    else if (lam == 0)
    {
        return 0;
    }
    else
    {
        return lcg_poisson_mult(seedp, lam);
    }
}

double lcg_gauss(int *seedp)
{
    double f, x1, x2, r2;
    do {
        x1 = 2.0*lcg_double(seedp) - 1.0;
        x2 = 2.0*lcg_double(seedp) - 1.0;
        r2 = x1*x1 + x2*x2;
    }
    while (r2 >= 1.0 || r2 == 0.0);
    /* Box-Muller transform */
    f = sqrt(-2.0*log(r2)/r2);
    return f*x2;
}

double lcg_normal(int *seedp, double loc, double scale)
{
    return loc + scale*lcg_gauss(seedp);
}

double lcg_uniform(int *seedp, double loc, double scale)
{
    return loc + scale*lcg_double(seedp);
}

int main()
{
    int i;
    double d;
    int seed;
    double sum = 0;

    for (int i = 0; i < 100; i++) {
      // printf("sum:%d i:%d counts:%d \n ", sum, i, counts[i]);
      printf("sum:%f \n ", sum);
      sum += 1;
    };
    printf("%d", 4);

    // seed = lcg_init_state();
    // printf("seed is %d\n", seed);
    //
    // for (i = 0; i < 10; i++) {
    //     seed = lcg_rand(seed);
    //     d = lcg_poisson(&seed, 30);
    //     /* d = ((double) seed) / MODULUS; */
    //     printf("poisson %d %f\n", seed, d);
    //     d = lcg_normal(&seed, 50, 2);
    //     printf("normal %d %f\n", seed, d);
    //     d = lcg_double(&seed) * 1;
    //     printf("double %d %f\n", seed, d);
    // }


    return 0;
}
