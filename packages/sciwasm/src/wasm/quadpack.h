#ifndef QUADPACK_H
#define QUADPACK_H

void dqagse(double(*fcn)(double* x), const double a, const double b,
            const double epsabs, const double epsrel, const int limit,
            double* result, double* abserr, int* neval, int* ier,
            double* alist, double* blist, double* rlist, double* elist,
            int* iord, int* last);

void dqagie(double(*fcn)(double* x), const double bound, const int inf,
            const double epsabs, const double epsrel, const int limit,
            double* result, double* abserr, int* neval, int* ier,
            double* alist, double* blist, double* rlist, double* elist,
            int* iord, int* last);

#endif
