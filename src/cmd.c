// hydromad: Hydrological Modelling and Analysis of Data
//
// Copyright (c) Felix Andrews <felix@nfrac.org>


#include <R_ext/Arith.h>

#include <R.h>

#ifndef min
#define min(a, b) ((a < b)?(a):(b))
#define max(a, b) ((a < b)?(b):(a))
#endif

#define my_isok(x) (!ISNA(x) & !ISNAN(x))


void
sma_cmd(double *P, double *E, int *n,
	double *d, double *g, double *e, double *M_0, 
	double *U, double *M, double *ET)
{
	int t;
	double Mf;
	double M_prev = *M_0;
	for (t = 0; t < *n; t++) {
		// rainfall reduces CMD (Mf)
		if (P[t] > 0) {
			if (M_prev < *d) {
				Mf = M_prev * exp(-P[t] / *d);
			} else if (M_prev < *d + P[t]) {
				Mf = *d * exp((-P[t] + M_prev - *d) / *d);
			} else {
				Mf = M_prev - P[t];
			}
		} else {
			Mf = M_prev;
		}
		// drainage (rainfall not accounted for in -dM)
		U[t] = max(0, P[t] - M_prev + Mf);
		// evapo-transpiration
		ET[t] = *e * E[t] * min(1, exp(2 * (1 - Mf / *g)));
		ET[t] = max(0, ET[t]);
		// mass balance
		M[t] = M_prev - P[t] + U[t] + ET[t];
		M[t] = max(0, M[t]);
		M_prev = M[t];
	}
}

