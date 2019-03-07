#include <math.h>
#include <stdio.h>

/* Based on Grigg and Hughes (2018), https://onlinelibrary.wiley.com/doi/full/10.1002/hyp.13282 */
/* justin.hughes@csiro.au */

void sma_gr8j(double *P, 
              double *E, 
			  double *LAI,
              int *n,
              double *x1,
              double *mem,
              double *x5,
              double *x6,
              double *x7,
			  double *x8,
			  double *avLAI,
              double *S_0,
              double *U, 
              double *S, 
              double *ET){
  
  int t;
  double Pn, En, Ps,  perc, sCheck, Sp;
  double S_prev = (*S_0) * (*x1);
  double Smax = *mem * *x1;
  
  
  for (t = 0; t < *n; t++) {
    Pn = fmax(P[t] - E[t], 0);
    En = fmax(E[t] - P[t], 0);
    
    
    // production
    Ps = 0;
    ET[t] = 0;
    
    
    if (Pn > 0) {
      // part of Pn fills the production store
      Sp = fmax(0, S_prev);
      Ps = ( *x1 * (1.0 - pow((Sp/ *x1), *x5)) * tanh(Pn / (*x1)) /
        (1.0 + (Sp/ *x1) * tanh(Pn / (*x1))) );
    }
    
    
    if (En > 0) {
      // actual evaporation
      sCheck = fmax(0, (Smax - *x1 + S_prev));
      ET[t] = En * *x6 * pow((sCheck/Smax), *x7) * pow(LAI[t]/ *avLAI, *x8);
    }
    
    // initial update production store
    S[t] = S_prev - ET[t] + Ps;
    Sp = fmax(0, S[t]);
    
    
    // percolation leakage
    perc = Sp * ( 1.0 - pow(1.0 + pow((4.0/9.0) * Sp / (*x1), 4.0), -0.25) );
    S[t] = S[t] - perc;
    U[t] = perc + (Pn - Ps);
    S_prev = S[t];
  }
  
  
}

