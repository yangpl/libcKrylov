/* CGNE method
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2020-2022 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "csolver.h"

/*< CGNE (Craig's Method), see Saad's book algorithm 8.5 >*/
void csolve_cgne(int n, complex *x, complex *b, cop_t Aop, cop_t Atop, int niter, double tol, int verb)
{
  int i, iter;
  double rsold, rsnew, ps, alpha, beta;
  complex *r, *p, *Ap;

  Ap = malloc(n*sizeof(complex));
  r = malloc(n*sizeof(complex));
  p = malloc(n*sizeof(complex));

  Aop(n, x, Ap);
  for(i=0; i<n; i++) r[i] = b[i]-Ap[i];
  Atop(n, r, p);
  rsold = creal(cdotprod(n, r, r));
  
  for(iter=0; iter<niter; iter++){
    if(verb) printf("iter=%d error=%e\n", iter, sqrt(rsold));
    ps = creal(cdotprod(n, p, p));
    alpha = rsold/ps;
    Aop(n, p, Ap);
    for(i=0; i<n; i++) {
      x[i] += alpha*p[i];
      r[i] -= alpha*Ap[i];
    }
    rsnew = creal(cdotprod(n, r, r));
    if(sqrt(rsnew)<tol*sqrt(rsold)) break;
    beta = rsnew/rsold;
    Atop(n, r, Ap);//Ap=At * r
    for(i=0; i<n; i++) p[i] = Ap[i] + beta*p[i];
    rsold = rsnew;
  }

  free(Ap);
  free(r);
  free(p);
 
}
