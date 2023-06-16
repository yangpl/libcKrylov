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

#include "solver.h"

/*< CGNE (Craig's Method), see Saad's book algorithm 8.5 >*/
void solve_cgne(int n, double *x, double *b, op_t Aop, op_t Atop, int niter, double tol, int verb)
{
  int i, iter;
  double rsold, rsnew, ps, alpha, beta;
  double *r, *p, *Ap;

  Ap = malloc(n*sizeof(double));
  r = malloc(n*sizeof(double));
  p = malloc(n*sizeof(double));

  Aop(n, x, Ap);
  for(i=0; i<n; i++) r[i] = b[i]-Ap[i];
  Atop(n, r, p);
  rsold = dotprod(n, r, r);
  
  for(iter=0; iter<niter; iter++){
    if(verb) printf("iter=%d error=%e\n", iter, sqrt(rsold));
    ps = dotprod(n, p, p);
    alpha = rsold/ps;
    Aop(n, p, Ap);
    for(i=0; i<n; i++) {
      x[i] += alpha*p[i];
      r[i] -= alpha*Ap[i];
    }
    rsnew = dotprod(n, r, r);
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
