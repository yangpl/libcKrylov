/* CGNR method
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2020-2022 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "solver.h"

//CGNR method
void solve_cgnr(int n, double *x, double *b, op_t Aop, op_t Atop, int niter, double tol, int verb)
{
  int i, iter;
  double *r, *z, *p, *w;
  double zsold, zsnew, ws, rs, rs0, alpha, beta;
  
  r = malloc(n*sizeof(double));
  z = malloc(n*sizeof(double));
  p = malloc(n*sizeof(double));
  w = malloc(n*sizeof(double));
  
  Aop(n, x, w);//w=A*x
  for(i=0; i<n; i++) r[i] = b[i] - w[i];
  Atop(n, r, z);//z=At*r
  for(i=0; i<n; i++) p[i] = z[i];
  zsold = dotprod(n, z, z);

  for(iter=0; iter<niter; iter++){
    rs = dotprod(n, r, r);
    if(verb) printf("iter=%d error=%e\n", iter, sqrt(rs));
    if(iter==0) rs0 = rs;
    if(sqrt(rs)<tol*sqrt(rs0)) break;
    
    Aop(n, p, w);//w=Ap
    ws = dotprod(n, w, w);
    alpha = zsold/ws;
    for(i=0; i<n; i++) {
      x[i] += alpha*p[i];
      r[i] -= alpha*w[i];
    }
    Atop(n, r, z);
    zsnew = dotprod(n, z, z);
    beta = zsnew/zsold;
    for(i=0; i<n; i++) p[i] = z[i] + beta*p[i];
    zsold = zsnew;
  }  

  free(r);
  free(z);
  free(w);
  free(p);
}
