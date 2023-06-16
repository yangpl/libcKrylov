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
#include <string.h>
#include <math.h>
#include <complex.h>

#include "csolver.h"

//CGNR method
void csolve_cgnr(int n, complex *x, complex *b, cop_t Aop, cop_t Atop, int niter, double tol, int verb)
{
  int i, iter;
  complex *r, *z, *p, *w;
  double zsold, zsnew, ws, rs, rs0, alpha, beta;
  
  r = malloc(n*sizeof(complex));
  z = malloc(n*sizeof(complex));
  p = malloc(n*sizeof(complex));
  w = malloc(n*sizeof(complex));
  
  Aop(n, x, w);//w=A*x
  for(i=0; i<n; i++) r[i] = b[i] - w[i];
  Atop(n, r, z);//z=At*r
  for(i=0; i<n; i++) p[i] = z[i];
  zsold = creal(cdotprod(n, z, z));

  for(iter=0; iter<niter; iter++){
    rs = creal(cdotprod(n, r, r));
    if(verb) printf("iter=%d error=%e\n", iter, sqrt(rs));
    if(iter==0) rs0 = rs;
    if(sqrt(rs)<tol*sqrt(rs0)) break;
    
    Aop(n, p, w);//w=Ap
    ws = creal(cdotprod(n, w, w));
    alpha = zsold/ws;
    for(i=0; i<n; i++) {
      x[i] += alpha*p[i];
      r[i] -= alpha*w[i];
    }
    Atop(n, r, z);
    zsnew = cdotprod(n, z, z);
    beta = zsnew/zsold;
    for(i=0; i<n; i++) p[i] = z[i] + beta*p[i];
    zsold = zsnew;
  }  

  free(r);
  free(z);
  free(w);
  free(p);
}
