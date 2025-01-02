/* - conjugate gradient method
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

#include "csolver.h"

//linear solver using conjugate gradient method
void csolve_cg(int n, complex *x, complex *b, cop_t Aop, int niter, double tol, int verb)
{
  int i, k;
  double rsold, rsnew, rs0, pAp, alp, beta;

  complex *r = malloc(n*sizeof(complex));
  complex *p = malloc(n*sizeof(complex));
  complex *Ap = malloc(n*sizeof(complex));

  Aop(n, x, Ap);//Ap=Ax0
  for(i=0; i<n; i++) {
    r[i] = b[i]-Ap[i];//r=b-Ax
    p[i] = r[i];
  }
  rsold = creal(cdotprod(n, r, r));
  rs0 = rsold;

  for(k=0; k<niter; k++){
    if(verb) printf("CG, k=%d rs=%e\n", k, rsold);

    Aop(n, p, Ap);//Ap=A*p
    pAp = creal(cdotprod(n, Ap, p));
    alp = rsold/pAp;

    for(i=0; i<n; i++){
      x[i] += alp*p[i];
      r[i] -= alp*Ap[i];
    }
    rsnew = creal(cdotprod(n, r, r));
    if(rsnew<tol*rs0) {
      if(verb) printf("converged at k=%d\n", k);
      break;
    }
    
    beta = rsnew/rsold;
    for(i=0; i<n; i++) p[i] = r[i] + beta*p[i];
    rsold = rsnew;
  }

  free(r);
  free(p);
  free(Ap);
}
