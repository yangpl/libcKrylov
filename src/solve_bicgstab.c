/* - BiCGStab (Bi-conjugate gradient with stabalization)
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

//linear solver using BiCGStab
void solve_bicgstab(int n, double *x, double *b, op_t Aop, int niter, double tol, int verb)
{
  int i, k;
  double rs0, rs, rho_old, rho_new, alpha, beta, omega;

  double *r = malloc(n*sizeof(double));
  double *r0 = malloc(n*sizeof(double));//rprime0
  double *p = malloc(n*sizeof(double));
  double *v = malloc(n*sizeof(double));
  double *s = malloc(n*sizeof(double));
  double *t = malloc(n*sizeof(double));
  
  Aop(n, x, v);//v=Ax
  for(i=0; i<n; i++) {
    r[i] = b[i]-v[i];//r=b-Ax
    r0[i] = r[i];
  }
  rho_old = 1.;
  alpha = 1.;
  omega = 1.;
  memset(p, 0, n*sizeof(double));
  memset(v, 0, n*sizeof(double));
  rs = dotprod(n, r, r);
  rs0 = rs;
  
  for(k=0; k<niter; k++){
    if(verb) printf("BiCGStab, k=%d rs=%e\n", k, rs);

    rho_new = dotprod(n, r0, r);
    beta = (rho_new/rho_old)*alpha/omega;
    for(i=0; i<n; i++) p[i] = r[i] + beta*(p[i]-omega*v[i]);

    Aop(n, p, v);
    alpha = rho_new/dotprod(n, r0, v);
    for(i=0; i<n; i++) s[i] = r[i] - alpha*v[i];

    Aop(n, s, t);
    omega = dotprod(n, t, s)/dotprod(n, t, t);

    for(i=0; i<n; i++){
      x[i] += alpha*p[i] + omega*s[i];
      r[i] = s[i] - omega*t[i];
    }
    rs = dotprod(n, r, r);
    if(rs<tol*rs0) {
      if(verb) printf("converged at k=%d\n", k);
      break;
    }
    rho_old = rho_new;
  }

  free(r);
  free(r0);
  free(p);
  free(v);
  free(s);
  free(t);
}
