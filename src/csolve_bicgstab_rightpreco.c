#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#include "csolver.h"

//Reference: Jie Chen, 2016 JSC, Right preconditioned/Flexible BiCGStab paper
void csolve_bicgstab_rightpreco(int n, complex *x, complex *b, cop_t Aop, cop_t invMop, int niter, double tol, int verb)
{
  int i, k;
  double rs0, rs;
  complex rho_old, rho_new, alpha, beta, omega;

  complex *r = malloc(n*sizeof(complex));
  complex *r0 = malloc(n*sizeof(complex));
  complex *p = malloc(n*sizeof(complex));
  complex *q = malloc(n*sizeof(complex));
  complex *v = malloc(n*sizeof(complex));
  complex *s = malloc(n*sizeof(complex));
  complex *t = malloc(n*sizeof(complex));
  
  Aop(n, x, v);//v=Ax
  for(i=0; i<n; i++) {
    r[i] = b[i]-v[i];//r=b-Ax
    p[i] = r[i];
    r0[i] = r[i];
  }
  rho_old = cdotprod(n, r, r0);
  rs = creal(cdotprod(n, r, r));
  rs0 = rs;
  
  for(k=0; k<niter; k++){
    if(verb) printf("bicgstab k=%d |r|=%e\n", k, sqrt(rs));

    invMop(n, p, q);//q=invM*p
    Aop(n, q, v);//v=Aq
    alpha = cdotprod(n, v, r0);
    alpha = rho_old/alpha;
    for(i=0; i<n; i++) s[i] = r[i] - alpha*v[i];
    
    invMop(n, s, r);//r=invM*s
    Aop(n, r, t);//t=Ar
    omega = cdotprod(n, s, t)/creal(cdotprod(n, t, t));

    for(i=0; i<n; i++){
      x[i] += alpha*q[i] + omega*r[i];
      r[i] = s[i] - omega*t[i];
    }
    rs = creal(cdotprod(n, r, r));
    if(sqrt(rs)<tol*sqrt(rs0)) {
      if(verb) printf("converged at k=%d\n", k);
      break;
    }

    rho_new = cdotprod(n, r, r0);
    beta = (rho_new/rho_old)*alpha/omega;
    for(i=0; i<n; i++) p[i] = r[i] + beta*(p[i]-omega*v[i]);

    rho_old = rho_new;
  }

  free(r);
  free(r0);
  free(p);
  free(v);
  free(s);
  free(t);
}
