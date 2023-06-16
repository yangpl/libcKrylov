#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "solver.h"


//Reference: Jie Chen, 2016 JSC, Right preconditioned/Flexible BiCGStab paper
void solve_bicgstab_rightpreco(int n, double *x, double *b, op_t Aop, op_t invMop, int niter, double tol, int verb)
{
  int i, k;
  double rs0, rs;
  double rho_old, rho_new, alpha, beta, omega;

  double *r = malloc(n*sizeof(double));
  double *r0 = malloc(n*sizeof(double));//rprime0
  double *p = malloc(n*sizeof(double));
  double *q = malloc(n*sizeof(double));
  double *v = malloc(n*sizeof(double));
  double *s = malloc(n*sizeof(double));
  double *t = malloc(n*sizeof(double));
  
  Aop(n, x, v);//v=Ax
  for(i=0; i<n; i++) {
    r[i] = b[i]-v[i];//r=b-Ax
    p[i] = r[i];
    r0[i] = r[i];
  }
  rho_old = dotprod(n, r0, r);
  rs = dotprod(n, r, r);
  rs0 = rs;
  
  for(k=0; k<niter; k++){
    if(verb) printf("iter=%d error=%e\n", k, sqrt(rs));

    invMop(n, p, q);//q=invM*p
    Aop(n, q, v);//v=Aq
    alpha = dotprod(n, r0, v);
    alpha = rho_old/alpha;
    for(i=0; i<n; i++) s[i] = r[i] - alpha*v[i];
    
    invMop(n, s, r);//r=invM*s
    Aop(n, r, t);//t=Ar
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

    rho_new = dotprod(n, r0, r);
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
