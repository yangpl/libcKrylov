#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#include "csolver.h"


//linear solver using BiCGStab, algorithm 7.7 in Saad book
void csolve_bicgstab(int n, complex *x, complex *b, op_t Aop, int niter, double tol, int verb)
{
  int i, k;
  double rs0, rs;
  complex rho_old, rho_new, alpha, beta, omega;

  complex *r = alloc1complex(n);
  complex *r0 = alloc1complex(n);//rprime0
  complex *p = alloc1complex(n);
  complex *v = alloc1complex(n);
  complex *s = alloc1complex(n);
  complex *t = alloc1complex(n);
  
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

    Aop(n, p, v);//v=Ap
    alpha = rho_old/cdotprod(n, v, r0);
    for(i=0; i<n; i++) s[i] = r[i] - alpha*v[i];

    Aop(n, s, t);//t=As
    omega = cdotprod(n, t, s)/creal(cdotprod(n, t, t));

    for(i=0; i<n; i++){
      x[i] += alpha*p[i] + omega*s[i];
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

  free1complex(r);
  free1complex(r0);
  free1complex(p);
  free1complex(v);
  free1complex(s);
  free1complex(t);
}

