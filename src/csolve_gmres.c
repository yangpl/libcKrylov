#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#include "csolver.h"


//GMRES without preconditioning
void csolve_gmres(int n, complex *x, complex *b, cop_t Aop, int niter, double tol, int m, int verb)
{
  int i, j, k, iter;

  complex *w = malloc(n*sizeof(complex));
  complex *r = malloc(n*sizeof(complex));
  complex *v = malloc(n*(m+1)*sizeof(complex));//n1=n, n2=m+1
  complex *h = malloc(m*(m+1)*sizeof(complex));//n1=m, n2=m+1
  complex *g = malloc((m+1)*sizeof(complex));
  complex *y = malloc((m+1)*sizeof(complex));
  complex *c = malloc((m+1)*sizeof(complex));
  double *s = malloc((m+1)*sizeof(double));
  double beta, tmp, r0;
  complex ss;
  
  for(iter=0; iter<niter; iter++){
    /*---------------------------------------------------*/
    Aop(n, x, w);//w=A*x
    for(i=0; i<n; i++) r[i] = b[i] -w[i];
    beta = sqrt(cdotprod(n, r, r));
    for(i=0; i<n; i++) v[i] = r[i]/beta;
    memset(g, 0, (m+1)*sizeof(complex));
    g[0] = beta;
    memset(h, 0, (m+1)*m*sizeof(complex));
    if(iter==0) r0 = beta;
    else if(beta<tol*r0) return;
    if(verb) printf("iter=%d |error|=%e\n", iter, beta);
 
    for(j=0; j<m; j++){
      Aop(n, &v[j*n], w);//r=Av;
      for(i=0; i<=j; i++) {
	h[i*m+j] = cdotprod(n, w, &v[i*n]);
	for(k=0; k<n; k++) w[k] -= h[i*m+j]*v[i*n+k];
      }
      h[(j+1)*m+j] = sqrt(creal(cdotprod(n, w, w)));

      if(cabs(h[(j+1)*m+j])==0.0) { m=j+1; break; }
      for(i=0; i<n; i++) v[(j+1)*n+i] = w[i]/h[(j+1)*m+j];

      //solve least-squares problem by QR factorization using Givens rotations
      //min \|g - H(1:m,1:m) y\|^2, g=(beta,0,...,0)
      //min \|G*g - G*Hy\|^2, G=G(i-1,i,theta)*...*G(2,3,theta)*G(1,2,theta)
      if(j>0){
	for(i=0; i<j; i++){
	  //apply G12, G23,..., G_{j-1,j} to the last column of H_{j,*}
	  ss = conj(c[i])*h[i*m+j] + s[i]*h[(i+1)*m+j];
	  h[(i+1)*m+j] = -s[i]*h[i*m+j] + c[i]*h[(i+1)*m+j];
	  h[i*m+j] = ss;
	}
      }

      //compute c=cos(theta) and s=sin(theta)
      tmp = sqrt(h[j*m+j]*conj(h[j*m+j]) + h[(j+1)*m+j]*conj(h[(j+1)*m+j]));
      s[j] = creal(h[(j+1)*m+j])/tmp;
      c[j] = h[j*m+j]/tmp;
      h[j*m+j] = conj(c[j])*h[j*m+j] + s[j]*h[(j+1)*m+j];
      h[(j+1)*m+j] = 0.;
      //g=G(j,j+1,theta)g with g[j+1]=0
      g[j+1] = -s[j]*g[j];
      g[j] = conj(c[j])*g[j];

      tmp = cabs(g[j+1]);
      if(tmp<tol*beta) { m=j+1; break; }
      //else printf("j=%d error=%e\n", j, tmp); 
    }

    //now, H becomes an upper triangule matrix, problem min\|g-Hy\|^2 is g=Hy
    //solve it by backward substitution, y = H(1:m,1:m)\g(1:m)
    y[m-1] = g[m-1]/h[(m-1)*m+m-1];
    for(i=m-2; i>=0; i--){
      y[i] = g[i];
      for(j=i+1; j<m; j++) y[i] -= h[i*m+j]*y[j];
      y[i] /= h[i*m+i];
    }

    //x=x0+Vm*y
    for(j=0; j<m; j++){
      for(i=0; i<n; i++){
	x[i] += y[j]*v[j*n+i];
      }
    }      
  }//end for iter
  
  free(r);
  free(w);
  free(v);
  free(h);
  free(g);
  free(y);
  free(c);
  free(s);
}
