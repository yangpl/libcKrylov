/* - GMRES with left and right preconditioning
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


//GMRES with right preconditioning
void solve_gmres_rightpreco(int n, double *x, double *b, op_t Aop, op_t invMop, int niter, double tol, int m, int verb)
{
  int i, j, k, iter;

  double *w = malloc(n*sizeof(double));
  double *r = malloc(n*sizeof(double));
  double *v = malloc(n*(m+1)*sizeof(double));
  double *h = malloc(m*(m+1)*sizeof(double));
  double *g = malloc((m+1)*sizeof(double));
  double *y = malloc((m+1)*sizeof(double));
  double *c = malloc((m+1)*sizeof(double));
  double *s = malloc((m+1)*sizeof(double));
  double eps = 1e-15;
  double beta, tmp;

  for(iter=0; iter<niter; iter++){
    /*---------------------------------------------------*/
    Aop(n, x, w);//w=A*x
    for(i=0; i<n; i++) r[i] = b[i]-w[i];
    beta = sqrt(dotprod(n, r, r));
    if(verb) printf("GMRES, iter=%d error=%e\n", iter, beta);
    for(i=0; i<n; i++) v[i] = r[i]/beta;
    memset(g, 0, (m+1)*sizeof(double));
    g[0] = beta;
    memset(h, 0, (m+1)*m*sizeof(double));
 
    for(j=0; j<m; j++){
      invMop(n, &v[j*n], r);//r=Minv*v_j
      Aop(n, r, w);//w=Ar;
      for(i=0; i<=j; i++) {
	h[i*m+j] = dotprod(n, &v[i*n], w);
	for(k=0; k<n; k++) w[k] -= h[i*m+j]*v[i*n+k];
      }
      h[(j+1)*m+j] = sqrt(dotprod(n, w, w));

      if(h[(j+1)*m+j]<eps) { m=j+1; break; }
      for(i=0; i<n; i++) v[(j+1)*m+i] = w[i]/h[(j+1)*m+j];

      //solve least-squares problem by QR factorization using Given rotations
      //min \|g - H(1:m,1:m) y\|^2, g=(beta,0,...,0)
      //min \|G*g - G*Hy\|^2, G=G(i-1,i,theta)*...*G(2,3,theta)*G(1,2,theta)
      if(j>0){
	for(i=0; i<j; i++){
	  //apply G12, G23,..., G_{j-1,j} to the last column of H_{j,*}
	  tmp = c[i]*h[i*m+j] + s[i]*h[(i+1)*m+j];
	  h[(i+1)*m+j] = -s[i]*h[i*m+j] + c[i]*h[(i+1)*m+j];
	  h[i*m+j] = tmp;
	}
      }

      //compute c=cos(theta) and s=sin(theta)
      if(fabs(h[j*m+j])>fabs(h[(j+1)*m+j])){
	tmp = h[(j+1)*m+j]/h[j*m+j];
	c[j] = 1./sqrt(1. + tmp*tmp);
	s[j] = c[j]*tmp;
      }else{
	tmp = h[j*m+j]/h[(j+1)*m+j];
	s[j] = 1./sqrt(1. + tmp*tmp);
	c[j] = s[j]*tmp;
      }
      h[j*m+j] = c[j]*h[j*m+j] + s[j]*h[(j+1)*m+j];
      h[(j+1)*m+j] = 0.;
      //g=G(j,j+1,theta)g with g[j+1]=0
      g[j+1] = -s[j]*g[j];
      g[j] = c[j]*g[j];

      tmp = fabs(g[j+1]);
      if(tmp<tol*beta) { m=j+1; break; }
      //else printf("j=%d error=%e\n", j, tmp); 
    }

    //now, H becomes an upper triangule matrix, problem min\|g-Hy\|^2 is g=Hy
    //solve it by backward substitution, y = H(1:m,1:m)\g(1:m)
    y[m-1] = g[m-1]/h[(m-1)*m + m-1];
    for(i=m-2; i>=0; i--){
      y[i] = g[i];
      for(j=i+1; j<m; j++) y[i] -= h[i*m+j]*y[j];
      y[i] /= h[i*m+i];
    }

    //r=Vm*y
    for(i=0; i<n; i++){
      r[i] = 0;
      for(j=0; j<m; j++){
	r[i] += y[j]*v[j*n+i];
      }
    }      
    invMop(n, r, w);//w=Minv*r
    for(i=0; i<n; i++) x[i] += w[i];//x=x0+w
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
