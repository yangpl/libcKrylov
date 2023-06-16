#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "solver.h"

double *r;
void A_toeplitz_init(double *r_)
{
  r = r_;
}

void A_toeplitz_apply(int n, double *x, double *y)
{
  int i, j;

  for(i=0; i<n; i++){
    y[i] = 0;
    for(j=0; j<n; j++) y[i] += r[i-j+n-1]*x[j];//Aij=r[i-j+n-1]
  }

}

void At_toeplitz_apply(int n, double *y, double *x)
{
  int i, j;

  for(j=0; j<n; j++){
    x[j] = 0;
    for(i=0; i<n; i++) x[j] += r[i-j+n-1]*y[i];//Aij=r[i-j+n-1]
  }
}

int main()
{
  int i;
  int n = 500;
  int method = 2;
  int verb = 1;
  int niter = 200;
  double tol = 1e-6;
  double *x, *b, *r;

  r = malloc((2*n-1)*sizeof(double));
  x = malloc(n*sizeof(double));
  b = malloc(n*sizeof(double));

  srand48(1001);
  for(i=0; i<2*n-1; i++) r[i] = drand48();
  for(i=0; i<n; i++) x[i] = drand48();

  A_toeplitz_init(r);
  A_toeplitz_apply(n, x, b);
  memset(x, 0, n*sizeof(double));
  
  if(method==1)
    solve_cgnr(n, x, b, A_toeplitz_apply, At_toeplitz_apply, niter, tol, verb);
  else if(method==2)
    solve_cgne(n, x, b, A_toeplitz_apply, At_toeplitz_apply, niter, tol, verb);
  
  //A_toeplitz_close();//not needed, nothing to do

  free(x);
  free(b);
  free(r);

}
