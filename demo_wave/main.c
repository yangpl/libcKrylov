#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "csolver.h"
#include "ccoo.h"


void Acoo_init(ccoo_t *A, float *v0, int n1, int n2, int ndamp, float dz, float dx, float omega);
void Acoo_close(ccoo_t *A);
void Acoo_apply(int n, complex *x, complex *y);

int main()
{
  int n1 = 201;
  int n2 = 201;
  float dz = 25;
  float dx = 25;
  int ndamp = 12;
  float freq = 10.;
  float omega = 2*3.1415926*freq;
  float *v0;
  ccoo_t *A;
  
  int verb = 1;
  int method = 1;
  int niter = 400;
  int isz = n1/2;
  int isx = n2/2;
  float tol = 1e-6;
  int restart = 10;
  int i, j;
  complex *x, *b;

  v0 = malloc(n1*n2*sizeof(float));
  for(j=0; j<n2; j++)
    for(i=0; i<n1; i++)
      v0[j*n1+i] = 2500;

  A = malloc(sizeof(ccoo_t));
  Acoo_init(A, v0, n1, n2, ndamp, dz, dx, omega);

  x = malloc(A->nrow*sizeof(complex));
  b = malloc(A->nrow*sizeof(complex));
  memset(x, 0, n1*n2*sizeof(complex));
  b[isz + n1*isx] = 1.;

  if(method==1)
    csolve_bicgstab(n1*n2, x, b, Acoo_apply, niter, tol, verb);
  else if(method==2)  
    csolve_gmres(n1*n2, x, b, Acoo_apply, niter, tol, restart, verb);

  FILE *fp = fopen("wave.bin", "wb");
  for(j=0; j<n2; j++){
    for(i=0; i<n1; i++){
      float tmp = creal(x[i + n1*j]);
      fwrite(&tmp, sizeof(float), 1, fp);
    }
  }
  fclose(fp);

  Acoo_close(A);
  free(A);
  free(v0);
}
