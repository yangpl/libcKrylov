#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "ccoo.h"

ccoo_t *A;

void Acoo_init(ccoo_t *A_, float *v0, int n1, int n2, int ndamp, float dz, float dx, float omega)
{
  int i, j, k;
  float xi1, xi1mh, xi1ph;
  float xi2, xi2mh, xi2ph;
  float tmp, damp1, damp2;
  float scale = 40.;
  static float PI = 3.141592653589793238462643 ;
  
  k = 0;
  for(j=0; j<n2; j++){
    for(i=0; i<n1; i++){
      k++;
      if(j-1>=0){
	k++;
      }
      if(j+1<n2){
	k++;
      }
      if(i-1>=0){
	k++;
      }
      if(i+1<n1){
	k++;
      }
    }
  }
  A = A_;
  A->nnz = k;
  A->nrow = n1*n2;
  A->ncol = n1*n2;
  A->row_ind = malloc(A->nnz*sizeof(int));
  A->col_ind = malloc(A->nnz*sizeof(int));
  A->val = malloc(A->nnz*sizeof(complex));
  
  k = 0;
  for(j=0; j<n2; j++){
    damp2 = 0;
    if(j<ndamp) damp2 = scale*cos(0.5*PI*(j+1)/(ndamp+1));
    if(j>=n2-ndamp) damp2 = scale*cos(0.5*PI*(n2-j)/(ndamp+1));
    xi2 = 1. - I*damp2/omega;

    damp2 = 0;
    if(j<ndamp) damp2 = scale*cos(0.5*PI*(j+0.5)/(ndamp+1));
    if(j>=n2-ndamp) damp2 = scale*cos(0.5*PI*(n2-j-0.5)/(ndamp+1));
    xi2mh = 1. - I*damp2/omega;
    
    damp2 = 0;
    if(j<ndamp) damp2 = scale*cos(0.5*PI*(j+1.5)/(ndamp+1));
    if(j>=n2-ndamp) damp2 = scale*cos(0.5*PI*(n2-j+0.5)/(ndamp+1));
    xi2ph = 1. - I*damp2/omega;
    
    for(i=0; i<n1; i++){
      damp1 = 0;
      if(i<ndamp) damp1 = scale*cos(0.5*PI*(i+1)/(ndamp+1));
      if(i>=n1-ndamp) damp1 = scale*cos(0.5*PI*(n1-i)/(ndamp+1));
      xi1 = 1. - I*damp1/omega;

      damp1 = 0;
      if(i<ndamp) damp1 = scale*cos(0.5*PI*(i+0.5)/(ndamp+1));
      if(i>=n1-ndamp) damp1 = scale*cos(0.5*PI*(n1-i-0.5)/(ndamp+1));
      xi1mh = 1. - I*damp1/omega;

      damp1 = 0;
      if(i<ndamp) damp1 = scale*cos(0.5*PI*(i+1.5)/(ndamp+1));
      if(i>=n1-ndamp) damp1 = scale*cos(0.5*PI*(n1-i+0.5)/(ndamp+1));
      xi1ph = 1. - I*damp1/omega;

      //---------------------------------------------------
      tmp = omega/v0[j*n1+i];
      A->val[k] = -tmp*tmp + (1./xi2mh + 1./xi2ph)/(xi2*dx*dx) + (1./xi1mh + 1./xi1ph)/(xi1*dz*dz);
      A->row_ind[k] = i + n1*j;
      A->col_ind[k] = i + n1*j;
      k++;
      if(j-1>=0){
	A->val[k] = -1./(xi2mh*xi2*dx*dx);
	A->row_ind[k] = i + n1*j;
	A->col_ind[k] = i + n1*(j-1);
	k++;
      }
      if(j+1<n2){
	A->val[k] = -1./(xi2ph*xi2*dx*dx);
	A->row_ind[k] = i + n1*j;
	A->col_ind[k] = i + n1*(j+1);
	k++;
      }
      if(i-1>=0){
	A->val[k] = -1./(xi1mh*xi1*dz*dz);
	A->row_ind[k] = i + n1*j;
	A->col_ind[k] = i-1 + n1*j;
	k++;
      }
      if(i+1<n1){
	A->val[k] = -1./(xi1ph*xi1*dz*dz);
	A->row_ind[k] = i + n1*j;
	A->col_ind[k] = i+1 + n1*j;
	k++;
      }
    }
  }

}

void Acoo_close(ccoo_t *A)
{
  free(A->val);
  free(A->row_ind);
  free(A->col_ind);
}

//compute y=Ax by applying operator A in coo format
void Acoo_apply(int n, complex *x, complex *y)
{
  int k;

  memset(y, 0, n*sizeof(complex));  
  /* y_i += a_{ij}*x_i where i=A->row[k], j=A->col[k], a_{ij}=A->val[k] */
  for(k=0; k<A->nnz; k++)
    y[A->row_ind[k]] += A->val[k]*x[A->col_ind[k]];
}
  
