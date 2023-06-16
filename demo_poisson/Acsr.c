#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct{
  int nrow;
  int ncol;
  int nnz;
  int *row_ptr;
  int *col_ind;
  double *val;
}csr_t;
csr_t *A;

double *diag;

//================================================================
//build matrix in coo format for laplace operator A=-Dx^2 - Dy^2
void Acsr_init(int nx, int ny, double dx, double dy)
{
  int i, j, k, row_ind;
  double _dx2, _dy2;
  
  //count the number of nonzeros in sparse banded matrix A
  k = 0;
  for(j=0; j<=ny; j++){
    for(i=0; i<=nx; i++){
      //the diagonal element A(i,j)
      k++;
      
      //the off-diagonal element
      if(j-1>=0){//element A(i,j-1)
	k++;	
      }
      if(i-1>=0){//element A(i-1,j)
	k++;	
      }
      if(i+1<=nx){//element A(i+1,j)
	k++;
      }
      if(j+1<=ny){//element A(i,j+1)
	k++;
      }      
    }
  }
  A = malloc(sizeof(csr_t));
  A->nnz = k;//number of non-zeros
  A->nrow = (nx+1)*(ny+1);
  A->ncol = (nx+1)*(ny+1);
  
  A->row_ptr = malloc((A->nrow+1)*sizeof(int));
  A->col_ind = malloc(A->nnz*sizeof(int));
  A->val = malloc(A->nnz*sizeof(double));
  diag = malloc(A->nrow*sizeof(double));

  _dx2 = 1./(dx*dx);
  _dy2 = 1./(dy*dy);
  k = 0;
  for(j=0; j<=ny; j++){
    for(i=0; i<=nx; i++){
      //the diagonal element A(i,j)
      A->val[k] = 2.*(_dx2 + _dy2);
      A->col_ind[k] = i + (nx+1)*j;
      k++;
      
      //the off-diagonal element
      if(i-1>=0){//element A(i-1,j)
	A->val[k] = -_dx2;
	A->col_ind[k] = (i-1) + (nx+1)*j;
	k++;	
      }
      if(i+1<=nx){//element A(i+1,j)
	A->val[k] = -_dx2;
	A->col_ind[k] = (i+1) + (nx+1)*j;
	k++;
      }
      if(j-1>=0){//element A(i,j-1)
	A->val[k] = -_dy2;
	A->col_ind[k] = i + (nx+1)*(j-1);
	k++;	
      }
      if(j+1<=ny){//element A(i,j+1)
	A->val[k] = -_dy2;
	A->col_ind[k] = i + (nx+1)*(j+1);
	k++;
      }
      
      row_ind = i + (nx+1)*j;
      A->row_ptr[row_ind+1] = k;
    }
  }
  
}


void Acsr_close()
{
  free(A->row_ptr);
  free(A->col_ind);
  free(A->val);
  free(A);
  free(diag);
}

//compute y=Ax by applying operator A in coo format
void Acsr_apply(int n, double *x, double *y)
{
  int i, j, k;

  for(i=0; i<A->nrow; i++){
    y[i] = 0;
    for(k=A->row_ptr[i]; k<A->row_ptr[i+1]; k++){
      j = A->col_ind[k];
      y[i] += A->val[k]*x[j];
    }
  }

}


/* //2. y=Aap^{-1}*x, LU-SGS for approximation of A^{-1} */
/* void sgs_apply(int n, double *x, double *y) */
/* { */
/*   int i, j, k; */
/*   double s; */
  
/*   y = x; */
/*   //------------------------------------------------------- */
/*   //xx=U^{-1}z */
/*   for(i=A->nrow-1; i>=0; i--){ */
/*     s = 0; */
/*     for(k=A->row_ptr[i]; k<A->row_ptr[i+1]; k++){ */
/*       j = A->col_ind[k]; */
/*       if(j==i) diag[i] = A->val[k]; */
/*       if(j>i) s += A->val[k]*y[j]; */
/*     } */
/*     y[i] = (x[i] - s)/diag[i]; */
/*   } */
/*   //xx=D*z; */
/*   for(i=0; i<A->nrow; i++) y[i] = diag[i]*x[i]; */
/*   //xx=L^{-1}z */
/*   for(i=0; i<A->nrow; i++){ */
/*     s = 0; */
/*     for(k=A->row_ptr[i]; k<A->row_ptr[i+1]; k++){ */
/*       j = A->col_ind[k]; */
/*       if(j<i) s += A->val[k]*y[j]; */
/*     } */
/*     y[i] = (x[i] - s)/diag[i]; */
/*   } */

/* } */
