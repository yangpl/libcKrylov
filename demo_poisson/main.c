/* Demo for conjugate gradient method with preconditioning
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2020-2022 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "solver.h"


void gmg_init(int nx, int ny, double dx, double dy,
	      int itermax_, int v1_, int v2_, int cycleopt_, float tol_, int lmax_);
void gmg_apply(int n, double *b, double *x);//solve Ax=b by multigrid
void gmg_close();


void Acsr_init(int nx, int ny, double dx, double dy);
void Acsr_apply(int n, double *x, double *y);
void Acsr_close();


int main(int argc, char *argv[])
{
  int nx, ny, method, niter, n, verb;
  double dx, dy, x2, y2, tol;
  int i, j, k;
  double *x, *xt, *b;
  FILE *fp;

  verb = 1;/* verbosity */
  nx = 128;/* dimension in x */
  ny = 128;/* dimension in y */
  method = 2;
  niter = 200;
  tol = 1e-8;
  
  dx = 1./nx;
  dy = 1./ny;
  n = (nx+1)*(ny+1);
  
  b = malloc(n*sizeof(double));
  x = malloc(n*sizeof(double));
  xt = malloc(n*sizeof(double));
  
  for(j=0; j<=ny; j++){
    y2 = j*dy;
    y2 = y2*y2;
    for(i=0; i<=nx; i++){
      x2 = i*dx;
      x2 = x2*x2;
      k = i + (nx+1)*j;
      //b[k] = 2.*((1.-6.*x2)*y2*(1-y2) + (1.-6.*y2)*x2*(1.-x2));
      xt[k] = (x2-x2*x2)*(y2-y2*y2);
    }
  }
  memset(x, 0, n*sizeof(double));//x=0 as initialization
  Acsr_init(nx, ny, dx, dy);
  Acsr_apply(n, xt, b);//b=A*xt


  int itermax = 10;/* maximum number of iterations */  
  int v1 = 2;/* number of pre-smoothing */
  int v2 = 2;/* number of post-smoothing */
  int cycleopt = 1;//1=v cycle; 2=f cycle
  int lmax = 0;//number of multigrid levels
  int nx_ = nx;
  while(nx_>1) {nx_=nx_>>1; lmax++;}

  
  if(method==0){//several geometrical multigrid is already an exact iterative solver
    gmg_init(nx, ny, dx, dy, itermax, v1, v2, cycleopt, tol, lmax);//preconditioner initialization
    gmg_apply(n, b, x);//direct solve Ax=b by multigrid
    gmg_close();
  }else if(method==1)
    solve_cg(n, x, b, Acsr_apply, niter, tol, verb);
  else if(method==2){//CG preconditioned by GMG, with only 1 V cycle per preconditioner
    gmg_init(nx, ny, dx, dy, 1, v1, v2, cycleopt, tol, lmax);//preconditioner initialization
    solve_pcg(n, x, b, Acsr_apply, gmg_apply, niter, tol, verb);
    gmg_close();
  }

  
  /* output true signal and reconstructed one */
  fp = fopen("true.bin", "wb");
  for(i=0; i<n; i++){
  float tmp = xt[i];
   fwrite(&tmp, sizeof(float), 1, fp);
   }
  fclose(fp);


  fp = fopen("error.bin", "wb");
  for(i=0; i<n; i++){
  float tmp = xt[i]-x[i];
   fwrite(&tmp, sizeof(float), 1, fp);
   }
  fclose(fp);

  Acsr_close();
  free(b);
  free(x);
  free(xt);

  return 0;
}
