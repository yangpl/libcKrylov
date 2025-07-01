/*
1. Compute initial residual:
    r₀ = b - A * x₀

2. Initialize:
    x = x₀
    u = r₀
    beta = ||u||₂
    u = u / beta

    v = Aᵗ * u
    alpha = ||v||₂
    v = v / alpha

    w = v
    phibar = beta
    rhobar = alpha

3. For k = 1 to max_it:
    a) Bi-diagonalization:
        u = A * v - alpha * u
        beta = ||u||₂
        u = u / beta

        v = Aᵗ * u - beta * v
        alpha = ||v||₂
        v = v / alpha

    b) Apply Givens rotation:
        rho = sqrt(rhobar² + beta²)
        c = rhobar / rho
        s = beta / rho
        theta = s * alpha
        rhobar = -c * alpha
        phi = c * phibar
        phibar = s * phibar

    c) Update solution:
        x = x + (phi / rho) * w
        w = v - (theta / rho) * w

    d) Check convergence:
        if |phibar| < tol * ||b||₂:  break

4. Return x
*/#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#include "csolver.h"

complex cdotprod(int n, complex *a, complex *b);

// LSQR algorithm
int csolve_lsqr(int n, complex *x, complex *b, cop_t A, cop_t AT, int niter, double tol, int verb)
{
  double alpha, beta, phibar, rhobar, normr;
  double rho, c, s, theta, phi;
  int iter, i;
  complex *u = malloc(n*sizeof(complex));
  complex *v = malloc(n*sizeof(complex));
  complex *w = malloc(n*sizeof(complex));
  complex *um = malloc(n*sizeof(complex));
  complex *vn = malloc(n*sizeof(complex));

  A(n, x, u);
  for(i=0; i<n; ++i) u[i] = b[i] - u[i];
  beta = sqrt(creal(cdotprod(n, u, u)));
  for(i=0; i<n; ++i) u[i] /= beta;
  
  AT(n, u, v);
  alpha = sqrt(creal(cdotprod(n, v, v)));
  for(i=0; i<n; ++i) v[i] /= alpha;
  memcpy(w, v, n*sizeof(complex));
  
  phibar = beta;
  rhobar = alpha;
  normr = beta;
  for(iter=0; iter<niter; ++iter){
    // A*v - alpha*u -> u
    A(n, v, um);
    for(i=0; i<n; ++i) u[i] = um[i] - alpha*u[i];
    beta = sqrt(creal(cdotprod(n, u, u)));//beta=||u||_2
    for(i=0; i<n; ++i) u[i] /= beta;

    // A^T*u - beta*v -> v
    AT(n, u, vn);
    for(i=0; i<n; ++i) v[i] = vn[i] - beta*v[i];
    alpha = sqrt(creal(cdotprod(n, v, v)));//alpha=||v||_2
    for(i=0; i<n; ++i) v[i] /= alpha;

    // Apply Givens rotation
    rho = sqrt(rhobar * rhobar + beta * beta);
    c = rhobar / rho;
    s = beta / rho;
    theta = s * alpha;
    rhobar = -c * alpha;
    phi = c * phibar;
    phibar = s * phibar;

    c = phi/rho;
    s = theta/rho;
    for(i = 0; i < n; ++i){
      x[i] = x[i] + c*w[i];
      w[i] = v[i] - s*w[i];
    }
    normr = fabs(phibar);
    if(verb) printf("Iter %d, |r| = %.6e\n", iter, normr);
    //printf("normr=%e, tol=%e, beta=%e, tol*beta=%e\n", normr, tol, beta, tol*beta);
    //if(normr<tol*beta) break;
  }

  free(u);
  free(v);
  free(w);
  free(um);
  free(vn);
    
  return iter;
}
/*
// Main function
int main() {
  complex x[N] = {1., 2., -1};
  complex b[M];
  matvec_A(M, N, x, b);

  memset(x, 0, N*sizeof(complex));
  int iters = lsqr(M, N, matvec_A, matvec_AT, b, x, 100, 1e-10, 1);

  printf("LSQR finished in %d iterations.\n", iters);
  printf("Solution x = [");
  for (int i = 0; i < N; ++i)  printf(" %.6f+ I*%.6f,", creal(x[i]), cimag(x[i]));
  printf("]\n");
    
  return 0;
}
*/
