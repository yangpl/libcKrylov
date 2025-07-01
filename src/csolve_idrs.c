/*
Input: 
    A     – matrix (can be matrix-free)
    b     – right-hand side vector
    x     – initial guess
    s     – dimension of shadow space (e.g., s = 4)
    tol   – convergence tolerance
    maxit – maximum number of iterations

Output: 
    x     – approximate solution to Ax = b

1. r = b - A * x                      // Compute initial residual
2. Initialize P[0:s-1] as random vectors of size n (shadow space basis)
3. Initialize G[0:s-1], U[0:s-1] as zero vectors

4. Compute residual norm = ||r||
5. while residual norm > tol and iteration < maxit:
       a. v = A * r
       b. omega = (vᵗ * r) / (vᵗ * v)     // Relaxation step
       c. x = x + omega * r              // Update solution
       d. r = r - omega * v              // Update residual

       e. for k = 0 to s-1:              // Bi-orthogonalize r against P
              alpha = (P[k]ᵗ * r)
              r = r - alpha * G[k]

       f. residual norm = ||r||
       g. Shift G and U arrays:
              for k = s-1 downto 1:
                   G[k] = G[k-1]
                   U[k] = U[k-1]
       h. G[0] = r
       i. U[0] = A * G[0]

       j. iteration += 1

6. return x
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#include "csolver.h"

complex cdotprod(int n, complex *a, complex *b);

// IDR(s) Solver
int csolve_idrs(int n, complex *x, complex *b, cop_t Aop, int s, int niter, double tol, int verb)
{
  int i, j, k, iter;
  complex omega, alpha;
  complex *r = calloc(n, sizeof(complex));
  complex *v = calloc(n, sizeof(complex));
  complex **G = malloc(s * sizeof(complex *));
  complex **U = malloc(s * sizeof(complex *));
  complex **P = malloc(s * sizeof(complex *));
  complex **M = malloc(s * sizeof(complex *));
  double res_norm;
  
  for(i = 0; i < s; i++) {
    G[i] = calloc(n, sizeof(complex));
    U[i] = calloc(n, sizeof(complex));
    P[i] = calloc(n, sizeof(complex));
    M[i] = calloc(s, sizeof(complex));
  }

  Aop(n, x, v);//v=Ax
  for(i = 0; i < n; i++) r[i] = b[i] - v[i];//r=b-Ax

  res_norm = sqrt(creal(cdotprod(n, r, r)));
  if(verb) printf("Iter = 0, ||r|| = %e\n", res_norm);

  // Random P vectors
  for(i = 0; i < s; i++) {
    for(j = 0; j < n; j++) P[i][j] = ((double)rand()/RAND_MAX) + ((double)rand()/RAND_MAX)*I;
  }

  iter = 0;
  while (res_norm > tol && iter < niter) {
    // omega step
    Aop(n, r, v);
    omega = cdotprod(n, v, r) / cdotprod(n, v, v);
    for(i = 0; i < n; i++){
      x[i] += omega * r[i];
      r[i] -= omega * v[i];
    }

    // Bi-orthogonalize r against P
    for(k = 0; k < s; k++) {
      alpha = cdotprod(n, P[k], r);//alpha=<P_k,r>
      for(i = 0; i < n; i++) r[i] -= alpha * G[k][i];
    }

    res_norm = sqrt(creal(cdotprod(n, r, r)));
    iter++;
    if(verb) printf("Iter = %d, ||r|| = %e\n", iter, res_norm);

    // Update shadow space
    for(k = s-1; k > 0; k--) {
      memcpy(G[k], G[k-1], n*sizeof(complex));
      memcpy(U[k], U[k-1], n*sizeof(complex));
    }
    memcpy(G[0], r, n*sizeof(complex));
    Aop(n, G[0], U[0]);
  }

  for(i = 0; i < s; i++) {
    free(G[i]);
    free(U[i]);
    free(P[i]);
    free(M[i]);
  }
  free(G); 
  free(U); 
  free(r); 
  free(v);
  return iter;
}


//IDR(s) with right preconditioning
int solve_idrs_rightpreco(int n, complex *x, complex *b, cop_t Aop, cop_t Minv, int s, int niter, double tol, int verb)
{
  int i, j, k, iter;
  complex omega, alpha;
  complex *y = calloc(n, sizeof(complex));
  complex *t = calloc(n, sizeof(complex));
  complex *r = calloc(n, sizeof(complex));
  complex *v = calloc(n, sizeof(complex));
  complex **G = malloc(s * sizeof(complex *));
  complex **U = malloc(s * sizeof(complex *));
  complex **P = malloc(s * sizeof(complex *));
  complex **M = malloc(s * sizeof(complex *));
  double res_norm;
  
  for(i = 0; i < s; i++) {
    G[i] = calloc(n, sizeof(complex));
    U[i] = calloc(n, sizeof(complex));
    P[i] = calloc(n, sizeof(complex));
    M[i] = calloc(s, sizeof(complex));
  }

  Aop(n, x, v);//v=Ax
  for(i = 0; i < n; i++) r[i] = b[i] - v[i];//r=b-Ax

  res_norm = sqrt(creal(cdotprod(n, r, r)));
  if(verb) printf("Iter = 0, ||r|| = %e\n", res_norm);

  // Random P vectors
  for(i = 0; i < s; i++) {
    for(j = 0; j < n; j++) P[i][j] = ((double)rand()/RAND_MAX) + ((double)rand()/RAND_MAX)*I;
  }

  iter = 0;
  while (res_norm > tol && iter < niter) {
    // omega step
    Minv(n, r, t);//M*t=r, t=M^{-1}r
    Aop(n, t, v);//v=At
    omega = cdotprod(n, v, r) / cdotprod(n, v, v);
    for(i = 0; i < n; i++){
      y[i] += omega * r[i];
      r[i] -= omega * v[i];
    }

    // Bi-orthogonalize r against P
    for(k = 0; k < s; k++) {
      alpha = cdotprod(n, P[k], r);//alpha=<P_k,r>
      for(i = 0; i < n; i++) r[i] -= alpha * G[k][i];
    }

    res_norm = sqrt(creal(cdotprod(n, r, r)));
    iter++;
    if(verb) printf("Iter = %d, ||r|| = %e\n", iter, res_norm);

    // Update shadow space
    for(k = s-1; k > 0; k--) {
      memcpy(G[k], G[k-1], n*sizeof(complex));
      memcpy(U[k], U[k-1], n*sizeof(complex));
    }
    memcpy(G[0], r, n*sizeof(complex));
    Minv(n, G[0], t);//t=M^{-1} G[0]
    Aop(n, t, U[0]);
  }
  Minv(n, y, x);//x=M^{-1} y
  
  for(i = 0; i < s; i++) {
    free(G[i]);
    free(U[i]);
    free(P[i]);
    free(M[i]);
  }
  free(G); 
  free(U); 
  free(r); 
  free(v);
  free(y);
  free(t);
  return iter;
}


/*
complex *A;
// y = A * x
void matvec(int n, complex *x, complex *y)
{
  int i, j;
  for(i = 0; i < n; i++) {
    y[i] = 0.0;
    for(j = 0; j < n; j++)  y[i] += A[i*n + j] * x[j];
  }
}


int main() {
  int n = 10;
  int s = 4;//dimension of shadow space
  double tol = 1e-8;
  int niter = 1000;
  int verb = 1;
  A = calloc(n*n, sizeof(complex));
  complex *b = calloc(n, sizeof(complex));
  complex *x = calloc(n, sizeof(complex));

  srand48((long int) time(NULL));//seed the random number generator
  for(int i = 0; i < n; i++) {
    A[i*n + i] = 1.0 + drand48();
    b[i] = 1.0 + I;
    x[i] = 0.0;
  }

  solve_idrs(n, x, b, matvec, s, niter, tol, verb);

  printf("Final solution x:\n");
  for(int i = 0; i < n; i++) printf("x[%d] = %.5f + %.5fi\n", i, creal(x[i]), cimag(x[i]));

  free(A);
  free(b); 
  free(x);
  return 0;
}

*/
