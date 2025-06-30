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
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Matrix size
#define M 3  // Rows
#define N 3  // Columns

// Define dense matrix A and vector b
double A_data[M][N] = {
  {1, 2, 6},
  {3, 4, 3},
  {5, 6, 2}
};

double b_data[M] = {5, 11, 17};

// Function pointer types for matrix-vector operations
typedef void (*Aop)(int m, int n, double *x, double *Ax);
typedef void (*Atop)(int m, int n, double *x, double *ATx);

// Dot product
double dot(int n, const double *a, const double *b) {
  double sum = 0.0;
  for (int i = 0; i < n; ++i) sum += a[i] * b[i];
  return sum;
}

// y = A x
void matvec_A(int m, int n, double *x, double *Ax) {
  for (int i = 0; i < m; ++i) {
    Ax[i] = 0.0;
    for (int j = 0; j < n; ++j) Ax[i] += A_data[i][j] * x[j];
  }
}

// y = A^T x
void matvec_AT(int m, int n, double *x, double *ATx) {
  for (int j = 0; j < n; ++j) {
    ATx[j] = 0.0;
    for (int i = 0; i < m; ++i) ATx[j] += A_data[i][j] * x[i];
  }
}

// LSQR algorithm
int lsqr(int m, int n, Aop A, Atop AT, const double *b, double *x, int niter, double tol, int verb)
{
  double alpha, beta, phibar, rhobar, normr;
  double rho, c, s, theta, phi;
  int iter, i;
  double *u = malloc(m*sizeof(double));
  double *v = malloc(n*sizeof(double));
  double *w = malloc(n*sizeof(double));
  double *um = malloc(m*sizeof(double));
  double *vn = malloc(n*sizeof(double));

  A(m, n, x, u);
  for(i=0; i<m; ++i) u[i] = b[i] - u[i];
  beta = sqrt(dot(m, u, u));
  for(i=0; i<m; ++i) u[i] /= beta;
  
  AT(m, n, u, v);
  alpha = sqrt(dot(n, v, v));
  for(i=0; i<n; ++i) v[i] /= alpha;
  memcpy(w, v, n*sizeof(double));
  
  phibar = beta;
  rhobar = alpha;
  normr = beta;
  for(iter=0; iter<niter; ++iter){
    // A*v - alpha*u -> u
    A(m, n, v, um);
    for(i=0; i<m; ++i) u[i] = um[i] - alpha*u[i];
    beta = sqrt(dot(m, u, u));//beta=||u||_2
    for(i=0; i<m; ++i) u[i] /= beta;

    // A^T*u - beta*v -> v
    AT(m, n, u, vn);
    for(i=0; i<n; ++i) v[i] = vn[i] - beta*v[i];
    alpha = sqrt(dot(n, v, v));//alpha=||v||_2
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
    if(normr<tol*beta) break;
  }

  free(u);
  free(v);
  free(w);
  free(um);
  free(vn);
    
  return iter;
}

// Main function
int main() {
  double x[N] = {0};

  int iters = lsqr(M, N, matvec_A, matvec_AT, b_data, x, 100, 1e-10, 1);

  printf("LSQR finished in %d iterations.\n", iters);
  printf("Solution x = [");
  for (int i = 0; i < N; ++i)  printf(" %.6f", x[i]);
  printf("]\n");
    
  return 0;
}
