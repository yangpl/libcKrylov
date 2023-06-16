#ifndef __csolver_h__
#define __csolver_h__

#include <complex.h>

typedef void (*cop_t)(int, complex*, complex*); //define type for linear operator

complex cdotprod(int n, complex *a, complex *b);

void csolve_cg(int n, complex *x, complex *b, cop_t Aop, int niter, double tol, int verb);

//linear solver using BiCGStab
void csolve_bicgstab(int n, complex *x, complex *b, cop_t Aop, int niter, double tol, int verb);

//Reference: Jie Chen, 2016 JSC, Right preconditioned/Flexible BiCGStab paper
void csolve_bicgstab_rightpreco(int n, complex *x, complex *b, cop_t Aop, cop_t invMop, int niter, double tol, int verb);

//GMRES without preconditioning
void csolve_gmres(int n, complex *x, complex *b, cop_t Aop, int niter, double tol, int m, int verb);

//GMRES with right preconditioning
void csolve_gmres_rightpreco(int n, complex *x, complex *b, cop_t Aop, cop_t invMop, int niter, double tol, int m, int verb);

//CGNR method
void csolve_cgnr(int n, complex *x, complex *b, cop_t Aop, cop_t Atop, int niter, double tol, int verb);

/*< CGNE (Craig's Method), see Saad's book algorithm 8.5 >*/
void csolve_cgne(int n, complex *x, complex *b, cop_t Aop, cop_t Atop, int niter, double tol, int verb);

#endif
