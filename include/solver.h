#ifndef __solver_h__
#define __solver_h__

typedef void (*op_t)(int, double*, double*); //define type for linear operator

double dotprod(int n, double *a, double *b);

//linear solver using conjugate gradient method
void solve_cg(int n, double *x, double *b, op_t Aop, int niter, double tol, int verb);

//linear solver using conjugate gradient method
void solve_pcg(int n, double *x, double *b, op_t Aop, op_t invMop, int niter, double tol, int verb);

//linear solver using BiCGStab
void solve_bicgstab(int n, double *x, double *b, op_t Aop, int niter, double tol, int verb);

//Reference: Jie Chen, 2016 JSC, Right preconditioned/Flexible BiCGStab paper
void solve_bicgstab_rightpreco(int n, double *x, double *b, op_t Aop, op_t invMop, int niter, double tol, int verb);

//GMRES without preconditioning
void solve_gmres(int n, double *x, double *b, op_t Aop, int niter, double tol, int m, int verb);

//GMRES with right preconditioning
void solve_gmres_rightpreco(int n, double *x, double *b, op_t Aop, op_t invMop, int niter, double tol, int m, int verb);

//CGNR method
void solve_cgnr(int n, double *x, double *b, op_t Aop, op_t Atop, int niter, double tol, int verb);

/*< CGNE (Craig's Method), see Saad's book algorithm 8.5 >*/
void solve_cgne(int n, double *x, double *b, op_t Aop, op_t Atop, int niter, double tol, int verb);


#endif
