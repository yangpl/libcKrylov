# libcKrylov
A generic implementation of linear solvers in Krylov space

Author: Pengliang Yang, Harbin Institute of Technology

E-mail: ypl.2100@gmail.com

Linear solvers in Krylov space
===================================================

* Conjugate gradient (CG) 

* BiCGStab (Bi-conjugate gradient with stabalization)

* BiCGStab with right preconditioning

* GMRES (generalized minimum residual)

* GMRES with right preconditioning

* CGNR

* CGNE

How to install and run
====================================================
1. first compile source code for Krylov solvers library

	cd src;
	
	make
	
   After this step, the implementation of these linear solvers will be compiled to
   form a library in /lib, named libcKrylov.a. To remove the old library and compile 
   for a new one, 'make clean' and 'make'.
   
   
2. Run the demo for frequency domain wave equation:

	cd demo_wave;
	
	make;
	
	./main

   Run the demo for Poisson equation:
   
	cd demo_poisson;
	
	make;
	
	./main
	
   Note that you are free to modify the parameters in main.c to test different algorithm options and the performance.
