# NonlocalEquations
Matlab code using a 1D collocation method for solving linear nonlocal equations.

The collocation method is provided in the work of Huang and Oberman: "Numerical Methods For The Fractional Laplacian: A Finite Difference-Quadrature Approach", SIAM J. Numer. Anal. 52, 3056 (2014) 

Here we've implemented the code for more general nonlocal kernels, and it currently supports solving Dirichlet boundary value problems. 

## Instructions

A detailed write up containing descriptions of the method is contained in the file, 1D_Collocation_Write_Up.pdf. We've provided this as a quick introduction to the ideas behind the Huang/Oberman paper so that the code can be quickly used.

To understand the code and file structure, we will focus on the Oberman example included in the repository. Within the System folder is contained the folder ObermanEx. This folder defines the equations we want to solve as well as code to initialize the equations. Specifically, CreateInitialData.m will construct a file, ObEx1.mat, in the DATA folder which defines the system to solve, it's specific parameter values, and the Dirichlet boundary value problem.

Once this .mat file is created, run the file Main.m. This will solve the equations defined in the .mat file and store the solution back into the .mat file.
