
1. The GSL Scientific library is needed for Robust fit. Minor modification is needed to make it faster.
   The maximum iteration is changed to 5 and convergence is not needed.
2. R "glmnet" library is needed for lasso regression. The long-term goal is to have R interface be 
   removed as the source code of "glmnet" is really in Fortran.

