# PyLRQP - Low Rank Quadratic Programming

This package provides a python version of the LowRankQP R package[1];
only the cholesky based method has been ported.

There is one function in this package:

    
    def LowRankQP(Vmat, dvec, Amat, bvec, uvec)

     Solves min α'Vα - α'd subject to Aα = b and 0 ≤ α ≤ u

      Vmat - p x p symmetric matrix
      dvec - p x 1 vector
      Amat - p x k matrix
      bvec - k x 1 vector
      uvec - p x 1 vector



[1] https://cran.r-project.org/web/packages/LowRankQP/index.html

