# coding: utf-8
import _LowRankQP
help(_LowRankQP.solve)
import np
import numpy as np
Vmat = np.array([7, 0.627, -1.187, 0.627, 0.056, -0.106, -1.187, -0.106,
  1.452])
dvec = np.array([-4.931, -0.442, 1.335])
Amat = np.array([1.0]*3)
bvec = np.array([1.0])
uvec =np.array([1.0]*3)
alpha = np.array([0.0]*6)
beta = np.array([0.0]*6)
zeta = np.array([0.0]*6)
xi = np.array([0.0]*6)
LowRankQP.solve(Vmat, dvec,Amat, bvec, uvec, alpha, beta, xi, zeta)
