# coding: utf-8
import _LowRankQP
from collections import namedtuple
from numpy import zeros

_result_type = namedtuple('pylrqp_result', 'alpha beta xi zeta')

def LowRankQP(Vmat, dvec, Amat, bvec, uvec):
    k,p = Amat.shape
    results = _result_type(
            zeros(p),
            zeros(k),
            zeros(p),
            zeros(p)
            )

    _LowRankQP.solve(Vmat, dvec, Amat, bvec, uvec, 
            results.alpha, results.beta, results.xi, results.zeta)
    return results
