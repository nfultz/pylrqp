#include "Python.h"
#include "numpy/arrayobject.h"

#include "LowRankQP.h"

static PyObject* solve(PyObject *dummy, PyObject *args) {

    PyArrayObject       *pao_Q, *pao_c, *pao_A, *pao_b, *pao_u, 
            *pao_alpha, *pao_beta, *pao_xi, *pao_zeta;

    double *Q, *c, *A, *b, *u, *alpha, *beta, *xi, *zeta;

    if(
        !PyArg_ParseTuple(args, "O!O!O!O!O!O!O!O!O!",
                                      &PyArray_Type, &pao_Q ,
                                      &PyArray_Type, &pao_c ,
                                      &PyArray_Type, &pao_A ,
                                      &PyArray_Type, &pao_b ,
                                      &PyArray_Type, &pao_u ,

                                      &PyArray_Type, &pao_alpha ,
                                      &PyArray_Type, &pao_beta  ,
                                      &PyArray_Type, &pao_xi    ,
                                      &PyArray_Type, &pao_zeta)
    )
        return NULL;
                         
    Q     = (double*) pao_Q->data;
    c     = (double*) pao_c->data;
    A     = (double*) pao_A->data;
    b     = (double*) pao_b->data;
    u     = (double*) pao_u->data;
    alpha = (double*) pao_alpha->data;
    beta  = (double*) pao_beta->data;
    xi    = (double*) pao_xi->data;
    zeta  = (double*) pao_zeta->data;

    int *n=pao_u->dimensions;
    int *p=pao_b->dimensions;


    int method[] = {2};
    int verbose[] = {0};
    int niter[] = {200};

    LowRankQP(n, n, p, method, verbose, niter, 
            Q, c, A, b, u, 
            alpha, beta, xi, zeta);


    return Py_BuildValue("");


}

static struct PyMethodDef methods[] = {
    {"solve", solve, METH_VARARGS, "Solves a low rank QP problem"},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_LowRankQP(void) {
    (void)Py_InitModule("_LowRankQP", methods);
    import_array();
}
