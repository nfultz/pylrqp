
/**** BLAS ****/
extern double dasum_(int*, double*, int*);
extern int idamax_( const int*     N, const double* DX, const int*     INCX );

extern void dscal_( const int*     N, const double* DA, const double* DX, const int*     INCX );

extern void daxpy_(const int *n, const double *alpha, const double *dx, const
        int *incx, double *dy, const int *incy);

extern void   dcopy_(const int *n, const double *dx, const int *incx,
                                      double *dy, const int *incy);
extern double ddot_(const int *n, const double *dx, const int *incx, const
        double *dy, const int *incy);


extern void   dgemv_(const char *trans, const int *m, const int *n, const
        double *alpha, const double *a, const int *lda, const double *x, const
        int *incx, const double *beta, double *y, const int *incy);

extern void   dgemm_(const char *transa, const char *transb, const int *m,
        const int *n, const int *k, const double *alpha, const double *a, const
        int *lda, const double *b, const int *ldb, const double *beta, double
        *c, const int *ldc);

extern void dtbsv_(const   char*   UPLO, const char*   TRANS, const char* DIAG,
        const int*     N, const int*     K, const double* A, const int*     LDA,
        const double* X, const int*     INCX);

extern void dgbmv_( const char*   TRANS, const int*     M, const int*     N,
        const int*     KL, const int*     KU, const double* ALPHA, const
        double*  A, const int*     LDA, const double*  X, const int*     INCX,
        const double* BETA, const double*  Y, const int*     INCY )     ;

/**** LAPACK ****/

extern void dpotrf_(const char* uplo, const int* n,
                 double* a, const int* lda, int* info);
extern void dpotrs_(const char* uplo, const int* n,
                 const int* nrhs, const double* a, const int* lda, double* b,
                 const int* ldb, int* info);
