
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#define MAX(A,B) ( (A) > (B) ? (A):(B))
#define MIN(A,B) ( (A) < (B) ? (A):(B))

#define CHOL 2

#define PRED 1
#define CORR 2
/**** BLAS ****/
extern double dasum_(int*, double*, int*);

extern void daxpy_(const int *n, const double *alpha,
       const double *dx, const int *incx,
                                                                    double *dy, const int *incy);
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

/**** LAPACK ****/

extern void dpotrf_(const char* uplo, const int* n,
                 double* a, const int* lda, int* info);
extern void dpotrs_(const char* uplo, const int* n,
                 const int* nrhs, const double* a, const int* lda, double* b,
                 const int* ldb, int* info);


/*****************************************************************************n

/* Global Variables */
#define EPSTERM   5.0E-11
#define EPSROUND  1.0E-8
#define EPSINIT   1.0E-2
#define EPSIPM    1.0E-2
#define EPSPERT   1.0E-14

/*****************************************************************************/

void PrintMatrix( char* name, double* vals, int* rows, int* cols )
{
    int i, j;
    printf("%8s = [\n", name);

    for (i=0;i<(*rows);i++)
    {
        for (j=0;j<(*cols);j++)
        {
            printf("%15.10e ", vals[i+j*(*rows)]);
        }
        printf(";\n");
    }
    printf("];\n");
}

/*****************************************************************************/

void VectorVectorCopy( double* lhs, double* rhs, int* n )
{
    int one=1;
    dcopy_( n, rhs, &one, lhs, &one );
    /* dcopy( n, rhs, &one, lhs, &one ); */
}

/*****************************************************************************/

void VectorVectorMult( double* alpha, double* x, double* y, int* n )
{
    int one=1;
    daxpy_( n, alpha, x, &one, y, &one );
    /* daxpy( n, alpha, x, &one, y, &one ); */
}

/*****************************************************************************/

double VectorVectorDot( double* x, double* y, int* n )
{
    int one=1;
    return ddot_( n, x, &one, y, &one );
    /* return ddot( n, x, &one, y, &one ); */
}

/*****************************************************************************/

void MatrixVectorMult( double* alpha, double* A, int trans, double* x, double*
beta, double* b, int* rows, int* cols )
{
    int one=1;
    if (trans)
    {
        dgemv_("T", rows, cols, alpha, A, rows, x, &one, beta, b, &one );
        /* dgemv('T', *rows, *cols, *alpha, A, *rows, x, one, *beta, b, one ); */
    }
    else
    {
        dgemv_("N", rows, cols, alpha, A, rows, x, &one, beta, b, &one );
        /* dgemv('N', *rows, *cols, *alpha, A, *rows, x, one, *beta, b, one ); */
    }
}

/*****************************************************************************/

void MatrixCholFactorize( double* A, int* n, int* info )
{
    dpotrf_( "L", n, A, n, info );
    /* dpotrf( 'L', *n, A, *n, info ); */
}

/*****************************************************************************/

void MatrixCholSolve( double* A, int* n, double* rhs, int *nrhs, int* info )
{
    dpotrs_("L", n, nrhs, A, n, rhs, n, info );
    /* dpotrs('L', *n, *nrhs, A, *n, rhs, *n, info ); */
}


/*****************************************************************************/

void MatrixMatrixCopy( double* lhs, double* rhs, int* rows, int* cols )
{
    int i;
    int len = (*rows)*(*cols);
    for (i=0;i<len;i++) lhs[i] = rhs[i];
}



/*****************************************************************************/

void LRQPHeader()
{
    printf("ITER  PRIM            DUAL            COMP            GAP           TERM\n");
}

/*****************************************************************************/

void LRQPInitPoint( int *n, int *m, int *p, double *Q, double *c, double *A,
    double *b, double *u, double *alpha, double* beta, double *xi, double *zeta,
    double *w, double *temp )
{
    int i;
    double pone =  1.0;
    double mone = -1.0;
    double zero =  0.0;
    for (i=0;i<(*n);i++) alpha[i] = u[i] > 1 ? EPSINIT : u[i] * EPSINIT;
    for (i=0;i<(*p);i++) beta[i]  = 0.0;
    MatrixVectorMult( &pone, Q, 1, alpha, &zero, w, n, m );
    for (i=0;i<(*n);i++) temp[i] += -w[i];
    VectorVectorMult( &mone, c, temp, n );
    for (i=0;i<(*n);i++)
    {
        xi[i]   = MAX(EPSINIT,temp[i]);
        zeta[i] = MAX(EPSINIT,xi[i]-temp[i]);
    }
}

/*****************************************************************************
*
*  MATLAB CODE FOR QpIpmCalcStats
*
*    w               = Q'*alpha;
*    UminusAlpha     = (u - alpha);
*    XiOnUminusAlpha = xi./UminusAlpha;
*    ZetaOnAlpha     = zeta./alpha;
*
*    if (n==m)
*        r1 = - c - w   - A*beta - xi + zeta;
*        quad = alpha'*w;
*    else
*        r1 = - c - Q*w   - A*beta - xi + zeta;
*        quad = w'*w;
*    end;
*
*    r2      = b - A'*alpha;
*    comp    = alpha'*zeta + UminusAlpha'*xi;
*    cTalpha = c'*alpha;
*    gap     = abs(quad + cTalpha + u'*xi + b'*beta);
*    term    = comp / (abs( 0.5*quad + cTalpha ) + 1.0);
*    t       = comp*(((1.0-mult+eps)/(10+mult))^2)/(2*n);
*    D       = spdiags( XiOnUminusAlpha + ZetaOnAlpha + 1.0e-14, 0, n, n );
*
******************************************************************************/

void LRQPCalcStats( int *n, int *m, int *p, double *Q, double *c, double *A,
    double *b, double *u, double *alpha, double* beta, double *xi, 
    double *zeta, double *dalpha, double* dbeta, double *dxi, double *dzeta,
    double *UminusAlpha, double *XiOnUminusAlpha, double *ZetaOnAlpha, 
    double* w, double *r1, double *r2, double *D, double *prim, double *dual, 
    double *comp, double *gap, double *term, double * mult, double *t)
{
    int i;
    int one = 1;
    double quad;
    double cTalpha;
    double temp;
    double pone =  1.0;
    double mone = -1.0;
    double zero =  0.0;
    MatrixVectorMult( &pone, Q, 1, alpha, &zero, w, n, m );

    VectorVectorCopy( UminusAlpha, u, n ); // UminusAlpha = u
    VectorVectorMult( &mone, alpha, UminusAlpha, n ); // UminusAlpha = UminusAlpha + -1 * alpha

    for (i=0;i<(*n);i++) XiOnUminusAlpha[i] = xi[i]/UminusAlpha[i];
    for (i=0;i<(*n);i++) ZetaOnAlpha[i] = zeta[i]/alpha[i];

        MatrixVectorMult( &mone, A, 0, beta, &zero, r1, n, p );
        VectorVectorMult( &mone, w,  r1, n );
        VectorVectorMult( &mone, c,  r1, n );
        VectorVectorMult( &mone, xi, r1, n );
        VectorVectorMult( &pone, zeta, r1, n );
        quad = VectorVectorDot( alpha, w, n );

        VectorVectorCopy( r2, b, p );
        MatrixVectorMult( &mone, A, 1, alpha, &pone, r2, n, p );
        *dual = dasum_(p, r2, &one ); //sum(abs(r2))

    *prim   = dasum_(n, r1, &one ); // sum(abs(r1))
    *comp   = VectorVectorDot( alpha, zeta, n ) + VectorVectorDot( UminusAlpha, xi, n );
    cTalpha = VectorVectorDot( c, alpha, n );
    *gap = fabs( quad + cTalpha + VectorVectorDot( u, xi, n ) + VectorVectorDot( b, beta, p ) );
    *term   = *comp / ( fabs( 0.5*quad + cTalpha) + 1.0);
    temp    = (1.0 - *mult + EPSIPM)/(10.0 + *mult);
    *t      = *comp*(temp*temp)/(2*(*n));
    for (i=0;i<(*n);i++) D[i] = XiOnUminusAlpha[i] + ZetaOnAlpha[i] + EPSPERT;
}

/*****************************************************************************/

void LRQPDisplay( int i, double *prim, double *dual, double *comp, double *gap,
                  double *term )
{
    printf("%3d %15.7e %15.7e %15.7e %15.7e %15.7e \n", i, *prim, *dual, *comp, *gap, *term );
}

/*****************************************************************************/

void LRQPSummary( int i, int niter, int method, int n, int m, double *prim, 
    double *dual, double *comp, double *gap, double *term )
{
    if (i==niter)
    {
        printf("LowRankQP FAILED TO CONVERGE\n");
        printf("    Try increasing niter, or rescaling problem.\n");
    }
    else
    {
        printf("LowRankQP CONVERGED IN %d ITERATIONS\n\n", i+1 );
        printf("    Primal Feasibility    = %15.7e\n", *prim);
        printf("    Dual Feasibility      = %15.7e\n", *dual);
        printf("    Complementarity Value = %15.7e\n", *comp);
        printf("    Duality Gap           = %15.7e\n", *gap);
        printf("    Termination Condition = %15.7e\n", *term);
    }
}

/******************************************************************************
*
*  MATLAB CODE FOR IpmSolve
*
*    if (strcmp(info.method,'FULL'))
*        sol = info.M \ rhs;
*    elseif (strcmp(info.method,'CHOL'))
*        z   = info.M'\rhs;
*        sol = info.M\z;
*    elseif (strcmp(info.method,'SMW'))
*        z   = D \ rhs;
*        t   = info.M'\(Q'*z);
*        t   = info.M \ t;
*        sol = z - D \ (Q*t);
*    elseif (strcmp(info.method,'PFCF'))
*        PfcfSolve()
*    end;
*
******************************************************************************/

void LRQPSolve( int *n, int *m, int *nrhs, int *method, double *Q, double *D,
    double *rhs, double *sol, double *M, int* pivN, double *buffMxP, double *P,
    double *Beta, double *Lambda )
{
    int i, j;
    int info    = 0;
    int one     = 1;
    double pone =  1.0;
    double mone = -1.0;
    double zero =  0.0;

    

    MatrixMatrixCopy( sol, rhs, n, nrhs );
        MatrixCholSolve( M, n, sol, nrhs, &info );
}

/******************************************************************************
*
*  MATLAB CODE FOR IpmSolve
*
*    R      = QpIpmSolve( PREDinfo, Q, D, A );
*    r3     = - zeta;    OR    r3 = (t - Dalpha.*Dzeta)./alpha - zeta;
*    r4     = - xi;      OR    r4 = (t + Dalpha.*Dxi )./UminusAlpha - xi;
*    r5     = r1 + r3 - r4
*    r      = QpIpmSolve( info, Q, D, r5 );
*    Dbeta  = (A'*R) \ (A'*r - r2);
*    Dalpha = r - R*Dbeta;
*    Dzeta  = r3 - ZetaOnAlpha.*Dalpha;
*    Dxi    = r4 + XiOnUminusAlpha.*Dalpha;
*
******************************************************************************/

void LRQPCalcDx( int *n, int *m, int *p, int *method, double *Q, double *c,
    double *A, double *b, double * u, double *alpha, double* beta, double *xi,
    double *zeta, double *dalpha, double* dbeta, double *dxi, double *dzeta,
    double *UminusAlpha, double *ZetaOnAlpha, double *XiOnUminusAlpha,
    double *buffMxP, double *buffMx1, double* buffPxP, double *buffPx1,
    int *pivN, double *R, double *r, double *r1, double* r2, double *r3,
    double *r4, double* r5, double *D, double *M, double *t, double *P,
    double *Beta, double *Lambda, double *LambdaTemp, double *T, int predcorr)
{

    int i, j;
    int    info = 0;
    int    one  = 1;
    double pone =  1.0;
    double mone = -1.0;
    double zero =  0.0;
    if (predcorr==PRED)
    {
        LRQPSolve( n, m, p, method, Q, D, A, R, M, pivN, buffMxP, P, Beta, Lambda );
    }
    for (i=0;i<(*n);i++)
    {
        r3[i] = - zeta[i];
        r4[i] = - xi[i];
    }
    if (predcorr==CORR)
    {
        for (i=0;i<(*n);i++) r3[i] += ( *t - (dalpha[i] * dzeta[i]) )/alpha[i];
        for (i=0;i<(*n);i++) r4[i] += ( *t + (dalpha[i] * dxi[i]) )/UminusAlpha[i];
    }
    for (i=0;i<(*n);i++) r5[i] = r1[i] + r3[i] - r4[i];

        LRQPSolve( n, m, &one, method, Q, D, r5, r, M, pivN, buffMx1, P, Beta, Lambda );
        VectorVectorCopy( buffPx1, r2, p );
        MatrixVectorMult( &pone, A, 1, r, &mone, buffPx1, n, p );
        
        dgemm_("T","N", p,p,n,&pone, A, n, R, n, &zero, buffPxP, p); // buffPxP = A' * R
        MatrixCholFactorize( buffPxP, p, &info );
        MatrixCholSolve( buffPxP, p, buffPx1, &one, &info);
        VectorVectorCopy( dbeta, buffPx1, p );
        VectorVectorCopy( dalpha, r , n);
        MatrixVectorMult( &mone, R, 0, dbeta, &pone, dalpha, n, p );
    
    for (i=0;i<(*n);i++) dzeta[i] = r3[i] - (ZetaOnAlpha[i] * dalpha[i]);
    for (i=0;i<(*n);i++) dxi[i]   = r4[i] + (XiOnUminusAlpha[i] * dalpha[i]);

}

/*****************************************************************************
*
*  MATLAB CODE FOR QpIpmStep
*
*    temp = [ 1.0 (-alpha./Dalpha)' (UminusAlpha./Dalpha)' (-xi./Dxi)'(-zeta./Dzeta)' ];
*    mult = 0.99*min( temp( find(temp>0.0) ) );
*    alpha = alpha + mult*Dalpha;
*    beta  = beta  + mult*Dbeta;
*    xi    = xi    + mult*Dxi;
*    zeta  = zeta  + mult*Dzeta;
*
*****************************************************************************/

void LRQPStep( int *n, int *p, double *alpha, double* beta, double *xi,
    double *zeta, double *dalpha, double* dbeta, double *dxi, double *dzeta,
    double *UminusAlpha, double *mult)
{
    int i;
    *mult= 1.0;
    for (i=0;i<(*n);i++)
    {
        if (dalpha[i]<0.0) *mult = MIN(*mult,(-alpha[i]/dalpha[i]));
        if (dalpha[i]>0.0) *mult = MIN(*mult,(UminusAlpha[i]/dalpha[i]));
        if (dxi[i]<0.0)    *mult = MIN(*mult,(-xi[i]/dxi[i]));
        if (dzeta[i]<0.0)  *mult = MIN(*mult,(-zeta[i]/dzeta[i]));
    }
    *mult *= 0.99;
    VectorVectorMult( mult, dalpha, alpha, n );
    VectorVectorMult( mult, dbeta,  beta,  p );
    VectorVectorMult( mult, dxi,    xi,    n );
    VectorVectorMult( mult, dzeta,  zeta,  n );
}

/******************************************************************************/

void LowRankQP( int *n, int *m, int *p, int* method, int* verbose, int* niter, 
    double *Q, double *c, double *A, double *b, double *u, double *alpha,
    double* beta, double *xi, double *zeta)
{
    int i;
    int info = 0;
    
    long int start;
    long int finish;

    /* Iteration Display variables */
    double mult = 0.0;
    double prim = 0.0;
    double dual = 0.0;
    double comp = 0.0;
    double gap  = 0.0;
    double term = 0.0;
    double t    = 0.0;

    /* Step direction vectors */
    double *dalpha = (double *) calloc( (*n), sizeof(double) );
    double *dbeta;
    double *dxi    = (double *) calloc( (*n), sizeof(double) );
    double *dzeta  = (double *) calloc( (*n), sizeof(double) );

    /* Some repeatedly occuring vectors */
    double *UminusAlpha     = (double *) calloc( *n, sizeof(double) );
    double *XiOnUminusAlpha = (double *) calloc( *n, sizeof(double) );
    double *ZetaOnAlpha     = (double *) calloc( *n, sizeof(double) );

    /* Some vectors used during calculations */
    double *w  = (double *) calloc( *m, sizeof(double) );
    double *r1 = (double *) calloc( *n, sizeof(double) );
    double *r2;
    double *r3 = (double *) calloc( *n, sizeof(double) );
    double *r4 = (double *) calloc( *n, sizeof(double) );
    double *r5 = (double *) calloc( *n, sizeof(double) );
    double *D  = (double *) calloc( *n, sizeof(double) );
    double *r  = (double *) calloc( *n, sizeof(double) );
    double *R;

    /* Various Buffers */
    double *buffMxP;
    double *buffPxP;
    double *buffPx1;

    double *M;
    int    *pivN;

    double *buffNxM;
    double *buffMx1;

    double *P;
    double *Beta;
    double *Lambda;
    double *LambdaTemp;
    double *T;

    /* Vectors to be created if p!=0 */
        dbeta   = (double *) calloc( (*p), sizeof(double) );
        r2      = (double *) calloc( (*p), sizeof(double) );
        R       = (double *) calloc( (*n)*(*p), sizeof(double) );
        buffMxP = (double *) calloc( (*m)*(*p), sizeof(double) );
        buffPxP = (double *) calloc( (*p)*(*p), sizeof(double) );
        buffPx1 = (double *) calloc( (*p), sizeof(double) );

        M    = (double *) calloc( (*n)*(*n), sizeof(double) );
        pivN = (int *) calloc( *n, sizeof(int) );

    /* Main Loop */
    if ( *verbose ) LRQPHeader();
    LRQPInitPoint( n, m, p, Q, c, A, b, u, alpha, beta, xi, zeta, w, r1 );

    for (i=0;i<(*niter);i++)
    {
        /* start = clock(); */
        LRQPCalcStats( n, m, p, Q, c, A, b, u, alpha, beta, xi, zeta, dalpha,
            dbeta, dxi, dzeta, UminusAlpha, XiOnUminusAlpha, ZetaOnAlpha, w, r1,
            r2, D, &prim, &dual, &comp, &gap, &term, &mult, &t );
        /* finish = clock();
        TimeIpmCalcStats += (finish - start)/CLK_TCK; */

        if ( *verbose ) LRQPDisplay( i+1, &prim, &dual, &comp, &gap, &term  );
        if ( term < EPSTERM ) break;

        //LRQPFactorize( n, m, method, Q, D, M, pivN, buffNxM, P, Beta, Lambda,
        //    LambdaTemp, T );
        MatrixMatrixCopy( M, Q, n, n );
        // diag(M) = Diag(M) + D
        for (int i=0;i<(*n);i++) M[i+i*(*n)] += D[i];
        MatrixCholFactorize( M, n, &info );

        LRQPCalcDx( n, m, p, method, Q, c, A, b, u, alpha, beta, xi, zeta,
            dalpha, dbeta, dxi, dzeta, UminusAlpha, ZetaOnAlpha, 
            XiOnUminusAlpha, buffMxP, buffMx1, buffPxP, buffPx1, pivN, R, r, r1,
            r2, r3, r4, r5, D, M, &t, P, Beta, Lambda, LambdaTemp, T, PRED);

        LRQPCalcDx( n, m, p, method, Q, c, A, b, u, alpha, beta, xi, zeta,
            dalpha, dbeta, dxi, dzeta, UminusAlpha, ZetaOnAlpha,
            XiOnUminusAlpha, buffMxP, buffMx1, buffPxP, buffPx1, pivN, R, r, r1,
            r2, r3, r4, r5, D, M, &t, P, Beta, Lambda, LambdaTemp, T, CORR);

        LRQPStep( n, p, alpha, beta, xi, zeta, dalpha, dbeta, dxi, dzeta,
            UminusAlpha, &mult );
    }
    
    LRQPSummary( i, *niter, *method, *n, *m, &prim, &dual, &comp, &gap, &term );

    /* Free Memory */
    free(dalpha);      free(dxi);             free(dzeta);
    free(UminusAlpha); free(XiOnUminusAlpha); free(ZetaOnAlpha);
    free(r1); free(r3); free(r4); free(r5);
    free(D);  free(w);  free(r);  
        
        free( dbeta );
        free( r2 );
        free( R );
        free( buffMxP );
        free( buffPxP );
        free( buffPx1 );

        free( M );
        free( pivN );
}

int main(){
  int n = 2, m = 2, p = 1, method = 2,  verbose=0, niter=200;
  
  double alpha[] = {0,0,0,0,0,0};
  double  beta[] = {0,0,0,0,0,0};
  double    xi[] = {0,0,0,0,0,0};
  double  zeta[] = {0,0,0,0,0,0};

  double Vmat[] = {20, -1.54, -1.54, 9.10};
  double dvec[] = {0.2,24.2};
  double Amat[] = {-1, 4};
  double bvec[] = {0.2};
  double uvec[] = {1,1};

  LowRankQP( &n, &m, &p, &method, &verbose, &niter,
              Vmat, dvec, Amat, bvec, uvec,
              alpha, beta, xi, zeta);

  printf("Sum of dvec: %f\n", (dasum_(&n, dvec, &p)));

  PrintMatrix("alpha", alpha, &n, &p); 
  PrintMatrix("beta", beta, &p, &p); 
  PrintMatrix("xi", xi, &n, &p); 
  PrintMatrix("zeta", zeta, &n, &p); 

}
