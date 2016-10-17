
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

void LRQPHeader()
{
    printf("ITER  PRIM            DUAL            COMP            GAP           TERM\n");
}

/*****************************************************************************/

void LRQPInitPoint( int *n, int *p, double *Q, double *c, double *A,
    double *b, double *u, double *alpha, double* beta, double *xi, double *zeta,
    double *w, double *temp )
{
    int i;
    int one = 1;
    double pone =  1.0;
    double mone = -1.0;
    double zero =  0.0;
    for (i=0;i<(*n);i++) alpha[i] = u[i] > 1 ? EPSINIT : u[i] * EPSINIT;
    for (i=0;i<(*p);i++) beta[i]  = 0.0;

    dgemv_("T", n, n, &pone, Q, n, alpha, &one, &zero, w, &one ); // w = Q' * alpha
    for (i=0;i<(*n);i++) temp[i] += -w[i];
    daxpy_(n, &mone, c, &one, temp, &one); //temp = temp - c
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

    dgemv_("T", n, m, &pone, Q, n, alpha, &one, &zero, w, &one ); // w = Q' * alpha


    dcopy_(n, u, &one, UminusAlpha, &one);  // UminusAlpha = u
    daxpy_(n, &mone, alpha, &one, UminusAlpha, &one); //// UminusAlpha = UminusAlpha + -1 * alpha

    for (i=0;i<(*n);i++) XiOnUminusAlpha[i] = xi[i]/UminusAlpha[i];
    for (i=0;i<(*n);i++) ZetaOnAlpha[i] = zeta[i]/alpha[i];


        dgemv_("N", n, p, &mone, A, n, beta, &one, &zero, r1, &one); // r1 = -A * beta
        daxpy_(n, &mone, w, &one, r1, &one ); // r1= r1 - w
        daxpy_(n, &mone, c, &one, r1, &one ); // r1= r1 - c
        daxpy_(n, &mone, xi, &one, r1, &one ); // r1= r1 - xi
        daxpy_(n, &pone, zeta, &one, r1, &one ); // r1= r1 + zeta

        quad = ddot_(n, alpha, &one, w, &one );


        dcopy_(p, b, &one, r2, &one); // r2 = b

        dgemv_("T", n, p, &mone, A, n, alpha, &one, &pone, r2, &one); // r2 = r2 - A' * alpha
        *dual = dasum_(p, r2, &one ); //sum(abs(r2))

    *prim   = dasum_(n, r1, &one ); // sum(abs(r1))

    *comp   = ddot_(n,  alpha, &one, zeta,  &one ) + ddot_(n, UminusAlpha, &one, xi,  &one );

    cTalpha = ddot_(n, c, &one, alpha, &one );

    *gap = fabs( quad + cTalpha + ddot_(n, u, &one, xi, &one ) + ddot_(p, b, &one, beta, &one ) );
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

void LRQPSummary( int i, int niter, double *prim, 
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
*        z   = info.M'\rhs;
*        sol = info.M\z;
*
*    dcopy_(&len, rhs, &one, sol, &one); // copy rhs to sol
*    dpotrs_("L", n, nrhs, M, n, sol, n, &info ); // sol = M^-1 * sol
******************************************************************************/

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

void LRQPCalcDx( int *n, int *p, double *Q, double *c,
    double *A, double *b, double * u, double *alpha, double* beta, double *xi,
    double *zeta, double *dalpha, double* dbeta, double *dxi, double *dzeta,
    double *UminusAlpha, double *ZetaOnAlpha, double *XiOnUminusAlpha,
    double* buffPxP, double *buffPx1,
    double *R, double *r, double *r1, double* r2, double *r3,
    double *r4, double* r5, double *D, double *M, double *t, int predcorr)
{

    int i, j;
    int    info = 0;
    int    one  = 1;
    double pone =  1.0;
    double mone = -1.0;
    double zero =  0.0;
    for (i=0;i<(*n);i++)
    {
        r3[i] -= zeta[i];
        r4[i] -= xi[i];
    }
    for (i=0;i<(*n);i++) r5[i] = r1[i] + r3[i] - r4[i];

        //LRQPSolve( n, m, &one, method, Q, D, r5, r, M, pivN, buffMx1, P, Beta, Lambda );
        dcopy_(n, r5, &one, r, &one); // copy r5 to r
        dpotrs_("L", n, &one, M, n, r, n, &info ); // r = M^-1 * r

        dcopy_(p, r2, &one, buffPx1, &one); // buffPx1 = r2

        dgemv_("T", n, p, &pone, A, n, r, &one, &mone, buffPx1, &one); // buffPx1 = A' * r - buffPx1
        
        dgemm_("T","N", p,p,n,&pone, A, n, R, n, &zero, buffPxP, p); // buffPxP = A' * R

        dpotrf_( "L", p, buffPxP, p, &info ); // Cholesky Factor of buffPxP

        dpotrs_("L", p, &one, buffPxP, p, buffPx1, p, &info ); // buffPx1 = buffPxP ^-1 * buffPx1


        dcopy_(p, buffPx1, &one, dbeta, &one); //dbeta = buffPx1

        dcopy_(n, r, &one, dalpha, &one); // dalpha = r

        dgemv_("N", n,p,&mone, R, n, dbeta, &one, &pone, dalpha, &one); // dalpha = dalpha - R*dbeta
    
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
    int one = 1;
    *mult= 1.0;
    for (i=0;i<(*n);i++)
    {
        if (dalpha[i]<0.0) *mult = MIN(*mult,(-alpha[i]/dalpha[i]));
        if (dalpha[i]>0.0) *mult = MIN(*mult,(UminusAlpha[i]/dalpha[i]));
        if (dxi[i]<0.0)    *mult = MIN(*mult,(-xi[i]/dxi[i]));
        if (dzeta[i]<0.0)  *mult = MIN(*mult,(-zeta[i]/dzeta[i]));
    }
    *mult *= 0.99;
    daxpy_(n, mult, dalpha, &one, alpha, &one);
    daxpy_(p, mult, dbeta, &one, beta, &one);
    daxpy_(n, mult, dxi, &one, xi, &one);
    daxpy_(n, mult, dzeta, &one, zeta, &one);
}

/******************************************************************************/

void LowRankQP( int *n, int *m, int *p, int* method, int* verbose, int* niter, 
    double *Q, double *c, double *A, double *b, double *u, double *alpha,
    double* beta, double *xi, double *zeta)
{
    int i;
    int info = 0;
    int n2 = *n * *n;
    int np = *n * *p;
    int one = 1;

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
    double *dbeta  = (double *) calloc( (*p), sizeof(double) );
    double *dxi    = (double *) calloc( (*n), sizeof(double) );
    double *dzeta  = (double *) calloc( (*n), sizeof(double) );

    /* Some repeatedly occuring vectors */
    double *UminusAlpha     = (double *) calloc( *n, sizeof(double) );
    double *XiOnUminusAlpha = (double *) calloc( *n, sizeof(double) );
    double *ZetaOnAlpha     = (double *) calloc( *n, sizeof(double) );

    /* Some vectors used during calculations */
    double *w  = (double *) calloc( *m, sizeof(double) );
    double *r1 = (double *) calloc( *n, sizeof(double) );
    double *r2 = (double *) calloc( *p, sizeof(double) );;
    double *r3 = (double *) calloc( *n, sizeof(double) );
    double *r4 = (double *) calloc( *n, sizeof(double) );
    double *r5 = (double *) calloc( *n, sizeof(double) );
    double *D  = (double *) calloc( *n, sizeof(double) );
    double *r  = (double *) calloc( *n, sizeof(double) );
    double *R  = (double *) calloc( (*n)*(*p), sizeof(double) );;

    /* Various Buffers */
    double *buffPxP = (double *) calloc( (*p)*(*p), sizeof(double) );;
    double *buffPx1 = (double *) calloc( (*p), sizeof(double) );;

    double *M = (double *) calloc( (*n)*(*n), sizeof(double) );;


    /* Vectors to be created if p!=0 */


    /* Main Loop */
    if ( *verbose ) LRQPHeader();
    LRQPInitPoint( n, p, Q, c, A, b, u, alpha, beta, xi, zeta, w, r1 );

    for (i=0;i<(*niter);i++)
    {
        LRQPCalcStats( n, m, p, Q, c, A, b, u, alpha, beta, xi, zeta, dalpha,
            dbeta, dxi, dzeta, UminusAlpha, XiOnUminusAlpha, ZetaOnAlpha, w, r1,
            r2, D, &prim, &dual, &comp, &gap, &term, &mult, &t );

        if ( *verbose ) LRQPDisplay( i+1, &prim, &dual, &comp, &gap, &term  );
        if ( term < EPSTERM ) break;

        //LRQPFactorize( n, m, method, Q, D, M, pivN, buffNxM, P, Beta, Lambda,
        //    LambdaTemp, T );
        //MatrixMatrixCopy( M, Q, n, n );
        dcopy_(&n2, Q, &one, M, &one); // copy Q into M
        
        for (int i=0;i<(*n);i++) M[i+i*(*n)] += D[i];// diag(M) = Diag(M) + D

        dpotrf_( "L", n, M, n, &info ); // Cholesky factor of M

        //LRQPSolve( n, m, p, method, Q, D, A, R, M, pivN, buffMxP, P, Beta, Lambda );
        dcopy_(&np, A, &one, R, &one); // copy A to R
        dpotrs_("L", n, p, M, n, R, n, &info ); // R = M^-1 * R
        for (i=0;i<(*n);i++) r3[i] = 0;
        for (i=0;i<(*n);i++) r4[i] = 0;
        
        LRQPCalcDx( n, p, Q, c, A, b, u, alpha, beta, xi, zeta,
            dalpha, dbeta, dxi, dzeta, UminusAlpha, ZetaOnAlpha, 
            XiOnUminusAlpha, buffPxP, buffPx1, R, r, r1,
            r2, r3, r4, r5, D, M, &t, PRED);

        for (i=0;i<(*n);i++) r3[i] = ( t - (dalpha[i] * dzeta[i]) )/alpha[i];
        for (i=0;i<(*n);i++) r4[i] = ( t + (dalpha[i] * dxi[i]) )/UminusAlpha[i];
        
        LRQPCalcDx( n, p, Q, c, A, b, u, alpha, beta, xi, zeta,
            dalpha, dbeta, dxi, dzeta, UminusAlpha, ZetaOnAlpha,
            XiOnUminusAlpha, buffPxP, buffPx1, R, r, r1,
            r2, r3, r4, r5, D, M, &t, CORR);

        LRQPStep( n, p, alpha, beta, xi, zeta, dalpha, dbeta, dxi, dzeta,
            UminusAlpha, &mult );
    }
    
    LRQPSummary( i, *niter, &prim, &dual, &comp, &gap, &term );

    /* Free Memory */
    free(dalpha);      free(dxi);             free(dzeta);
    free(UminusAlpha); free(XiOnUminusAlpha); free(ZetaOnAlpha);
    free(r1); free(r3); free(r4); free(r5);
    free(D);  free(w);  free(r);
        
        free( dbeta );
        free( r2 );
        free( R );
        free( buffPxP );
        free( buffPx1 );

        free( M );
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
