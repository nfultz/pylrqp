
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#define MAX(A,B) ( (A) > (B) ? (A):(B))
#define MIN(A,B) ( (A) < (B) ? (A):(B))

#define CHOL 2

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

/*******************/

struct IterVars {
    double mult;
    double prim;
    double dual;
    double comp;
    double gap ;
    double term;
    double t   ;
};

struct sol {
    double* alpha;
    double* beta;
    double* xi;
    double* zeta;
};


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
    double *b, double *u, struct sol *solution, double *w)
{
    int i;
    int one = 1;
    int izero = 0;
    double pone =  1.0;
    double mone = -1.0;
    double zero =  0.0;

    for (i=0;i<(*n);i++) solution->alpha[i] = u[i] > 1 ? EPSINIT : u[i] * EPSINIT;

    dcopy_(p, &zero, &izero, solution->beta, &one); // beta = 0

    dgemv_("T", n, n, &mone, Q, n, solution->alpha, &one, &zero, w, &one ); // w = - Q' * alpha
    daxpy_(n, &mone, c, &one, w, &one); //w = w - c
    for (i=0;i<(*n);i++)
    {
        solution->xi[i]   = MAX(EPSINIT,w[i]);
        solution->zeta[i] = MAX(EPSINIT,solution->xi[i]-w[i]);
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
*    w = - c - w   - A*beta - xi + zeta;
*    quad = alpha'*w;
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
    double *b, double *u,
    struct sol *solution,
    double *UminusAlpha, double *XiOnUminusAlpha, double *ZetaOnAlpha, 
    double* w, double *r2, double *D, struct IterVars *ivars)
{
    int i;
    int one = 1;
    double quad;
    double cTalpha;
    double temp;
    double pone =  1.0;
    double mone = -1.0;
    double zero =  0.0;

    dgemv_("T", n, m, &pone, Q, n, solution->alpha, &one, &zero, w, &one ); // w = Q' * alpha
    quad = ddot_(n, solution->alpha, &one, w, &one );


    dcopy_(n, u, &one, UminusAlpha, &one);  // UminusAlpha = u
    daxpy_(n, &mone, solution->alpha, &one, UminusAlpha, &one); //// UminusAlpha = UminusAlpha + -1 * alpha

    for (i=0;i<(*n);i++) XiOnUminusAlpha[i] = solution->xi[i]/UminusAlpha[i];
    for (i=0;i<(*n);i++) ZetaOnAlpha[i] = solution->zeta[i]/solution->alpha[i];


        dgemv_("N", n, p, &mone, A, n, solution->beta, &one, &mone, w, &one); // w = -A * beta -w
        daxpy_(n, &mone, c, &one, w, &one ); // w= w - c
        daxpy_(n, &mone, solution->xi, &one, w, &one ); // w= w - xi
        daxpy_(n, &pone, solution->zeta, &one, w, &one ); // w= w + zeta



        dcopy_(p, b, &one, r2, &one); // r2 = b

        dgemv_("T", n, p, &mone, A, n, solution->alpha, &one, &pone, r2, &one); // r2 = r2 - A' * alpha

    ivars->dual = dasum_(p, r2, &one ); //sum(abs(r2))

    ivars->prim   = dasum_(n, w, &one ); // sum(abs(w))

    ivars->comp   = ddot_(n,  solution->alpha, &one, solution->zeta,  &one ) + ddot_(n, UminusAlpha, &one, solution->xi,  &one );

    cTalpha = ddot_(n, c, &one, solution->alpha, &one );

    ivars->gap = fabs( quad + cTalpha + ddot_(n, u, &one, solution->xi, &one ) + ddot_(p, b, &one, solution->beta, &one ) );
    ivars->term   = ivars->comp / ( fabs( 0.5*quad + cTalpha) + 1.0);
    temp    = (1.0 - ivars->mult + EPSIPM)/(10.0 + ivars->mult);
    ivars->t      = ivars->comp*(temp*temp)/(2*(*n));
    for (i=0;i<(*n);i++) D[i] = XiOnUminusAlpha[i] + ZetaOnAlpha[i] + EPSPERT;
}

/*****************************************************************************/

void LRQPDisplay( int i, struct IterVars *ivars)
{
    printf("%3d %15.7e %15.7e %15.7e %15.7e %15.7e \n", i,
            ivars->prim,
            ivars->dual,
            ivars->comp,
            ivars->gap,
            ivars->term );
}

/*****************************************************************************/

void LRQPSummary( int i, int niter, struct IterVars *ivars)
{
    if (i==niter)
    {
        printf("LowRankQP FAILED TO CONVERGE\n");
        printf("    Try increasing niter, or rescaling problem.\n");
    }
    else
    {
        printf("LowRankQP CONVERGED IN %d ITERATIONS\n\n", i+1 );
        printf("    Primal Feasibility    = %15.7e\n", ivars->prim);
        printf("    Dual Feasibility      = %15.7e\n", ivars->dual);
        printf("    Complementarity Value = %15.7e\n", ivars->comp);
        printf("    Duality Gap           = %15.7e\n", ivars->gap);
        printf("    Termination Condition = %15.7e\n", ivars->term);
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
*    r     = w + r3 - r4
*    r      = QpIpmSolve( info, Q, D, r );
*    Dbeta  = (A'*R) \ (A'*r - r2);
*    Dalpha = r - R*Dbeta;
*    Dzeta  = r3 - ZetaOnAlpha.*Dalpha;
*    Dxi    = r4 + XiOnUminusAlpha.*Dalpha;
*
******************************************************************************/

void LRQPCalcDx( int *n, int *p, double *Q, double *c,
    double *A, double *b, double * u,
    struct sol *solution,
    struct sol *step,
    double *UminusAlpha, double *ZetaOnAlpha, double *XiOnUminusAlpha,
    double* buffPxP, double *buffPx1,
    double *R, double *r, double *w, double* r2, double *r3,
    double *r4, double *M, double *t)
{

    int i, j;
    int    info = 0;
    int    one  = 1;
    double pone =  1.0;
    double mone = -1.0;
    double zero =  0.0;


    daxpy_(n, &mone, solution->zeta, &one, r3, &one); //r3 -= zeta
    daxpy_(n, &mone, solution->xi,   &one, r4, &one); // r4 -= xi

    for (i=0;i<(*n);i++) r[i] = w[i] + r3[i] - r4[i];

        //LRQPSolve( n, m, &one, method, Q, D, r5, r, M, pivN, buffMx1, P, Beta, Lambda );
        //dcopy_(n, r5, &one, r, &one); // copy r5 to r
        dpotrs_("L", n, &one, M, n, r, n, &info ); // r = M^-1 * r

        dcopy_(p, r2, &one, buffPx1, &one); // buffPx1 = r2

        dgemv_("T", n, p, &pone, A, n, r, &one, &mone, buffPx1, &one); // buffPx1 = A' * r - buffPx1
        
        dgemm_("T","N", p,p,n,&pone, A, n, R, n, &zero, buffPxP, p); // buffPxP = A' * R

        dpotrf_( "L", p, buffPxP, p, &info ); // Cholesky Factor of buffPxP

        dpotrs_("L", p, &one, buffPxP, p, buffPx1, p, &info ); // buffPx1 = buffPxP ^-1 * buffPx1


        dcopy_(p, buffPx1, &one, step->beta, &one); //dbeta = buffPx1

        dcopy_(n, r, &one, step->alpha, &one); // dalpha = r

        dgemv_("N", n,p,&mone, R, n, step->beta, &one, &pone, step->alpha, &one); // dalpha = dalpha - R*dbeta
    
    for (i=0;i<(*n);i++) step->zeta[i] = r3[i] - (ZetaOnAlpha[i] * step->alpha[i]);
    for (i=0;i<(*n);i++) step->xi[i]   = r4[i] + (XiOnUminusAlpha[i] * step->alpha[i]);

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

void LRQPStep( int *n, int *p, struct sol *solution,
        struct sol *step,
    double *UminusAlpha, double *mult)
{
    int i;
    int one = 1;
    *mult= 1.0;
    for (i=0;i<(*n);i++)
    {
        if (step->alpha[i]<0.0) *mult = MIN(*mult,(-solution->alpha[i]/step->alpha[i]));
        if (step->alpha[i]>0.0) *mult = MIN(*mult,(UminusAlpha[i]/step->alpha[i]));
        if (step->xi[i]<0.0)    *mult = MIN(*mult,(-solution->xi[i]/step->xi[i]));
        if (step->zeta[i]<0.0)  *mult = MIN(*mult,(-solution->zeta[i]/step->zeta[i]));
    }
    *mult *= 0.99;
    daxpy_(n, mult, step->alpha, &one, solution->alpha, &one);
    daxpy_(p, mult, step->beta,  &one, solution->beta,  &one);
    daxpy_(n, mult, step->xi,    &one, solution->xi,    &one);
    daxpy_(n, mult, step->zeta,  &one, solution->zeta,  &one);
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
    int np1 = *n + 1;
    int one = 1;
    int izero = 0;
    double pone =  1.0;
    double zero =  0.0;

    /* Iteration Display variables */
    struct IterVars ivars = {
      .mult = 0.0,
      .prim = 0.0,
      .dual = 0.0,
      .comp = 0.0,
      .gap  = 0.0,
      .term = 0.0,
      .t    = 0.0
    };

    struct sol solution = { alpha, beta, xi, zeta };
    struct sol step;

    /* Step direction vectors */
    step.alpha = (double *) calloc( (*n), sizeof(double) );
    step.beta  = (double *) calloc( (*p), sizeof(double) );
    step.xi    = (double *) calloc( (*n), sizeof(double) );
    step.zeta  = (double *) calloc( (*n), sizeof(double) );

    /* Some repeatedly occuring vectors */
    double *UminusAlpha     = (double *) calloc( *n, sizeof(double) );
    double *XiOnUminusAlpha = (double *) calloc( *n, sizeof(double) );
    double *ZetaOnAlpha     = (double *) calloc( *n, sizeof(double) );

    /* Some vectors used during calculations */
    double *w  = (double *) calloc( *m, sizeof(double) );
    double *r2 = (double *) calloc( *p, sizeof(double) );;
    double *r3 = (double *) calloc( *n, sizeof(double) );
    double *r4 = (double *) calloc( *n, sizeof(double) );
    double *D  = (double *) calloc( *n, sizeof(double) );
    double *r  = (double *) calloc( *n, sizeof(double) );
    double *R  = (double *) calloc( (*n)*(*p), sizeof(double) );;

    /* Various Buffers */
    double *buffPxP = (double *) calloc( (*p)*(*p), sizeof(double) );;
    double *buffPx1 = (double *) calloc( (*p), sizeof(double) );;

    double *M = (double *) calloc( (*n)*(*n), sizeof(double) );;


    /* Main Loop */
    if ( *verbose ) LRQPHeader();
    LRQPInitPoint( n, p, Q, c, A, b, u, &solution, w);

    for (i=0;i<(*niter);i++)
    {
        LRQPCalcStats( n, m, p, Q, c, A, b, u, &solution,
            UminusAlpha, XiOnUminusAlpha, ZetaOnAlpha, w,
            r2, D, &ivars);

        if ( *verbose ) LRQPDisplay( i+1, &ivars);
        if ( ivars.term < EPSTERM ) break;

        //LRQPFactorize( n, m, method, Q, D, M, pivN, buffNxM, P, Beta, Lambda,
        //    LambdaTemp, T );
        //MatrixMatrixCopy( M, Q, n, n );
        dcopy_(&n2, Q, &one, M, &one); // copy Q into M
        
        daxpy_(n,&pone,D, &one, M, &np1); // diag(M) = Diag(M) + D

        dpotrf_( "L", n, M, n, &info ); // Cholesky factor of M

        //LRQPSolve( n, m, p, method, Q, D, A, R, M, pivN, buffMxP, P, Beta, Lambda );
        dcopy_(&np, A, &one, R, &one); // copy A to R
        dpotrs_("L", n, p, M, n, R, n, &info ); // R = M^-1 * R

        dcopy_(n,&zero, &izero, r3, &one); // r3 = 0
        dcopy_(n,&zero, &izero, r4, &one); // r4 = 0
        
        LRQPCalcDx( n, p, Q, c, A, b, u, &solution,
            &step, UminusAlpha, ZetaOnAlpha, 
            XiOnUminusAlpha, buffPxP, buffPx1, R, r, w,
            r2, r3, r4, M, &ivars.t);

        for (i=0;i<(*n);i++) r3[i] = ( ivars.t - (step.alpha[i] * step.zeta[i]) )/solution.alpha[i];
        for (i=0;i<(*n);i++) r4[i] = ( ivars.t + (step.alpha[i] * step.xi[i]) )/UminusAlpha[i];
        
        LRQPCalcDx( n, p, Q, c, A, b, u, &solution,
            &step, UminusAlpha, ZetaOnAlpha,
            XiOnUminusAlpha, buffPxP, buffPx1, R, r, w,
            r2, r3, r4, M, &ivars.t);

        LRQPStep( n, p, &solution, &step, UminusAlpha, &ivars.mult );
    }
    
    LRQPSummary( i, *niter, &ivars);

    /* Free Memory */
    free(step.alpha);      free(step.xi);             free(step.zeta);
    free(UminusAlpha); free(XiOnUminusAlpha); free(ZetaOnAlpha);
    free(r3); free(r4);
    free(D);  free(w);  free(r);
        
        free( step.beta );
        free( r2 );
        free( R );
        free( buffPxP );
        free( buffPx1 );

        free( M );
}

int main(){
  int n = 3, m = 3, p = 1, method = 2,  verbose=0, niter=200;
  
  double alpha[] = {0,0,0,0,0,0};
  double  beta[] = {0,0,0,0,0,0};
  double    xi[] = {0,0,0,0,0,0};
  double  zeta[] = {0,0,0,0,0,0};

  double Vmat[] = {7, 0.627, -1.187, 0.627, 0.056, -0.106, -1.187, -0.106, 
  1.452};
  double dvec[] = {-4.931, -0.442, 1.335};
  double Amat[] = {1, 1, 1};
  double bvec[] = {1};
  double uvec[] = {1,1,1};

  LowRankQP( &n, &m, &p, &method, &verbose, &niter,
              Vmat, dvec, Amat, bvec, uvec,
              alpha, beta, xi, zeta);

  printf("Sum of dvec: %f\n", (dasum_(&n, dvec, &p)));

  PrintMatrix("alpha", alpha, &n, &p);
  PrintMatrix("beta", beta, &p, &p);
  PrintMatrix("xi", xi, &n, &p);
  PrintMatrix("zeta", zeta, &n, &p);

}
