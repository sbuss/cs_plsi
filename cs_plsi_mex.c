#include "cs_mex.h"
#include <string.h>

#ifndef CS_INT
#define CS_INT int
#endif

double *cs_mex_get_double_mat (CS_INT m, CS_INT n, const mxArray *X);
double *cs_mex_put_double_mat (CS_INT m, CS_INT n, const double *b, mxArray **X);

/* get a MATLAB dense column vector */
double *cs_mex_get_double_mat (CS_INT m, CS_INT n, const mxArray *X)
{
    cs_mex_check (0, m, n, 0, 0, 1, X) ;
    return (mxGetPr (X)) ;
}

/* return a double vector to MATLAB */
double *cs_dl_mex_put_double_mat (CS_INT m, CS_INT n, const double *b, mxArray **X)
{
    /* TODO: This might not work as expected */
    double *x ;
    CS_INT j;
    *X = mxCreateDoubleMatrix (m, n, mxREAL) ;      /* create x */
    x = mxGetPr (*X) ;
    for (j = 0; j < m*n; j++) {
        x [j] = b [j] ;       /* copy x = b */
    }
    return (x) ;
}

void mexFunction (
    int nargout,
    mxArray *pargout[],
    int nargin,
    const mxArray *pargin[]
)
{
    /* cplsi MATLAB interface
     * 
     * Required arguments:
     *   - cs X, the data matrix, should be m-by-n, where m is # users, n is # articles
     *   - int k, the number of Zs to find
     *
     * Optional arguments: (just kidding, they're all required!)
     *   - double B, the tempering factor, beta
     *   - pZ, pWgZ, pDgZ 
     *
     * Returns:
     *   - mxArray pZ, k-by-1
     *   - mxArray pWgZ, n-by-k
     *   - mxArray pDgZ, k-by-m
     *   - double logLike
     */

    cs_dl Xmatrix, *X;
    CS_INT m,n,K;
    double B;
    double *pZ,    *pWgZ,    *pDgZ;
    double *pZout, *pWgZout, *pDgZout;
    double logLike;
    bool debug = 0;
     
    if(nargout != 4 || nargin < 6) {
        /* Currently don't support not passing in pZ, pWgZ, pDgZ */
        mexErrMsgTxt("Usage: [pZo, pWgZo, pDgZo, ll] = cplsi(X,k,beta,pZ,pWgZ,pDgZ)");
    }

    /* Get Data */
    X = cs_dl_mex_get_sparse(&Xmatrix, 0, 1, pargin[0]);   /* Get X */
    m = X->m; n = X->n;
    if (debug) printf("m: %d, n: %d\n", m, n);

    K = ((int)mxGetScalar(pargin[1]));                  /* Get k */
    if (debug) printf("K is %d\n" , K);

    if(nargin > 2) {
        B = ((double)mxGetScalar(pargin[2]));           /* Get beta */
    } else {
        B = 1.0;
    }
    if (debug) printf("B is %f\n", B);
    
    pZ = cs_dl_mex_get_double(K,pargin[3]);                /* Get pZ */
    if (debug) printf("Got pZ\n");
    pWgZ = cs_mex_get_double_mat(K,m,pargin[4]);        /* Get pWgZ */
    if (debug) printf("Got pWgZ\n");
    pDgZ = cs_mex_get_double_mat(K,n,pargin[5]);        /* Get pDgZ */
    if (debug) printf("Got pDgZ\n");
    
    /* TODO: Generate initial pZ, pWgZ, pDgZ if not passed in */

    /* Allocate space for the "out" variables */
    /*pZout = cs_calloc(K, sizeof(double));
    memcpy(pZout,pZ,K*sizeof(double));
    if (debug) printf("Allocated pZout\n");
    pWgZout = cs_calloc(K*m, sizeof(double));
    memcpy(pWgZout,pWgZ,K*m*sizeof(double));
    if (debug) printf("Allocated pWgZout\n");
    pDgZout = cs_calloc(K*n, sizeof(double));
    memcpy(pDgZout,pDgZ,K*n*sizeof(double));
    if (debug) printf("Allocated pDgZout\n");*/
    /* ignore this, trigger recompile */
    if (debug) printf("Putting pZout\n");
    pZout = cs_dl_mex_put_double(K, pZ, &(pargout[0]));
    if (debug) printf("Putting pWgZout\n");
    pWgZout = cs_dl_mex_put_double_mat(K, m, pWgZ, &(pargout[1]));
    if (debug) printf("Putting pDgZout\n");
    pDgZout = cs_dl_mex_put_double_mat(K, n, pDgZ, &(pargout[2]));
    /* Call cplsi, and put results in the "out" variables */
    
    CS_INT one = 4;
    CS_INT two = 872622;
    printf("%ld\n", one * two);
    if (debug) printf("Calling cs_plsi\n");
    
    logLike = 0;
    logLike = cs_plsi(X, K, B, pZ, pWgZ, pDgZ, pZout, pWgZout, pDgZout);

    /* Return the "out" variables to matlab */
    if (debug) printf("Putting loglike\n");
    cs_dl_mex_put_double(1, &logLike, &(pargout[3]));
    if (debug) printf("Returning\n");

    return;
}
