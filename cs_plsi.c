#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "cs.h"

#ifndef CS_INT
#define CS_INT int
#endif

bool converged(double llast, double ll);
bool closeEnough(double one, double two);
double step(const cs *X,        const CS_INT K,          const double B, 
            const double *pZin, const double *pWgZin, const double *pDgZin, 
            double *pZout,      double *pWgZout,      double *pDgZout);


/* Test for convergance */
bool converged(double llast, double ll) {
    printf("  %f, %e\n", ll, fabs(ll-llast));
    return fabs(ll - llast) < 1e-5;
}

/* Round off errors means we have to do ~= 0 */
bool closeEnough(double one, double two) {
    return abs(one - two) < 1e-50;
}

/* CS_PLSI , based on Probabilistic Latent Semantic Indexing, as described in
Hoffman 1999, on a sparse data matrix X for K many latent variables.

 W := vector of N many words {w_1, w_2, ..., w_N}
 D := vector of M many documents {d_1, d_2, ..., d_M}
 X := the (word, document) co-occurrence matrix:
      {x_nm} forall (n, m), where x_nm is 1 or 0
 Z := the 'true topic' matrix:
      {z_nm} forall (n, m), where each z_nm has 1-of-K representation. So, I'm
      denoting p(z_nmk = 1) as p(z_k | w_n, d_m)

 note that p(Z) is a Kx1 vector [p(z_1), p(z_2), ... p(z_K)] where p(z_k) is a
 scalar, while p(Z | W, D) is an NxM matrix where each element p(Z | W, D)_nm
 is a Kx1 vector p(Z | w_n, d_m).

 TODO: Fix this description of p(Z | W, D)
 p(Z | W, D) is realized as a NxMxK matrix where p(Z | W, D)_nmk = p(z_k| w_n,
 d_m), following the definition of p(z_k | w_n, d_m) given above.
 
 Initial values for p(Z), p(W | Z), and p(D | Z) are optional. If omitted, 
 they are randomly initialized.

 required args:
	cs X	-- NxM data matrix over N words {w_n} and M docs {d_m}, where
		   X_nm indicates the number of occurrences of word w_n in
		   document d_m.
	int K	-- # latent variable states to solve for

        double *pZout   -|
        double *pWgZout -+-- Allocated space to put the ML results
        double *pDgZout -|

 optional args:
	double b 	-- "control paramater" beta, see eq. (9)
	double *pZ	-- p(Z)		-- Kx1 vec; pZ(k) = p(z_k)
	double *pWgZ	-- p(W | Z)	-- NxK mat; pWgZ[n, k] = p(w_n | z_k)
	double *pDgZ	-- p(D | Z)	-- KxM mat; pDgZ[k, m] = p(d_m | z_k)
 returns:
        double *logLike -- The log likelihood of the solution
*/
double cs_plsi(const cs *X,        const CS_INT K,          const double B, 
               const double *pZin, const double *pWgZin, const double *pDgZin, 
               double *pZout,      double *pWgZout,      double *pDgZout) {
    /* It is worth noting that pWgZ should be transposed, and be passed in 
    as a kxn matrix. pDgZ should be normal, a kxm matrix. This is done for
    cache efficiency. */
    bool debug = 0;

    CS_INT n = X->n; /* Number of columns in X */
    CS_INT m = X->m; /* Number of rows in X */
    CS_INT k = 0;
    
    double ll, llast = 0;
    if(debug) printf("About to malloc\n");
    double *pZt = cs_malloc(K, sizeof(double));
    if(debug) printf("Malloc'd pZt\n");
    double *pWgZt = cs_malloc(m*K, sizeof(double));
    if(debug) printf("Malloc'd pWgZt\n");
    double *pDgZt = cs_malloc(n*K, sizeof(double));
    if(debug) printf("Malloc'd pDgZt\n");
    if(debug) printf("About to memcpy\n");
    memcpy(pZt, pZin, K * sizeof(double));
    if(debug) printf("Memcpy'd pZ\n");
    memcpy(pWgZt, pWgZin, K * m * sizeof(double));
    if(debug) printf("Memcpy'd pWgZ\n");
    memcpy(pDgZt, pDgZin, K * n * sizeof(double));
    if(debug) printf("Memcpy'd pDgZ\n");
    ll = step(X, K, B, pZt, pWgZt, pDgZt, pZout, pWgZout, pDgZout);
    CS_INT iter = 0;

    if(debug) printf("Got ll back, testing for convergance\n");

    while(!converged(llast, ll)) {
        if(debug) printf("About to memcpy\n");
        memcpy(pZt, pZout, K * sizeof(double));
        if(debug) printf("Memcpy'd pZ\n");
        memcpy(pWgZt, pWgZout, K * m * sizeof(double));
        if(debug) printf("Memcpy'd pWgZ\n");
        memcpy(pDgZt, pDgZout, K * n * sizeof(double));
        if(debug) printf("Memcpy'd pDgZ\n");

        for ( k = 0; k < K; k++) {
            pZout[k] = 0.0;
        }
        for ( k = 0; k < K*m; k++) {
            pWgZout[k] = 0.0;
        }
        for( k = 0; k < K*n; k++) {
            pDgZout[k] = 0.0;
        }
        llast = ll;
        ll = step(X, K, B, pZt, pWgZt, pDgZt, pZout, pWgZout, pDgZout);
        iter++;
        /*if(iter % 100 == 0) {*/
            printf("iter = %d\n", iter);
        /*}*/
    }
    cs_free(pZt);
    cs_free(pWgZt);
    cs_free(pDgZt);
    printf("iter = %d\n", iter);
    return ll;

}

double step(const cs *X,        const CS_INT K,          const double B, 
            const double *pZin, const double *pWgZin, const double *pDgZin, 
            double *pZout,      double *pWgZout,      double *pDgZout) {
    bool debug = 0;
    CS_INT p;
    CS_INT j;
    CS_INT k;
    CS_INT nn, mm;
    if(debug) printf("About to calloc\n");
    double *s = cs_calloc(K, sizeof(double));
    if(debug) printf("Calloc'd\n");
    double st;
    double pZnorm = 0;
    double ll = 0;
    double lls;
    double sum;

    
    double *Xx = X->x; /* non-zero values in X */
    CS_INT *Xp = X->p; /* # elements in column */
    CS_INT *Xi = X->i; /* row indices of X */
    CS_INT n = X->n; /* Number of columns in X */
    CS_INT m = X->m; /* Number of rows in X */
    /* As a data structure reference, the matrix :
        A = [ 4.5  0  3.2  0  ]
            [ 3.1 2.9  0  0.9 ]
            [  0  1.7 3.0  0  ]
            [ 3.5 0.4  0  1.0 ]

    Would be stored as:
        int p []    = { 0,             3,             6,        8,       10 };
        int i []    = { 0,   1,   3,   1,   2,   3,   0,   2,   1,   3   };
        double x [] = { 4.5, 3.1, 3.5, 2.9, 1.7, 0.4, 3.2, 3.0, 0.9, 1.0 };
    */
    if(debug) printf("About to enter loop\n");

    for ( j = 0; j < n; j++) { /* Column Index */
        for ( p = Xp[j]; p < Xp[j+1]; p++) { /* Index into Xx for X(n,m) */
            sum = 0;
            nn = j;
            mm = Xi[p];
            for ( k = 0; k < K; k++) {
                st = pow(pWgZin[k + K * mm] * pDgZin[k + K * nn], B) * pZin[k];
                sum += st;
                s[k] = st;
            }
            ll += log(sum) * Xx[p];
            for ( k = 0; k < K; k++) {
                st = (s[k] /* Xx[p]*/) / sum;
                pWgZout[k + K * mm] += st;
                pDgZout[k + K * nn] += st;
                pZout[k] += st;
                pZnorm += st;
            }
        }
    }

    if(debug) printf("About to normalize\n");

    CS_INT c = 0;
    /* Normalize pWgZ and pDgZ */
    for ( j = 0; j < m; j++) {
        c = K*j;
        for (k = 0; k < K; k++) {
            pWgZout[k + c] /= pZout[k];
        }
    }
    if(debug) printf("Normalized pWgZout\n");
    for ( j = 0; j < n; j++) {
        c = K*j;
        for (k = 0; k < K; k++) {
            pDgZout[k + c] /= pZout[k];
        }
    }
    if(debug) printf("Normalized pDgZout\n");
    
    /* Normalize pZ */
    for ( k = 0; k < K; k++) {
        pZout[k] /= pZnorm;
    }
    if (debug) {
        printf("pZ:\n");
        for ( k = 0; k < K; k++) {
            printf("%d, %f\n", k, pZout[k]);
        }
    }

    cs_free(s);
    if(debug) printf("Returning\n");

    return ll; 
}


