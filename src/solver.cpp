/*
 *  solver.cpp
 *  smoke
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "solver.h"
#include "utility.h"

const char *solver_name[] = { "Gauss-Seidel", "Conjugate Gradient", "Multigrid", NULL };
	
// Clamped Fetch
static double x_ref( double **x, int i, int j, int n ) {
	i = min(max(0,i),n-1);
	j = min(max(0,j),n-1);
	return x[i][j];
}

// Ans = Ax
static void compute_Ax( double **x, double **ans, int n ) {
	double h2 = 1.0/(n*n);
	for( int i=0; i<n; i++ ) for( int j=0; j<n; j++ ) {
		ans[i][j] = (x_ref(x,i+1,j,n)+x_ref(x,i-1,j,n)+x_ref(x,i,j+1,n)+x_ref(x,i,j-1,n)-4.0*x[i][j])/h2;
	}
}

// Gauss-Seidel Iteration
static void gaussseidel( double **x, double **b, int n, int t ) {
	double h2 = 1.0/(n*n);
	for( int k=0; k<t; k++ ) {
		for( int i=0; i<n; i++ ) for( int j=0; j<n; j++ ) {
			x[i][j] = (x_ref(x,i+1,j,n)+x_ref(x,i-1,j,n)+x_ref(x,i,j+1,n)+x_ref(x,i,j-1,n)-h2*b[i][j]) / 4.0;
		}
	}
}

// ans = x^T * x
static double product( double **x, double **y, int n ) {
	double ans = 0.0;
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			ans += x[i][j]*y[i][j];
		}
	}
	return ans;
}

// x = 0
static void clear( double **x, int n ) {
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			x[i][j] = 0.0;
		}
	}
}

// x <= y
static void copy( double **x, double **y, int n ) {
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			x[i][j] = y[i][j];
		}
	}
}
				 
// Ans = x + a*y
static void op( double **x, double **y, double **ans, double a, int n ) {
	static double **tmp = alloc2D(n);
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			tmp[i][j] = x[i][j]+a*y[i][j];
		}
	}
	copy(ans,tmp,n);
}

static void smooth( double **x, double **b, int n, int t ) {
	// Smooth Using Gaus-Seidel Method
	gaussseidel( x, b, n, t );
}

// r = b - Ax
static void residual( double **x, double **b, double **r, int n ) {
	compute_Ax(x,r,n);
	op( b, r, r, -1.0, n );
}

// Shrink the image
static void shrink( double **fine, double **coarse, int fn ) {
	for( int i=0; i<fn/2; i++ ) {
		for( int j=0; j<fn/2; j++ ) {
			// TODO: Interpolate Smoothly.
			coarse[i][j] = (fine[2*i][2*j]+fine[2*i+1][2*j]+fine[2*i][2*j+1]+fine[2*i+1][2*j+1]) / 4.0;
		}
	}
}

// Expand the image
static void expand( double **coarse, double **fine, int fn ) {
	for( int i=0; i<fn; i++ ) {
		for( int j=0; j<fn; j++ ) {
			// TODO: Interpolate Smoothly
			fine[i][j] = coarse[i/2][j/2];
		}
	}
}

// (V-Cycle Only) Multigrid Method
// TODO: Might Be Better Implement Full Multigrid Method
#define MAX_LAYER		8
static void mgv( double **x, double **b, int n, int recr=0 ) {
	
	// Memory Saving Part
	static double **fine_r[MAX_LAYER];
	static double **fine_e[MAX_LAYER];
	static double **coarse_r[MAX_LAYER];
	static double **coarse_e[MAX_LAYER];
	static double initialized = false;
	if( ! initialized  ) {
		for( int n=0; n<MAX_LAYER; n++ ) {
			fine_r[n] = NULL;
			fine_e[n] = NULL;
			coarse_r[n] = NULL;
			coarse_e[n] = NULL;
		}
		initialized = true;
	}
	
	if( ! fine_r[recr] ) fine_r[recr] = alloc2D(n);
	if( ! fine_e[recr] ) fine_e[recr] = alloc2D(n);
	if( ! coarse_r[recr] ) coarse_r[recr] = alloc2D(n/2);
	if( ! coarse_e[recr] ) coarse_e[recr] = alloc2D(n/2);
	
	clear(fine_r[recr],n);
	clear(fine_e[recr],n);
	clear(coarse_r[recr],n/2);
	clear(coarse_e[recr],n/2);
	
///////////// Beginning of V-Cycle
	
	// Pre-smoothing
	smooth( x, b, n, 4 );
	
	// Compute Residual
	residual( x, b, fine_r[recr], n );
	
	// Restrict
	shrink(fine_r[recr],coarse_r[recr],n);

	if( n <= 2 ) {
		// TODO: Should Be Solved Exactly
		smooth( coarse_e[recr], coarse_r[recr], n/2, 10 );
	} else {
		// Recursively Call Itself
		mgv(coarse_e[recr],coarse_r[recr],n/2,recr+1);
	}
	
	// Interpolate
	expand(coarse_e[recr],fine_e[recr],n);
	
	// Correct ( x = x + e )
	op( x, fine_e[recr], x, 1.0, n );
	
	// Post-smoothing
	smooth( x, b, n, 4 );
}

static void conjGrad( double **x, double **b, int n ) {
	// Pre-allocate Memory
	static double **r = alloc2D(n);
	static double **p = alloc2D(n);
	static double **Ap = alloc2D(n);
	clear(r,n);
	clear(p,n);
	clear(Ap,n);
	
	residual( x, b, r, n );					// r = b-Ax
	copy( p, r, n );						// p = r
	for( int k=0; k<n*n; k++ ) {
		compute_Ax( p, Ap, n );				// Ap
		double pAp = product( p, Ap, n );	// p^T * Ap
		double rr1 = product( r, r, n );	// r^T * r
		double a;
		if( pAp ) { 
			a = rr1/pAp;					// a = r^T * r / p^T * Ap
		} else break;
		op( x, p, x, a, n );				// x = x + a*p
		op( r, Ap, r, -a, n );				// r = r - a*Ap
		double rr2 = product( r, r, n );	// r1^T * r1
		if( rr2/n < 1.0e-8 ) break;
		if( rr1 ) {
			double b = rr2/rr1;
			op( r, p, p, b, n );			// p = r + b*p
		}
	}
}

double solver::solve( int method, int numiter, double **x, double **b, int n ) {
	static double **r = alloc2D(n);
	clear(r,n);
	
	switch(method) {
		case 0:
			// Gaus-Seidel
			gaussseidel(x,b,n,numiter);
			break;
		case 1:
			// Conjugate Gradient Method
			conjGrad(x,b,n);
			break;
		case 2:
			// Multigrid Method
			mgv(x,b,n);
			smooth( x, b, n, 8 );
			break;
	}
	residual( x, b, r, n );
	return sqrt(product( r, r, n ))/(n*n);
}