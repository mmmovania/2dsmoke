/*
 *  solver.h
 *  smoke
 *
 *  Created by Ryoichi Ando on 2/7/11.
 *
 */

// Method:
// 0: Gauss-Seldel Method
// 1: Conjugate Gradient Method
// 2: Multigrid V-Cycle

extern const char *solver_name[];

namespace solver {
	// Solve Ax = b
	// RETURN: Residual
	
	// NOTICE: A is a Nullspace Matrix
	double solve( int method, int numiter, double **x, double **b, int n );
}
