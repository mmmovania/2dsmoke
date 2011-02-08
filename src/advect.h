/*
 *  advect.h
 *  smoke
 *
 *  Created by Ryoichi Ando on 2/7/11.
 *
 */

// Method:
// 0: Upwind
// 1: WENO5
// 2: QUICK
// 3: Semi-Lagrangian
// 4: MacCormack

// Interp:
// 0: Linear Interpolation
// 1: Cubic Interpolation

// u:
// Staggered Velocity Field

// n:
// Size of Velocity Field Grid Size

// cn:
// Size of Concentration Grid Size

// dt:
// Timestep Stride

extern const char *advection_name[];
extern const char *interp_name[];

namespace advect {
	void advect( int method, int interp, double ***u, double **c, int n, int cn, double dt );
}