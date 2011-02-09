/*
 *  advect.h
 *  smoke
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
// 1: Cubic Spline Interpolation
// 2: Monotonic Cubic Interpolation

// Integrator:
// 0: Euler
// 1: Modified Euler
// 2: Runge-Kutta

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
extern const char *integrator_name[];

namespace advect {
	void advect( int method, int interp, int integrator, double ***u, double **c, int n, int cn, double dt );
}