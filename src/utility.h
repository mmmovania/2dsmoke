/*
 *  utility.h
 *  smoke
 *
 */

#define max(i,j) (i>j?i:j)
#define min(i,j) (i>j?j:i)

#define FOR_EVERY_X_FLOW(N)	for( int xn=0; xn<(N+1)*N; xn++ ) { int i=xn%(N+1); int j=xn/(N+1);
#define FOR_EVERY_Y_FLOW(N)	for( int yn=0; yn<(N+1)*N; yn++ ) { int i=yn%N; int j=yn/N;
#define FOR_EVERY_CELL(N)	for( int ci=0; ci<N*N; ci++ ) { int i=ci%N; int j=ci/N;
#define END_FOR }

#ifdef _OPENMP
#include <omp.h>
#define OPENMP_FOR		_Pragma("omp parallel for" )
#define OPENMP_SECTION  _Pragma("omp section" )
#define OPENMP_BEGIN	_Pragma("omp parallel" ) {
#define OPENMP_END		}
#define OPENMP_FOR_P	_Pragma("omp for" )
#else
#define OPENMP_FOR
#define OPENMP_SECTION
#define OPENMP_BEGIN
#define OPENMP_END
#define OPENMP_FOR_P
#endif

double **alloc2D( int n );
void free2D( double **ptr );
void copy2D( double **dst, double **src, int n );
void op2D( double **dst, double **src1, double **src2, double a, double b, int n ); // dst = a*src1 + b*src2