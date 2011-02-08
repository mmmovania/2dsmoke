/*
 *  utility.h
 *  smoke
 *
 */

#define max(i,j) (i>j?i:j)
#define min(i,j) (i>j?j:i)

#define FOR_EVERY_X_FLOW(N)	for( int i=0; i<=N; i++ ) for( int j=0; j<N; j++ ) {
#define FOR_EVERY_Y_FLOW(N)	for( int i=0; i<N; i++ ) for( int j=0; j<=N; j++ ) {
#define FOR_EVERY_CELL(N)		for( int i=0; i<N; i++ ) for( int j=0; j<N; j++ ) {
#define END_FOR }

double **alloc2D( int n );
void free2D( double **ptr );