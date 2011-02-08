/*
 *  utility.cpp
 *  smoke
 *
 *  Created by Ryoichi Ando on 2/7/11.
 *
 */

#include "utility.h"
#include <stdio.h>
#include <stdlib.h>

double ** alloc2D( int n ) {
	double **ptr = new double *[n+1];
	for( int i=0; i<n; i++ ) {
		ptr[i] = new double[n+1];
		for( int j=0; j<n+1; j++ ) ptr[i][j] = 0.0;
	}
	ptr[n] = NULL;
	return ptr;
}

void free2D( double **ptr ) {
	for( int i=0; ptr[i]; i++ ) delete [] ptr[i];
	delete [] ptr;
}