/*
 *  utility.cpp
 *  smoke
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

void copy2D( double **dst, double **src, int n ) {
	FOR_EVERY_CELL(n) {
		dst[i][j] = src[i][j];
	} END_FOR
}

void op2D( double **dst, double **src1, double **src2, double a, double b, int n ) {
	FOR_EVERY_CELL(n) {
		dst[i][j] = a*src1[i][j]+b*src2[i][j];
	} END_FOR
}
