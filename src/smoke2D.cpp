/*
 *  smoke2D.cpp
 *  smoke
 *
 */

#include "smoke2D.h"
#include "solver.h"
#include "utility.h"
#include "advect.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if defined(__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <sys/time.h>
#elif defined(WIN32)
#include "glut.h"
#include <windows.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#include <sys/time.h>
#endif

int	N;		// Fluid Grid Size
int	M;		// Smoke Grid Size

#define		DT		0.1				// Watch for a CFL Limit in case of Derivative Advection

#define NUM_ITER	500

static int solver_num = 2;
static int advection_num = 3;
static int interp_num = 1;

static double ***u = NULL;		// Access Bracket u[DIM][X][Y] ( Staggered Grid )
static double **c = NULL;		// Equivalent to c[N][N]
static double **p = NULL;		// Equivalent to p[N][N]
static double **d = NULL;		// Equivalent to d[N][N]
static double **vort = NULL;	// Equivalent to vort[N][N]

static double residual = 0.0;
static unsigned long solverTime = 0;
static unsigned long advectTime = 0;
static unsigned long simTime = 0;

static bool show_velocity = true;
static bool show_pressure = true;
static bool dragging = false;

void smoke2D::init( int gsize ) {
	N = gsize;
	M = gsize*2;
		
	// Allocate Variables
	if( ! p ) p = alloc2D(N);	
	if( ! d ) d = alloc2D(N);
	if( ! c ) c = alloc2D(M);
	if( ! vort ) vort = alloc2D(N);
	if( ! u ) {
		u = new double **[3];
		u[0] = alloc2D(N+1);
		u[1] = alloc2D(N+1);
	}
	
	// Clear Variables
	FOR_EVERY_X_FLOW(N) {
		u[0][i][j] = 0.0;
	} END_FOR
	
	FOR_EVERY_Y_FLOW(N) {
		u[1][i][j] = 0.0;
	} END_FOR
	
	FOR_EVERY_CELL(M) {
		c[i][j] = 0.0;
	} END_FOR
	
	// Turn On Blending
	glEnable(GL_BLEND);
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
}

void smoke2D::reshape( int w, int h ) {
	double margin = 0.0;
	glViewport(0, 0, w, h);
	glLoadIdentity();
	glOrtho(-margin,1.0+margin,-margin,1.0+margin,-1.0,1.0);
}

static unsigned long getMicroseconds()
{
#if defined(_WIN32)
	LARGE_INTEGER nFreq, Time;
	QueryPerformanceFrequency(&nFreq);
	QueryPerformanceCounter(&Time);
	return (double)Time.QuadPart / nFreq.QuadPart * 1000000;
#else
	struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec*1000000 + tv.tv_usec;
#endif
}

static unsigned long tickTime() {
	static unsigned long prevTime = getMicroseconds();
	unsigned long curTime = getMicroseconds();
	unsigned long res = curTime-prevTime;
	prevTime = curTime;
	return res;
}

void raw_drawBitmapString( const char *string)
{
	while (*string) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *string++);
}

static int cnt = 0;
static void drawBitmapString( const char *string )
{
	if( cnt++ ) raw_drawBitmapString( " + " );
	raw_drawBitmapString( string );
}

static void comp_divergence() {
	double h = 1.0/N;
	FOR_EVERY_CELL(N) {
		double div = (u[0][i+1][j]-u[0][i][j]) + (u[1][i][j+1]-u[1][i][j]);
		d[i][j] = div/h;
	} END_FOR
}

static void enforce_boundary() {
	FOR_EVERY_X_FLOW(N) {
		if( i==0 || i==N ) u[0][i][j] = 0.0;
	} END_FOR
	
	FOR_EVERY_Y_FLOW(N) {
		if( j==0 || j==N ) u[1][i][j] = 0.0;
	} END_FOR
}

static void compute_pressure() {
	// Clear Pressure
	FOR_EVERY_CELL(N) {
		p[i][j] = 0.0;
	} END_FOR
	
	tickTime();
	// Solve Ap = d ( p = Pressure, d = Divergence )
	residual = solver::solve( solver_num, NUM_ITER, p, d, N );
	solverTime = tickTime();
}

static void subtract_pressure() {
	double h = 1.0/N;
	FOR_EVERY_X_FLOW(N) {
		if( i>0 && i<N ) u[0][i][j] -= (p[i][j]-p[i-1][j])/h;
	} END_FOR
	
	FOR_EVERY_Y_FLOW(N) {
		if( j>0 && j<N ) u[1][i][j] -= (p[i][j]-p[i][j-1])/h;
	} END_FOR
}

static void decrease_concentration() {
	FOR_EVERY_CELL(N) {
		c[i][j] *= 0.997;
	} END_FOR
}

static void advection() {
	tickTime();
	advect::advect(advection_num,interp_num,u,c,N,M,DT);
	advectTime = tickTime();
}

static void vorticityConfinement() {
	double h = 1.0/N;
	double e = 0.1;
	static double ** vcAdd[2] = {alloc2D(N),alloc2D(N)};
	
	// Compute Vorticty
	FOR_EVERY_CELL(N) {
		vort[i][j] = 0.5*(u[1][i+1][j]-u[0][i][j] - (u[0][i][j+1]-u[0][i][j]))/h;
	} END_FOR
	
	FOR_EVERY_CELL(N) {
		if( i==0 || i==N-1 || j==0 || j==N-1 ) continue;
		double w = vort[i][j];
		double n[2] = { (fabs(vort[i+1][j])-fabs(vort[i-1][j]))*0.5/h, (fabs(vort[i][j+1])-fabs(vort[i][j-1]))*0.5/h };
		double len = hypot( n[0], n[1] );
		vcAdd[0][i][j] = 0.0;
		vcAdd[1][i][j] = 0.0;
		if( len )
		{
			double NL[2] = { n[0]/len, n[1]/len };
			double Nw[2] = { NL[1]*w, -NL[0]*w };
			vcAdd[0][i][j] = DT*e*h*Nw[0];
			vcAdd[1][i][j] = DT*e*h*Nw[1];
		}
	} END_FOR
	
	FOR_EVERY_X_FLOW(N) {
		if( i>0 && i<N-1 ) u[0][i][j] += 0.5*vcAdd[0][i][j]+0.5*vcAdd[0][i-1][j];
	} END_FOR
	
	FOR_EVERY_Y_FLOW(N) {
		if( j>0 && j<N-1 ) u[1][i][j] += 0.5*vcAdd[1][i][j]+0.5*vcAdd[1][i][j-1];
	} END_FOR
}

static void computeStep() {
	
	unsigned long startTime = getMicroseconds();
	
	enforce_boundary();
	comp_divergence();
	compute_pressure();
	subtract_pressure();
	decrease_concentration();
	advection();
	//vorticityConfinement();
	
	simTime = getMicroseconds()-startTime;
}

void smoke2D::display() {
	
	// Simulate One Step
	computeStep();
	
	// Draw Concentration
#if 1
	FOR_EVERY_CELL(M) {
		if( i == M-1 || j == M-1 ) continue;
		double h = 1.0/M;
		double p[2] = {i*h+h/2.0,j*h+h/2.0};
		double color[3] = { 0.4, 0.6, 1.0 };
		double ex = show_velocity && dragging ? 0.3 : 1.0;
		glBegin(GL_QUADS);
		glColor4d(color[0],color[1],color[2],c[i][j]*ex);
		glVertex2d(p[0],p[1]);
		glColor4d(color[0],color[1],color[2],c[i+1][j]*ex);
		glVertex2d(p[0]+h,p[1]);
		glColor4d(color[0],color[1],color[2],c[i+1][j+1]*ex);
		glVertex2d(p[0]+h,p[1]+h);
		glColor4d(color[0],color[1],color[2],c[i][j+1]*ex);
		glVertex2d(p[0],p[1]+h);
		glEnd();
	} END_FOR
#endif
	
	if( dragging && show_pressure ) {
#if 0
		// Draw Divergence
		FOR_EVERY_CELL(N) {
			double div = 10.0*d[i][j];
			glColor4d(div>0,0.0,div<0,fabs(div));
			
			double h = 1.0/N;
			double p[2] = {i*h,j*h};
			glBegin(GL_QUADS);
			glVertex2d(p[0],p[1]);
			glVertex2d(p[0]+h,p[1]);
			glVertex2d(p[0]+h,p[1]+h);
			glVertex2d(p[0],p[1]+h);
			glEnd();
		} END_FOR
#elif 1
		// Draw Pressure
		double minv = 1.0e8;
		double maxv = -1.0e8;
		FOR_EVERY_CELL(N) {
			if( p[i][j]<minv ) minv = p[i][j];
			if( p[i][j]>maxv ) maxv = p[i][j];
		} END_FOR
		 
		FOR_EVERY_CELL(N) {
			double press = 5000.0*p[i][j];//(p[i][j]-minv)/(maxv-minv) - 1.0;
			glColor4d(press>0,0.0,press<0,fabs(press));
			
			double h = 1.0/N;
			double p[2] = {i*h,j*h};
			glBegin(GL_QUADS);
			glVertex2d(p[0],p[1]);
			glVertex2d(p[0]+h,p[1]);
			glVertex2d(p[0]+h,p[1]+h);
			glVertex2d(p[0],p[1]+h);
			glEnd();
		} END_FOR
#elif 0
		// Draw Vorticy
		FOR_EVERY_CELL(N) {
			double vt = 0.5*vort[i][j];
			glColor4d(vt>0,0.0,vt<0,fabs(vt));
			
			double h = 1.0/N;
			double p[2] = {i*h,j*h};
			glBegin(GL_QUADS);
			glVertex2d(p[0],p[1]);
			glVertex2d(p[0]+h,p[1]);
			glVertex2d(p[0]+h,p[1]+h);
			glVertex2d(p[0],p[1]+h);
			glEnd();
		} END_FOR
#endif
	}
	
#if 0
	// Draw Vertical Grid
	glColor4d(1.0,1.0,1.0,1.0);
	glLineWidth(1.0);
	for( int i=0; i<N+1; i++ ) {
		double h = 1.0/N;
		glBegin(GL_LINES);
		glVertex2d(h*i,0.0);
		glVertex2d(h*i,1.0);
		glEnd();
	}
	// Draw horizontal Grid
	for( int j=0; j<N+1; j++ ) {
		double h = 1.0/N;
		glBegin(GL_LINES);
		glVertex2d(0.0,h*j);
		glVertex2d(1.0,h*j);
		glEnd();
	}
#endif
	
	// Draw X flow
#if 0
	glColor4d(0.0,0.0,1.0,1.0);
	FOR_EVERY_X_FLOW {
		double h = 1.0/N;
		double p[2] = {i*h,j*h+h/2.0};
		glBegin(GL_LINES);
		glVertex2d(p[0],p[1]);
		glVertex2d(p[0]+DT*u[0][i][j],p[1]);
		glEnd();
	} END_FOR
	
	// Draw Y Flow
	glColor4d(1.0,0.0,0.0,1.0);
	FOR_EVERY_Y_FLOW {
		double h = 1.0/N;
		double p[2] = {i*h+h/2.0,j*h};
		glBegin(GL_LINES);
		glVertex2d(p[0],p[1]);
		glVertex2d(p[0],p[1]+DT*u[1][i][j]);
		glEnd();
	} END_FOR
#endif

	if( show_velocity && dragging ) {
		// Draw Cell Center Flow
		glColor4d(1.0,1.0,0.0,0.8);
		FOR_EVERY_CELL(N) {
			double h = 1.0/N;
			double p[2] = {i*h+h/2.0,j*h+h/2.0};
			double v[2] = {0.5*u[0][i][j]+0.5*u[0][i+1][j],0.5*u[1][i][j]+0.5*u[1][i][j+1]};
			double s = 8.0;
			glBegin(GL_LINES);
			glVertex2d(p[0],p[1]);
			glVertex2d(p[0]+s*DT*v[0],p[1]+s*DT*v[1]);
			glEnd();
		} END_FOR
	}
	
	// Display Method Text
	glColor4f(1.0,1.0,1.0,1.0);
	glRasterPos2d(0.04, 0.065);
	char tmp[64];
	cnt = 0; // Reset Message Counter
	sprintf( tmp, "%s (Time=%.2fms, Residual=%.2e)", solver_name[solver_num], solverTime/(double)1000, residual );
	drawBitmapString(tmp);
	
	glRasterPos2d(0.04, 0.03);
	if( advection_num > 2 )
		sprintf( tmp, "%s (Time=%.2fms, Interp=%s)", advection_name[advection_num], advectTime/(double)1000, interp_name[interp_num] );
	else 
		sprintf( tmp, "%s (Time=%.2fms)", advection_name[advection_num], advectTime/(double)1000  );
	
	cnt = 0;
	drawBitmapString(tmp);
	
	glRasterPos2d(0.04, 0.95);
	sprintf( tmp, "SimTime/Frame=%.2fms", simTime/(double)1000 );
	raw_drawBitmapString(tmp);
	
	glRasterPos2d(0.04, 0.9);
	raw_drawBitmapString("Press \"a\" to switch advection scheme");
	
	glRasterPos2d(0.04, 0.87);
	raw_drawBitmapString("Press \"s\" to switch pressure solver");
	
	glRasterPos2d(0.04, 0.84);
	raw_drawBitmapString("Press \"i\" to switch Semi-Lagrangian interpolation method");
	
	glRasterPos2d(0.04, 0.81);
	raw_drawBitmapString("Press \"p\" to toggle pressure view");
	
	glRasterPos2d(0.04, 0.78);
	raw_drawBitmapString("Press \"v\" to toggle velocity view");
	
	glRasterPos2d(0.04, 0.75);
	raw_drawBitmapString("Press \"c\" to clear all");
}

void smoke2D::keyDown( unsigned char key ) {
	switch(key) {
		case 'c':
			init(N);
			break;
		case 'v':
			show_velocity = ! show_velocity;
			break;
		case 'p':
			show_pressure = ! show_pressure;
			break;
		case 's':
			solver_num ++;
			if( ! solver_name[solver_num] ) solver_num = 0;
			break;
		case 'a':
			advection_num ++;
			if( ! advection_name[advection_num] ) advection_num = 0;
			break;
		case 'i':
			interp_num ++;
			if( ! interp_name[interp_num] ) interp_num = 0;
			break;
		case '\e':
			exit(0);
			break;
	}
}

void smoke2D::mouse( double x, double y, int state ) {
	if( state == 1 ) {
		dragging = true;
		motion( x, y, 0.0, 0.0 );
	} else {
		dragging = false;
	}
}

void smoke2D::motion( double x, double y, double dx, double dy ) {
	int i = min(N-1,max(0,x*N));
	int j = min(N-1,max(0,y*N));
	double m = 1.0;
	u[0][i][j] = u[0][i+1][j] = min(m/N/DT,max(-m/N/DT,m*N*dx));
	u[1][i][j] = u[1][i][j+1] = min(m/N/DT,max(-m/N/DT,m*N*dy));
	
	i = min(M-1,max(0,x*M));
	j = min(M-1,max(0,y*M));
	int w = M/N;
	if( i>w && i<M-w-1 && j>w && j < M-w-1 ) {
		for( int ii = -w; ii <= w; ii++ ) {
			for( int jj = -w; jj <= w; jj++ ) {
				if( hypot(ii,jj) <= w ) {
					c[i+ii][j+jj] = 2.0;
				}
			}
		}
	}
}












