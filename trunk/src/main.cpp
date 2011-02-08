#include "smoke2D.h"
#include <stdio.h>
#include <stdlib.h>
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

int win_x = 500;
int win_y = 500;
int prev_x, prev_y;
int mstat = 1;

static void display(void) {
	glClear(GL_COLOR_BUFFER_BIT);
	smoke2D::display();
	glutSwapBuffers();
}

static void init( int gsize ) {
	glClearColor(0.0, 0.0, 0.0, 1.0);
	smoke2D::init(gsize);
}

static void reshape(int w, int h) {
	win_x = w;
	win_y = h;
	smoke2D::reshape(w,h);
}

static void idle ( void ) {
	glutPostRedisplay ();
}

static void keyboard( unsigned char key, int x, int y ) {
	smoke2D::keyDown(key);
}

static void mouse ( int button, int state, int x, int y ) {
	prev_x = x;
	prev_y = y;
	mstat = state;
	smoke2D::mouse( x/(GLdouble)win_x, 1.0 - y/(GLdouble)win_y, ! state );
}

static void motion ( int x, int y ) {
	if( mstat == 0 ) {
		smoke2D::motion( x/(GLdouble)win_x, 1.0 - y/(GLdouble)win_y,
						 (x-prev_x)/(GLdouble)win_x, -(y-prev_y)/(GLdouble)win_y );
	}
	prev_x = x;
	prev_y = y;
}


int main (int argc, char * argv[]) {
	
#ifdef _OPENMP
	int grid_size = 128;
#else
	int grid_size = 64;
#endif
	
	if( argc == 2  ) {
		sscanf( argv[1], "%d", &grid_size );
	}
	
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GL_DOUBLE);
	glutInitWindowPosition ( 100, 100 );
	glutInitWindowSize ( win_x, win_y );
	glutCreateWindow(argv[0]);
	glutIdleFunc(idle);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouse);
	glutMotionFunc (motion);
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	init(grid_size);
	glutMainLoop();
	return 0;
}
