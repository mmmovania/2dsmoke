/*
 *  smoke2D.h
 *  smoke
 *
 *  Created by mscp on 2/7/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

namespace smoke2D {
	void init( int gsize );
	void reshape( int w, int h );
	void display();
	void mouse( double x, double y, int state );
	void motion( double x, double y, double dx, double dy );
	void keyDown( unsigned char key );
}