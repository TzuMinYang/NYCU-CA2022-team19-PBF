#include "GLWindow2D.h"
#include <Utils/Utils.h>
#include <MathLib/MathLib.h>
#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

/**
	Default constructor
*/
GLWindow2D::GLWindow2D( int posX, int posY, int sizeX, int sizeY ) : GLWindow(posX, posY, sizeX, sizeY) {
	minX = 0;
	maxX = (double)viewportWidth / MAX(viewportWidth, viewportHeight);
	minY = 0;
	maxY = (double)viewportHeight / MAX(viewportWidth, viewportHeight);
}

GLWindow2D::GLWindow2D() : GLWindow(){
	minX = 0;
	maxX = (double)viewportWidth / MAX(viewportWidth, viewportHeight);
	minY = 0;
	maxY = (double)viewportHeight / MAX(viewportWidth, viewportHeight);
}

// sets up the window for drawing
void GLWindow2D::preDraw() {
	glPushAttrib(GL_ENABLE_BIT | GL_TRANSFORM_BIT | GL_VIEWPORT_BIT | GL_SCISSOR_BIT | GL_POINT_BIT | GL_LINE_BIT | GL_TRANSFORM_BIT);

	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_DEPTH_TEST);

	glColor4d(1, 1, 1, 1);
	glPointSize(1);
	glLineWidth(1);

	glViewport(viewportX, viewportY, viewportWidth,viewportHeight);
	glScissor(viewportX, viewportY, viewportWidth, viewportHeight);
	glEnable(GL_SCISSOR_TEST);

	// Save projection matrix and sets it to a simple orthogonal projection
    glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(minX, maxX, minY, maxY);
	
	// Save model-view matrix and sets it to identity
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	// Draw on a clean background
	glClear(GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	glColor4d(bgColorR, bgColorG, bgColorB, bgColorA);

	glBegin(GL_QUADS);
	glVertex2d(minX, minY);
	glVertex2d(maxX, minY);
	glVertex2d(maxX, maxY);
	glVertex2d(minX, maxY);
	glEnd();

	glLineWidth(5);
	glColor3d(0, 0, 0);
	glBegin(GL_LINE_LOOP);
		glVertex2d(minX, minY);
		glVertex2d(minX, maxY);
		glVertex2d(maxX, maxY);
		glVertex2d(maxX, minY);
	glEnd();

	glLineWidth(1);

	// Draw a bit of a border as a visual indication that the window is or is not active...
	if (isActive() || isSelected()) {
		glColor3d(1 - bgColorR, 1 - bgColorG, 1 - bgColorB);
		double s = MIN(maxX, maxY) * 0.02;
		glBegin(GL_LINES);
			glVertex2d(s, 5*s);
			glVertex2d(s, s);
			glVertex2d(s, s);
			glVertex2d(5*s, s);

			glVertex2d(s, maxY-5*s);
			glVertex2d(s, maxY-s);
			glVertex2d(s, maxY-s);
			glVertex2d(5*s, maxY-s);

			glVertex2d(maxX - s, 5*s);
			glVertex2d(maxX - s, s);
			glVertex2d(maxX - s, s);
			glVertex2d(maxX - 5*s, s);

			glVertex2d(maxX - s, maxY - 5 * s);
			glVertex2d(maxX - s, maxY - s);
			glVertex2d(maxX - s, maxY - s);
			glVertex2d(maxX - 5*s, maxY - s);
		glEnd();
	}

}

// clean up
void GLWindow2D::postDraw() {
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
    glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glPopAttrib();

}

void GLWindow2D::draw() {
	preDraw();

//	glColor4d(bgColorR, bgColorG, bgColorB, bgColorA);
//	glBegin(GL_QUADS);
//	glVertex2d(minX, minY);
//	glVertex2d(maxX, minY);
//	glVertex2d(maxX, maxY);
//	glVertex2d(minX, maxY);
//	glEnd();

	postDraw();
}
