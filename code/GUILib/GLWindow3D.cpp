#include "GLApplication.h"

#include "GLWindow3D.h"
#include <Utils/Utils.h>
#include "GLUtils.h"

#include "GLContentManager.h"
#include <GUILib/GLMesh.h>
#include "GLCamera.h"
#include "GLTrackingCamera.h"

/**
	Default constructor
*/
GLWindow3D::GLWindow3D(int posX, int posY, int sizeX, int sizeY) : GLWindow(posX, posY, sizeX, sizeY) {
	init();
}

GLWindow3D::GLWindow3D() : GLWindow() {
	init();
}

void GLWindow3D::draw() {
	preDraw();
	drawScene();
	drawBorders();
	postDraw();
}

void GLWindow3D::init() {
	camera = new GLTrackingCamera();
}

void GLWindow3D::setupLights() {

	GLfloat bright[] = { 0.8f, 0.8f, 0.8f, 1.0f };
	GLfloat mediumbright[] = { 0.3f, 0.3f, 0.3f, 1.0f };

	glLightfv(GL_LIGHT1, GL_DIFFUSE, bright);
	glLightfv(GL_LIGHT2, GL_DIFFUSE, mediumbright);

	GLfloat light0_position[] = { 0.0f, 10000.0f, 10000.0f, 0.0f };
	GLfloat light0_direction[] = { 0.0f, -10000.0f, -10000.0f, 0.0f };

	GLfloat light1_position[] = { 0.0f, 10000.0f, -10000.0f, 0.0f };
	GLfloat light1_direction[] = { 0.0f, -10000.0f, 10000.0f, 0.0f };

	GLfloat light2_position[] = { 0.0f, -10000.0f, 0.0f, 0.0f };
	GLfloat light2_direction[] = { 0.0f, 10000.0f, -0.0f, 0.0f };

	glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
	glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
	glLightfv(GL_LIGHT2, GL_POSITION, light2_position);

	glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, light0_direction);
	glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, light1_direction);
	glLightfv(GL_LIGHT2, GL_SPOT_DIRECTION, light2_direction);

	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHT2);
}

void GLWindow3D::preDraw() {
	glPushAttrib(GL_ENABLE_BIT | GL_TRANSFORM_BIT | GL_VIEWPORT_BIT | GL_SCISSOR_BIT | GL_POINT_BIT | GL_LINE_BIT | GL_TRANSFORM_BIT);

	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);

	// set viewport
	glViewport(viewportX, viewportY, viewportWidth, viewportHeight);

	glMatrixMode(GL_PROJECTION);						// Select The Projection Matrix
	glPushMatrix();
	glLoadIdentity();									// Reset The Projection Matrix
	gluPerspective(45.0, (double)viewportWidth / viewportHeight, 0.05, 150.0);
	glMatrixMode(GL_MODELVIEW);							// Select The Modelview Matrix
	glPushMatrix();
	glLoadIdentity();									// Reset The Modelview Matrix

	glEnable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);

	if (camera)
		camera->applyCameraTransformations();

	setupLights();
}

void GLWindow3D::drawBorders() {
	//relative coordinates normalized to 0..1, with a square aspect ratio (i.e. largest dimension goes to 1, smaller dimension is relative to larger one)
	double minX = 0, minY = 0, maxX = (double)viewportWidth / MAX(viewportWidth, viewportHeight), maxY = (double)viewportHeight / MAX(viewportWidth, viewportHeight);

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	// Save projection matrix and sets it to a simple orthogonal projection
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(minX, maxX, minY, maxY);

	// Save model-view matrix and sets it to identity
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	// Draw borders
	glLineWidth(5);
	glColor3d(0, 0, 0);
	if (selected) glColor3d(1, .3, .3);
	glBegin(GL_LINE_LOOP);
	glVertex2d(minX, minY);
	glVertex2d(minX, maxY);
	glVertex2d(maxX, maxY);
	glVertex2d(maxX, minY);
	glEnd();

	glLineWidth(1);

	// Draw a bit of a border as a visual indication that the window is or is not active...
	if (isActive()) {
		glLineWidth(8);
		glColor3d(0, 0, 1);
		double s = MIN(maxX, maxY) * 0.02;
		glBegin(GL_LINES);
		glVertex2d(s, 5 * s);
		glVertex2d(s, s);
		glVertex2d(s, s);
		glVertex2d(5 * s, s);

		glVertex2d(s, maxY - 5 * s);
		glVertex2d(s, maxY - s);
		glVertex2d(s, maxY - s);
		glVertex2d(5 * s, maxY - s);

		glVertex2d(maxX - s, 5 * s);
		glVertex2d(maxX - s, s);
		glVertex2d(maxX - s, s);
		glVertex2d(maxX - 5 * s, s);

		glVertex2d(maxX - s, maxY - 5 * s);
		glVertex2d(maxX - s, maxY - s);
		glVertex2d(maxX - s, maxY - s);
		glVertex2d(maxX - 5 * s, maxY - s);
		glEnd();
	}

	glEnable(GL_DEPTH_TEST);

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
}

// clean up
void GLWindow3D::postDraw() {
	glMatrixMode(GL_MODELVIEW);							// Select The Modelview Matrix
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);						// Select The Projection Matrix
	glPopMatrix();
	glPopAttrib();
}

//all these methods should returns true if the event is processed, false otherwise...

//triggered when mouse moves
bool GLWindow3D::onMouseMoveEvent(double xPos, double yPos) {
	if (camera)
		if (camera->onMouseMoveEvent(xPos, yPos) == true)
			return true;
	return false;
}

//triggered when mouse buttons are pressed
bool GLWindow3D::onMouseButtonEvent(int button, int action, int mods, double xPos, double yPos) {
	if (camera) {
		if(camera->onMouseButtonEvent(button, action, mods, xPos, yPos))
			return true;
	}

	return false;
}

//triggered when using the mouse wheel
bool GLWindow3D::onMouseWheelScrollEvent(double xOffset, double yOffset) {
	if (camera)
		if (camera->onMouseWheelScrollEvent(xOffset, yOffset))
			return true;
	return false;
}

