#include "GLWindow.h"
#include "GLApplication.h"
#include <Utils/Utils.h>
#include <MathLib/MathLib.h>
#include "GLShaderMaterial.h"

/**
	Default constructor
*/
GLWindow::GLWindow( int posX, int posY, int sizeX, int sizeY ) {
	setViewportParameters(posX, posY, sizeX, sizeY);
}

GLWindow::GLWindow() {
	setViewportParameters(0,0,100,100);
}

//sets viewport parameters
void GLWindow::setViewportParameters(int posX, int posY, int sizeX, int sizeY){
	this->viewportX = posX;
	this->viewportY = posY;
	this->viewportWidth = sizeX;
	this->viewportHeight = sizeY;
}
double GLWindow::getViewportXFromWindowX(double wX) {
	return wX - viewportX;
}

// NOTE: Viewport origin is bottom-left of the viewport, mouse origin is top-left of the containing window
double GLWindow::getViewportYFromWindowY(double wY) {
	return GLApplication::getGLAppInstance()->getMainWindowHeight() - wY - viewportY;
}

// converts viewport coordinates to a relative position within the sub window
double GLWindow::getRelativeXFromViewportX(double vpX) {
	return (vpX) / (MAX(viewportWidth, viewportHeight));
}

double GLWindow::getRelativeYFromViewportY(double vpY) {
	return (vpY) / (MAX(viewportWidth, viewportHeight));
}

double GLWindow::getViewportXFromRelativeX(double relX) {
	return relX*MAX(viewportWidth, viewportHeight);
}

double GLWindow::getViewportYFromRelativeY(double relY) {
	return relY*MAX(viewportWidth, viewportHeight);
}

bool GLWindow::mouseIsWithinWindow(double mouseX, double mouseY){
	double vX = getViewportXFromWindowX(mouseX);
	double vY = getViewportYFromWindowY(mouseY);

	//Logger::consolePrint("mouseX: %lf, mouseY: %lf, vX: %lf, vY: %lf\n", mouseX, mouseY, vX, vY);

	return vX >= 0 && vY >= 0 && vX <= viewportWidth && vY <= viewportHeight;
}


GLTexture * GLWindow::getGroundTexture() {
	return GLContentManager::getTexture("../data/textures/ground_TileLight2.bmp");
}

//draws a volumetric design environment
void GLWindow::drawDesignEnvironmentBox() {
	checkOGLErrors();
	double sizeX = 10 * designEnvironmentScale;
	double sizeY = 5 * designEnvironmentScale;
	double sizeZ = 10 * designEnvironmentScale;

	glDisable(GL_LIGHTING);

	GLShaderMaterial* shaderMaterial = GLContentManager::getShaderMaterial("radialGradient");
	shaderMaterial->apply();

	//don't draw the sides if they are in front of the camera...
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glBegin(GL_QUADS);

	//we will only draw the walls if they are not in the way of the camera

	// right side

	glTexCoord2d(sizeY, sizeZ);   glVertex3d(sizeX, sizeY, sizeZ);
	glTexCoord2d(sizeY, -sizeZ);  glVertex3d(sizeX, sizeY, -sizeZ);
	glTexCoord2d(-sizeY, -sizeZ); glVertex3d(sizeX, -sizeY, -sizeZ);
	glTexCoord2d(-sizeY, sizeZ);  glVertex3d(sizeX, -sizeY, sizeZ);


	// left side

	glTexCoord2d(sizeY, sizeZ);   glVertex3d(-sizeX, sizeY, sizeZ);
	glTexCoord2d(-sizeY, sizeZ);  glVertex3d(-sizeX, -sizeY, sizeZ);
	glTexCoord2d(-sizeY, -sizeZ); glVertex3d(-sizeX, -sizeY, -sizeZ);
	glTexCoord2d(sizeY, -sizeZ);  glVertex3d(-sizeX, sizeY, -sizeZ);


	// top side

	glTexCoord2d(sizeX, sizeZ);   glVertex3d(sizeX, sizeY, sizeZ);
	glTexCoord2d(-sizeX, sizeZ);  glVertex3d(-sizeX, sizeY, sizeZ);
	glTexCoord2d(-sizeX, -sizeZ); glVertex3d(-sizeX, sizeY, -sizeZ);
	glTexCoord2d(sizeX, -sizeZ);  glVertex3d(sizeX, sizeY, -sizeZ);


	// bottom side

	glTexCoord2d(sizeX, sizeZ);   glVertex3d(sizeX, -sizeY, sizeZ);
	glTexCoord2d(sizeX, -sizeZ);  glVertex3d(sizeX, -sizeY, -sizeZ);
	glTexCoord2d(-sizeX, -sizeZ); glVertex3d(-sizeX, -sizeY, -sizeZ);
	glTexCoord2d(-sizeX, sizeZ);  glVertex3d(-sizeX, -sizeY, sizeZ);


	// far side

	glTexCoord2d(sizeX, sizeY);   glVertex3d(sizeX, sizeY, sizeZ);
	glTexCoord2d(sizeX, -sizeY);  glVertex3d(sizeX, -sizeY, sizeZ);
	glTexCoord2d(-sizeX, -sizeY); glVertex3d(-sizeX, -sizeY, sizeZ);
	glTexCoord2d(-sizeX, sizeY);  glVertex3d(-sizeX, sizeY, sizeZ);


	// near side

	glTexCoord2d(sizeX, sizeY);   glVertex3d(sizeX, sizeY, -sizeZ);
	glTexCoord2d(-sizeX, sizeY);  glVertex3d(-sizeX, sizeY, -sizeZ);
	glTexCoord2d(-sizeX, -sizeY); glVertex3d(-sizeX, -sizeY, -sizeZ);
	glTexCoord2d(sizeX, -sizeY);  glVertex3d(sizeX, -sizeY, -sizeZ);


	glEnd();
	glDisable(GL_CULL_FACE);

	shaderMaterial->end();
	glDisable(GL_TEXTURE_2D);
}

void GLWindow::drawGround() {
	V3D t1, t2;
	Globals::groundPlane.n.getOrthogonalVectors(t1, t2);
	GLTexture* texture = getGroundTexture();
	if (texture) {
		glEnable(GL_TEXTURE_2D);
		texture->activate();
	}
	glNormal3d(Globals::groundPlane.n[0], Globals::groundPlane.n[1], Globals::groundPlane.n[2]);
	P3D p;
	double size = 100;
	glBegin(GL_QUADS);
	p = Globals::groundPlane.p + t1 * size + t2 * size;
	glTexCoord2d(size, size);
	glVertex3d(p[0], p[1], p[2]);

	p = Globals::groundPlane.p + t1 * -size + t2 * size;
	glTexCoord2d(-size, size);
	glVertex3d(p[0], p[1], p[2]);

	p = Globals::groundPlane.p + t1 * -size + t2 * -size;
	glTexCoord2d(-size, -size);
	glVertex3d(p[0], p[1], p[2]);

	p = Globals::groundPlane.p + t1 * size + t2 * -size;
	glTexCoord2d(size, -size);
	glVertex3d(p[0], p[1], p[2]);
	glEnd();
	glDisable(GL_TEXTURE_2D);
}

