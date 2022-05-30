#pragma once
#define GLFW_INCLUDE_GLU
#include "GLWindow.h"
#include <Utils/Logger.h>
#include <Utils/Timer.h>
#include "GLWindow2D.h"
#include "GLCamera.h"
#include <GUILib/GLMesh.h>
#include <GLFW/glfw3.h>

/**
* Used for any window that needs to be displayed in the OpenGL window
*/
class GLWindow3D : public GLWindow {
protected:
	//have a flag that shows if this window is active or not, as indicated by a mouse hovering over it
	double bgColorR = 1, bgColorG = 1, bgColorB = 1, bgColorA = 0.5;

	// interactive camera
	GLCamera* camera;

	// sets up the window for drawing
	virtual void preDraw();
	virtual void drawBorders();

	virtual void setupLights();

	// clean up
	virtual void postDraw();

	void init();
public:

	// constructors
	GLWindow3D(int x, int y, int w, int h);
	GLWindow3D();

	// destructor
	virtual ~GLWindow3D() {};

	// Draw the App scene - camera transformations, lighting, shadows, reflections, etc apply to everything drawn by this method
	virtual void draw();


	//input callbacks...

	//all these methods should returns true if the event is processed, false otherwise...
	//triggered when mouse buttons are pressed
	virtual bool onMouseButtonEvent(int button, int action, int mods, double xPos, double yPos);
	//triggered when mouse moves
	virtual bool onMouseMoveEvent(double xPos, double yPos);

	//triggered when using the mouse wheel
	virtual bool onMouseWheelScrollEvent(double xOffset, double yOffset);

	virtual void drawScene() {};
	GLCamera* getCamera() { return camera; }
};

