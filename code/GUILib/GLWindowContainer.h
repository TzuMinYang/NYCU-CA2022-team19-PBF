#pragma once
#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>
#include "GLWindow.h"
#include <Utils/Logger.h>
#include "GLWindow3D.h"
#include "GLCamera.h"
#include "GLApplication.h"
#include <GUILib/GLMesh.h>

class GLWindowContainer : public GLWindow2D {
protected:
	double startX = 0;
	double startY = 0;
public:
	int nRows, nCols;

	int width;
	int height;
	
	// sub windows
	DynamicArray<GLWindow*> subWindows;

	// constructor
	GLWindowContainer();
	GLWindowContainer(int nRows, int nCols, int posX, int posY, int sizeX, int sizeY);

	// destructor
	virtual ~GLWindowContainer(void);

	virtual void draw();

	void addSubWindow(GLWindow* subWindow);

	//triggered when mouse buttons are pressed
	virtual bool onMouseButtonEvent(int button, int action, int mods, double xPos, double yPos);
	//triggered when mouse moves
	virtual bool onMouseMoveEvent(double xPos, double yPos);
	//triggered when using the mouse wheel
	virtual bool onMouseWheelScrollEvent(double xOffset, double yOffset);
	//any time a physical key is pressed, this event will trigger. Useful for reading off special keys...
	virtual bool onKeyEvent(int key, int action, int mods);
	//this one gets triggered on UNICODE characters only...
	virtual bool onCharacterPressedEvent(int key, int mods);
};

