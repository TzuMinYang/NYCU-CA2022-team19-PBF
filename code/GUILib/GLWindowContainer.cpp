#include "GLUtils.h"

#include "GLWindowContainer.h"
#include "GLWindow3D.h"
#include <Utils/Utils.h>

#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

#include "GLContentManager.h"
#include "GLWindowContainer.h"
#include "GLApplication.h"
#include <GUILib/GLMesh.h>
#include "GLCamera.h"
#include "GLTrackingCamera.h"
#include <iostream>

/**
Default constructor
*/
GLWindowContainer::GLWindowContainer(int nRows, int nCols, int posX, int posY, int sizeX, int sizeY) : GLWindow2D(posX, posY, sizeX, sizeY) {
	this->nRows = nRows;
	this->nCols = nCols;

	width = this->viewportWidth / nCols;
	height = this->viewportHeight / nRows;

	bgColorR = 0.6;
	bgColorG = 0.6;
	bgColorB = 0.6;
	bgColorA = 1.0;
}

GLWindowContainer::GLWindowContainer() : GLWindow2D() {
	nRows = 1;
	nCols = 5;
	width = this->viewportWidth / nCols;
	height = this->viewportHeight / nRows;
}

GLWindowContainer::~GLWindowContainer(void) {
	subWindows.clear();
}

void GLWindowContainer::addSubWindow(GLWindow* subWindow) {
	subWindow->setViewportParameters(viewportX + (subWindows.size() % nCols) * width, viewportY + (subWindows.size() / nCols) * height, width, height);
	subWindows.push_back(subWindow);
}

bool GLWindowContainer::onMouseButtonEvent(int button, int action, int mods, double xPos, double yPos){
	//update the mouse state for this window and all sub-windows...
	for (uint i = 0; i < subWindows.size(); i++) {
		if (action == GLFW_PRESS && !subWindows[i]->isActive())
			subWindows[i]->selected = false;
		subWindows[i]->updateSelectedFlag(button, action);
		subWindows[i]->updateClickedFlag(button, action);
	}

	if (!mouseIsWithinWindow(xPos, yPos)) {
		return false;
	}

	for (uint i = 0; i < subWindows.size(); i++) {
		if (subWindows[i]->isActive()) {
			if (subWindows[i]->onMouseButtonEvent(button, action, mods, xPos, yPos))
				return true;
		}
	}

	return true;
}

bool GLWindowContainer::onMouseMoveEvent(double xPos, double yPos) {
	for (uint i = 0; i < subWindows.size(); i++) {
		if (subWindows[i]->isClicked())
			if (subWindows[i]->onMouseMoveEvent(xPos, yPos) == true)
				return true;
	}
	return false;
}

void GLWindowContainer::draw() {
	preDraw();

	glScissor(viewportX, viewportY, viewportWidth, viewportHeight);
	glEnable(GL_SCISSOR_TEST);

	// because GLWindow2D disables depth test
	glEnable(GL_DEPTH_TEST);

	// draw subwindows and their content.
	for (uint i = 0; i < subWindows.size(); i++)
		subWindows[i]->draw();

	glDisable(GL_SCISSOR_TEST);

	postDraw();

}

bool GLWindowContainer::onMouseWheelScrollEvent(double xOffset, double yOffset) {
	for (uint i = 0; i < subWindows.size(); i++) {
		if (subWindows[i]->isActive()) {
			if (subWindows[i]->onMouseWheelScrollEvent(xOffset, yOffset))
				return true;
		}
	}

	if (isActive()) {
		for (uint i = 0; i < subWindows.size(); i++) 
			subWindows[i]->setViewportParameters((int)startX + viewportX + (i % nCols) * width, (int)startY + viewportY + (i / nCols) * height, width, height);
		return true;
	}
	return false;
}

bool GLWindowContainer::onKeyEvent(int key, int action, int mods){
	if ((key == GLFW_KEY_LEFT_SHIFT && action != GLFW_RELEASE) || mods & GLFW_MOD_SHIFT)
		Logger::consolePrint("%d %d %d\n", key, action, mods);

	return false;
}

bool GLWindowContainer::onCharacterPressedEvent(int key, int mods){
	return false;
}

