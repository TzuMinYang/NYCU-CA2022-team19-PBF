#include "GlobalMouseState.h"
#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

double GlobalMouseState::lastMouseX = 0, GlobalMouseState::lastMouseY = 0;
double GlobalMouseState::mouseMoveX = 0, GlobalMouseState::mouseMoveY = 0;
bool GlobalMouseState::rButtonPressed = false, GlobalMouseState::lButtonPressed = false, GlobalMouseState::mButtonPressed = false;

void GlobalMouseState::updateMouseState(double xPos, double yPos, int button, int action, int mods) {
	lastMouseX = xPos;
	lastMouseY = yPos;

	if (button == GLFW_MOUSE_BUTTON_LEFT)
		lButtonPressed = action != GLFW_RELEASE;
	if (button == GLFW_MOUSE_BUTTON_MIDDLE)
		mButtonPressed = action != GLFW_RELEASE;
	if (button == GLFW_MOUSE_BUTTON_RIGHT)
		rButtonPressed = action != GLFW_RELEASE;
}

void GlobalMouseState::updateMouseMove(double xPos, double yPos) {
	mouseMoveX = lastMouseX - xPos;
	mouseMoveY = -lastMouseY + yPos;
	lastMouseX = xPos; lastMouseY = yPos;
}

