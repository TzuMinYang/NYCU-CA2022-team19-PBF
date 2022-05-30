#pragma once

#include <MathLib/V3D.h>
#include <MathLib/Ray.h>

/**
  * Keep track of what the mouse is doing, accessible by any window or any app
  */
class GlobalMouseState {
public:
	static void updateMouseState(double xPos, double yPos, int button, int action, int mods);
	static void updateMouseMove(double xPos, double yPos);

	//keep track of the last mouse position
	static double lastMouseX, lastMouseY;
	static double mouseMoveX, mouseMoveY;
	static bool rButtonPressed, lButtonPressed, mButtonPressed;

	GlobalMouseState() {}

	~GlobalMouseState() {}

};

