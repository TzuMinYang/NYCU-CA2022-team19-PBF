#pragma once

#include <MathLib/V3D.h>
#include <MathLib/Ray.h>

/**
  * Used for any type of widget, window, etc, that expects input from mouse, keyboards, and so on...
  */
class InteractiveWidget {
protected:

	InteractiveWidget();
	virtual ~InteractiveWidget();
public:
	//draw function
	virtual void draw() = 0;

	//all these methods should returns true if the event is processed, false otherwise...
	//any time a physical key is pressed, this event will trigger. Useful for reading off special keys...
	virtual bool onKeyEvent(int key, int action, int mods) = 0;
	//this one gets triggered on UNICODE characters only...
	virtual bool onCharacterPressedEvent(int key, int mods) = 0;
	//triggered when mouse buttons are pressed
	virtual bool onMouseButtonEvent(int button, int action, int mods, double xPos, double yPos) = 0;
	//triggered when mouse moves
	virtual bool onMouseMoveEvent(double xPos, double yPos) = 0;
	//triggered when using the mouse wheel
	virtual bool onMouseWheelScrollEvent(double xOffset, double yOffset) = 0;

	Ray getRayFromScreenCoords(double x, double y);

};

