#pragma once

#include "TransformWidget.h"

#include <MathLib/V3D.h>
#include <MathLib/Quaternion.h>
#include <MathLib/Ray.h>
#include <Utils/Logger.h>


/**
This class draws a widget to translate and rotate a point. It also handles mouse
interactions.
*/
class RotateWidget : public TransformWidget {
public:
	// position in world
	P3D pos;

public:
	RotateWidget();
	virtual ~RotateWidget();

	/**
	Draws the widget
	*/
	virtual void draw() = 0;
	virtual bool isPicked() = 0;

	//all these methods should returns true if the event is processed, false otherwise...
	//any time a physical key is pressed, this event will trigger. Useful for reading off special keys...
	virtual bool onKeyEvent(int key, int action, int mods) { return false; }
	//this one gets triggered on UNICODE characters only...
	virtual bool onCharacterPressedEvent(int key, int mods) { return false; }
	//triggered when mouse buttons are pressed
	virtual bool onMouseButtonEvent(int button, int action, int mods, double xPos, double yPos) {
		return active;
	}
	//triggered when mouse moves
	virtual bool onMouseMoveEvent(double xPos, double yPos) { return false; }
	//triggered when using the mouse wheel
	virtual bool onMouseWheelScrollEvent(double xOffset, double yOffset) { return false; }

	virtual void setOrientation(const Quaternion& q) = 0;

	virtual Quaternion getOrientation() = 0;
};

