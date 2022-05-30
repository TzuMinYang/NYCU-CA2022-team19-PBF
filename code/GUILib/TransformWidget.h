#pragma once

#include "InteractiveWidget.h"

#include <MathLib/V3D.h>
#include <MathLib/Quaternion.h>
#include <MathLib/Ray.h>

enum TransformWidgetAxes {
	AXIS_X = 1,
	AXIS_Y = 2,
	AXIS_Z = 4
};

/**
	This class draws a widget to translate and rotate a point. It also handles mouse
	interactions.
*/
class TransformWidget : public InteractiveWidget {
public:
	// Widget scale
	double scale = 0.2;
	// Flag indicating if the widget is visible. If false, it is also unpickable
	bool visible = true;
	// Flag indicating if the widget is active, i.e., if the mouse is over it
	bool active = false;
	// Flag indicating if the widget is being dragged
	bool dragging = false;
public:
	TransformWidget();
	virtual ~TransformWidget();

	/**
		Draws the widget
	*/
	virtual void draw() = 0;
	virtual bool isPicked() = 0;

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
};

