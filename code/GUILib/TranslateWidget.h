#pragma once

#include "TransformWidget.h"

#include <MathLib/V3D.h>
#include <MathLib/Quaternion.h>
#include <MathLib/Ray.h>

class TranslateWidget;

class TranslateWidgetAxis{
public:
	V3D tAxis;
	//pointer to the translate widget that this axis belongs to
	TranslateWidget* pWidget = NULL;

	// Last picked axis
	bool isPicked = false;
	// Offset from widget origin to point where axis was picked - expressed in local coordinates of the widget
	V3D pickingOffset;

	void draw();

	void pickWith(const Ray& ray);

	TranslateWidgetAxis(const V3D& axis, TranslateWidget* parentWidget) {
		tAxis = axis;
		pWidget = parentWidget;
	}
};


/**
	This class draws a widget to translate and rotate a point. It also handles mouse
	interactions.
*/
class TranslateWidget : public TransformWidget {
public:
	// Translation
	P3D pos;
	// Orientation
	Quaternion orientation;

	// Axes to be controlled
	DynamicArray<TranslateWidgetAxis> tAxes;

	// This stores the applied translation or rotation
	P3D pickedPoint;

public:
	TranslateWidget(uint axes = AXIS_X | AXIS_Y | AXIS_Z);
	virtual ~TranslateWidget();

	/**
		Draws the widget
	*/
	virtual void draw();
	virtual bool isPicked();


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
	virtual bool onMouseMoveEvent(double xPos, double yPos);
	//triggered when using the mouse wheel
	virtual bool onMouseWheelScrollEvent(double xOffset, double yOffset) { return false; }
};

