#pragma once

#include "RotateWidget.h"

#include <MathLib/V3D.h>
#include <MathLib/Quaternion.h>
#include <MathLib/Ray.h>
#include <Utils/Logger.h>


class RotateWidgetV2;

class RotateWidgetV2Torus {
public:
	//pointer to the translate widget that this axis belongs to
	RotateWidgetV2* pWidget = NULL;

	double lastPickedAngle = 0;
	V3D zeroRotationAxis;
	V3D rotationAxis;

	// Last picked axis
	bool isPicked = false;

	void draw();

	void pickWith(const Ray& ray);
	double getAngleOnTorusForRay(const Ray& ray);

	RotateWidgetV2Torus(RotateWidgetV2* parentWidget, const V3D& axis) {
		pWidget = parentWidget;
		this->rotationAxis = axis;
		V3D t1, t2; rotationAxis.getOrthogonalVectors(t1, t2);
		zeroRotationAxis = t1;
	}
};

/**
This class draws a widget to translate and rotate a point. It also handles mouse
interactions.
*/
class RotateWidgetV2 : public RotateWidget {
public:
	//rotation. The rotation widgets remain fixed (world coords), but the overal orientation changes
	Quaternion orientation;
	// sub-widgets that actually do the work - the axis widget, together with the position, define the plane of the torus
	DynamicArray<RotateWidgetV2Torus> rAxes;
public:
	RotateWidgetV2(uint axes = AXIS_X | AXIS_Y | AXIS_Z);
	virtual ~RotateWidgetV2();

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

	void setOrientation(const Quaternion& q) {
		orientation = q;
	}

	Quaternion getOrientation() {
		return orientation;
	}
};

