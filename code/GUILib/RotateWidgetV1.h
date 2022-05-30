#pragma once

#include "RotateWidget.h"

#include <MathLib/V3D.h>
#include <MathLib/Quaternion.h>
#include <MathLib/Ray.h>
#include <Utils/Logger.h>


class RotateWidgetV1;

class RotateWidgetV1Axis {
public:
	//pointer to the translate widget that this axis belongs to
	RotateWidgetV1* pWidget = NULL;
	//the rotation axis, in local-coords, is the z-axis. 
	V3D localAxis = V3D(0,0,1);
	//We will keep track of this axis by storing its world rotation
	Quaternion rot;

	V3D getAxis() {
		return rot.rotate(localAxis);
	}

	void setRotationToMatch(const V3D& newAxis) {
		V3D rotAxis = newAxis.cross(localAxis);
		if (IS_ZERO(rotAxis.length()))
			rotAxis = V3D(1, 0, 0);
		else
			rotAxis.toUnit();

//		Logger::consolePrint("%lf %lf %lf\n", rotAxis[0], rotAxis[1], rotAxis[2]);

		double rotAngle = newAxis.angleWith(localAxis, rotAxis);
		rot = getRotationQuaternion(-rotAngle, rotAxis);
		assert(IS_ZERO((newAxis - getAxis()).length()));
	}

	// Last picked axis
	bool isPicked = false;

	void draw();

	void pickWith(const Ray& ray);

	RotateWidgetV1Axis(RotateWidgetV1* parentWidget) {
		pWidget = parentWidget;
	}
};


class RotateWidgetV1Torus {
public:
	//pointer to the translate widget that this axis belongs to
	RotateWidgetV1* pWidget = NULL;

	//we will keep track of the twist angle
	double twistAngle = 0;

	double lastPickedAngle = 0;
	V3D zeroAxis = V3D(1,0,0);

	// Last picked axis
	bool isPicked = false;

	void draw();

	void pickWith(const Ray& ray);
	double getAngleOnTorusForRay(const Ray& ray);


	RotateWidgetV1Torus(RotateWidgetV1* parentWidget) {
		pWidget = parentWidget;
	}
};

/**
This class draws a widget to translate and rotate a point. It also handles mouse
interactions.
*/
class RotateWidgetV1 : public RotateWidget {
public:
	// sub-widgets that actually do the work - the axis widget, together with the position, define the plane of the torus
	RotateWidgetV1Axis* axis;
	RotateWidgetV1Torus* torus;

public:
	RotateWidgetV1();
	virtual ~RotateWidgetV1();

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
//		computeEulerAnglesFromQuaternion(q, V3D(0, 0, 1), V3D(0, 1, 0), V3D(0, 0, 1), torus->twistAngle, axis->phi, axis->theta);
		//we could compute euler angles as above, but that will result in problems from singularities. Do it in a few different steps then...
		axis->setRotationToMatch(q.rotate(axis->localAxis));
		torus->twistAngle = (axis->rot.getComplexConjugate() * q).getRotationAngle(axis->localAxis);
	}

	Quaternion getOrientation() {
		//the orientation is given by first a twist, then a rotation that aligns tAxis to the z-axis.
		return axis->rot * getRotationQuaternion(torus->twistAngle, axis->localAxis);
	}
};

