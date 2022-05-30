#pragma once

#include <GUILib/GLCamera.h>

class GLTrackingCamera : public GLCamera{
private:

public:
	//explicitely keep track of the rotations about the "up" and "horizontal right" directions to avoid singularities in decomposition of large rotations
	double rotAboutUpAxis = 0, rotAboutRightAxis = 0;
	//the location of the camera is determined as a function of its orientation, its target, and how far along the target it is
	double camDistance;
	//this is the target of the camera, in world coordinates
	P3D camTarget;
	//we need to keep track of the view direction and up axis
	V3D camViewDirection, camUpAxis;

public:
	GLTrackingCamera(double distToTarget = -4, const V3D& cameraViewDirection = V3D(0, 0, -1), const V3D& cameraUpAxis = V3D(0, 1, 0)); // y-up
//	GLTrackingCamera(double distToTarget = -4, const V3D& cameraViewDirection = V3D(-1, 0, 0), const V3D& cameraUpAxis = V3D(0, 0, 1)); // z-up
	
	virtual ~GLTrackingCamera(void);

	virtual Quaternion getRotationToOpenGLCoordinateSystem() {
		Quaternion q;
		q.setRotationFrom(camUpAxis.cross(-camViewDirection), camUpAxis, -camViewDirection);
		return q.getComplexConjugate();
	}

	virtual P3D getCameraPosition();
	virtual P3D getCameraTarget();
	virtual V3D getWorldUpAxis();
	virtual Quaternion getCameraRotation();
	virtual void setCameraTarget(const P3D& p);
	virtual void setRotations(const V3D& r);

	//all these methods should returns true if the event is processed, false otherwise...
	//any time a physical key is pressed, this event will trigger. Useful for reading off special keys...
	virtual bool onKeyEvent(int key, int action, int mods) { return false; }
	//this one gets triggered on UNICODE characters only...
	virtual bool onCharacterPressedEvent(int key, int mods) { return false; }
	//triggered when mouse buttons are pressed
	virtual bool onMouseButtonEvent(int button, int action, int mods, double xPos, double yPos);
	//triggered when mouse moves
	virtual bool onMouseMoveEvent(double xPos, double yPos);
	//triggered when using the mouse wheel
	virtual bool onMouseWheelScrollEvent(double xOffset, double yOffset);

};


