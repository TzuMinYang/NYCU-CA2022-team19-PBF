#pragma once

#include <MathLib/P3D.h>
#include <MathLib/V3D.h>
#include <MathLib/Quaternion.h>
#include <GUILib/InteractiveWidget.h>

/*======================================================*
 * Interface to open GL camera class.
 *======================================================*/
class GLCamera : public InteractiveWidget {
public:

public:
	GLCamera(void);
	virtual ~GLCamera(void);

	//this method is used to apply the transofmations 
	void applyCameraTransformations();

	virtual P3D getCameraPosition() = 0;
	virtual P3D getCameraTarget() = 0;
	virtual void setCameraTarget(const P3D& p) = 0;
	virtual V3D getWorldUpAxis() = 0;
	virtual Quaternion getCameraRotation() = 0;

	virtual Quaternion getRotationToOpenGLCoordinateSystem() {
		//openGL uses a y-up, right-handed coordinate system. If the camera is used to set up a different coordinate system, then this is the rotation that transforms things to openGL's way of doing things
		return Quaternion();
	}

	//draw function
	virtual void draw() {}

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


