#include "GLTrackingCamera.h"
#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>
#include <Utils/Logger.h>
#include "GlobalMouseState.h"

//TODO: fixed camera target (with translation allowed only about vertical axis)

GLTrackingCamera::GLTrackingCamera(double distToTarget, const V3D& camViewDirection, const V3D& camUpAxis){
	assert(IS_ZERO(camViewDirection.dot(camUpAxis)));
	camDistance = distToTarget;
	camTarget = P3D(0,0,0);
	this->camViewDirection = camViewDirection;
	this->camUpAxis = camUpAxis;
}

GLTrackingCamera::~GLTrackingCamera(void){
}

P3D GLTrackingCamera::getCameraPosition(){
	Quaternion orientation = getCameraRotation();
	//we know we are looking down the z-axis of the camera coordinate frame. So we need to compute what that axis is in world coordinates.
	V3D viewDirWorld = orientation.rotate(camViewDirection);
	//this is the camera's up vector, expressed in world coordinates 
	//we know where we are looking, so we can compute the position of the camera in world coordinates.
	return camTarget + viewDirWorld * camDistance;
}

Quaternion GLTrackingCamera::getCameraRotation(){
//build the rotation of the camera...
	return getRotationQuaternion(rotAboutUpAxis, camUpAxis) * getRotationQuaternion(rotAboutRightAxis, camUpAxis.cross(camViewDirection));
}

V3D GLTrackingCamera::getWorldUpAxis(){
	//this is the camera's up vector, expressed in world coordinates 
	return getCameraRotation().rotate(camUpAxis);
}

P3D GLTrackingCamera::getCameraTarget(){
	return camTarget;
}

void GLTrackingCamera::setCameraTarget(const P3D& p){
	camTarget = p;
}

void GLTrackingCamera::setRotations(const V3D& r){
//	rotations = r;
}

//triggered when mouse buttons are pressed
bool GLTrackingCamera::onMouseButtonEvent(int button, int action, int mods, double xPos, double yPos) {
	return false;
}

//triggered when mouse moves
bool GLTrackingCamera::onMouseMoveEvent(double xPos, double yPos) {
	V3D screenCameraOffset(GlobalMouseState::mouseMoveX, GlobalMouseState::mouseMoveY);
	//translate (x,y) offset of the mouse to a world coordinates offset for the camera camTarget
	V3D worldOffset = getCameraRotation().rotate(getRotationToOpenGLCoordinateSystem().inverseRotate(screenCameraOffset));

	if (GlobalMouseState::rButtonPressed) {
		//TODO: the multiplier should probably correspond to a function that maps screen displacements to displacements in the viewport...
		camTarget -= worldOffset * (0.0005 * camDistance);
		return true;
	}

	if (GlobalMouseState::lButtonPressed) {
		rotAboutUpAxis += screenCameraOffset[0] * 0.01;
		rotAboutRightAxis += screenCameraOffset[1] * 0.003;
		return true;
	}

	if (GlobalMouseState::mButtonPressed) {
		camTarget += getCameraRotation().rotate(camViewDirection) * screenCameraOffset[1] * 0.01;
		return true;
	}

	return false;
}

//triggered when using the mouse wheel
bool GLTrackingCamera::onMouseWheelScrollEvent(double xOffset, double yOffset) {
	camDistance *= 1.0 - yOffset * 0.05;
	if (camDistance > -0.1) camDistance = -0.1;
	return true;
}

