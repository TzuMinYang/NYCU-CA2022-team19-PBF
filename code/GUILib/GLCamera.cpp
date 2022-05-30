#include "GLCamera.h"
#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

GLCamera::GLCamera(void){
}

GLCamera::~GLCamera(void){
}

void GLCamera::applyCameraTransformations(){
	V3D up = getWorldUpAxis();
	P3D camPos = getCameraPosition();
	P3D target = getCameraTarget();
	gluLookAt(camPos[0], camPos[1], camPos[2], target[0], target[1], target[2], up[0], up[1], up[2]);
}


