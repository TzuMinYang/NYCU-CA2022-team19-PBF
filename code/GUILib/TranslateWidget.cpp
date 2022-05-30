#include "GLUtils.h"
#include "TranslateWidget.h"
#include <MathLib/MathLib.h>
#include <Utils/Logger.h>
#include "GlobalMouseState.h"

TranslateWidget::TranslateWidget(uint axes){
	if (axes & AXIS_X)
		tAxes.push_back(TranslateWidgetAxis(V3D(1, 0, 0), this));
	if (axes & AXIS_Y)
		tAxes.push_back(TranslateWidgetAxis(V3D(0, 1, 0), this));
	if (axes & AXIS_Z)
		tAxes.push_back(TranslateWidgetAxis(V3D(0, 0, 1), this));
}

TranslateWidget::~TranslateWidget(){
}

bool TranslateWidget::isPicked() {
	for (uint i = 0; i < tAxes.size(); i++) {
		if (tAxes[i].isPicked)
			return true;
	}

	return false;
}


void TranslateWidget::draw(){
	if (visible) {
		glPushAttrib(GL_ENABLE_BIT | GL_LIGHTING_BIT | GL_CURRENT_BIT);
		glEnable(GL_LIGHTING);
		glDisable(GL_TEXTURE_2D);
//		glDisable(GL_DEPTH_TEST);
		glPushMatrix();

		//to keep the size of the arrows roughly the same on screen, mess around with the scale...
		P3D camPos = getRayFromScreenCoords(0, 0).origin;
		double z = V3D(camPos, pos).length();
		scale = 50 * z / 1000;

		// Get and set the transformation matrix
		glTranslated(pos[0], pos[1], pos[2]);
		V3D rotAxis = orientation.v.unit();
		glRotated(DEG(orientation.getRotationAngle(rotAxis)), rotAxis[0], rotAxis[1], rotAxis[2]);

		int selCount = 0;
		for (auto it = tAxes.begin(); it != tAxes.end(); ++it) {
			if (it->isPicked)
				selCount += 1;
		}

		if (selCount <= 1) {
			glColor3d(0.8, 0.8, 0);
			drawSphere(P3D(0, 0, 0), scale * 0.2);
		}
		else {
			glColor3d(1, 1, 0);
			drawSphere(P3D(0, 0, 0), scale * 0.25);
		}

		for (auto it = tAxes.begin(); it != tAxes.end(); ++it)
			it->draw();

		glPopMatrix();
		glPopAttrib();

		glColor3d(1, 1, 1); 
	}
}

//triggered when mouse moves
bool TranslateWidget::onMouseMoveEvent(double xPos, double yPos) {

	if (visible == false) return false;
	Ray mouseRay = getRayFromScreenCoords(xPos, yPos);

	if (GlobalMouseState::lButtonPressed) {
		bool axisPicked = false;
		P3D newPickedPoint; mouseRay.getDistanceToPoint(pos, &newPickedPoint);
		V3D overallOffset(pickedPoint, newPickedPoint);
		V3D axisProjectedOffset;

		for (auto it = tAxes.begin(); it != tAxes.end(); ++it){
			if (it->isPicked) {
				axisPicked = true;
				axisProjectedOffset += orientation.rotate(it->tAxis) * overallOffset.getComponentAlong(orientation.rotate(it->tAxis));
			}
		}
		if (axisPicked) {
			pos += axisProjectedOffset;
			pickedPoint = newPickedPoint;
			return true;
		}
	}
	else {
		active = false;
		for (auto it = tAxes.begin(); it != tAxes.end(); ++it){
			it->pickWith(mouseRay);
			if (it->isPicked)
				active = true;
		}
		mouseRay.getDistanceToPoint(pos, &pickedPoint);
		if (active) return true;
	}
	return false;
}

void TranslateWidgetAxis::pickWith(const Ray& ray) {
	P3D closestPtToRay;
	// Get the distance between the ray and the axis
	double dist = ray.getDistanceToSegment(pWidget->pos, pWidget->pos + pWidget->orientation * tAxis * pWidget->scale, &closestPtToRay);

	if (dist < pWidget->scale / 10) {
		isPicked = true;
		pickingOffset = pWidget->orientation.inverseRotate(V3D(pWidget->pos, closestPtToRay));
	} else
		isPicked = false;
}

void TranslateWidgetAxis::draw() {
	double alpha = isPicked ? 3 : 1;
	glColor3d(tAxis[0] / 3 * alpha, tAxis[1] / 3 * alpha, tAxis[2] / 3 * alpha);
	drawArrow(P3D(), P3D() + tAxis*pWidget->scale * (1+(alpha-1)/20.0), pWidget->scale * (1 + (alpha - 1) / 20.0) / 10);
}
