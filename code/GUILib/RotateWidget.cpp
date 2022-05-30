#include "GLUtils.h"
#include "RotateWidget.h"
#include <MathLib/MathLib.h>
#include <Utils/Logger.h>
#include "GlobalMouseState.h"

RotateWidget::RotateWidget() {
	pos = P3D(0, 0, 0);
}

RotateWidget::~RotateWidget() {
}
