#include <GUILib/GLUtils.h>
#include "PBFApp.h"
#include "Constants.h"

MassSpringApp::MassSpringApp() {
	setWindowTitle("Assignment 2: Position Based Fluids");
	TwAddSeparator(mainMenuBar, "sep2", "");

	TwAddVarRW(mainMenuBar, "Draw Particles", TW_TYPE_BOOLCPP, &ParticleSystem::drawParticles, "");
	TwAddVarRW(mainMenuBar, "Enable Gravity", TW_TYPE_BOOLCPP, &ParticleSystem::enableGravity, "");

	showGroundPlane = false;
	showDesignEnvironmentBox = true;
	showReflections = false;

	particleSystem = ParticleSystemLoader::loadFromOBJ("../meshes/bunny670.obj");

	pickedParticle = -1;
}

MassSpringApp::~MassSpringApp(void) {
	delete particleSystem;
}

const double PICK_DISTANCE = 0.1;

//triggered when mouse moves
bool MassSpringApp::onMouseMoveEvent(double xPos, double yPos) {
	if (pickedParticle > -1) {
		double modelViewMatrix[16];
		double projMatrix[16];
		GLint viewport[4];

		glGetDoublev(GL_MODELVIEW_MATRIX, modelViewMatrix);
		glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
		glGetIntegerv(GL_VIEWPORT, viewport);

		// Find the depth of the original point, after projection to screen
		double xScreen, yScreen, zScreen;

		gluProject(pickedPosition.at(0), pickedPosition.at(1), pickedPosition.at(2),
			modelViewMatrix, projMatrix, viewport,
			&xScreen, &yScreen, &zScreen);

		// Use the same depth to determine the altered position from screen coordinates

		double xWorld, yWorld, zWorld;

		gluUnProject(xPos, viewport[3] - yPos, zScreen,
			modelViewMatrix, projMatrix, viewport,
			&xWorld, &yWorld, &zWorld);

		pickedPosition = P3D(xWorld, yWorld, zWorld);
		particleSystem->setPosition(pickedParticle, pickedPosition);

		return false;
	}
	if (GLApplication::onMouseMoveEvent(xPos, yPos) == true) return true;

	return false;
}

bool MassSpringApp::pickParticle(double screenX, double screenY) {
	Ray mouseRay = camera->getRayFromScreenCoords(screenX, screenY);
	int particlePicked = -1;
	P3D closestPosition;
	double closestDistance = PICK_DISTANCE;

	for (int i = 0; i < particleSystem->particleCount(); i++) {
		P3D x = particleSystem->getPositionOf(i);
		double distance = mouseRay.getDistanceToPoint(x, NULL);
		if (distance < closestDistance) {
			closestDistance = distance;
			closestPosition = x;
			particlePicked = i;
		}
	}

	if (particlePicked > -1) {
		pickedParticle = particlePicked;
		pickedPosition = closestPosition;
		return true;
	}
	return false;
}

//triggered when mouse buttons are pressed
bool MassSpringApp::onMouseButtonEvent(int button, int action, int mods, double xPos, double yPos) {

	if (button == 0) {
		// Left mouse
		if (action == 1) {
			// Down
			pickParticle(xPos, yPos);
		}
		else {
			// Up
			pickedParticle = -1;
		}
	}

	if (GLApplication::onMouseButtonEvent(button, action, mods, xPos, yPos)) {
		return true;
	}

	return false;
}

//triggered when using the mouse wheel
bool MassSpringApp::onMouseWheelScrollEvent(double xOffset, double yOffset) {
	if (GLApplication::onMouseWheelScrollEvent(xOffset, yOffset)) return true;

	return false;
}

bool MassSpringApp::onKeyEvent(int key, int action, int mods) {
	if (GLApplication::onKeyEvent(key, action, mods)) return true;

	if (action == 1) {
		if (key == 'G') {
			if (ParticleSystem::enableGravity) {
				Logger::consolePrint("Disabled gravity");
				ParticleSystem::enableGravity = false;
			}
			else {
				Logger::consolePrint("Enabled gravity");
				ParticleSystem::enableGravity = true;
			}
		}
	}

	return false;
}

bool MassSpringApp::onCharacterPressedEvent(int key, int mods) {
	if (GLApplication::onCharacterPressedEvent(key, mods)) return true;

	return false;
}


void MassSpringApp::loadFile(const char* fName) {
	Logger::consolePrint("Loading file \'%s\'...\n", fName);
	std::string fileName;
	fileName.assign(fName);

	std::string fNameExt = fileName.substr(fileName.find_last_of('.') + 1);
}

void MassSpringApp::saveFile(const char* fName) {
	Logger::consolePrint("SAVE FILE: Do not know what to do with file \'%s\'\n", fName);
}


// Run the App tasks
void MassSpringApp::process() {
	// Take enough steps so that we are always running in (close to) real time
	int numSteps = (int)((1. / 30.) / DELTA_T);
	if (numSteps < 1) numSteps = 1;
	for (int i = 0; i < numSteps; i++) {
		particleSystem->integrate_PBF(DELTA_T);
		if (pickedParticle > -1) {
			particleSystem->setPosition(pickedParticle, pickedPosition);
		}
	}
}

// Draw the App scene - camera transformations, lighting, shadows, reflections, etc apply to everything drawn by this method
void MassSpringApp::drawScene() {
	particleSystem->drawParticleSystem();
}

// This is the wild west of drawing - things that want to ignore depth buffer, camera transformations, etc. Not pretty, quite hacky, but flexible. Individual apps should be careful with implementing this method. It always gets called right at the end of the draw function
void MassSpringApp::drawAuxiliarySceneInfo() {

}

// Restart the application.
void MassSpringApp::restart() {

}

bool MassSpringApp::processCommandLine(const std::string& cmdLine) {

	istringstream iss(cmdLine);

	string command, argument;

	if ((iss >> command >> argument)) {
		if (command == "load") {
			if (argument.length() < 5) {
				return true;
			}
			string extension = argument.substr(argument.length() - 4, argument.length() - 1);
			if (extension == ".obj") {
				ParticleSystem* ps = ParticleSystemLoader::loadFromOBJ(argument);
				if (ps) {
					delete particleSystem;
					particleSystem = ps;
					ParticleSystem::enableGravity = false;
				}
				else {
					Logger::consolePrint("Failed to load particle system from %s", argument.c_str());
				}
			}
			else {
				Logger::consolePrint("Extension not recognized: %s", extension.c_str());
			}
			return true;
		}
		else if (GLApplication::processCommandLine(cmdLine)) return true;
	}

	if (GLApplication::processCommandLine(cmdLine)) return true;

	return false;
}

