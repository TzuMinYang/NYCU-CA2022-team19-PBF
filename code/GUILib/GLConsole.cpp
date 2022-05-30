#include "GLUtils.h"

#include "GLConsole.h"
#include <Utils/Utils.h>

#include "GLContentManager.h"
#include "FreeType.h"
#include "GLApplication.h"

/**
	Default constructor
*/
GLConsole::GLConsole( int posX, int posY, int sizeX, int sizeY ) : GLWindow2D(posX, posY, sizeX, sizeY) {
}

GLConsole::GLConsole() : GLWindow2D(){
}

void GLConsole::draw() {
	preDraw();

	if (blinker.timeEllapsed() > 0.5) {
		blinker.restart();
		showCursor = !showCursor;
	}

	FreeTypeFont* font = GLContentManager::getFont("../data/fonts/arial.ttf 14");
	//we want to get the last line to end at the bottom of the console window, so figure out where the first one should start then...
	double startHeight = font->getLineHeight() * (Logger::consoleOutput.size() + 1.5);

	glColor3d(1-bgColorR, 1-bgColorG, 1-bgColorB);
	font->print(getViewportXFromRelativeX(0.0)+20, getViewportYFromRelativeY(0) + startHeight, Logger::consoleOutput);

	if (isActive() || isSelected())
		glColor3d(1 - bgColorR, 1 - bgColorG, bgColorB);

	if (showCursor)
		font->print(getViewportXFromRelativeX(0.0) + 20, getViewportYFromRelativeY(0) + font->getLineHeight(), ">> " + inputLine + "_");
	else
		font->print(getViewportXFromRelativeX(0.0) + 20, getViewportYFromRelativeY(0) + font->getLineHeight(), ">> " + inputLine);

	postDraw();
	glDisable(GL_TEXTURE_2D);
}

//all these methods should returns true if the event is processed, false otherwise...
//any time a physical key is pressed, this event will trigger. Useful for reading off special keys...
bool GLConsole::onKeyEvent(int key, int action, int mods) {
	if (key == GLFW_KEY_BACKSPACE && action == GLFW_PRESS) {
		if (inputLine.size() > 0)
			inputLine.erase(inputLine.end()-1);
		return true;
	}

	if (key == GLFW_KEY_ENTER && action == GLFW_PRESS) {
		if (inputLine.length() > 0)
			if (!GLApplication::getGLAppInstance()->processCommandLine(inputLine))
				Logger::consolePrint("%s: Command not recognized...\n", inputLine.c_str());

		inputLine.clear();

		return true;
	}

	return false;
}

//this one gets triggered on UNICODE characters only...
bool GLConsole::onCharacterPressedEvent(int key, int mods) {
	if (isActive() || isSelected()) {
		inputLine += key;
		return true;
	}
	return false;
}
