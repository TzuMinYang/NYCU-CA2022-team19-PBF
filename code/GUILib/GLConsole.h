#pragma once

#include "GLWindow2D.h"
#include <Utils/Timer.h>
#include <string>

/**
  * Used for any window that needs to be displayed in the OpenGL window
  */
class GLConsole : public GLWindow2D {
protected:
	std::string inputLine;
	Timer blinker;
	bool showCursor=false;
public:

	GLConsole( int posX, int posY, int sizeX, int sizeY );
	GLConsole();

	virtual ~GLConsole() {}

	virtual void draw();

	//all these methods should returns true if the event is processed, false otherwise...
	//any time a physical key is pressed, this event will trigger. Useful for reading off special keys...
	virtual bool onKeyEvent(int key, int action, int mods);
	//this one gets triggered on UNICODE characters only...
	virtual bool onCharacterPressedEvent(int key, int mods);

};

