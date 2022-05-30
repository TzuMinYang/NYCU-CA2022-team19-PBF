#pragma once
#ifdef _WIN32
#include <include/glew.h>
#else
#include <GL/glew.h>
#endif


#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>
#include <AntTweakBar.h>
#include "GLWindow.h"
#include <Utils/Logger.h>
#include <Utils/Timer.h>
#include "GLWindow2D.h"
#include "GLConsole.h"
#include "GLCamera.h"
#include "GLContentManager.h"
#include "GlobalMouseState.h"

#pragma warning( disable : 4005)

/**
 * All Apps will instantiate this base class...
 */
class GLApplication : public GLWindow {
public:
	static GLApplication *getGLAppInstance();
	static void setGLAppInstance(GLApplication* instance);

public:
	// Pointer to menu bar
	TwBar *mainMenuBar;
	// Pointer to glfw window
	GLFWwindow* glfwWindow;

protected:
	//keep a timer here to see how long it's been since the last redraw
	Timer fpsTimer;
	// Console window
	GLConsole* consoleWindow;
	// interactive camera
	GLCamera* camera;

	float cameraRot[4] = { 0, 0, 0, 1 };

	bool appIsRunning = false;
	bool waitForFrameRate = true;
	//this is the desired frame rate, specified in FPS
	double desiredFrameRate = 30;
	double animationSpeedupFactor = 1.0;

	// The size of the window
	int mainWindowWidth, mainWindowHeight;

	// Background color
	float bgColor[3] = { 0.5f, 0.5f, 0.5f };

	// Sets up various settings for OpenGL
	void setupOpenGL();

	// Set up lights
	virtual void setupLights();

	// Main draw function - draws app content, shadows, reflections, etc...
	virtual void draw();
	
	// Draw information regarding the frame rate and performance
	void drawFPS(double timeSinceLastUpdate, double percentageOfTimeSpentProcessing);
	void init(int x, int y, int w, int h);

protected:
	bool saveScreenshots = false;
	bool showMenus = true;
	bool showFPS = true;
	bool showConsole = true;
	bool showGroundPlane = true;
	bool showDesignEnvironmentBox = false;
	
	bool showReflections = true, showShadows = false;


public:
	// constructors
	GLApplication(int x, int y, int w, int h);
	GLApplication();

	// destructor
	virtual ~GLApplication(void);

	// Run the App tasks
	virtual void process();
	// Draw the App scene - camera transformations, lighting, shadows, reflections, etc apply to everything drawn by this method
	virtual void drawScene();

	// draws reflections on the ground plane.
	void drawReflections();

	// This is the wild west of drawing - things that want to ignore depth buffer, camera transformations, etc. Not pretty, quite hacky, but flexible. Individual apps should be careful with implementing this method. It always gets called right at the end of the draw function
	virtual void drawAuxiliarySceneInfo();
	// Restart the application.
	virtual void restart();

	// Run the main loop
	virtual void runMainLoop();

	// adjusts the window size
	void setWindowTitle(char* windowTitle);

	//adjusts the size of the window
	virtual void adjustWindowSize(int w, int h);

	int getMainWindowWidth() {
		return mainWindowWidth;
	}

	int getMainWindowHeight() {
		return mainWindowHeight;
	}


	//input callbacks...

	//all these methods should returns true if the event is processed, false otherwise...
	//any time a physical key is pressed, this event will trigger. Useful for reading off special keys...
	virtual bool onKeyEvent(int key, int action, int mods);
	//this one gets triggered on UNICODE characters only...
	virtual bool onCharacterPressedEvent(int key, int mods);
	//triggered when mouse buttons are pressed
	virtual bool onMouseButtonEvent(int button, int action, int mods, double xPos, double yPos);
	//triggered when mouse moves
	virtual bool onMouseMoveEvent(double xPos, double yPos);
	//triggered when using the mouse wheel
	virtual bool onMouseWheelScrollEvent(double xOffset, double yOffset);

	virtual bool processCommandLine(const std::string& cmdLine);
	
	virtual void saveFile(const char* fName);
	virtual void loadFile(const char* fName);
	virtual void loadFiles(int n, const char** fNames);
};


// Callback function called by GLFW when window size changes
void windowResizeEventCallback(GLFWwindow* window, int width, int height);
// Callback function called by GLFW when key is pressed
void keyEventCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
void charEventCallback(GLFWwindow* window, unsigned int key, int mods);
void mouseButtonEventCallback(GLFWwindow* window, int button, int action, int mods);
void mouseMoveEventCallback(GLFWwindow* window, double xPos, double yPos);
void mouseScrollEventCallback(GLFWwindow* window, double xOffset, double yOffset);
void fileDropEventCallback(GLFWwindow* window, int n, const char** fNames);

void TW_CALL menuBarButtonEvent(void* clientData);



