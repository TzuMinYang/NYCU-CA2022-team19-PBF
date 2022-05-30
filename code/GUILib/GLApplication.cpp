#ifdef _WIN32
#include <include/glew.h>
#else
#include <GL/glew.h>
#endif

#include "GLApplication.h"
#include "GLUtils.h"
#include <Utils/Utils.h>
#include "GLTrackingCamera.h"
#include "GLContentManager.h"
#include "GLTexture.h"
#include "GLShaderMaterial.h"

GLApplication* glAppInstance = NULL;

GLApplication::GLApplication(int x, int y, int w, int h){
	GLApplication::setGLAppInstance(this);

	if (!glfwInit()) {
		// An error occured
		Logger::print("GLFW initialization failed\n");
		exit(0);
	}

	init(x, y, w, h);

	const GLFWvidmode * mode = glfwGetVideoMode(glfwGetPrimaryMonitor());

	x = 0;
	y = 0;
	w = mode->width;
	h = mode->height;
}

GLApplication::GLApplication() {
	GLApplication::setGLAppInstance(this);

	if (!glfwInit()) {
		// An error occured
		Logger::print("GLFW initialization failed\n");
		exit(0);
	}

	const GLFWvidmode * mode = glfwGetVideoMode(glfwGetPrimaryMonitor());

	int borderLeft = 10;
	int borderTop = 45;
	int borderRight = 10;
	int borderBottom = 70;

	init(borderLeft, borderTop, mode->width-borderLeft - borderRight, mode->height - borderTop - borderBottom);
}

void GLApplication::init(int x, int y, int w, int h) {
	mainWindowWidth = w;
	mainWindowHeight = h;

	glfwWindowHint(GLFW_SAMPLES, 4);
	/* Create a windowed mode window and its OpenGL context */
	glfwWindow = glfwCreateWindow(w, h, "", NULL, NULL);
	if (!glfwWindow) {
		Logger::print("Could not initialize GLFW window\n");
		glfwTerminate();
		exit(0);
	}
	setWindowTitle("Simulation And Control Playground");
	glfwSetWindowPos(glfwWindow, x, y);

	/* Make the window's context current */
	glfwMakeContextCurrent(glfwWindow);

	// Initialize AntTweakBar
//	TwDefine(" GLOBAL fontscaling=1.2 ");
	TwInit(TW_OPENGL, NULL);

	// Create a tweak bar
	mainMenuBar = TwNewBar("MainMenuBar");

	int mainMenuWidth = 260;

//	TwDefine(" MainMenuBar help='Play, Pause, etc.' position='0 0' size='260 800' valueswidth=100"); // set up the bar...
	TwDefine((" MainMenuBar help='Play, Pause, etc.' position='0 0' size='" + std::to_string(mainMenuWidth) + " " + std::to_string(h) + "' valueswidth=100").c_str());
	TwDefine(" GLOBAL fontsize=3 ");
	
//	TwAddVarRW(mainMenuBar, "bgColor", TW_TYPE_COLOR3F, &bgColor, " label='Background color' ");
	TwAddVarRW(mainMenuBar, "WaitForRefresh", TW_TYPE_BOOLCPP, &waitForFrameRate, " label='Wait For Framerate' group='Viz'");
	TwAddVarRW(mainMenuBar, "Menu Bars", TW_TYPE_BOOLCPP, &showMenus, " label='Show Menu Bars' key=m group='Viz'");
	TwAddVarRW(mainMenuBar, "Console", TW_TYPE_BOOLCPP, &showConsole, " label='Show Console' key=c group='Viz'");
	TwAddVarRW(mainMenuBar, "FPS", TW_TYPE_BOOLCPP, &showFPS, " label='Show FPS' group='Viz'");
	//TwAddVarRW(mainMenuBar, "Ground", TW_TYPE_BOOLCPP, &showGroundPlane, " label='Show Ground Plane' group='Viz'");
	//TwAddVarRW(mainMenuBar, "Reflections", TW_TYPE_BOOLCPP, &showReflections, " label='Show Reflections' group='Viz'");

	TwAddVarRW(mainMenuBar, "DesignEnvironmentBox", TW_TYPE_BOOLCPP, &showDesignEnvironmentBox, " label='Show Design Environment' group='Viz'");
	TwAddVarRW(mainMenuBar, "Record Screen Shots", TW_TYPE_BOOLCPP, &saveScreenshots, " label='Record Screenshots'");

	TwAddSeparator(mainMenuBar, "sep1", "");
	TwAddButton(mainMenuBar, "Restart", menuBarButtonEvent, (void*)"Restart", "");
	TwAddVarRW(mainMenuBar, "Play", TW_TYPE_BOOLCPP, &appIsRunning, " label='Play/Pause'");
	TwAddButton(mainMenuBar, "Step", menuBarButtonEvent, (void*)"Step", "");

	//setup callback functions...
	glfwSetWindowSizeCallback(glfwWindow, windowResizeEventCallback);	// Window resize
	glfwSetKeyCallback(glfwWindow, keyEventCallback);					// key event - for any physical key being pressed, repeated, or released
	glfwSetCharModsCallback(glfwWindow, charEventCallback);				// another type of key event - on text (UNICODE) key press or repeated
	glfwSetMouseButtonCallback(glfwWindow, mouseButtonEventCallback);	// mouse button click
	glfwSetCursorPosCallback(glfwWindow, mouseMoveEventCallback);		// mouse move event
	glfwSetScrollCallback(glfwWindow, mouseScrollEventCallback);		// mouse scroll
	glfwSetDropCallback(glfwWindow, fileDropEventCallback);				// file drop events

	adjustWindowSize(w, h);
	setupOpenGL();

	consoleWindow = new GLConsole((int)(mainMenuWidth*1.14), 0, (int)(mainWindowWidth - mainMenuWidth*1.14), 280);

	camera = new GLTrackingCamera();

	glewInit();
	GLContentManager::addShaderProgram("radialGradientShader", "../data/shaders/radialGradient/radialGradient.vert", "../data/shaders/radialGradient/radialGradient.frag");

	GLShaderMaterial* material = new GLShaderMaterial();
	material->readFromFile("../data/shaders/radialGradient/radialGradient.mat");
	GLContentManager::addShaderMaterial(material->getMaterialName().c_str(), material);
}

//triggered when mouse moves
bool GLApplication::onMouseMoveEvent(double xPos, double yPos) {
	if (showMenus)
		if (TwEventMousePosGLFW((int)xPos, (int)yPos)) return true;
	if (showConsole)
		consoleWindow->onMouseMoveEvent(xPos, yPos);
	if (camera)
		camera->onMouseMoveEvent(xPos, yPos);
	return false;
}

//triggered when mouse buttons are pressed
bool GLApplication::onMouseButtonEvent(int button, int action, int mods, double xPos, double yPos) {
	if (showMenus)
		if (TwEventMouseButtonGLFW(button, action))
		return true;
	if (showConsole)
		consoleWindow->onMouseButtonEvent(button, action, mods, xPos, yPos);
	if (camera)
		if (camera->onMouseButtonEvent(button, action, mods, xPos, yPos))
		return true;

	return false;
}

//triggered when using the mouse wheel
bool GLApplication::onMouseWheelScrollEvent(double xOffset, double yOffset) {
	if (showMenus)
		if (TwEventMouseWheelGLFW((int)yOffset)) 
			return true;
	if (camera)
		if (camera->onMouseWheelScrollEvent(xOffset, yOffset))
			return true;
	return false;
}

bool GLApplication::onKeyEvent(int key, int action, int mods) {
	//listen to key presses even if the menu is not visible...
	if (key == GLFW_KEY_ESCAPE) {
		glfwSetWindowShouldClose(glfwWindow, GL_TRUE);
		return true;
	}
	//if the console is active and listening, it means we're typing in it, so let it handle it first
	if (consoleWindow->onKeyEvent(key, action, mods))
		return true;

	if (TwEventKeyGLFW(key, action))
		return true;

	if (camera)
		if (camera->onKeyEvent(key, action, mods))
			return true;

	return false;
}

bool GLApplication::onCharacterPressedEvent(int key, int mods) {
	//if the console is active and listening, it means we're typing in it, so let it handle it first
	if (consoleWindow->onCharacterPressedEvent(key, mods))
		return true;

	//we listen to the keys even if the menu is not visible...
	if (TwEventCharGLFW(key, GLFW_PRESS)) 
		return true;

	if (camera)
		if (camera->onCharacterPressedEvent(key, mods))
			return true;

	if (key == ' ') {
		appIsRunning = !appIsRunning;
	}

	return false;
}

void GLApplication::loadFiles(int nFiles, const char** fNames) {
	for (int i = 0;i < nFiles;i++)
		GLApplication::getGLAppInstance()->loadFile(fNames[i]);
}

void GLApplication::loadFile(const char* fName) {
	Logger::consolePrint("LOAD FILE: Do not know what to do with file \'%s\'\n", fName);
}

void GLApplication::saveFile(const char* fName) {
	Logger::consolePrint("SAVE FILE: Do not know what to do with file \'%s\'\n", fName);
}

void GLApplication::setWindowTitle(char* windowTitle) {
	glfwSetWindowTitle(glfwWindow, windowTitle);
}

GLApplication::~GLApplication(void){
	// Terminate AntTweakBar and GLFW
	TwTerminate();
	glfwTerminate();
}

void GLApplication::runMainLoop() {
	// Main loop (repeated while window is not closed and [ESC] is not pressed)
	while (!glfwWindowShouldClose(glfwWindow)){
		double timeSpentProcessing = 0;
		fpsTimer.restart();
		if (appIsRunning)
			process();
		timeSpentProcessing = fpsTimer.timeEllapsed();

		draw();

		//wait until the required ammount of time has passed (respect the desired FPS requirement)
		if (waitForFrameRate)
			while (fpsTimer.timeEllapsed()< 1.0 / desiredFrameRate);

		if (showFPS) 
			drawFPS(fpsTimer.timeEllapsed(), timeSpentProcessing / fpsTimer.timeEllapsed());

		//print a screenshot if needed
		if (saveScreenshots) {
			static int screenShotNumber = 0;
			char fName[100];
			sprintf(fName, "..\\screenShots\\ss%05d.bmp", screenShotNumber);
			screenShotNumber++;
			int viewportSettings[4];//x, y, w, h
			glGetIntegerv(GL_VIEWPORT, viewportSettings);
			saveScreenShot(fName, 0, 0, mainWindowWidth, mainWindowHeight);
		}

		/* Swap front and back buffers */
		glfwSwapBuffers(glfwWindow);

		/* Poll for and process events */
		glfwPollEvents();
	}

}

/**
* draw information regarding the frame rate and performance
*/
void GLApplication::drawFPS(double timeSinceLastUpdate, double percentageOfTimeSpentProcessing) {
	//in order to avoid writing the a new fps rate at every redraw, we will instead compute averages and display those. We need: the total ellapsedTime since the last update
	static double ellapsedTime = 0;
	//the frame rate that is currently being displayed
	static double oldFrameRate = 0;
	//the number of frames that we've counted so far
	static int nrFramesSinceUpdate = 0;
	//the total time spent processing
	static double processingTime = 0;
	//and the processing time that we are displaying
	static double oldPerformanceRate = 0;

	ellapsedTime += timeSinceLastUpdate;
	processingTime += percentageOfTimeSpentProcessing;
	nrFramesSinceUpdate++;

	//only change the numbers that we display about 3 times per second
	if (ellapsedTime >= 1 / 3.0) {
		oldFrameRate = nrFramesSinceUpdate / ellapsedTime;
		oldPerformanceRate = processingTime / nrFramesSinceUpdate;
		nrFramesSinceUpdate = 0;
		ellapsedTime = 0;
		processingTime = 0;
	}

	glPushMatrix();
	glLoadIdentity();
	glTranslatef(0.0f, 0.0f, -1.0f);

	glColor3d(1.0, 1.0, 1.0);
	glprint(viewportWidth - 400, viewportHeight-15, "FPS: %7.2lf (processing: %7.2lf %%)\n", oldFrameRate, 100 * oldPerformanceRate);

	glPopMatrix();
}


//this method is used to set up the lights (position, direction, etc), relative to the camera position
void GLApplication::setupLights(){
	GLfloat bright[] = { 0.8f, 0.8f, 0.8f, 1.0f }; 
	GLfloat mediumbright[] = { 0.3f, 0.3f, 0.3f, 1.0f };

	glLightfv(GL_LIGHT1, GL_DIFFUSE, bright);
	glLightfv(GL_LIGHT2, GL_DIFFUSE, mediumbright);

	GLfloat light0_position[] = { 0.0f, 10000.0f, 10000.0f, 0.0f };
	GLfloat light0_direction[] = { 0.0f, -10000.0f, -10000.0f, 0.0f };

	GLfloat light1_position[] = { 0.0f, 10000.0f, -10000.0f, 0.0f };
	GLfloat light1_direction[] = { 0.0f, -10000.0f, 10000.0f, 0.0f };

	GLfloat light2_position[] = { 0.0f, -10000.0f, 0.0f, 0.0f };
	GLfloat light2_direction[] = { 0.0f, 10000.0f, -0.0f, 0.0f };

	glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
	glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
	glLightfv(GL_LIGHT2, GL_POSITION, light2_position);

	glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, light0_direction);
	glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, light1_direction);
	glLightfv(GL_LIGHT2, GL_SPOT_DIRECTION, light2_direction);

	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHT2);

	glDisable(GL_LIGHT3);
	glDisable(GL_LIGHT4);

}

// Run the App tasks
void GLApplication::process() {
	//do the work here...
}

// Draw the App scene - camera transformations, lighting, shadows, reflections, etc apply to everything drawn by this method
void GLApplication::drawScene() {
	int pass, numPass;
	bool wireframe = false;
	unsigned char cubeColor[] = { 255, 0, 0, 128 }; // Model color (32bits RGBA)

	glMatrixMode(GL_MODELVIEW);

//	static double turn = 0;
//	turn += 0.001;
//	glRotated(360.0*turn, 0.4, 1, 0.2);

	glBegin(GL_LINES);
	glColor3d(1, 0, 0);
	glVertex3d(0, 0, 0);
	glVertex3d(1, 0, 0);

	glColor3d(0, 1, 0);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 1, 0);

	glColor3d(0, 0, 1);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 0, 1);
	glEnd();

	glTranslated(-0.5, -0.5, -0.5);

	// Set color and draw model
	glColor4ubv(cubeColor);

	glLineWidth(3.0);

	if (wireframe){
		glDisable(GL_CULL_FACE);
		glDisable(GL_LIGHTING);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		numPass = 1;
	}
	else{
		glEnable(GL_CULL_FACE);
		glEnable(GL_LIGHTING);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		numPass = 2;
	}

	for (pass = 0; pass<numPass; ++pass){
		// Since the material could be transparent, we draw the convex model in 2 passes:
		// first its back faces, and second its front faces.
		glCullFace((pass == 0) ? GL_FRONT : GL_BACK);

		// Draw the model (a cube)
		glBegin(GL_QUADS);
		glNormal3f(0, 0, -1); glVertex3f(0, 0, 0); glVertex3f(0, 1, 0); glVertex3f(1, 1, 0); glVertex3f(1, 0, 0); // front face
		glNormal3f(0, 0, +1); glVertex3f(0, 0, 1); glVertex3f(1, 0, 1); glVertex3f(1, 1, 1); glVertex3f(0, 1, 1); // back face
		glNormal3f(-1, 0, 0); glVertex3f(0, 0, 0); glVertex3f(0, 0, 1); glVertex3f(0, 1, 1); glVertex3f(0, 1, 0); // left face
		glNormal3f(+1, 0, 0); glVertex3f(1, 0, 0); glVertex3f(1, 1, 0); glVertex3f(1, 1, 1); glVertex3f(1, 0, 1); // right face
		glNormal3f(0, -1, 0); glVertex3f(0, 0, 0); glVertex3f(1, 0, 0); glVertex3f(1, 0, 1); glVertex3f(0, 0, 1); // bottom face
		glNormal3f(0, +1, 0); glVertex3f(0, 1, 0); glVertex3f(0, 1, 1); glVertex3f(1, 1, 1); glVertex3f(1, 1, 0); // top face
		glEnd();
	}

	glLineWidth(1.0);
	glDisable(GL_LIGHTING);
}

// This is the wild west of drawing - things that want to ignore depth buffer, camera transformations, etc. Not pretty, quite hacky, but flexible. Individual apps should be careful with implementing this method. It always gets called right at the end of the draw function
void GLApplication::drawAuxiliarySceneInfo() {

}

// Restart the application.
void GLApplication::restart() {

}

/**
 *  Initializes the openGL settings
 */
void GLApplication::setupOpenGL(){
//	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
	glEnable(GL_LINE_SMOOTH);

	glShadeModel(GL_SMOOTH);
	glClearDepth(1.0f);									// Depth Buffer Setup
	glEnable(GL_DEPTH_TEST);							// Enables Depth Testing
	glDepthFunc(GL_LEQUAL);								// The Type Of Depth Testing To Do
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);	// Really Nice Perspective Calculations
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA ,GL_ONE_MINUS_SRC_ALPHA);
}

/**
 *	This method is used to draw the scene.
 */
void GLApplication::draw(){
	//clear the screen
	glClearColor(bgColor[0], bgColor[1], bgColor[2], 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	//setup the viewport and the perspective transformation matrix
	glViewport(viewportX, viewportY, viewportWidth, viewportHeight);
	glMatrixMode(GL_PROJECTION);						// Select The Projection Matrix
	glLoadIdentity();									// Reset The Projection Matrix
	gluPerspective(45.0, (double)viewportWidth / viewportHeight, 0.05, 150.0);
	glMatrixMode(GL_MODELVIEW);							// Select The Modelview Matrix
	glLoadIdentity();									// Reset The Modelview Matrix

	if (camera){
//		if (followCameraTarget){
//			V3D upAxis = camera->getCameraOrientation().inverseRotate(camera->getWorldUpAxis());
//			double verticalValue = camera->getCameraTarget().getComponentAlong(upAxis);
//			P3D cameraTarget = getCameraTarget();
//			cameraTarget.setComponentAlong(upAxis, verticalValue);
//			setCameraTarget(cameraTarget);
//		}
		camera->applyCameraTransformations();
	}else
		gluLookAt(-1, 0, 3, 0, 0, 0, 0, 1, 0);

	setupLights();

	if (showDesignEnvironmentBox) {
		glColor3d(1, 1, 1);
		glDisable(GL_LIGHTING);
		drawDesignEnvironmentBox();
	}

	if (showGroundPlane && !showReflections){
		glColor3d(1,1,1);
		glDisable(GL_LIGHTING);
		drawGround();
	}

	if (showReflections && showGroundPlane) drawReflections();
//	if (showShadows && showGroundPlane)	drawShadows();

	drawScene();

	if (showMenus) {
		//update the orientation of the camera orientation visualized in the menu...
		if (camera) {
			//transformation goes from world to camera to openGL coordinate system
			Quaternion glRelativeCamRot = camera->getRotationToOpenGLCoordinateSystem() * camera->getCameraRotation().getComplexConjugate();

			cameraRot[0] = (float)glRelativeCamRot.v[0];
			cameraRot[1] = (float)glRelativeCamRot.v[1];
			cameraRot[2] = (float)glRelativeCamRot.v[2];
			cameraRot[3] = (float)glRelativeCamRot.s;
		}

		TwDraw();
		//if the menu is not set up properly, it will have problems drawing and mess other things up... so get the opengl error here...

	}
	if (showConsole) consoleWindow->draw();

	drawAuxiliarySceneInfo();
}


/**
* draws reflections on the ground plane.
*/
void GLApplication::drawReflections() {
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glClear(GL_STENCIL_BUFFER_BIT);

	glEnable(GL_STENCIL_TEST); //Enable using the stencil buffer
	glColorMask(0, 0, 0, 0); //Disable drawing colors to the screen
	glDisable(GL_DEPTH_TEST); //Disable depth testing
	glStencilFunc(GL_ALWAYS, 1, 1); //Make the stencil test always pass
									//Make pixels in the stencil buffer be set to 1 when the stencil test passes
	glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
	drawGround();

	glColorMask(1, 1, 1, 1); //Enable drawing colors to the screen
	glEnable(GL_DEPTH_TEST); //Enable depth testing
							 //Make the stencil test pass only when the pixel is 1 in the stencil buffer
	glStencilFunc(GL_EQUAL, 1, 1);
	glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP); //Make the stencil buffer not change
	glPushMatrix();

	V3D scale(1, 1, 1);
	scale -= Globals::groundPlane.n * 2;

	glScalef((float)scale[0], (float)scale[1], (float)scale[2]);

	glEnable(GL_CLIP_PLANE0);

	double plane[4];
	Globals::groundPlane.getCartesianEquationCoefficients(plane[0], plane[1], plane[2], plane[3]);
	glClipPlane(GL_CLIP_PLANE0, plane);
	drawScene();
	glPopMatrix();

	glDisable(GL_CLIP_PLANE0);

	glDisable(GL_STENCIL_TEST);

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);
	glDisable(GL_LIGHTING);
	glColor4f(1, 1, 1, 0.85f);
	drawGround();
}




//adjusts the window size
void GLApplication::adjustWindowSize(int w, int h) {
	// Send the window size to AntTweakBar
	TwWindowSize(w, h);

	setViewportParameters(0, 0, w, h);
}

GLApplication* GLApplication::getGLAppInstance() {
	return glAppInstance;
}

void GLApplication::setGLAppInstance(GLApplication* instance) {
	if (glAppInstance == NULL)
		glAppInstance = instance;
}

bool GLApplication::processCommandLine(const std::string& cmdLine) {
	std::vector<std::string> lines;
	getCharSeparatedStringList(cmdLine.c_str(), lines, ' ');

	if (strcmp(lines[0].c_str(), "clear") == 0) {
		Logger::consoleOutput.clear();
		return true;
	}

	if (strcmp(lines[0].c_str(), "save") == 0 && cmdLine.length() > strlen("save")) {
		saveFile(cmdLine.c_str() + strlen("save") + 1);
		return true;
	}

	if (strcmp(lines[0].c_str(), "load") == 0 && cmdLine.length() > strlen("load")) {
		loadFile(cmdLine.c_str() + strlen("load") + 1);
		return true;
	}

	return false;
}

// Callback function called by GLFW when window size changes
void windowResizeEventCallback(GLFWwindow* window, int width, int height) {
	GLApplication::getGLAppInstance()->adjustWindowSize(width, height);
}

// Callback function called by GLFW when key is pressed
void keyEventCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	//	Logger::print("key press: %d (%c), scancode: %d (%c), key event: %d, mods: %d\n", key, key, scancode, scancode, action, mods);
	GLApplication::getGLAppInstance()->onKeyEvent(key, action, mods);
}

void charEventCallback(GLFWwindow* window, unsigned int key, int mods) {
	//	Logger::print("char event: %d (%c), mods: %d\n", key, key, mods);
	GLApplication::getGLAppInstance()->onCharacterPressedEvent(key, mods);
}

void mouseButtonEventCallback(GLFWwindow* window, int button, int action, int mods) {
	double xPos, yPos;
	glfwGetCursorPos(window, &xPos, &yPos);

	GlobalMouseState::updateMouseState(xPos, yPos, button, action, mods);

	//Logger::print("mouse button event: button %d, action: %d, mods: %d, xPos: %lf, yPos: %lf\n", button, action, mods, xPos, yPos);
	GLApplication::getGLAppInstance()->onMouseButtonEvent(button, action, mods, xPos, yPos);
}

void mouseMoveEventCallback(GLFWwindow* window, double xPos, double yPos) {
	GlobalMouseState::updateMouseMove(xPos, yPos);

	//	Logger::print("mouse move event: xPos: %lf, yPos: %lf\n", xPos, yPos);
	GLApplication::getGLAppInstance()->onMouseMoveEvent(xPos, yPos);
}

void mouseScrollEventCallback(GLFWwindow* window, double xOffset, double yOffset) {
	//	Logger::print("mouse scroll event: xOffset: %lf, yOffset: %lf\n", xOffset, yOffset);
	GLApplication::getGLAppInstance()->onMouseWheelScrollEvent(xOffset, yOffset);
}

void fileDropEventCallback(GLFWwindow* window, int n, const char** fNames) {
	GLApplication::getGLAppInstance()->loadFiles(n, fNames);
}

void TW_CALL menuBarButtonEvent(void* clientData) {
	if (strcmp((char*)clientData, "Step") == 0)
		GLApplication::getGLAppInstance()->process();
	if (strcmp((char*)clientData, "Restart") == 0)
		GLApplication::getGLAppInstance()->restart();
}
