#include <stdio.h>
#include <stdlib.h>
#include <../GUILib/GLApplication.h>
#include "PBFApp.h"

int main(void){

	GLApplication* theApp;
	theApp = new MassSpringApp();

	theApp->runMainLoop();
	delete theApp;
	return 0;
}


