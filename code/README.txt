Windows:

First, determine the full path to the libs file in this repository,
for example:

    C:\Users\christoy\cs15467_codebase\libs

Append this to your environment PATH variable. You can do this via
Control Panel > System > Advanced system settings > Environment Variables.
Find "Path" in your system variables, select it, and press edit. Go to
the very end of the value, add a semicolon, and paste in the path.
(Do not replace the entire value.)

Now download and install Visual Studio 2015 Community Edition from Dreamspark.
Once done, open SCP.sln. Visual Studio will prompt you to install some
dependencies, so install those, and open SCP.sln again.

Note: You may have to right click on the Assignment2 project in the Solution Explorer
and select "Set as startup project"

After this, you should be able to compile and run the code.

Linux:

This project has the following dependencies:

  - freetype
  - GLUT
  - GLEW
  - GLFW
  - AntTweakBar

The first three can be installed from most package managers; if using apt, the
recommended packages are libfreetype6-dev, glut3-dev, libglew-dev.

Version 3 of GLFW is required. At the moment, apt only has version 2, so
you will need to download glfw from the website (http://www.glfw.org/) and
build and install it from source. AntTweakBar should also be downloaded
and built and installed from source (http://anttweakbar.sourceforge.net/doc/).

