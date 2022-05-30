#pragma once

#include <Utils/Utils.h>

//FreeType Headers
#include <ft2build.h>
#include FT_FREETYPE_H

#include <vector>
#include <string>


//This holds all of the information related to any freetype font that we want to create.  
class FreeTypeFont {
public:
	FreeTypeFont();
	~FreeTypeFont();

	float height;			// the height of the font.
	uint * textures;	// texture id's 
	uint list_base;		// display list id

	//The init function will create a font of
	//of the height h from the file fname.
	void init(const char * fname);

	//Free all the resources assosiated with the font.
	void clean();

	double getLineHeight() {
		return	height / .63f;	/*about 1.5* that of text*/
	}

	//the x and y coordinates are relative to the current viewport
	void print(double x, double y, const char *text);
	void print(double x, double y, const std::string& line);
	void print(double x, double y, const std::vector<std::string>& lines);
};
