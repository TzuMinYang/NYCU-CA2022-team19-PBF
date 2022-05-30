//Include our header file.
#include "FreeType.h"
#include <ftglyph.h>
#include <ftoutln.h>
#include <fttrigon.h>
#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

#include <Utils/Utils.h>

///This function gets the first power of 2 >= the
///int that we pass it.
inline int next_p2(int a) {
	int rval = 1;
	while (rval < a) rval <<= 1;
	return rval;
}

///Create a display list coresponding to the give character.
void make_dlist(FT_Face face, char ch, uint list_base, uint * tex_base) {
	//The first thing we do is get FreeType to render our character into a bitmap.

	//Load the Glyph for our character.
	if (FT_Load_Glyph(face, FT_Get_Char_Index(face, ch), FT_LOAD_DEFAULT))
		throwError("FT_Load_Glyph failed");

	//Move the face's glyph into a Glyph object.
	FT_Glyph glyph;
	if (FT_Get_Glyph(face->glyph, &glyph))
		throwError("FT_Get_Glyph failed");

	//Convert the glyph to a bitmap.
	FT_Glyph_To_Bitmap(&glyph, ft_render_mode_normal, 0, 1);
	FT_BitmapGlyph bitmap_glyph = (FT_BitmapGlyph)glyph;

	//This reference will make accessing the bitmap easier
	FT_Bitmap& bitmap = bitmap_glyph->bitmap;

	//Use our helper function to get the widths of
	//the bitmap data that we will need in order to create
	//our texture.
	int width = next_p2(bitmap.width);
	int height = next_p2(bitmap.rows);

	//Allocate memory for the texture data.
	unsigned char* expanded_data = new unsigned char[2 * width * height];

	//Here we fill in the data for the expanded bitmap.
	//Notice that we are using two channel bitmap (one for
	//luminocity and one for alpha), but we assign
	//both luminocity and alpha to the value that we
	//find in the FreeType bitmap. 
	//We use the ?: operator so that value which we use
	//will be 0 if we are in the padding zone, and whatever
	//is the the Freetype bitmap otherwise.
	for (uint j = 0; j < (uint)height;j++) {
		for (uint i = 0; i < (uint)width; i++) {
			expanded_data[2 * (i + j*width)] = 255;
			expanded_data[2 * (i + j*width) + 1] =
				(i >= bitmap.width || j >= bitmap.rows) ?
				0 : bitmap.buffer[i + bitmap.width*j];
		}
	}

	//Now we just setup some texture paramaters.
	glBindTexture(GL_TEXTURE_2D, tex_base[ch]);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	//Here we actually create the texture itself, notice
	//that we are using GL_LUMINANCE_ALPHA to indicate that
	//we are using 2 channel data.
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height,
		0, GL_LUMINANCE_ALPHA, GL_UNSIGNED_BYTE, expanded_data);

	//With the texture created, we don't need to expanded data anymore
	delete[] expanded_data;

	//So now we can create the display list
	glNewList(list_base + ch, GL_COMPILE);

	glBindTexture(GL_TEXTURE_2D, tex_base[ch]);

	glPushMatrix();

	//first we need to move over a little so that
	//the character has the right amount of space
	//between it and the one before it.
	glTranslatef((float)bitmap_glyph->left, 0, 0);

	//Now we move down a little in the case that the
	//bitmap extends past the bottom of the line 
	//(this is only true for characters like 'g' or 'y'.
	glTranslatef(0, (float)bitmap_glyph->top - bitmap.rows, 0);

	//Now we need to account for the fact that many of
	//our textures are filled with empty padding space.
	//We figure what portion of the texture is used by 
	//the actual character and store that information in 
	//the x and y variables, then when we draw the
	//quad, we will only reference the parts of the texture
	//that we contain the character itself.
	float	x = (float)bitmap.width / (float)width,
		y = (float)bitmap.rows / (float)height;

	//Here we draw the texturemaped quads.
	//The bitmap that we got from FreeType was not 
	//oriented quite like we would like it to be,
	//so we need to link the texture to the quad
	//so that the result will be properly aligned.

	glBegin(GL_QUADS);
	glTexCoord2d(0, 0); glVertex2f(0, (float)bitmap.rows);
	glTexCoord2d(0, y); glVertex2f(0, 0);
	glTexCoord2d(x, y); glVertex2f((float)bitmap.width, 0);
	glTexCoord2d(x, 0); glVertex2f((float)bitmap.width, (float)bitmap.rows);
	glEnd();
	glPopMatrix();
	glTranslatef((float)(face->glyph->advance.x >> 6), 0, 0);

	//Finnish the display list
	glEndList();

	FT_Done_Glyph(glyph);
}

FreeTypeFont::FreeTypeFont() {

}

FreeTypeFont::~FreeTypeFont() {
	clean();
}

void FreeTypeFont::init(const char * fname) {
	std::vector<std::string> lines;
	getCharSeparatedStringList(fname, lines, ' ');
	int h = 12;

	if (lines.size() > 1)
		sscanf(lines[1].c_str(), "%d", &h);

	//Allocate some memory to store the texture ids.
	textures = new GLuint[128];

	this->height = (float)h;

	//Create and initilize a freetype font library.
	FT_Library library;
	if (FT_Init_FreeType(&library))
		throwError("FT_Init_FreeType failed");

	//The object in which Freetype holds information on a given
	//font is called a "face".
	FT_Face face;

	//This is where we load in the font information from the file.
	//Of all the places where the code might die, this is the most likely,
	//as FT_New_Face will die if the font file does not exist or is somehow broken.
	if (FT_New_Face(library, lines[0].c_str(), 0, &face))
		throwError("FT_New_Face failed (there is probably a problem with your font file)");

	//For some twisted reason, Freetype measures font size
	//in terms of 1/64ths of pixels.  Thus, to make a font
	//h pixels high, we need to request a size of h*64.
	//(h << 6 is just a prettier way of writting h*64)
	FT_Set_Char_Size(face, h << 6, h << 6, 96, 96);

	//Here we ask opengl to allocate resources for
	//all the textures and displays lists which we
	//are about to create.  
	list_base = glGenLists(128);
	glGenTextures(128, textures);

	//This is where we actually create each of the fonts display lists.
	for (unsigned char i = 0;i < 128;i++)
		make_dlist(face, i, list_base, textures);

	//We don't need the face information now that the display
	//lists have been created, so we free the assosiated resources.
	FT_Done_Face(face);

	//Ditto for the library.
	FT_Done_FreeType(library);
}

void FreeTypeFont::clean() {
	glDeleteLists(list_base, 128);
	glDeleteTextures(128, textures);
	delete[] textures;
}

///Print function that works with free type fonts
void FreeTypeFont::print(double x, double y, const char* text) {
	std::vector<std::string> lines;
	getCharSeparatedStringList(text, lines);

	print(x, y, lines);
}

///Print function that works with free type fonts
void FreeTypeFont::print(double x, double y, const std::vector<std::string>& lines) {
	double h = getLineHeight();

	for (uint i = 0;i < lines.size();i++)
		print(x, y - h*i, lines[i]);
}

///Print function that works with free type fonts
void FreeTypeFont::print(double x, double y, const std::string& line) {
	// We want a coordinate system where things coresponding to window pixels.
	glPushAttrib(GL_TRANSFORM_BIT);
	GLint	viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0, viewport[2], 0, viewport[3]);
	glPopAttrib();

	uint font = list_base;

	glPushAttrib(GL_LIST_BIT | GL_CURRENT_BIT | GL_ENABLE_BIT | GL_TRANSFORM_BIT);

	glMatrixMode(GL_MODELVIEW);
	glDisable(GL_LIGHTING);
	glEnable(GL_TEXTURE_2D);
	glDisable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);

	glListBase(font);

	float modelview_matrix[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, modelview_matrix);

	glPushMatrix();

	glLoadIdentity();
	glTranslated(x, y, 0);
	glMultMatrixf(modelview_matrix);
	glCallLists(line.length(), GL_UNSIGNED_BYTE, line.c_str());

	glPopMatrix();

	glPopAttrib();

	glPushAttrib(GL_TRANSFORM_BIT);
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glPopAttrib();
}


