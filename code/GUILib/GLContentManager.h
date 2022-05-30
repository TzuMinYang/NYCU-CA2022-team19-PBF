#pragma once

#include <MathLib/MathLib.h>
#include <map>
#include <string>

class FreeTypeFont;
class GLTexture;
class GLShader;
class GLShaderProgram;
class GLMesh;
class GLShaderMaterial;

/**
	This class is used to manage OpenGL resources, so that they only need to be loaded once.
*/
class GLContentManager{
private:
	static std::map<std::string, FreeTypeFont *> fontMap;
	static std::map<std::string, GLTexture *> texMap;
	static std::map<std::string, GLMesh *> meshMap;
	static std::map<std::string, GLShader *> shaderMap;					//keeps a list of fragment and vertex shaders
	static std::map<std::string, GLShaderProgram *> shaderProgramMap;	//keeps a list of shader programs, which are frag/vert shader pairs + knows how to pass in parameters to them
	static std::map<std::string, GLShaderMaterial *> shaderMaterialMap;	//the shader materials keep pointers to shader programs and provide appropriate parameters for all of its variables

public:
	static FreeTypeFont *getFont(const char* fontName);
	static GLTexture *getTexture(const char *filename);
	static GLMesh *getGLMesh(const char *filename);
	static GLShader *getShader(const char *filename, uint type);
	static GLShaderProgram *getShaderProgram(const char *name);
	static GLShaderMaterial *getShaderMaterial(const char *name);
	static void addShaderProgram(const char *name, const char *vsFilename, const char *fsFilename);
	//returns true if the material is added, false if it is already in the list, and thus ignored
	static bool addShaderMaterial(const char* name, GLShaderMaterial* shaderMaterial);

	static void destroyAllContent();
};


