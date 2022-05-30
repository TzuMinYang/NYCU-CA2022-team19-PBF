#ifdef _WIN32
#include <include/glew.h>
#else
#include <GL/glew.h>
#endif

#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

#include "GLContentManager.h"
#include "FreeType.h"
#include "GLTexture.h"
#include "GLShader.h"
#include "GLShaderMaterial.h"
#include "GLMesh.h"
#include "OBJReader.h"

std::map<std::string, FreeTypeFont *> GLContentManager::fontMap;
std::map<std::string, GLTexture *> GLContentManager::texMap;
std::map<std::string, GLShader *> GLContentManager::shaderMap;
std::map<std::string, GLShaderProgram *> GLContentManager::shaderProgramMap;
std::map<std::string, GLShaderMaterial *> GLContentManager::shaderMaterialMap;
std::map<std::string, GLMesh *> GLContentManager::meshMap;

void GLContentManager::destroyAllContent() {
	//should destory all content here...
	for (auto it = fontMap.begin(); it != fontMap.end(); ++it)
		delete it->second;

	for (auto it = texMap.begin(); it != texMap.end(); ++it)
		delete it->second;

	for (auto it = shaderMap.begin(); it != shaderMap.end(); ++it)
		delete it->second;

	for (auto it = shaderProgramMap.begin(); it != shaderProgramMap.end(); ++it)
		delete it->second;

	for (auto it = meshMap.begin(); it != meshMap.end(); ++it)
		delete it->second;
}

FreeTypeFont* GLContentManager::getFont(const char* fontName) {
	auto it = fontMap.find(fontName);
	if (it == fontMap.end()){
		FreeTypeFont *font = new FreeTypeFont();
		font->init(fontName);
		fontMap[fontName] = font;
		return font;
	}
	else
		return it->second;
}


GLTexture *GLContentManager::getTexture(const char *fileName){
	auto it = texMap.find(fileName);
	if (it == texMap.end()){
		GLTexture *tex = new GLTexture((char *)fileName);
		texMap[fileName] = tex;
		return tex;
	}
	else
		return it->second;
}

GLShader *GLContentManager::getShader(const char *fileName, uint type){
	auto it = shaderMap.find(fileName);
	if (it == shaderMap.end()){
		GLShader *shader = new GLShader();
		shader->loadFromFile(fileName, type);
		shaderMap[fileName] = shader;
		return shader;
	}
	else
		return it->second;
}

GLMesh *GLContentManager::getGLMesh(const char *fileName) {
	auto it = meshMap.find(fileName);
	if (it == meshMap.end()) {
		GLMesh *mesh = OBJReader::loadOBJFile(fileName);
		meshMap[fileName] = mesh;
		return mesh;
	}
	else
		return it->second;
}

GLShaderProgram *GLContentManager::getShaderProgram(const char *name){
	auto it = shaderProgramMap.find(name);
	if (it == shaderProgramMap.end()){
		throwError("Shader program not found");
		return NULL;
	}
	else
		return it->second;
}

GLShaderMaterial *GLContentManager::getShaderMaterial(const char *name) {
	auto it = shaderMaterialMap.find(name);
	if (it == shaderMaterialMap.end()) {
//		throwError("Shader material not found");
		return NULL;
	}
	else
		return it->second;
}

//returns true if the material is added, false if it is already in the list, and thus ignored
bool GLContentManager::addShaderMaterial(const char* name, GLShaderMaterial* shaderMaterial) {
	auto it = shaderMaterialMap.find(name);
	if (it == shaderMaterialMap.end()) {
		shaderMaterialMap[name] = shaderMaterial;
		return true;
	}
	return false;
}

void GLContentManager::addShaderProgram(const char *name, const char *vsFilename, const char *fsFilename){
	auto it = shaderProgramMap.find(name);
	if (it == shaderProgramMap.end()){
		GLShader *vs = getShader(vsFilename, GL_VERTEX_SHADER);
		GLShader *fs = getShader(fsFilename, GL_FRAGMENT_SHADER);

		GLShaderProgram *prog = new GLShaderProgram();
		prog->load(vs, fs);
		shaderProgramMap[name] = prog;
	}
}

