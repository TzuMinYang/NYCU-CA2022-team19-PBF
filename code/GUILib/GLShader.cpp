#ifdef _WIN32
#include <include/glew.h>
#else
#include <GL/glew.h>
#endif

#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

#include "GLShader.h"
#include <Utils/Utils.h>
#include <fstream>
#include <Utils/Logger.h>


GLShader::GLShader() : id(-1) {
}

GLShader::~GLShader() {
	//unload();
}

void GLShader::loadFromFile(const char *fileName, uint type) {
	FILE *f = fopen(fileName, "r");
	if (f == NULL)
		throwError("Shader file not found");

	fseek(f, 0, SEEK_END);
	long fsize = ftell(f);
	rewind(f);

	char *code = new char[fsize + 1];
	size_t sizeRead = fread(code, 1, fsize, f);
	code[sizeRead] = '\0';

	fclose(f);
	
	create(code, type);
	
	delete[] code;
}

void GLShader::create(const char *code, uint type) {
	id = glCreateShader(type);
	glShaderSource(id, 1, &code, NULL);
	glCompileShader(id);

	int val;
	glGetShaderiv(id, GL_COMPILE_STATUS, &val);
	if (val == GL_FALSE || glGetError() != 0){
		char msg[512];
		int len;
		glGetShaderInfoLog(id, 512, &len, msg);
//		Logger::consolePrint(msg);
		throwError((std::string("Error while compiling shader: ") + msg).c_str());
	}
}

void GLShader::unload(){
	if (id >= 0)
		glDeleteShader(id);
	id = -1;
}

GLShaderProgram::GLShaderProgram()
	:	programID(-1),
		vertexShader(NULL),
		fragmentShader(NULL){
}

GLShaderProgram::~GLShaderProgram(){
//	unload();
}

void GLShaderProgram::load(GLShader *vertexShader, GLShader *fragmentShader) {
	programID = glCreateProgram();
	this->vertexShader = vertexShader;
	if (vertexShader)
		glAttachShader(programID, vertexShader->id);
	this->fragmentShader = fragmentShader;
	if (fragmentShader)
		glAttachShader(programID, fragmentShader->id);

	glLinkProgram(programID);

	GLint val;
	glGetProgramiv(programID, GL_LINK_STATUS, &val);
	if (val == GL_FALSE || glGetError() != 0) {
		char msg[512];
		int len;
		glGetProgramInfoLog(programID, 512, &len, msg);
		Logger::consolePrint(msg);
		throwError("Error while compiling shader, see log");
	}
}

void GLShaderProgram::unload(){
	if (vertexShader)
		glDetachShader(programID, vertexShader->id);
	if (fragmentShader)
		glDetachShader(programID, fragmentShader->id);

	if (programID >= 0)
		glDeleteProgram(programID);

	programID = -1;
	vertexShader = NULL;
	fragmentShader = NULL;
}

void GLShaderProgram::bind(){
	glUseProgram(programID);
}

void GLShaderProgram::unbind(){
//	Logger::consolePrint("Before use program (0): %s \n", gluErrorString(glGetError()));
	glUseProgram(0);
//	Logger::consolePrint("After use program (0): %s \n", gluErrorString(glGetError()));
}

void GLShaderProgram::setUniform(const char *name, int value){
	bind();
	GLint loc = glGetUniformLocation(programID, name);
	if (loc >= 0)
		glUniform1i(loc, value);
	assert(!glGetError());
}

void GLShaderProgram::setUniform(const char *name, float value){
	bind();
	GLint loc = glGetUniformLocation(programID, name);
	if (loc >= 0)
		glUniform1f(loc, value);
//	Logger::consolePrint("Error for set uniform value for %s: %s\n", name, gluErrorString(glGetError()));
	assert(!glGetError());
}

void GLShaderProgram::setUniform(const char *name, float x, float y){
	bind();
	GLint loc = glGetUniformLocation(programID, name);
	if (loc >= 0)
		glUniform2f(loc, x, y);
	assert(!glGetError());
}

void GLShaderProgram::setUniform(const char *name, float x, float y, float z){
	bind();
	GLint loc = glGetUniformLocation(programID, name);
	if (loc >= 0)
		glUniform3f(loc, x, y, z);
	assert(!glGetError());
}

void GLShaderProgram::setUniform(const char *name, float x, float y, float z, float w){
	bind();
	GLint loc = glGetUniformLocation(programID, name);
	if (loc >= 0)
		glUniform4f(loc, x, y, z, w);
	assert(!glGetError());
}

void GLShaderProgram::setUniformSampler(const char *name, int channel){
	bind();
	GLint loc = glGetUniformLocation(programID, name);
	if (loc >= 0)
		glUniform1i(loc, channel);
	assert(!glGetError());
}

