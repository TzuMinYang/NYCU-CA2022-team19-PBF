#pragma once

#include <MathLib/MathLib.h>

/**
	This class is used to load and store vertex- and fragment shaders.
*/
class GLShader{
	friend class GLShaderProgram;

private:
	int id;

public:
	GLShader();
	~GLShader();

	/*
		Load the shader code from a file.
		'type' can be either GL_VERTEX_SHADER or GL_FRAGMENT_SHADER.
	*/
	void loadFromFile(const char *fileName, uint type);
	void create(const char *code, uint type);
	void unload();
};

/**
	This class is used to load and bind shader programs.
*/
class GLShaderProgram{
private:
	int programID;
	GLShader *vertexShader;
	GLShader *fragmentShader;

public:
	GLShaderProgram();
	~GLShaderProgram();

	void load(GLShader *vertexShader, GLShader *fragmentShader);
	void unload();

	void bind();
	static void unbind();

	void setUniform(const char *name, int value);
	void setUniform(const char *name, float value);
	void setUniform(const char *name, float x, float y);
	void setUniform(const char *name, float x, float y, float z);
	void setUniform(const char *name, float x, float y, float z, float w);
	void setUniformSampler(const char *name, int channel);
};
