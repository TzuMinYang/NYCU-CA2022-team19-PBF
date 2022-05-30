#pragma once

#include "GLTexture.h"
#include "GLShader.h"
#include <map>
#include <vector>
#include <string>
#include <MathLib/P3D.h>

struct GLShaderMaterialFloatParam {
	GLShaderMaterialFloatParam() {}
	GLShaderMaterialFloatParam(const GLShaderMaterialFloatParam &o) : dim(o.dim) { val[0]=o.val[0]; val[1]=o.val[1]; val[2]=o.val[2]; val[3]=o.val[3]; }

	GLShaderMaterialFloatParam(float val) : dim(1) { this->val[0] = val; }
	GLShaderMaterialFloatParam(float x, float y) : dim(2) { val[0]=x; val[1]=y; }
	GLShaderMaterialFloatParam(float x, float y, float z) : dim(3) { val[0]=x; val[1]=y; val[2]=z; }
	GLShaderMaterialFloatParam(float x, float y, float z, float w) : dim(4) { val[0]=x; val[1]=y; val[2]=z; val[3]=w; }

	float val[4];
	int dim;
};

struct GLShaderMaterialLightSetting {
	GLShaderMaterialLightSetting() : hasSetting(false) {}

	void readFromFile(FILE *f);

	bool hasSetting;
	std::vector<P3D> lightPositions;
};

/**
	This class is used to store a material, composed of a shader program and shader parameters, textures, etc.
*/
class GLShaderMaterial {
private:
	std::string materialName;
	GLShaderProgram *shaderProgram;
	std::map<std::string, GLTexture *> texParamMap;
	std::map<std::string, GLShaderMaterialFloatParam> floatParamMap;
	GLShaderMaterialLightSetting lightSetting;
public:
	//base colour of the material (white by default)
	double r = 1, g = 1, b = 1, a = 1;
	void setColor(double r, double g, double b, double a) {
		this->r = r;
		this->g = g;
		this->b = b;
		this->a = a;
	}
public:
	GLShaderMaterial();
	~GLShaderMaterial();

	void readFromFile(FILE *f);
	void readFromFile(const char* fName);

	void setShaderProgram(GLShaderProgram *shaderProgram);
	bool hasShaderProgram() const { return shaderProgram != NULL; }

	void setTextureParam(const char *name, GLTexture *tex);
	void setFloatParam(const char *name, float val);
	void setFloatParam(const char *name, float x, float y);
	void setFloatParam(const char *name, float x, float y, float z);
	void setFloatParam(const char *name, float x, float y, float z, float w);

	const GLShaderMaterialLightSetting &getLightSetting() const { return lightSetting; }
	void setLightSetting(const GLShaderMaterialLightSetting &lightSetting) { this->lightSetting = lightSetting; }

	void apply();
	void end();

	bool hasTextureParam() const { return !texParamMap.empty(); }

	std::string getMaterialName() {
		return materialName;
	}
};
