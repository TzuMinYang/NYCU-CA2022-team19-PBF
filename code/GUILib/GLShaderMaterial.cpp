//TODO: have a global option of loading or not loading textures/shaders, probably here in GLContentManager. For quick debugging, since loading these things can take a long time...

#ifdef _WIN32
#include <include/glew.h>
#else
#include <GL/glew.h>
#endif

#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

#include "GLShaderMaterial.h"
#include <Utils/Utils.h>
#include <GUILib/GLContentManager.h>

enum {
	SHM_SHADER = 1,
	SHM_END,
	SHM_FLOAT,
	SHM_FLOAT2,
	SHM_FLOAT3,
	SHM_FLOAT4,
	SHM_TEXTURE,
	SHM_LIGHT_SETTING,
	SHM_LIGHT_SETTING_END,
	SHM_LIGHT_POSITION,
	SHM_COLOR,
	SHM_NAME,
};

KeyWord smKeywords[] = {
	{ "name", SHM_NAME},
	{ "shader", SHM_SHADER },
	{ "float2", SHM_FLOAT2 },
	{ "float3", SHM_FLOAT3 },
	{ "float4", SHM_FLOAT4 },
	{ "float", SHM_FLOAT },
	{ "color", SHM_COLOR },
	{ "texture", SHM_TEXTURE },
	{ "lights", SHM_LIGHT_SETTING },
	{ "/lights", SHM_LIGHT_SETTING_END },
	{ "lightPosition", SHM_LIGHT_POSITION },
	{ "/material", SHM_END }
};

void GLShaderMaterialLightSetting::readFromFile(FILE *f) {
	hasSetting = true;
	lightPositions.clear();

	char buffer[200];
	while (!feof(f)) {
		fgets(buffer, 200, f);
		if (strlen(buffer) > 195)
			throwError("The input file contains a line that is longer than ~200 characters - not allowed");
		char *line = lTrim(buffer);
		int lineType = getLineType(line, smKeywords, sizeof(smKeywords) / sizeof(smKeywords[0]));
		switch (lineType) {
		case SHM_LIGHT_POSITION: {
				double x, y, z;
				if (sscanf(line, "%lf %lf %lf", &x, &y, &z) != 3)
					throwError("Incorrect light position in material description - 3 float arguments expected");
				lightPositions.push_back(P3D(x,y,z));
			}
			break;
		case SHM_LIGHT_SETTING_END:
			return;
		default:
			break;
		}
	}
}

GLShaderMaterial::GLShaderMaterial()
	:	shaderProgram(NULL) {
}

GLShaderMaterial::~GLShaderMaterial(){
}

void GLShaderMaterial::readFromFile(const char* fName) {
	FILE* f = fopen(fName, "r");
	readFromFile(f);
	fclose(f);
}

void GLShaderMaterial::readFromFile(FILE *f){
	char buffer[200];
	char nameBuffer[200];

	while (!feof(f)){
		fgets(buffer, 200, f);
		if (strlen(buffer) > 195)
			throwError("The input file contains a line that is longer than ~200 characters - not allowed");
		char *line = lTrim(buffer);
		int lineType = getLineType(line, smKeywords, sizeof(smKeywords) / sizeof(smKeywords[0]));
		switch (lineType){
		case SHM_NAME:
			sscanf(line, "%s", nameBuffer);
			materialName = std::string() + nameBuffer;
			break;
		case SHM_SHADER:
			sscanf(line, "%s", nameBuffer);
			shaderProgram = GLContentManager::getShaderProgram(nameBuffer);
			if (shaderProgram == NULL)
				Logger::consolePrint("Shader program not loaded: \'%s\'", nameBuffer);
			break;
		case SHM_END:
			// Done
			return;
		case SHM_FLOAT:{
				double x;
				if (sscanf(line, "%s %lf", nameBuffer, &x) != 2)
					throwError("Incorrect material description - 1 float argument expected");
				setFloatParam(nameBuffer, (float)x);
			}
			break;
		case SHM_FLOAT2:{
				double x, y;
				if (sscanf(line, "%s %lf %lf", nameBuffer, &x, &y) != 3)
					throwError("Incorrect material description - 2 float arguments expected");
				setFloatParam(nameBuffer, (float)x, (float)y);
			}
			break;
		case SHM_FLOAT3:{
				double x, y, z;
				if (sscanf(line, "%s %lf %lf %lf", nameBuffer, &x, &y, &z) != 4)
					throwError("Incorrect material description - 3 float arguments expected");
				setFloatParam(nameBuffer, (float)x, (float)y, (float)z);
			}
			break;
		case SHM_FLOAT4:{
				double x, y, z, w;
				if (sscanf(line, "%s %lf %lf %lf %lf", nameBuffer, &x, &y, &z, &w) != 5)
					throwError("Incorrect material description - 5 float arguments expected");
				setFloatParam(nameBuffer, (float)x, (float)y, (float)z, (float)w);
			}
			break;
		case SHM_TEXTURE:{
				char idBuffer[200];
				sscanf(line, "%s %s", idBuffer, nameBuffer);
				GLTexture *tex = GLContentManager::getTexture(nameBuffer);
				setTextureParam(idBuffer, tex);
			}
			break;
		case SHM_LIGHT_SETTING:
			lightSetting.readFromFile(f);
			break;
		case SHM_COLOR:
			if (sscanf(line, "%lf %lf %lf %lf", &r, &g, &b, &a) != 5)
				throwError("Incorrect color description - 4 float arguments expected");
			break;
		default:
			break;
			// throwError("Incorrect material description: \'%s\' - unexpected line.", buffer);
		}
	}
}

void GLShaderMaterial::setShaderProgram(GLShaderProgram *shaderProgram){
	this->shaderProgram = shaderProgram;
}

void GLShaderMaterial::setTextureParam(const char *name, GLTexture *tex){
	texParamMap[name] = tex;
}

void GLShaderMaterial::setFloatParam(const char *name, float val){
	floatParamMap[name] = GLShaderMaterialFloatParam(val);
}

void GLShaderMaterial::setFloatParam(const char *name, float x, float y){
	floatParamMap[name] = GLShaderMaterialFloatParam(x, y);
}

void GLShaderMaterial::setFloatParam(const char *name, float x, float y, float z){
	floatParamMap[name] = GLShaderMaterialFloatParam(x, y, z);
}

void GLShaderMaterial::setFloatParam(const char *name, float x, float y, float z, float w){
	floatParamMap[name] = GLShaderMaterialFloatParam(x, y, z, w);
}

void GLShaderMaterial::apply(){
	if (lightSetting.hasSetting){
		glPushAttrib(GL_LIGHTING_BIT);

		for (uint i=0; i<lightSetting.lightPositions.size(); i++){
			GLfloat pos[] = {(GLfloat)lightSetting.lightPositions[i][0], (GLfloat)lightSetting.lightPositions[i][1], (GLfloat)lightSetting.lightPositions[i][2], 1};
			glLightfv(GL_LIGHT0 + i, GL_POSITION, pos);
		}
	}

	glColor4d(r, g, b, a);

	if (shaderProgram == NULL){
		GLShaderProgram::unbind();
		glDisable(GL_TEXTURE_2D);
		return;
	}

	for (auto it = floatParamMap.begin(); it!=floatParamMap.end(); ++it){
		switch (it->second.dim){
		case 1:
			shaderProgram->setUniform(it->first.c_str(), it->second.val[0]);
			break;
		case 2:
			shaderProgram->setUniform(it->first.c_str(), it->second.val[0], it->second.val[1]);
			break;
		case 3:
			shaderProgram->setUniform(it->first.c_str(), it->second.val[0], it->second.val[1], it->second.val[2]);
			break;
		case 4:
			shaderProgram->setUniform(it->first.c_str(), it->second.val[0], it->second.val[1], it->second.val[2], it->second.val[3]);
			break;
		default:
			break;
		}
	}
	
	if (!texParamMap.empty()){
		int nTextures = (int)texParamMap.size();
		glEnable(GL_TEXTURE_2D);

		auto it = texParamMap.begin();
		for (int i=nTextures-1; i>=0; i--, it++){
			if (it->second != NULL){
				shaderProgram->setUniformSampler(it->first.c_str(), i);
				glActiveTexture(GL_TEXTURE0 + i);
				it->second->activate();
			}
		}

		glActiveTexture(GL_TEXTURE0);
	}

	shaderProgram->bind();
}

void GLShaderMaterial::end(){
	if (lightSetting.hasSetting)
		glPopAttrib();
	GLShaderProgram::unbind();
}

