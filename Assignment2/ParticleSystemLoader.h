#pragma once

#include "GUILib/GLMesh.h"
#include "ParticleSystem.h"
#include <string>

using namespace std;

class ParticleSystemLoader {
public:
	static ParticleSystem* loadFromMSS(string filename);
	static ParticleSystem* loadFromOBJ(string filename);
};