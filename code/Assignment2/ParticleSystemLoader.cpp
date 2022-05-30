#include "ParticleSystemLoader.h"
#include "GUILib/OBJReader.h"
#include <iostream>
#include <fstream>
#include "Utils/Logger.h"
#include <sstream>
#include <unordered_set>
#include "Constants.h"

struct SpringData {
	int p1;
	int p2;
	double stiffness;
};

struct ZeroLengthSpringData {
	int p;
	P3D x0;
	double stiffness;
};

struct pair_hash {
	inline std::size_t operator()(const std::pair<int, int> & v) const {
		return v.first * 31 + v.second;
	}
};

ParticleSystem* ParticleSystemLoader::loadFromOBJ(string filename) {
	GLMesh* mesh = 0;
	try {
		mesh = OBJReader::loadOBJFile(filename.c_str());
	}

	catch (char* e) {
		Logger::consolePrint("%s", e);
		return 0;
	}
	
	vector<ParticleInit> ps;
	unordered_set<pair<int, int>, pair_hash> edges;

	// First obtain vertices; these will become particles
	for (int i = 0; i < mesh->vertexCount; i++) {
		int base = 3 * i;
		ParticleInit p;
		p.position = P3D(mesh->vertexList[base], mesh->vertexList[base + 1], mesh->vertexList[base + 2]);
		p.velocity = V3D(0, 0, 0);
		p.mass = 1;
		ps.push_back(p);
	}

	// Logger::consolePrint("%d edges", edges.size());

	ParticleSystem* system = new ParticleSystem(ps);
	return system;
}