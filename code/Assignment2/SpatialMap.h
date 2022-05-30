#pragma once

#include "Particle.h"
#include <unordered_map>

struct IntTriple {
	int x, y, z;

	bool operator==(const IntTriple &other) const {
		return (x == other.x && y == other.y && z == other.z);
	}

	IntTriple addOffset(int dx, int dy, int dz) {
		IntTriple t2;
		t2.x = x + dx;
		t2.y = y + dy;
		t2.z = z + dz;
		return t2;
	}
};

inline std::size_t triple_hash(const IntTriple &t)  {
	return 31 * t.x + 73 * t.y + t.z;
}

typedef std::unordered_map<IntTriple, std::vector<int>, decltype(&triple_hash)> MapType;

class SpatialMap {
private:
	MapType particleMap;
	double bucketSize;

public:
	SpatialMap(double h);
	void clear();
	IntTriple indexOfPosition(Particle p);
	void add(int i, Particle p);
	void findNeighbors(Particle &p, std::vector<Particle> &particles);
};