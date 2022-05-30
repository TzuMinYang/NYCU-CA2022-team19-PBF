#include "SpatialMap.h"
#include "Utils/Logger.h"
#include "Constants.h"

using namespace std;

SpatialMap::SpatialMap(double h) :
	particleMap(100, triple_hash)
{
	bucketSize = h;
}

void SpatialMap::clear() {
	particleMap.clear();
}

IntTriple SpatialMap::indexOfPosition(Particle p) {
	IntTriple index;
	index.x = (int)std::floor(p.x_star[0] / bucketSize);
	index.y = (int)std::floor(p.x_star[1] / bucketSize);
	index.z = (int)std::floor(p.x_star[2] / bucketSize);
	return index;
}

void SpatialMap::add(int i, Particle p) {
	IntTriple index = indexOfPosition(p);
	
	if (particleMap.count(index) > 0) {
		//Then the bucket already exists, so just add this to the bin
		particleMap[index].push_back(i);
	}
	// Otherwise we need to create the bucket
	else {
		particleMap[index] = vector<int>();
		particleMap[index].push_back(i);
	}
}

bool closeEnough(Particle p1, Particle p2) {
	V3D diff = p1.x_star - p2.x_star;
	if (abs(diff[0]) >= KERNEL_H || abs(diff[1]) >= KERNEL_H || abs(diff[2]) >= KERNEL_H) {
		return false;
	}
	return (diff.norm() < KERNEL_H);
}

void SpatialMap::findNeighbors(Particle &p_i, vector<Particle> &particles) {
	IntTriple index = indexOfPosition(p_i);
	p_i.neighbors.clear();

	// Look in all neighboring grid cells
	for (int dx = -1; dx <= 1; dx++) {
		for (int dy = -1; dy <= 1; dy++) {
			for (int dz = -1; dz <= 1; dz++) {
				IntTriple neighborIndex = index.addOffset(dx, dy, dz);
				if (particleMap.count(neighborIndex) == 0) continue;
				for (auto &j : particleMap[neighborIndex]) {
					if (closeEnough(p_i, particles[j])) {
						p_i.neighbors.push_back(j);
					}
				}
			}
		}
	}
}