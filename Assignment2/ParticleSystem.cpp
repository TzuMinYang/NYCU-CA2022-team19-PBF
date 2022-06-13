#ifdef _WIN32
#include <include/glew.h>
#else
#include <GL/glew.h>
#endif

#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

#include "ParticleSystem.h"
#include "GUILib/OBJReader.h"
#include "Utils/Logger.h"
#include "Constants.h"
#include <math.h>

#include<iostream>
using namespace std;

GLuint makeBoxDisplayList();


ParticleSystem::ParticleSystem(vector<ParticleInit>& initialParticles)
	: particleMap(KERNEL_H)
{
	int numParticles = initialParticles.size();
	Logger::consolePrint("Created particle system with %d particles", numParticles);
	drawParticles = true;
	count = 0;

	// Create all particles from initial data
	for (auto ip : initialParticles) {
		Particle p;
		p.x_i = ip.position;
		p.v_i = ip.velocity;
		p.x_star = p.x_i;
		p.neighbors.clear();
		particles.push_back(p);
	}

	// Create floor and walls
	CollisionPlane floor(P3D(0, -1, 0), V3D(0, 1, 0));
	CollisionPlane left_wall(P3D(-1, 0, 0), V3D(1, 0, 0));
	CollisionPlane right_wall(P3D(1, 0, 0), V3D(-1, 0, 0));
	CollisionPlane front_wall(P3D(0, 0, -1), V3D(0, 0, 1));
	CollisionPlane back_wall(P3D(0, 0, 1), V3D(0, 0, -1));
	CollisionPlane ceiling(P3D(0, 2, 0), V3D(0, -1, 0));

	planes.push_back(floor);
	planes.push_back(left_wall);
	planes.push_back(right_wall);
	planes.push_back(front_wall);
	planes.push_back(back_wall);
	planes.push_back(ceiling);

	// Arrays to be passed to OpenGL
	positionArray = vector<double>(numParticles * 3);
	pointsIndexArray = vector<unsigned int>(numParticles);

	for (int i = 0; i < numParticles; i++) {
		pointsIndexArray[i] = i;
	}

	boxList = makeBoxDisplayList();
}

ParticleSystem::~ParticleSystem() {
	if (boxList >= 0) {
		glDeleteLists(boxList, 1);
	}
}

bool ParticleSystem::drawParticles = true;
bool ParticleSystem::enableGravity = true;

P3D ParticleSystem::getPositionOf(int i) {
	return particles[i].x_i;
}

// Set the position of particle i.
void ParticleSystem::setPosition(int i, P3D x) {
	particles[i].x_i = x;
	particles[i].x_star = x;
}

// Set the velocity of particle i.
void ParticleSystem::setVelocity(int i, V3D v) {
	particles[i].v_i = v;
}

int ParticleSystem::particleCount() {
	return particles.size();
}

// Gravitational constant.
const V3D GRAVITY = V3D(0, -9.8, 0);

// Applies external forces to particles in the system.
// This is currently limited to just gravity.
void ParticleSystem::applyForces(double delta) {
	if (enableGravity) {
		// Assume all particles have unit mass to simplify calculations.
		for (auto& p : particles) {
			p.v_i += (GRAVITY * delta);
		}
	}
}

// Integrate one time step.
void ParticleSystem::integrate_PBF(double delta) {
	applyForces(delta);
	// Predict positions for this timestep.
	for (auto& p : particles) {
		p.x_star = p.x_i + (p.v_i * delta);
	}

	// Find neighbors for all particles.
	particleMap.clear();
	for (int i = 0; i < particles.size(); i++) {
		particleMap.add(i, particles[i]);
	}

	// hash table
	for (auto& p_i : particles) {
		particleMap.findNeighbors(p_i, particles);
	}

	// brute_force_findNeighbors
	/*for (int i = 0; i < particles.size(); ++i) {
		particleMap.brute_force_findNeighbors(particles[i], i, particles);
	}*/

	// TODO: implement the solver loop.
	int iteration = 0;
	while (iteration < SOLVER_ITERATIONS) {
		// calculate lambda
		for (int i = 0; i < particles.size(); i++) {
			particles[i].density = computeDensity(i);
			particles[i].lambda_i = computeLambda(i);
		}

		// calculate delta position
		for (auto& p : particles) {
			p.delta_p = computeDeltaP(p);
			// collision detection
			for (auto& cp : planes)
				//p = cp.handleCollision_particle(p);
				p.x_star = cp.handleCollision(p.x_star);
		}

		// update position
		for (auto& p : particles) {
			p.x_star += p.delta_p;
		}
		iteration++;
	}

	for (auto& p : particles) {
		// TODO: edit this loop to apply vorticity and viscosity.
		p.v_i = (p.x_star - p.x_i) / delta;

		// apply vorticity
		p.v_i += (computeVorticity(p) * delta);
		// apply viscosity
		p.v_i += computeViscosity(p);
		p.x_i = p.x_star;
	}
}
// functions added by a
double ParticleSystem::computeDensity(int i) {
	double density = 0.0;
	// neighbor for optimization
	for (auto& j : particles[i].neighbors) {
		// mass = 1 -> ignore mass
		density += poly6WKernel(particles[i].x_star - particles[j].x_star, KERNEL_H);
	}
	return density;
}

// functions added by a
double  ParticleSystem::computeLambda(int i) {
	double constraint = (particles[i].density / REST_DENSITY) - 1.0;
	double constraintGradientSum = 0.0;
	for (auto& k : particles[i].neighbors) {
		// calculate constraint gradient
		V3D constraintGradient = V3D();
		if (k == i) {
			for (auto& j : particles[i].neighbors) {
				constraintGradient += spikyWKernel(particles[i].x_star - particles[j].x_star, KERNEL_H, false);
			}
		}
		else {
			constraintGradient = -spikyWKernel(particles[i].x_star - particles[k].x_star, KERNEL_H, true);
		}
		constraintGradient /= REST_DENSITY;
		constraintGradientSum += constraintGradient.length2();
	}

	return -constraint / (constraintGradientSum + CFM_EPSILON);
}

// functions added by a
V3D ParticleSystem::computeDeltaP(Particle i) {
	V3D deltaP = V3D();
	for (auto& j : i.neighbors) {
		deltaP += (spikyWKernel(i.x_star - particles[j].x_star, KERNEL_H, false) * (i.lambda_i + particles[j].lambda_i + computeSurfaceTension(i, particles[j])));
	}

	deltaP /= REST_DENSITY;
	return deltaP;
}

// functions added by s
V3D ParticleSystem::computeVorticity(Particle i) {
	i.vorticity_W = V3D();
	for (auto& j : i.neighbors) {
		i.vorticity_W += (particles[j].v_i - i.v_i).cross(spikyWKernel(i.x_star - particles[j].x_star, KERNEL_H, false));
	}

	// not sure: how to get vorticity_N ?
	V3D mass_center = (i.x_i + i.x_star) / 2;
	i.vorticity_N = (mass_center - i.x_star) / (mass_center - i.x_star).norm();

	return VORTICITY_EPSILON * (i.vorticity_N.cross(i.vorticity_W));
}

// functions added by a
double ParticleSystem::poly6WKernel(V3D r, double h) {

	if (r.length() > h || r.length() < 0) {
		return 0.0;
	}
	return 315.0 / (64.0 * PI * pow(h, 9)) * pow(pow(h, 2) - r.length2(), 3);
}

// functions added by a
V3D ParticleSystem::spikyWKernel(V3D r, double h, bool dir) {
	if (r.length() > h || r.length() < 0) {
		return V3D();
	}
	V3D direction = dir ? r.unit() : -r.unit();
	return -45.0 / (PI * pow(h, 6)) * pow(h - r.length(), 2) * direction;
}

// functions added by a
double ParticleSystem::computeSurfaceTension(Particle i, Particle j) {
	//return 0.0001;
	return -TENSILE_K * pow(poly6WKernel(i.x_star - j.x_star, KERNEL_H) / poly6WKernel(V3D(TENSILE_DELTA_Q, 0.0, 0.0), KERNEL_H), TENSILE_N);
}

// functions added by a
V3D ParticleSystem::computeViscosity(Particle i) {
	V3D delta_v = V3D();

	for (auto& j : i.neighbors) {
		delta_v += (particles[j].v_i - i.v_i) * poly6WKernel(i.x_star - particles[j].x_star, KERNEL_H);
	}

	return VISCOSITY_C * delta_v;
}

// Code for drawing the particle system is below here.

GLuint makeBoxDisplayList() {

	GLuint index = glGenLists(1);

	glNewList(index, GL_COMPILE);
	glLineWidth(3);
	glColor3d(0, 0, 0);
	glBegin(GL_LINES);
	glVertex3d(-1, -1, -1);
	glVertex3d(-1, -1, 1);

	glVertex3d(-1, -1, 1);
	glVertex3d(1, -1, 1);

	glVertex3d(1, -1, 1);
	glVertex3d(1, -1, -1);

	glVertex3d(1, -1, -1);
	glVertex3d(-1, -1, -1);

	glVertex3d(-1, -1, -1);
	glVertex3d(-1, 2, -1);

	glVertex3d(-1, -1, 1);
	glVertex3d(-1, 2, 1);

	glVertex3d(1, -1, 1);
	glVertex3d(1, 2, 1);

	glVertex3d(1, -1, -1);
	glVertex3d(1, 2, -1);

	glVertex3d(-1, 2, -1);
	glVertex3d(-1, 2, 1);

	glVertex3d(-1, 2, 1);
	glVertex3d(1, 2, 1);

	glVertex3d(1, 2, 1);
	glVertex3d(1, 2, -1);

	glVertex3d(1, 2, -1);
	glVertex3d(-1, 2, -1);
	glEnd();

	glEndList();
	return index;
}

void ParticleSystem::drawParticleSystem() {

	int numParticles = particles.size();
	int i = 0;

	glCallList(boxList);

	// Copy particle positions into array
	positionArray.clear();
	pointsIndexArray.clear();
	for (auto& p : particles) {
		positionArray.push_back(p.x_i[0]);
		positionArray.push_back(p.x_i[1]);
		positionArray.push_back(p.x_i[2]);
		pointsIndexArray.push_back(i);
		i++;
	}

	if (drawParticles && numParticles > 0) {
		// Draw all particles as blue dots
		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(3, GL_DOUBLE, 0, &(positionArray.front()));

		glColor4d(0.2, 0.2, 0.8, 1);
		glPointSize(15);
		glDrawElements(GL_POINTS, numParticles, GL_UNSIGNED_INT, &(pointsIndexArray.front()));

		glDisableClientState(GL_VERTEX_ARRAY);
	}

}
