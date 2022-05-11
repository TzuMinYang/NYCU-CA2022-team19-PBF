#include "integrator.h"

#include "configs.h"

void ExplicitEuler::integrate(const std::vector<Particles *> &particles, std::function<void(void)>) const {
  // TODO: Integrate velocity and acceleration
  //   1. Integrate velocity.
  //   2. Integrate acceleration.
  //   3. You should not compute position using acceleration. Since some part only update velocity. (e.g. impulse)
  // Note:
  //   1. You don't need the simulation function in explicit euler.
  //   2. You should do this first because it is very simple. Then you can chech your collision is correct or not.
  //   3. This can be done in 5 lines. (Hint: You can add / multiply all particles at once since it is a large matrix.)
  for (const auto &p : particles) {
    p->position() += p->velocity() * deltaTime;
    p->velocity() += p->acceleration() * deltaTime;
  }
}

void ImplicitEuler::integrate(const std::vector<Particles *> &particles,
                              std::function<void(void)> simulateOneStep) const {
  // TODO: Integrate velocity and acceleration
  //   1. Backup original particles' data.
  //   2. Integrate velocity and acceleration using explicit euler to get Xn+1.
  //   3. Compute refined Xn+1 using (1.) and (2.).
  // Note:
  //   1. Use simulateOneStep with modified position and velocity to get Xn+1.
  std::vector<Particles> current;
  for (const auto &p : particles) {
    current.emplace_back(*p);
    p->position() += p->velocity() * deltaTime;
    p->velocity() += p->acceleration() * deltaTime;
  }
  simulateOneStep();
  for (size_t i = 0; i < current.size(); ++i) {
    particles[i]->position() = current[i].position() + particles[i]->velocity() * deltaTime;
    particles[i]->velocity() = current[i].velocity() + particles[i]->acceleration() * deltaTime;
  }
}

void MidpointEuler::integrate(const std::vector<Particles *> &particles,
                              std::function<void(void)> simulateOneStep) const {
  // TODO: Integrate velocity and acceleration
  //   1. Backup original particles' data.
  //   2. Integrate velocity and acceleration using explicit euler to get Xn+1.
  //   3. Compute refined Xn+1 using (1.) and (2.).
  // Note:
  //   1. Use simulateOneStep with modified position and velocity to get Xn+1.
  std::vector<Particles> current;
  for (const auto &p : particles) {
    current.emplace_back(*p);
    p->position() += p->velocity() * deltaTime * 0.5f;
    p->velocity() += p->acceleration() * deltaTime * 0.5f;
  }
  simulateOneStep();
  for (size_t i = 0; i < current.size(); ++i) {
    particles[i]->position() = current[i].position() + particles[i]->velocity() * deltaTime;
    particles[i]->velocity() = current[i].velocity() + particles[i]->acceleration() * deltaTime;
  }
}

void RungeKuttaFourth::integrate(const std::vector<Particles *> &particles,
                                 std::function<void(void)> simulateOneStep) const {
  // TODO: Integrate velocity and acceleration
  //   1. Backup original particles' data.
  //   2. Compute k1, k2, k3, k4
  //   3. Compute refined Xn+1 using (1.) and (2.).
  // Note:
  //   1. Use simulateOneStep with modified position and velocity to get Xn+1.
  std::vector<Particles> current;
  for (const auto &p : particles) current.emplace_back(*p);
  // k1
  std::vector<Particles> k1(current);
  for (size_t i = 0; i < current.size(); ++i) {
    k1[i].position() = particles[i]->velocity() * deltaTime;
    k1[i].velocity() = particles[i]->acceleration() * deltaTime;
    particles[i]->position() = current[i].position() + k1[i].position() * 0.5f;
    particles[i]->velocity() = current[i].velocity() + k1[i].velocity() * 0.5f;
  }
  simulateOneStep();
  // k2
  std::vector<Particles> k2(current);
  for (size_t i = 0; i < current.size(); ++i) {
    k2[i].position() = particles[i]->velocity() * deltaTime;
    k2[i].velocity() = particles[i]->acceleration() * deltaTime;
    particles[i]->position() = current[i].position() + k2[i].position() * 0.5f;
    particles[i]->velocity() = current[i].velocity() + k2[i].velocity() * 0.5f;
  }
  simulateOneStep();
  // k3
  std::vector<Particles> k3(current);
  for (size_t i = 0; i < current.size(); ++i) {
    k3[i].position() = particles[i]->velocity() * deltaTime;
    k3[i].velocity() = particles[i]->acceleration() * deltaTime;
    particles[i]->position() = current[i].position() + k3[i].position();
    particles[i]->velocity() = current[i].velocity() + k3[i].velocity();
  }
  simulateOneStep();
  // k4
  std::vector<Particles> k4(current);
  for (size_t i = 0; i < current.size(); ++i) {
    k4[i].position() = particles[i]->velocity() * deltaTime;
    k4[i].velocity() = particles[i]->acceleration() * deltaTime;

    particles[i]->position() =
        current[i].position() +
        (k1[i].position() + 2.0f * k2[i].position() + 2.0f * k3[i].position() + k4[i].position()) / 6.0f;

    particles[i]->velocity() =
        current[i].velocity() +
        (k1[i].velocity() + 2.0f * k2[i].velocity() + 2.0f * k3[i].velocity() + k4[i].velocity()) / 6.0f;
  }
}
