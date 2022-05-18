#include "shape.h"

#include <Eigen/Dense>
#include "configs.h"
#include "integrator.h"

using Eigen::Matrix4f;

// Gravitational constant.
const Eigen::Vector4f GRAVITY(0, -9.8f, 0, 0);

Shape::Shape(int size, float mass_) noexcept :
    _particles(size, mass_), modelMatrix(Matrix4f::Identity()), normalMatrix(Matrix4f::Identity()) {}

void Shape::setModelMatrix(const Eigen::Ref<const Eigen::Matrix4f>& _modelMatrix) {
  modelMatrix = _modelMatrix;
  normalMatrix.topLeftCorner<3, 3>() = _modelMatrix.topLeftCorner<3, 3>().inverse().transpose();
}

void Shape::computeExternalForce(double delta) {
  /* comment by a
  for (size_t i = 0; i < _particles.mass().size(); ++i) {
    if (_particles.mass(i) == 0.0f) {
      _particles.acceleration(i).setZero();
    } else {
      _particles.acceleration(i) = Eigen::Vector4f(0, -9.8f, 0, 0);
      _particles.acceleration(i) -= _particles.velocity(i) * viscousCoef * _particles.inverseMass(i);
    }
  }*/
  for (size_t i = 0; i < _particles.mass().size(); ++i) {
    _particles.velocity(i) += (GRAVITY * delta);
  }
}

void Shape::integrate_PBF(double delta) { computeExternalForce(delta); }
