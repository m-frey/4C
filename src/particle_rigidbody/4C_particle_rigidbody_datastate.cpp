/*---------------------------------------------------------------------------*/
/*! \file
\brief rigid body data state container
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_particle_rigidbody_datastate.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
void ParticleRigidBody::RigidBodyDataState::init()
{
  // nothing to do
}

void ParticleRigidBody::RigidBodyDataState::setup()
{
  // nothing to do
}

void ParticleRigidBody::RigidBodyDataState::allocate_stored_states(const int numrigidbodies)
{
  mass_.resize(numrigidbodies, 0.0);
  inertia_.resize(numrigidbodies, std::vector<double>(6, 0.0));
  position_.resize(numrigidbodies, std::vector<double>(3, 0.0));
  velocity_.resize(numrigidbodies, std::vector<double>(3, 0.0));
  acceleration_.resize(numrigidbodies, std::vector<double>(3, 0.0));
  rotation_.resize(numrigidbodies, std::vector<double>(4, 0.0));
  angularvelocity_.resize(numrigidbodies, std::vector<double>(3, 0.0));
  angularacceleration_.resize(numrigidbodies, std::vector<double>(3, 0.0));
  force_.resize(numrigidbodies, std::vector<double>(3, 0.0));
  torque_.resize(numrigidbodies, std::vector<double>(3, 0.0));
}

FOUR_C_NAMESPACE_CLOSE
