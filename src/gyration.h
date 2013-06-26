#ifndef GYRATION_H
#define GYRATION_H

#include <vector>
#include "typedefs.h"
#include "particlesystem.h"

tensor getgytensor(const ParticleSystem&, const std::vector<int>&);

#endif
