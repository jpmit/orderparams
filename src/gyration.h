#ifndef GYRATION_H
#define GYRATION_H

#include <vector>
#include "boost/multi_array.hpp"
#include "diagonalize.h"
#include "particle.h"
#include "box.h"

boost::multi_array<double,2> gytensor(const std::vector<Particle>&);
std::vector<Particle> replicate(const std::vector<Particle>&);

#endif
