#ifndef QLMFUNCTIONS_H
#define QLMFUNCTIONS_H

#include <vector>
#include <complex>
#include "particle.h"
#include "box.h"
#include "boost/multi_array.hpp"

void qlmtildes(boost::multi_array<std::complex<double>,2>&,
					const std::vector<int>&, const int);
boost::multi_array<std::complex<double>,2 > qlms(const std::vector<Particle>&,
															 const Box& simbox,
															 std::vector<int>& numneigh,
															 std::vector<std::vector<int> >&,
															 const int);
std::vector<int> xtalpars(const boost::multi_array<std::complex<double>,2 >&,
								  const std::vector<int>&,
								  const std::vector<std::vector<int> >&,
								  const int, const int,
								  const double, const int);
double Qpars(const boost::multi_array<std::complex<double>,2 >&,
				 const std::vector<int>&, const int);


#endif
