#ifndef QLMFUNCTIONS_H
#define QLMFUNCTIONS_H

#include <vector>
#include <complex>
#include "particle.h"
#include "box.h"
#include "boost/multi_array.hpp"

boost::multi_array<std::complex<double>,2 > qlmtildes(const boost::multi_array<std::complex<double>,2>&,
																		const std::vector<int>&, const int);
boost::multi_array<std::complex<double>,2 > qlmbars(const boost::multi_array<std::complex<double>,2>&,
																	 const std::vector<std::vector<int> >&, const int);
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
std::vector<double> qls(const boost::multi_array<std::complex<double>,2 >&);
std::vector<std::complex<double> > wls(const boost::multi_array<std::complex<double>,2 >&);

#endif
