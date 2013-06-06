#ifndef QDATA_H
#define QDATA_H

#include <vector>
#include "particle.h"
#include "box.h"
#include "boost/multi_array.hpp"

// QData is a class to store data used for the Steindhardt bond order
// parameters. The class exists since we don't want to recompute
// spherical harmonics many times.

class QData
{
public:
	  QData(std::vector<Particle> allps, Box& sbox,
			  int nsur, int nlin, double linval, int lv);

	  // member functions for returning OPs 
	  
	  int getNCluster();
	  double getQCluster();
	  double getQGlobal();
	  double getClusterShape();

	  // ql, \bar{ql}, wl, \bar{wl} for each particle i
	  std::vector<double> ql;
	  std::vector<double> qlbar;
	  std::vector<double> wl;
	  std::vector<double> wlbar;
	  
private:
	  double getqcluster(const boost::multi_array<std::complex<double>,2>&,
								const std::vector<int>&,
								const std::vector<int>&);
	  std::vector<int> getxpars(const boost::multi_array<std::complex<double>,2>&,
										 const std::vector<int>&,
										 const std::vector<std::vector<int> >&);

	  /* we store all parameters passed in constructor */
	  
	  std::vector<Particle> allpars; // particle positions
	  Box simbox;
	  int nsurf;      // num surface particles (don't count as part of the cluster)
	  int nlinks;     // min num xtal links required for particle to be xtal
	  double linkval; // min value for link to be considered xtal
	  int lval;       // spherical harmonic num, usually 4 or 6

	  /* parameters stored to return order parameters */
	  
	  std::vector<int> xps; // indexes of particles identified as xtal
	  std::vector<int> cnums; // indexes of pars in largest cluster
	  double qcluster; // Q value of cluster
	  double qglobal; // Q value of entire system (excluding surface particles)


};

#endif
