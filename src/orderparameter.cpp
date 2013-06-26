#include <vector>
#include <cmath>
#include <cstring>
#include <complex>
#include <iostream>
#include "boost/multi_array.hpp"
#include "orderparameter.h"
#include "box.h"
#include "particle.h"
#include "constants.h"
#include "opfunctions.h"
#include "conncomponents.h"
#include "qlmfunctions.h"
#include "gyration.h"
#include "utility.h"
#include "readwrite.h"

using std::vector;
using std::complex;
using std::norm;

typedef boost::multi_array<double,2> tensor;
typedef boost::multi_array<complex<double>,2> array2d;
typedef boost:: adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> graph;

/* Size of largest cluster in the system identified by local bond order parameter.
 */

double NCluster::operator() (const vector<Particle>& particles, const Box& simbox)
{
	  vector<Particle>::size_type npar = particles.size();
	  
  	  vector<int> numneigh(npar,0); // number of neighbours for each particle
  	  vector<vector<int> > lneigh; // vector of neighbour particle numbers for each par
	  lneigh.resize(npar);

	  // get \tilde{qlm}[i] for every particle i
	  array2d qlmt = qlms(particles, simbox, numneigh, lneigh, lval);
	  qlmtildes(qlmt, numneigh, lval);

	  // get crystal particle numbers according to link threshold and min number of links
	  vector<int> xps = xtalpars(qlmt, numneigh, lneigh, nsurf, nlinks, linkval, lval);

	  // get the graph of xtal particles, with each particle a vertice and each link an edge
	  graph xgraph = getxgraph(particles, xps, simbox);

	  // use the graph to compute largest cluster
	  return bopxbulk(xgraph);
}

/* Q value of largest cluster in the system (can be Q6, Q4 etc. depending on lval).
 */

double QCluster::operator() (const vector<Particle>& particles, const Box& simbox)
{
	  vector<Particle>::size_type npar = particles.size();
	  
  	  vector<int> numneigh(npar,0); // number of neighbours for each particle
  	  vector<vector<int> > lneigh; // vector of neighbour particle numbers for each par
	  lneigh.resize(npar);

	  // get both \bar{qlm}[i] and \tilde{qlm}[i] for every particle i
	  array2d qlmb = qlms(particles, simbox, numneigh, lneigh, lval);
	  array2d qlmt = qlmb;
	  qlmtildes(qlmt, numneigh, lval);

	  // get crystal particle numbers according to link threshold and min number of links
	  vector<int> xps = xtalpars(qlmt, numneigh, lneigh, nsurf, nlinks, linkval, lval);	  

	  // get the graph of xtal particles, with each particle a vertice and each link an edge
	  graph xgraph = getxgraph(particles, xps, simbox);

	  // figure out which particles are in the largest cluster
	  vector<int> clusternums = largestcomponent(xgraph);
	  reindex(clusternums, xps);

	  return Qpars(qlmb, clusternums, lval);
}

/* Gyration tensor of largest cluster (exact value to return currently in progress).
 */

double GTensor::operator() (const vector<Particle>& particles, const Box& simbox)
{
	  vector<Particle>::size_type npar = particles.size();
	  
  	  vector<int> numneigh(npar,0); // number of neighbours for each particle
  	  vector<vector<int> > lneigh; // vector of neighbour particle numbers for each par
	  lneigh.resize(npar);

	  // get \tilde{qlm}[i] for every particle i
	  array2d qlmt = qlms(particles, simbox, numneigh, lneigh, lval);
	  qlmtildes(qlmt, numneigh, lval);

	  // get crystal particle numbers according to link threshold and min number of links
	  vector<int> xps = xtalpars(qlmt, numneigh, lneigh, nsurf, nlinks, linkval, lval);	  	  

	  // get the graph of xtal particles, with each particle a vertice and each link an edge
	  graph xgraph = getxgraph(particles, xps, simbox);

	  // figure out which particles are in the largest cluster
	  vector<int> clusternums = largestcomponent(xgraph);

	  // build up vector of particles which are largest cluster only
	  vector<Particle>::size_type ncl = clusternums.size();
	  vector<Particle> clusterpars(ncl);
	  vector<Particle>::size_type i;
	  for (i = 0; i != ncl; ++i) {
			 // this is a bit confusing:
			 // clusternums[i] gives an index into xps, which
			 // in turn gives an index into particles
			 clusterpars[i].pos[0] = particles[xps[clusternums[i]]].pos[0];
			 clusterpars[i].pos[1] = particles[xps[clusternums[i]]].pos[1];
			 clusterpars[i].pos[2] = particles[xps[clusternums[i]]].pos[2];
			 clusterpars[i].symbol = 'S';
	  }

	  // take away periodic bcs
	  vector<Particle> newclusterpars = xtalposnoperiodic(clusterpars, simbox);
	  writexyz(newclusterpars, "testrep.xyz");

	  // get radius of gyration tensor
	  tensor gyt = gytensor(newclusterpars);
	  std::cout << "Gyration tensor:" << std::endl;
	  for (int j = 0; j != 3; ++j) {
			 std::cout << gyt[j][0] << " " << gyt[j][1] << " "
						  << gyt[j][2] << std::endl;
	  }

	  // diagonalise complete tensor
	  double fullg[] = {gyt[0][0], gyt[0][1], gyt[0][2],
							  gyt[1][0], gyt[1][1], gyt[1][2],
							  gyt[2][0], gyt[2][1], gyt[2][2]};
	  double fullres[9];
	  double fulleig[3];
	  diagonalize(fullg, 3, fullres, fulleig);
	  //printeig(fullres, fulleig, 3);

	  // diagonalise 2*2 tensor (x-y)
	  double g[] = {gyt[0][0], gyt[0][1],
						 gyt[1][0], gyt[1][1]};
	  double res[4];
	  double eig[2];
	  diagonalize(g, 2, res, eig);
	  //printeig(res, eig, 2);
	  
	  return eig[0]*eig[0] + eig[1]*eig[1];
}
