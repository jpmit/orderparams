#include <iostream>
#include <vector>
#include <complex>
#include "boost/multi_array.hpp"
#include "qdata.h"
#include "box.h"
#include "particle.h"
#include "qlmfunctions.h"
#include "constants.h"
#include "conncomponents.h"
#include "orderparameter.h"
#include "utilityfunctions.h"
#include "gyration.h"

using std::vector;
using std::complex;
using std::norm;

typedef boost::multi_array<complex<double>,2> array2d;
typedef boost::multi_array<double,2> tensor;
typedef boost:: adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> graph;

QData::QData(vector<Particle> allps, Box& sbox, int nsur, int nlin, double linval, int lv)
	  : allpars(allps),
		 simbox(sbox),
		 nsurf(nsur),
		 nlinks(nlin),
		 linkval(linval),
		 lval(lv)
{
	  /* the constructor is used to compute xps, clusternums, and qcluster.
		  These are stored in object */

	  vector<Particle>::size_type npar = allpars.size();
  	  vector<int> numneigh(npar,0); // num neighbours for each particle
  	  vector<vector<int> > lneigh; // vector of neighbour particle nums for each par
	  lneigh.resize(npar);

	  /* get both \bar{qlm}[i] and \tilde{qlm}[i] for every particle i.
		  \bar{qlm}[i] are used to compute big Q
		  \tidle{qlm}[i] are the normalised vectors used to compute links
		  for each particle, and hence to decide which particles are xtal.
		  The xtal particles are then subjected to a cluster analysis, which
		  gives the size of the largest cluster.
	  */
		  
	  array2d qlmb = qlms(allpars, simbox, numneigh, lneigh, lval);
	  array2d qlmt = qlmb;
	  qlmbars(qlmb, numneigh, lval);	  
	  qlmtildes(qlmt, numneigh, lval);

	  // xtal particle nums according to link threshold and min number of links
	  xps = xtalpars(qlmt, numneigh, lneigh, nsurf, nlinks, linkval, lval);

	  // graph of xtal particles, with each particle a vertex and each link an edge
	  graph xgraph = getxgraph(allpars, xps, simbox);

	  // indexes into xps of particles that are in the largest cluster
	  cnums = largestcomponent(xgraph);
	  reindex(cnums, xps);

	  // store global q and q of the cluster in object 
	  qcluster = Qpars(qlmb, cnums, lval);
	  vector<int> pnums = range(nsurf, allpars.size());
	  qglobal = Qpars(qlmb, pnums, lval);
}

/* Number of particles in largest cluster. */

int QData::getNCluster()
{
	  return cnums.size();
}

/* Q value of largest cluster. */

double QData::getQCluster()
{
	  return qcluster;
}

/* Global Q value i.e. for whole system. */

double QData::getQGlobal()
{
	  return qglobal;
}

/* 'Shape' of largest cluster. */

double QData::getClusterShape()
{
	  /* get radius of gyration tensor of largest cluster. */
	  vector<Particle> cpars;
	  int ncl = cnums.size();
	  cpars.resize(ncl);
	  for (vector<Particle>::size_type i = 0;
			 i != ncl; ++i) {
			 cpars[i] = allpars[cnums[i]];
	  }
	  cpars = xtalposnoperiodic(cpars, simbox);
	  tensor gyt = gytensor(cpars);

	  /* diagonalise (x-y) part of tensor */
	  double g[] = {gyt[0][0], gyt[0][1],
						 gyt[1][0], gyt[1][1]};
	  double res[4];
	  double eig[2];
	  diagonalize(g, 2, res, eig);

	  /*
	  for (int j = 0; j != 3; ++j) {
			 std::cout << gyt[j][0] << " " << gyt[j][1] << " "
						  << gyt[j][2] << std::endl;
	  }
	  */
	  
	  return eig[0]*eig[0] + eig[1]*eig[1];
}
