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
typedef boost:: adjacency_list <boost::vecS, boost::vecS,
										  boost::undirectedS> graph;

// QData stores qlm in its various forms for each particle in the system
// Usually l=6 or 4 but the code will accept any value
//
// The first thing to do is to compute the 2l + 1 dimensional vector
// qlm(i) (m = -l,..,0,..,l) for each particle i in the system.
// This is defined as:
// qlm(i) = 1/N_b(i) \sum_{j=1}^{N_b(i)} Ylm(r_ij)
// Here N_b(i) is the number of neighbours belonging to particle i.
// Two particles are considered neighbours if they are within some cut
// off distance r_cut (This needs to be specified).
// The sum is over all neighbouring particles.
// The functions Ylm are the spherical harmonics, and r_ij is a vector
// from particle i to particle j.
//
// We can use the information qlm(i) in a large number of ways:
//
// i)   To compute the Ql value (e.g. Q6) of a group of particles.
//      (usually all the particles in the system)
//      This is usually called a 'global order parameter'
//      This is defined as:
//      Ql = sqrt( 4pi/(2l + 1) * sum_{m=-l}^{l} |<qlm(i)>|^2 )
//      Where <..> denotes an average over all particles, and |..|
//      denotes absolute value.
// 
// ii)  To compute the ql value of a single particle.
//      This is the same as the formula for Ql, except the < > is
//      removed:
//      ql(i) = sqrt( 4pi/(2l + 1) * sum_{m=-l}^{l} |qlm(i)|^2 )
//      See Lechner and Dellago JCP 129 (2008), equation (3).
//
// iii) To compute the Wl value (e.g. W6) of a group of particles
//      See Steinhart, Nelson, Ronchetti PRB 28 784 (1983)
//      equation (1.5) for the formula.
//
// iv)  To compute the wl value of a single particle
//      See Lechner and Dellago JCP 129 (2008), equation (4).
//      The analogy between i) and ii) is the same as between
//      iii) and iv)
//
// v)   We can normalise the qlm(i).
//      I call the normalised versions \tilde{qlm(i)} or qlmtilde/qlmt
//      in the code below.
//      Then, we can take the dot product (qlmtilde(i) dot qlmtilde(j))
//      between two neighbouring particles i and j
//      Sometimes this dot product is referred to as S_ij in the
//      literature.
//      if S_ij > 0.65 , particles i and j are said to form a 'link'
//      The threshold value 0.65 is arbitrary, others use values
//      between 0.5 and 0.8, the threshold can be specified by the
//      parameter 'linval' below.
//      If particle i has at least 6 links (others use 5-8), it is said
//      to be in a crystalline environment. Again the threshold can
//      be specified by using the parameter 'nlin' in the code.
//      So using this criterion, we can go through every particle in
//      the system, and say whether or not it is in a crystalline
//      (xtal) environment.
//      The size of the largest cluster for which all particles are
//      in a crystalline environment is the familiar Frenkel/ ten Wolde
//      order parameter (which I call N_cl in my papers).
//
// vi)  Next,we can compute the average qlm(i)
//
//

QData::QData(vector<Particle> allps, Box& sbox, int nsur, int nlin, double linval, int lv)
	  : allpars(allps),
		 simbox(sbox),
		 nsurf(nsur),
		 nlinks(nlin),
		 linkval(linval),
		 lval(lv)
{
	  // the constructor is used to compute xps, clusternums, and qcluster.
	  // These are stored in object 

	  vector<Particle>::size_type npar = allpars.size();
  	  vector<int> numneigh(npar,0); // num neighbours for each particle
  	  vector<vector<int> > lneigh; // vector of neighbour particle
	                               // nums for each par
	  lneigh.resize(npar);

	  array2d qlm = qlms(allpars, simbox, numneigh, lneigh, lval);
	  array2d qlmt = qlmtildes(qlm, numneigh, lval);
	  // Lechner dellago eq 6
	  array2d qlmb = qlmbars(qlm, lneigh, lval);

	  // get qls and wls
	  ql = qls(qlm);
	  wl = wls(qlm);	  
	  // lechner dellago eq 5
	  qlbar = qls(qlmb);
	  wlbar = wls(qlmb);

	  // xtal particle nums according to link threshold and min number
	  // of links
	  xps = xtalpars(qlmt, numneigh, lneigh, nsurf, nlinks, linkval, lval);

	  // graph of xtal particles, with each particle a vertex and each
	  // link an edge
	  graph xgraph = getxgraph(allpars, xps, simbox);

	  // indexes into xps of particles that are in the largest cluster
	  cnums = largestcomponent(xgraph);
	  reindex(cnums, xps);

	  // store global q and q of the cluster in object 
	  qcluster = Qpars(qlm, cnums, lval);
	  vector<int> pnums = range(nsurf, allpars.size());
	  qglobal = Qpars(qlm, pnums, lval);
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
