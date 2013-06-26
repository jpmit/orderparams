#include <iostream>
#include <vector>
#include "boost/multi_array.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include "particle.h"
#include "box.h"
#include "conncomponents.h"
#include "utility.h"

/* gyration.cpp - functions to compute the gyration tensor.
	Be warned, periodic boundary conditions need to be removed before computing.
	Here this is achieved via the function replicate.
*/

using std::vector;

typedef boost::multi_array<double,2> tensor;
typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> graph;

/* Center of mass of particles.
 */

vector<double> cofmass(const vector<Particle>& particles)
{
	  vector<double> cm(3,0.0);
	  
	  vector<Particle>::size_type i;
	  vector<Particle>::size_type npar = particles.size();
	  for (i = 0; i != npar; ++i) {
			 cm[0] += particles[i].pos[0];
			 cm[1] += particles[i].pos[1];
			 cm[2] += particles[i].pos[2];
	  }
	  cm[0] /= npar;
	  cm[1] /= npar;
	  cm[2] /= npar;

	  //std::cout << cm[0] << " " << cm[1] << " " << cm[2] << std::endl;
	  
	  return cm;
}

/* Radius of gyration tensor (not diagonalised).
 */

tensor gytensor(const vector<Particle>& particles)
{
	  vector<double> cmass = cofmass(particles);
	  tensor gyt(boost::extents[3][3]);
	  std::fill(gyt.origin(), gyt.origin() + gyt.size(), 0.0);	  

	  double rcm[3];
	  vector<Particle>::size_type i;
	  vector<Particle>::size_type npar = particles.size();
	  for (i = 0; i != npar; ++i) {
			 // distance of particle from center of mass
			 rcm[0] = particles[i].pos[0] - cmass[0];
			 rcm[1] = particles[i].pos[1] - cmass[1];
			 rcm[2] = particles[i].pos[2] - cmass[2];
			 // contribution to gyration tensor (top half)
			 for (int j = 0; j != 3; ++j) {
					for (int k = j; k != 3; ++k) {
						  gyt[j][k] += rcm[j]*rcm[k];
					}
			 }
	  }
	  // compute remaining entries (since symmetric) and normalise
	  gyt[1][0] = gyt[0][1];
	  gyt[2][0] = gyt[0][2];
	  gyt[2][1] = gyt[1][2];
	  for (int j = 0; j != 3; ++j) {
			 for (int k = 0; k != 3; ++k) {
					gyt[j][k] /= npar;
			 }
	  }
	  return gyt;
}

/* Take positions of particles in largest cluster.  
	Return particles in the largest cluster, but without periodic BCS.
	The trick here is to replicate the system in x and y directions,
	then to find the largest cluster in this large system
	(note the system is assumed not to be periodic in z)
*/

vector<Particle> xtalposnoperiodic(const vector<Particle>& cpars, const Box& simbox)
{
	  const vector<Particle>::size_type ncl = cpars.size();
	  const vector<Particle> repcpars = replicate(cpars, simbox);

	  // hack to ignore periodic bcs
	  Box bigbox = simbox;
	  bigbox.lboxx *= 20;
	  bigbox.lboxy *= 20;

	  graph xgraph = getxgraph(repcpars, range(0, repcpars.size()), bigbox);
	  vector<int> cluspars = largestcomponent(xgraph);

	  // create vector of cluster positions
	  vector<Particle> ret;
	  ret.resize(ncl);
	  for (vector<Particle>::size_type i = 0; i != ncl; ++i) {
			 ret[i] = repcpars[cluspars[i]];
	  }

	  return ret;
}

/* Replicate particles in x and y directions.
 */

vector<Particle> replicate(const vector<Particle>& pars, const Box& simbox)
{
	  const vector<Particle>::size_type npar = pars.size();
	  vector<Particle> newpars;
	  newpars.resize(9 * npar);

	  for (vector<Particle>::size_type i = 0; i != npar; ++i) {
			 for (int j = 0; j != 9; ++j)
					newpars[i + j*npar] = pars[i];
			 // top
			 newpars[i + npar].pos[1] += simbox.lboxy;
			 // top right
			 newpars[i + 2*npar].pos[0] += simbox.lboxx;
			 newpars[i + 2*npar].pos[1] += simbox.lboxy;			 
			 // right
			 newpars[i + 3*npar].pos[0] += simbox.lboxx;
			 // bottom right
			 newpars[i + 4*npar].pos[0] += simbox.lboxx;
			 newpars[i + 4*npar].pos[1] -= simbox.lboxy;			 			 
			 // bottom
			 newpars[i + 5*npar].pos[1] -= simbox.lboxy;			 			 
			 // bottom left
			 newpars[i + 6*npar].pos[0] -= simbox.lboxx;
			 newpars[i + 6*npar].pos[1] -= simbox.lboxy;			 			 
			 // left
			 newpars[i + 7*npar].pos[0] -= simbox.lboxx;
			 // top left
			 newpars[i + 8*npar].pos[0] -= simbox.lboxx;
			 newpars[i + 8*npar].pos[1] += simbox.lboxy;
	  }

	  return newpars;
}

/* uncomment and compile with conncomponents.cpp to test gyration tensor
int main()
{

	  Particle p1,p2;
	  p1.pos[0] = 0;
	  p1.pos[1] = 0;
	  p1.pos[2] = 0;
	  p2.pos[0] = 1;
	  p2.pos[1] = 0;
	  p2.pos[2] = 0;	  
			 
	  vector<Particle> ps;
	  ps.push_back(p1);
	  ps.push_back(p2);

	  tensor g = gytensor(ps);

	  for (int i = 0; i != 3; ++i) {
			 for (int j = 0; j != 3; ++j) {
					std::cout << g[i][j] << " ";
			 }
			 std::cout << std::endl;
	  }
}
*/
