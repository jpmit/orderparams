#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include "particle.h"
#include "box.h"

using std::vector;
using std::cout;
using std::endl;
typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> graph;

/* Create graph with crystal pars as nodes and vertices between neighbours.
 */

graph getxgraph(const vector<Particle>& particles, const vector<int>& xpars, const Box& simbox)
{
	  graph G; // graph containing crystal particles as nodes
	  vector<int>::size_type nxtal = xpars.size();
	  vector<int>::size_type i,j;
	  double sep;
	  
	  for (i = 0; i != nxtal; ++i) 
			 for (j = i + 1; j != nxtal; ++j)
					if (simbox.isneigh(particles[xpars[i]],particles[xpars[j]],sep))
						  add_edge(i, j, G);

	  return G;
}

/* Bopxbulk is the number of particles in largest cluster (largest connected component).
 */

int bopxbulk(const graph& g)
{
	  // compute number of connected components
	  vector<int> component(num_vertices(g));
	  int num = connected_components(g, &component[0]);

	  // sort each particle into one of the connected components
	  vector<int> ncomp(num,0);
	  for (vector<int>::size_type i = 0; i != component.size(); ++i)
			 ++ncomp[component[i]];

	  // max element returns iterator so dereference to give integer
	  return *max_element(ncomp.begin(),ncomp.end());
}

/* Return vector of ints containing nodes (particle nums) of largest connected component.
 */

vector<int> largestcomponent(const graph& G)
{
	  // compute number of connected components
	  vector<int> component(num_vertices(G));
	  int num = connected_components(G, &component[0]);

	  // sort each particle into one of the connected components
	  vector<int> ncomp(num,0);
	  for (vector<int>::size_type i = 0; i != component.size(); ++i)
			 ++ncomp[component[i]];

	  int maxcomp = distance(ncomp.begin(),max_element(ncomp.begin(),ncomp.end()));
	  
	  vector<int> ret;
	  for (vector<int>::size_type i = 0; i != component.size(); ++i) {
			 if (component[i] == maxcomp)
					ret.push_back(i);
	  }

	  return ret;
}
