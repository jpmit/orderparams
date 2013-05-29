#ifndef CONNCOMPONENTS_H
#define CONNCOMPONENTS_H

#include <boost/graph/adjacency_list.hpp>
#include "particle.h"
#include "box.h"

boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS> getxgraph(const std::vector<Particle>&,
																									 const std::vector<int>&, const Box&);
int bopxbulk(const boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS>&);
std::vector<int> largestcomponent(const boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS>&);

#endif
