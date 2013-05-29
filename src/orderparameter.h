#ifndef ORDERPARAMETER_H
#define ORDERPARAMETER_H

#include <vector>
#include "particle.h"
#include "box.h"
#include "boost/multi_array.hpp"

/* OrderParameter is an abstract base class.
	The virtual operator() ensures that any instance of a derived
	class must implement this method, which is an evaluation of the
	order parameter.
*/

class OrderParameter
{
public:
	  virtual double operator() (const std::vector<Particle>&, const Box& simbox) = 0;
};

/* NCluster is a class for computing the largest crystal cluster.
 */

class NCluster : public OrderParameter
{
public:
     NCluster(int nsur, int nlin, double linval, int lv):
	  nsurf(nsur), nlinks(nlin), linkval(linval), lval(lv){}
	  double operator() (const std::vector<Particle>&, const Box&);
	  std::vector<int> xpars(const boost::multi_array<std::complex<double>,2>&, const std::vector<int>&,
									 const std::vector<std::vector<int> >&);
protected:
	  int nsurf; // number of surface particles (these don't count as part of the cluster)
	  int nlinks; // min number of xtal links required for particle to be xtal
	  double linkval; // value for link to be considered xtal
	  int lval; // e.g. lval = 6 for 6 fold symmetry, 4 for 4 fold
};

/* QCluster computes the Q value of the largest cluster in the system.
 */

class QCluster : public NCluster
{
public:
     QCluster(int nsur, int nlink, double linval, int lv):
	  NCluster(nsur, nlink, linval, lv){}
	  double operator() (const std::vector<Particle>&, const Box&);
};

/* GTensor computes the gyration tensor of the largest cluster in the system (in progress).
 */

class GTensor : public NCluster
{
public:
     GTensor(int nsur, int nlink, double linval, int lv):
	  NCluster(nsur, nlink, linval, lv){}
	  double operator() (const std::vector<Particle>&, const Box&);
};

/* QGlobal computers the Q value of the entire system (ignoring surface particles).
 */

class QGlobal : public QCluster
{
public:
     QGlobal(int nsur, int lv):
	  QCluster(nsur, 0, 0.0, lv){}
};

#endif
