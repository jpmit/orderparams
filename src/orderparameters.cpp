#include <vector>
#include "constants.h"
#include "qlmfunctions.h"
#include "qdata.h"

using std::vector;

// Size of cluster according to Lechner Dellago (LD) method
int csizeld(const vector<int>& ldcnums)
{
	  return ldcnums.size();
}

// Size of cluster according to Ten-Wolde Frenkel (TF) method
int csizetf(const vector<int>& tfcnums)
{
	  return tfcnums.size();
}

// Average Q value of a group of particles, e.g. those in cluster
double qavgroup(const QData& qdata, const vector<int>& pnums)
{
	  return Qpars(qdata.qlm, pnums, qdata.lval);
}


