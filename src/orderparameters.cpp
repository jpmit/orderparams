#include <vector>
#include "constants.h"

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
