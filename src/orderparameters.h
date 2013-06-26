#ifndef ORDERPARAMETERS_H
#define ORDERPARAMETERS_H

#include <iostream>
#include "constants.h"
#include "qdata.h"

int csizeld(const std::vector<int>&);
int csizetf(const std::vector<int>&);
double qavgroup(const QData&, const std::vector<int>&);

template <class etype>
double parfrac(const std::vector<etype>& pclass, const std::vector<int>& cnums,
					const etype plabel)
{
	  int num = 0;
	  //std::cout << cnums.size() << std::endl;
	  
	  for (std::vector<int>::size_type i = 0; i != cnums.size(); ++i) {
			 //std::cout << cnums[i] << " " << pclass[i] << std::endl;
			 if (pclass[cnums[i]] == plabel) {
					++num;
			 }
	  }

	  return static_cast<double>(num) / cnums.size();
}


#endif
