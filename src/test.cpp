#include <iostream>
#include <vector>
#include "readwrite.h"
#include "box.h"
#include "orderparameter.h"

using std::cout;
using std::endl;
using std::vector;

int main()
{
	  vector<Particle> allpars = readxyz("fcctest.xyz");

	  // for hcp test
	  //Box simbox(22.4492409662,21.3857742696,9.16486424,1.5,true);

	  // for fcc test
	  Box simbox(22.4492409662,21.3857742696,8.24837778,1.5,true);
	  NCluster ncluster(0,6,0.65,6);
//	  allpars.resize(1320);
	  QCluster q6cluster(0,6,0.65,6);
	  QCluster q4cluster(0,6,0.65,4);
	  QGlobal q6(0,6);
	  QGlobal q4(0,4);
	  std::cout << ncluster(allpars,simbox)
					<< " Q6cluster " << q6cluster(allpars,simbox)
		         << " Q4cluster " << q4cluster(allpars,simbox)
			      << " Q6 " << q6(allpars,simbox)
					<< " Q4 " << q4(allpars,simbox)
					<< std::endl;
}
