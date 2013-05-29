#include <iostream>
#include <vector>
#include <map>
#include <string>
#include "readwrite.h"
#include "box.h"
#include "orderparameter.h"
#include "gyration.h"
#include "qdata.h"

using std::cout;
using std::endl;
using std::vector;
using std::map;
using std::string;

/* Compute the values of all order parameters for a given configuration
 */

int main(int argc, char* argv[])
{
	  if (argc != 2) {
			 cout << "Syntax: " << argv[0] << " filename.xyz" << endl;
			 return 1;
	  }

	  // read particle positions from xyz file passed on command line
	  vector<Particle> allpars = readxyz(argv[1]);

	  // read nsurf, lboxx, lboxy, lboxz, zperiodic from params.out file
	  map<string, string> params = readparams("params.out");
	  int nsurf = atoi(params["nparsurf"].c_str());
	  double lboxx = atof(params["lboxx"].c_str());
	  double lboxy = atof(params["lboxy"].c_str());
	  double lboxz = atof(params["lboxz"].c_str());
	  map<string, bool> bmap;
	  bmap["True"] = true;
	  bmap["False"] = false;	 
	  bool zperiodic = bmap[params["zperiodic"]];
	  cout << lboxx << " " << lboxy << " " << lboxz << " " << nsurf << " " << zperiodic
			 << endl;
	  
	  // arguments are lboxx, lboxy, lboxz, neighboursep, periodic?
	  Box simbox(22.4492409662,21.3857742696,14.2055395823,1.5,false);

	  // test OPS using both interfaces
	  int nlin = 6;
	  double linval = 0.65;

	  // 1st interface
/*	  
	  NCluster ncluster(nsurf,nlin,linval,6);
	  QCluster q6cluster(nsurf,nlin,linval,6);
	  QCluster q4cluster(nsurf,nlin,linval,4);
	  QGlobal q6(nsurf,6);
	  QGlobal q4(nsurf,4);
	  GTensor gyt(nsurf,nlin,linval,6);
	  cout << "Ncluster " << ncluster(allpars,simbox)
			 << " Q6cluster " << q6cluster(allpars,simbox)
			 << " Q4cluster " << q4cluster(allpars,simbox)
			 << " Q6 " << q6(allpars,simbox)
			 << " Q4 " << q4(allpars,simbox)
			 << " Shape " << gyt(allpars,simbox)
			 << std::endl;
*/
	  
	  // second interface
	  QData q6data(allpars,simbox,nsurf,nlin,linval,6);
	  //QData q4data(allpars,simbox,nsurf,nlin,linval,4);
	  cout << "filename: " << argv[1] << std::endl;
	  cout << "Nclusterq6|Q6cluster|Shapeq6" << std::endl;
	  cout << q6data.getNCluster() << " "
			 << q6data.getQCluster() << " "
			 << q6data.getClusterShape() << std::endl;
}
