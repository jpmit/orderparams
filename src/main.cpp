#include <iostream>
#include <vector>
#include <map>
#include <iomanip>
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
	  if (argc != 2 and argc != 3) {
			 cout << "Syntax: " << argv[0] << " filename.xyz" << endl
					<< "Or:     " << argv[0] << " filename.xyz paramfile"
					<< endl;
			 return 1;
	  }

	  // read particle positions from xyz file passed on command line
	  vector<Particle> allpars = readxyz(argv[1]);

	  // read nsurf, lboxx, lboxy, lboxz, zperiodic from params.out file
	  // if we didn't give a name for the parameter file, default to params.out
	  string pfile;
	  if (argc == 2) {
			 pfile = "params.out";
	  }
	  else {
			 pfile = argv[2];
	  }
					
	  map<string, string> params = readparams(pfile);
	  int nsurf = atoi(params["nparsurf"].c_str());
	  double lboxx = atof(params["lboxx"].c_str());
	  double lboxy = atof(params["lboxy"].c_str());
	  double lboxz = atof(params["lboxz"].c_str());
	  map<string, bool> bmap;
	  bmap["True"] = true;
	  bmap["False"] = false;	 
	  bool zperiodic = bmap[params["zperiodic"]];
//	  cout << lboxx << " " << lboxy << " " << lboxz << " " << nsurf << " " << zperiodic
//			 << endl;
	  
	  // arguments are lboxx, lboxy, lboxz, neighboursep, periodic?
	  Box simbox(lboxx,lboxy,lboxz,1.35,zperiodic);

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
	  QData q4data(allpars,simbox,nsurf,nlin,linval,4);

	  // print q4, q4bar, w4, w4bar, q6, q6bar, w6, w6bar
	  // for all particles
	  cout << "# q4 q4bar w4 w4bar q6 q6bar w6 w6bar" << endl;
	  cout << std::setprecision(6) << std::fixed;
	  for (int i = 0; i != allpars.size(); ++i) {
			 cout	<< q4data.ql[i] << " " << q4data.qlbar[i] << " "
					<< q4data.wl[i] << " " << q4data.wlbar[i] << " "
					<< q6data.ql[i] << " " << q6data.qlbar[i] << " "
					<< q6data.wl[i] << " " << q6data.wlbar[i] << endl;
	  }
	  
	  //cout << "filename: " << argv[1] << std::endl;
	  //cout << "Nclusterq6|Q6cluster|Shapeq6" << std::endl;
	  //cout << q6data.getNCluster() << " "
	  //	 << q6data.getQCluster() << " "
	  //		 << q6data.getClusterShape() << std::endl;
}
