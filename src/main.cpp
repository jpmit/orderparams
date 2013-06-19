#include <iostream>
#include <iomanip>
#include <string>
#include "particlesystem.h"
#include "orderparameter.h"
#include "gyration.h"
#include "qdata.h"

using std::cout;
using std::endl;
using std::string;

// Compute the values of all order parameters for a given configuration

int main(int argc, char* argv[])
{
	  if (argc != 2) {
			 cout << "Syntax: " << argv[0] << " paramfile" << endl;
			 return 1;
	  }

	  // get name of input file
	  string pfile = argv[1];
	  
	  // create the particle system from name of input file
	  // input file must contain the following fields:
	  // filename   - name of xyz file to read positions from
	  // lboxx      - x dimension of simulation box
	  // lboxy      - y dimension     ""
	  // lboxz      - z dimension     ""
	  // nsep       - neighbour separation in units of sigma
	  // zperiodic  - box periodic or not (either "True" or "False")
	  // nparsurf   - number of surface particles
	  // q6link     - threshold for Sij to be considered a link
	  // q6numlinks - number of links a particle needs to be xtal
	  ParticleSystem psystem(pfile);

	  // create the cluster data
	  // warning: at the moment the number of links, and the
	  // threshold value for a link is the same for both l=4 and l=6
	  // (psystem.linval and psystem.nlinks respectively)
	  QData q6data(psystem.allpars,psystem.simbox,psystem.nsurf,
						psystem.nlinks,psystem.linval,6);
	  QData q4data(psystem.allpars,psystem.simbox,psystem.nsurf,
						psystem.nlinks,psystem.linval,4);

	  // print q4, q4bar, w4, w4bar, q6, q6bar, w6, w6bar
	  // for all particles
	  // and find which are hcp, fcc etc.
	  int nbcc, nhcp, nfcc, nliq;
	  nbcc = nhcp = nfcc = nliq = 0;
	  
	  cout << "# q4 q4bar w4 w4bar q6 q6bar w6 w6bar" << endl;
	  cout << std::setprecision(6) << std::fixed;
	  for (int i = 0; i != psystem.allpars.size(); ++i) {
			 cout	<< q4data.ql[i] << " " << q4data.qlbar[i] << " "
					<< q4data.wl[i] << " " << q4data.wlbar[i] << " "
					<< q6data.ql[i] << " " << q6data.qlbar[i] << " "
					<< q6data.wl[i] << " " << q6data.wlbar[i] << endl;

			 // based on qs and ws, work out if fcc, hcp, bcc, or other
			 // need to optimize these thresholds later!
			 if (i >= psystem.nsurf) {
					if (q6data.qlbar[i] > 0.3) {
						  // the particle is a solid particle
						  if (q6data.wlbar[i] > 0.0) {
								 // the particle is BCC
								 nbcc += 1;
//								 allpars[i].symbol = 'Y';
						  }
						  else {
								 if (q4data.wlbar[i] > 0.0) {
										// the particle is HCP
										nhcp += 1;
//										allpars[i].symbol = 'P';
								 }
								 else {
										// the particle is FCC
										nfcc += 1;
//										allpars[i].symbol = 'S';
								 }
						  }
					}
					else {
						  // the particle is liquid
						  nliq += 1;
//						  allpars[i].symbol = 'N';
					}
			 }
	  }

	  // print nhcp, nfcc, nbcc, nliq
	  cout << "nbcc nhcp nfcc nliq" << endl << nbcc << " " << nhcp
			 << " " << nfcc << " " << nliq << endl;

	  // save file
//	  if (printsymbols) {
//			 writexyz(allpars, savefile);
//	  }
}

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
