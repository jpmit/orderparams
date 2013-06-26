#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include "particlesystem.h"
#include "orderparameters.h"
#include "qdata.h"
#include "constants.h"
#include "utility.h"

using std::cout;
using std::endl;
using std::string;

// Compute the values of all order parameters for a given configuration
// Currently there are 45 different order parameters output, but a lot
// of these are 'duplicated' because there are two different clusters
// (see below).

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

	  // compute the qlm data
	  // warning: at the moment the number of links, and the
	  // threshold value for a link is the same for both l=4 and l=6
	  // (psystem.linval and psystem.nlinks respectively)
	  QData q6data(psystem, 6);
	  QData q4data(psystem, 4);

	  
	  // from q6data and q4 data, classify each particle as bcc, hcp etc.
	  // using Lechner Dellago approach.
	  vector<LDCLASS> ldclass = classifyparticlesld(psystem, q4data,
																	q6data);

	  // from q6 data only, classify each particle as either crystalline
	  // or liquid, using TenWolde Frenkel approach
	  vector<TFCLASS> tfclass = classifyparticlestf(psystem, q6data);

	  // indices into particle vector (psystem.allpars) of those
	  // particles in the ten-Wolde Frenkel largest cluster and
	  // those in the Lechner Dellago cluster.
	  vector<int> tfcnums = largestclustertf(psystem, tfclass);
	  vector<int> ldcnums = largestclusterld(psystem, ldclass);

	  // compute each order parameter in turn and print to stdout.
	  // See orderparams.cpp for these functions.

	  //////////////////////////////////////////////////////////////////
	  // The following order parameters are associated in some way with
	  // properties of the largest cluster.  There are two approaches to
	  // determining this cluster, which I call Lecher Dellage (LD) and
	  // ten-Wolde Frenkel (TF), and thus two different clusters.  All
	  // of the OPs are computed for both clusters.
	  /////////////////////////////////////////////////////////////////
	  
	  // Size of cluster by LD method
	  cout << "N_ld " << csizeld(ldcnums) << endl;

	  // Size of cluster by TF method
	  cout << "N_tf " << csizetf(tfcnums) << endl;
	  
	  // fraction of bcc pars in LD cluster
	  cout << "n_bccLD " << parfrac(ldclass, ldcnums, BCC) << endl;
	  
	  // fraction of bcc pars in TF cluster	  
	  cout << "n_bccTF " << parfrac(ldclass, tfcnums, BCC) << endl;
	  
	  // fraction of fcc pars in LD cluster
	  cout << "n_fccLD " << parfrac(ldclass, ldcnums, FCC) << endl;
	  
	  // fraction of fcc pars in TF cluster	  
	  cout << "n_fccTF " << parfrac(ldclass, tfcnums, FCC) << endl;
	  
	  // fraction of hcp pars in LD cluster
	  cout << "n_hcpLD " << parfrac(ldclass, ldcnums, HCP) << endl;
	  
	  // fraction of hcp pars in TF cluster	  
	  cout << "n_hcpTF " << parfrac(ldclass, tfcnums, HCP) << endl;
	  
	  // fraction of icos pars in LD cluster
	  cout << "n_icosLD " << parfrac(ldclass, ldcnums, ICOS) << endl;
	  
	  // // fraction of icos pars in TF cluster	  
	  cout << "n_icosTF " << parfrac(ldclass, tfcnums, ICOS) << endl;
	 
	  // average Q6 of LD cluster
	  cout << "Q6clusLD " << qavgroup(q6data, ldcnums) << endl;

	  // average Q6 of TF cluster
	  cout << "Q6clusTF " << qavgroup(q6data, tfcnums) << endl;
	  	 
	  // average Q4 of LD cluster
	  cout << "Q4clusLD " << qavgroup(q4data, ldcnums) << endl;

	  // average Q4 of TF cluster
	  cout << "Q4clusTF " << qavgroup(q4data, tfcnums) << endl;

	  // // number of liquid like particles with at least one neighbour in
	  // // LD cluster.  Note that we could pass either q6data.lneigh or
	  // // q4data.lneigh, since these are identical 
	  // cout << "N_s " << nliquidneighld(ldclass, ldcnums, q6data.lneigh);

	  // // same as above but for TF cluster
	  // cout << "N_s " << nliquidneightf(tfclass, tfcnums, q6data.lneigh);	  

	  // // total number of connections for all liquid-like particles with
	  // // at least one neighbour in cluster for LD cluster
	  // cout << "N_l " << nconnectionsliqld(ldcnums, q6data, ldclass);

	  // // same as above but for TF cluster
	  // cout << "N_l " << nconnectionsliqtf(tfcnums, q6data, tfclass);

	  // average q6 of liquid-like particles with at least one neighbour
	  // in cluster for LD cluster

	  // same as above but for TF cluster

	  // average q4 of liquid-like particles with at least one neighbour
	  // in cluster for LD cluster

	  // same as above but for LD cluster

	  // potential energy (per particle) for LD cluster

	  // potential energy (per particle) for TF cluster

	  // largest eigenvalue of gyration tensor for LD cluster

	  // largest eigenvalue of gyration tensor for TF cluster

	  // middle eigenvalue of gyration tensor for LD cluster

	  // middle eigenvalue of gyration tensor for TF cluster

	  // smallest eigenvalue of gyration tensor for LD cluster

	  // smallest eigenvalue of gyration tensor for TF cluster

	  // (3,3) element of non-diagonalized gyration tensor for LD
	  // cluster

	  // (3,3) element of non-diagonalized gyration tensor for TF
	  // cluster

	  // largest eigenvalue of top-diagonalised gyration tensor for LD
	  // cluster

	  // largest eigenvalue of top-diagonalised gyration tensor for TF
	  // cluster

	  // smallest eigenvalue of top-diagonalised gyration tensor for LD
	  // cluster

	  // smallest eigenvalue of top-diagonalised gyration tensor for TF
	  // cluster

     //////////////////////////////////////////////////////////////////
	  // These order parameter are 'global' i.e. for the entire system
	  // (that is, no mention of a cluster of any kind!).
	  // Note that we exclude surface particles from the calculations.
	  //////////////////////////////////////////////////////////////////

	  // indexes of all particles (minus surface particles)
	  vector<int> pindices = range(psystem.nsurf, psystem.allpars.size());
	  
	  // fraction of bcc particles in entire system
	  cout << "s_bcc " << parfrac(ldclass, pindices, BCC) << endl;
	  
	  // fraction of fcc particles in entire system
	  cout << "s_fcc " << parfrac(ldclass, pindices, FCC) << endl;
	  
	  // fraction of hcp particles in entire system
	  cout << "s_hcp " << parfrac(ldclass, pindices, HCP) << endl;

	  // fraction of icosahedral particles in entire system
	  cout << "s_icos " << parfrac(ldclass, pindices, ICOS) << endl;
	  
	  // potential energy per (moving) particle of entire system

	  // Average q6 of all particles in system
	  cout << "Q6 " << qavgroup(q6data, pindices) << endl;

	  // average q4 of all particles in system
	  cout << "Q4 " << qavgroup(q4data, pindices) << endl;	  
	  
}
