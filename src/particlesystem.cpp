#include <map>
#include <string>
#include <cstdlib>
#include "particlesystem.h"
#include "readwrite.h"
#include "box.h"

using std::map;
using std::string;

ParticleSystem::ParticleSystem(string pfile)
{

	  // read parameters from specified file
	  map<string, string> params = readparams(pfile);

	  // particle positions from xyz file
	  allpars = readxyz(params["filename"].c_str());

	  // get box parameters and use to create box
	  double lboxx = atof(params["lboxx"].c_str());
	  double lboxy = atof(params["lboxy"].c_str());
	  double lboxz = atof(params["lboxz"].c_str());
	  nsep = atof(params["stillsep"].c_str());
	  map<string, bool> bmap;
	  bmap["True"] = true;
	  bmap["False"] = false;	 
	  bool zperiodic = bmap[params["zperiodic"]];
	  simbox = Box(lboxx,lboxy,lboxz,nsep,zperiodic);

	  // number of surface particles
	  nsurf = atoi(params["nparsurf"].c_str());

	  // warning: at the moment these are used for both l=4 and l=6
	  linval = atof(params["q6link"].c_str());
	  nlinks = atoi(params["q6numlinks"].c_str());
	  
}

