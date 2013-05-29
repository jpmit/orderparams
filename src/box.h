#ifndef BOX_H
#define BOX_H

#include <cmath>
#include <vector>
#include "particle.h"
#include "box.h"

/* Simulation box; the member functions handle the periodic bcs.
 */

class Box
{
public:
     Box(double lx, double ly, double lz, double ns = 1.5, bool pz = false) :
	  lboxx(lx), lboxy(ly), lboxz(lz), nsep(ns), nsepsq(ns*ns), periodicz(pz){}
	  inline void sep(const Particle& p1, const Particle& p2, double* s) const;
	  inline double sepsq(const Particle& p1, const Particle& p2) const;
	  inline bool isneigh(const Particle& p1, const Particle& p2, double& r2) const;
	  inline bool isneigh(double* s, double&r2) const;
	  inline bool posvalid(double* pos) const;
	  inline bool getvalidifnot(double* pos) const;
	  friend std::vector<Particle>
			 initpositions(int npar, const Box& simbox, double rcinit);
	  friend std::vector<Particle>
			 replicate(const std::vector<Particle>&, const Box&);
	  friend std::vector<Particle> xtalposnoperiodic(const std::vector<Particle>&,
																	 const Box&);
private:
	  double lboxx;
	  double lboxy;
	  double lboxz;
	  double nsep;
	  double nsepsq;
	  bool periodicz;
};

/* Is pos = {x, y, z} valid?  If so make pos modulo periodic bcs.
 */

inline bool Box::posvalid(double* pos) const
{
	  // handle z boundary
	  if (periodicz) {
			 if (pos[2] < 0.0) 
					pos[2] = pos[2] + lboxz;
			 else if (pos[2] > lboxz)
					pos[2] = pos[2] - lboxz;
	  }
	  else { // z not periodic
			 if (pos[2] > lboxz || pos[2] < 0.0)
					return false;
	  }

	  // handle x and y boundaries
	  if (pos[0] < 0.0)
			 pos[0] = pos[0] + lboxx;
	  else if (pos[0] > lboxx)
			 pos[0] = pos[0] - lboxx;
	  if (pos[1] < 0.0)
			 pos[1] = pos[1] + lboxy;
	  else if (pos[1] > lboxy)
			 pos[1] = pos[1] - lboxy;

	  return true;
}

/* square of separation between two particles.
 */

inline double Box::sepsq(const Particle& p1, const Particle& p2) const
{
	  double s[3];
	  sep(p1,p2,s);
	  return (s[0]*s[0] + s[1]*s[1] + s[2]*s[2]);
}

/* separation between two particles modulo periodic bcs.
 */

inline void Box::sep(const Particle& p1, const Particle& p2, double* s) const
{
	  double sepx,sepy,sepz;
	  sepx = p1.pos[0] - p2.pos[0];
	  sepy = p1.pos[1] - p2.pos[1];
	  sepz = p1.pos[2] - p2.pos[2];
	  
	  if (sepx > 0.5*lboxx) {
			 sepx = sepx - lboxx;
	  }
	  else if (sepx < -0.5*lboxx) {
			 sepx = sepx + lboxx;
	  }
	  if (sepy > 0.5*lboxy) {
			 sepy = sepy - lboxy;
	  }
	  else if (sepy < -0.5*lboxy) {
			 sepy = sepy + lboxy;
	  }
	  if (periodicz) {					
					if (sepz > 0.5*lboxz) {
						  sepz = sepz - lboxz;
					}
					else if (sepz < -0.5*lboxz) {
						  sepz = sepz + lboxz;
					}
			 }

	  s[0] = sepx;
	  s[1] = sepy;
	  s[2] = sepz;
	  return;
}

/* Are p1 and p2 neighbours?  Also return square separation modulo periodic bcs.
 */

inline bool Box::isneigh(const Particle& p1, const Particle& p2, double& rsq) const
{
	  double s[3];

	  // compute separation between particles (store in s)
	  sep(p1,p2,s);

	  if ((std::abs(s[0]) < nsep) && (std::abs(s[1]) < nsep)
			&& (std::abs(s[2]) < nsep)) {
			 rsq = s[0]*s[0] + s[1]*s[1] + s[2]*s[2];
			 if (rsq < nsepsq)
					return true;
	  }
	  return false;
}

/* Same as above but first argument points to separation array.
 */

inline bool Box::isneigh(double *s, double& rsq) const
{
	  if ((std::abs(s[0]) < nsep) && (std::abs(s[1]) < nsep)
			&& (std::abs(s[2]) < nsep)) {
			 rsq = s[0]*s[0] + s[1]*s[1] + s[2]*s[2];
			 if (rsq < nsepsq)
					return true;
	  }
	  return false;
}

#endif
