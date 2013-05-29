#include "boost/multi_array.hpp"
#include <complex>
#include "constants.h"
#include "particle.h"
#include "box.h"
#include "opfunctions.h"

using std::complex;
using std::vector;

typedef boost::multi_array<complex<double>,2> array2d;
	  
/* Return a vector whose elements are crystalline particles as idenfified by the order parameter.
 */

vector<int> xtalpars(const array2d& qlmt, const vector<int>& numneigh,
							const vector<vector<int> >& lneigh, const int nsurf,
							const int nlinks, const double linkval,
							const int lval)
{
	  vector<int> xtals;
	  array2d::index npar = qlmt.shape()[0];
	  int nlin,k;
	  double linval;

	  // compute dot product \tilde{qlm}(i).\tilde{qlm}(j) for each neighbour pair,
	  // whenever this is greater than the threshold (linkval), we call this a crystal link
	  for (array2d::index i = nsurf; i != npar; ++i) {
			 // only interested in particle if it has at least nlinks neighbours
			 if (numneigh[i] >= nlinks) {
					nlin = 0;
					for (int j = 0; j != numneigh[i]; ++j) {
						  k = lneigh[i][j];
						  linval = 0.0;
						  for (int m = 0; m != 2*lval + 1; ++m)
								 linval += qlmt[i][m].real()*qlmt[k][m].real() +
										     qlmt[i][m].imag()*qlmt[k][m].imag();
						  if (linval >= linkval)
								 nlin = nlin + 1;
					}
					// if particle has >=nlinks crystal links, it is in a crystal environment
					if (nlin >= nlinks) {
						  xtals.push_back(i);
					}
			 }
	  }
	  return xtals;
}

/* Get Q of all particles in pnums, which gives indexes into particles.
	This can be used to get Q global, or Q cluster, depending on pnums.
*/

double Qpars(const array2d& qlmb, const vector<int>& pnums, const int lval)
{
	  vector<complex<double> > qlmaverage(2*lval + 1,0.0);

	  for (array2d::index i = 0; i != pnums.size(); ++i) {
			 for (int m = 0; m != 2*lval + 1; ++m) {
					qlmaverage[m] += qlmb[pnums[i]][m];
			 }
	  }
	  
	  for (int m = 0; m != 2*lval + 1; ++m) {
			 qlmaverage[m] = qlmaverage[m]/ (static_cast<double>(pnums.size()));
	  }

	  double qvalue = 0.0;
	  for (int m = 0; m != 2*lval + 1; ++m) {
			 qvalue += norm(qlmaverage[m]);
	  }

	  qvalue = sqrt(qvalue*(4.0*PI/(2*lval + 1)));
	  return qvalue;
}

/* Convert matrix of Nb(i)*\bar{qlm}(i) to matrix of \tilde{qlm}(i) (in Amanda notation)
 */

void qlmtildes(array2d& qlm, const vector<int>& numneigh, const int lval)
{
	  // normalise each of rows in the matrix, this gives qlmtilde
	  int npar = qlm.shape()[0];	  
	  for (int i = 0; i != npar; ++i) {
			 if (numneigh[i] >= 1) {
					double qnorm = 0.0;
					for (int k = 0; k != 2*lval + 1; ++k) 
						  qnorm = qnorm + norm(qlm[i][k]);
					qnorm = sqrt(qnorm);
					for (int k = 0; k != 2*lval + 1; ++k) 
						  qlm[i][k] = qlm[i][k]/qnorm;
			 }
	  }
}

/* Convert matrix of Nb(i)*\bar{qlm}(i) to matrix of \bar{qlm}(i) (in Amanda notation)
 */

void qlmbars(array2d& qlm, const vector<int>& numneigh, const int lval)
{
	  // divide each row of the matrix by the number of neighbours, this gives qlmbar
	  int npar = qlm.shape()[0];
	  for (int i = 0; i != npar; ++i) {
			 if (numneigh[i] >= 1) {
					for (int k = 0; k != 2*lval + 1; ++k) 
						  qlm[i][k] = qlm[i][k]/((double) numneigh[i]);
			 }
	  }
}

/* Return matrix of Nb(i)*\bar{qlm}(i) (in Amanda notation).
 */

array2d qlms(const vector<Particle>& particles, const Box& simbox,
				 vector<int>& numneigh, vector<vector<int> >& lneigh, const int lval)
{
	  vector<Particle>::size_type npar = particles.size();
	  
     // 2d array of complex numbers to store qlm for each particle
	  array2d qlm(boost::extents[npar][2*lval + 1]);
	  std::fill( qlm.origin(), qlm.origin() + qlm.size(), 0.0 );
	  
	  double r2,r,costheta,phi,rh;
	  double sep[3];
	  int k,m;
	  vector<Particle>::size_type i,j;
	  
	  for (i = 0; i != npar; ++i) {
			 for (j = 0; j != npar; ++j) {
					if (i != j) {
						  simbox.sep(particles[i],particles[j], sep);
						  if (simbox.isneigh(sep,r2)) {
								 // particles i and j are neighbours
								 ++numneigh[i];
								 lneigh[i].push_back(j);

								 // compute angles cos(theta) and phi in spherical coords
								 r = sqrt(r2);
								 costheta = sep[2]/r;
								 rh = sqrt(sep[0]*sep[0] + sep[1]*sep[1]);
								 if ((sep[0] == 0.0) && (sep[1] == 0.0)) {
										phi = 0.0;
								 }
								 else if (sep[1] > 0.0) {
										phi = acos(sep[0]/rh);
								 }
								 else {
										phi = 2.0*PI - acos(sep[0]/rh);
								 }

								 // compute contribution of particle j to q6 of particle i
								 for (k = 0; k != 2*lval + 1; ++k) {
										m = -lval + k;
										// spherical harmonic
										qlm[i][k] += ylm(lval, m, costheta, phi);
								 }
						  }
					}
			 }
	  
			 // Now we have an array of Nb(i)*\bar{qlm}(i) in Amanda notation
			 // Note that we either want to convert this to \bar{qlm}(i),
			 // which is needed for computing global OPS e.g Q_6^cl
			 // OR we want return normalized values \tidle{qlm} in Amanda notation
			 // which is needed for computing local order parameters
	  }

	  return qlm;
}
