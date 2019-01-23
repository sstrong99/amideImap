#ifndef CALCFREQ_H
#define CALCFREQ_H

#include "traj.h"
#include "charges.h"
#include "vecManip.h"
#include "groFile.h"

#include <cmath>
#include <vector>
#include <cstdio>
#include <algorithm>

class CalcFreq {
 public:
  CalcFreq(const int nchrom,const Charges &chg, const GroFile &gro,
	const vector<int> &grpSt, const int *grpInd, const int nCutGrp);
  ~CalcFreq();
  void compute(const Traj &traj, const vector<int> &inds);
  void print(FILE *fFreq, FILE *fDip);

  float getDip(const int ii,const int jj) const { return dip[ii][jj]; }
  float getHam(const int ii)              const { return ham[ii];     }
  int   getTS()                           const { return ts;          }

 private:
  const int    nchrom;
  const int    nham;
  const int    natoms;

  const string *type;
  const bool   *backbone;
  const int    *chain;
  const int    *resnumAll;
  const float  *q;

  const vector<int> grpSt;   //starting atom of each cutoff group
  const int        *grpInd;  //group ind of each atom
  const int         nCutGrp; //number of groups for cutoff (eg: cg, res, or atoms)

  rvec box;

  int ts;

  float *ham;      //hamiltonian in upper triangular form
  rvec  *dip;      //transition dipole moments

  uint calcAngles(const int atomI, const int resI, const rvec *x,vector<int> &angleID, vector<float> &phi, vector<float> &psi);
  uint search(const int st, const int resI, const string &whichtype);
  float calcDihedral(const rvec &, const rvec &,
		     const rvec &, const rvec &);
  static inline void normalize(rvec &v);
  static float interp2(const float &, const float &, const float z[]);
  void calcCOG(const rvec *x,rvec *cog);
  vector<uint> getExcludes(const uint nnC[3]);

  //**************************************************************************
  //MAP parameters
  //TODO: could read these from file

  //efield cutoffs
  static constexpr float cut     = 2.0*A0INV;      //20A
  static constexpr float cut2    = cut*cut;
  static constexpr float cutExt  = cut+0.3*A0INV;  //add 3A //TODO: 2.5 ok?
  static constexpr float cutExt2 = cutExt*cutExt;

  //freq map function
  static inline float map(const float &Ec, const float &En)
  { return 1618 + 7729*Ec - 3576*En; };

  //transition dipole map
  static void calcTD(const rvec &rCO, const rvec &rCN, rvec &out);
  static void calcTDpos(const rvec &,const rvec &,const rvec &,rvec &);
  static constexpr float coupConst = 2584.615;         //cm^-1 * a0^3
  static constexpr float tdMag     = 2.73;             //D/A/amu^1/2
  static constexpr float tdAngle = 10.0*(PI/180.0);

  //NN frequency shifts, copied from Jansen's AmideImap program
  static constexpr float dTheta = 30.0;
  static constexpr int   nTheta = 13;
  static constexpr float NtermShift[169]={ 0.404,-1.205, 0.789,13.541,22.541,28.548,50.551,28.548,22.541,13.541, 0.789,-1.205, 0.404,-1.887,-1.141, 5.762,14.643,24.000,34.820,24.308,17.676,13.092, 4.738,-2.997,-2.603,-1.887,-3.387,-0.124, 9.386,14.208,22.015,30.925,14.221,11.718, 4.889, 1.667,-3.653,-4.087,-3.387,-0.069, 8.421,13.625,15.852,21.411,21.014,14.828,-0.075,-0.872, 6.685, 0.784,-1.270,-0.069, 6.132,14.360,15.442,17.498,21.387,17.706,-2.303,-8.094,-4.944, 4.026,-7.002, 2.536, 6.132,10.251, 3.635, 9.426,13.395,13.176,-0.075,-7.130,-7.167,-3.443, 4.315, 3.121, 5.735,10.251,-7.081, 7.816, 6.365, 5.061, 1.329,-4.588,-5.550,-4.588, 1.329, 5.061, 6.365, 7.816,-7.081,10.251, 5.735, 3.121, 4.315,-3.443,-7.167,-7.130,-0.075,13.176,13.395, 9.426, 3.635,10.251, 6.132, 2.536,-7.002, 4.026,-4.944,-8.094,-2.303,17.706,21.387,17.498,15.442,14.360, 6.132,-0.069,-1.270, 0.784, 6.685,-0.872,-0.075,14.828,21.014,21.411,15.852,13.625, 8.421,-0.069,-3.387,-4.087,-3.653, 1.667, 4.889,11.718,14.221,30.925,22.015,14.208, 9.386,-0.124,-3.387,-1.887,-2.603,-2.997, 4.738,13.092,17.676,24.308,34.820,24.000,14.643, 5.762,-1.141,-1.887, 0.404,-1.205, 0.789,13.541,22.541,28.548,50.551,28.548,22.541,13.541, 0.789,-1.205, 0.404};
static constexpr float CtermShift[169]={-13.410,-9.262,-0.184, 8.535,17.984,30.352,48.440,30.352,17.984, 8.535,-0.184,-9.262,-13.410,-10.729,-3.802, 3.846,14.153,26.874,41.926,34.992,16.279, 5.857, 1.667,-5.281,-11.611,-10.729,-4.913, 0.205, 4.514,16.736,31.155,35.848,19.171, 5.586,-5.799,-8.260,-8.499,-8.982,-4.913,-0.695,-0.631, 6.910,16.212,24.955,20.718,10.336,-4.799,-20.484,-22.949,-10.223,-4.079,-0.695, 3.455, 2.482, 7.735,15.208,13.495, 9.584,-4.401,-22.259,-30.507,-27.853, 2.842, 1.872, 3.455,10.288,16.772,16.805,13.984, 7.758,-6.108,-18.810,-25.888,-28.228,-13.221, 5.214, 9.661,10.288, 1.788,15.052,12.962, 6.891,-11.364,-22.800,-22.201,-22.800,-11.364, 6.891,12.962,15.052, 1.788,10.288, 9.661, 5.214,-13.221,-28.228,-25.888,-18.810,-6.108, 7.758,13.984,16.805,16.772,10.288, 3.455, 1.872, 2.842,-27.853,-30.507,-22.259,-4.401, 9.584,13.495,15.208, 7.735, 2.482, 3.455,-0.695,-4.079,-10.223,-22.949,-20.484,-4.799,10.336,20.718,24.955,16.212, 6.910,-0.631,-0.695,-4.913,-8.982,-8.499,-8.260,-5.799, 5.586,19.171,35.848,31.155,16.736, 4.514, 0.205,-4.913,-10.729,-11.611,-5.281, 1.667, 5.857,16.279,34.992,41.926,26.874,14.153, 3.846,-3.802,-10.729,-13.410,-9.262,-0.184, 8.535,17.984,30.352,48.440,30.352,17.984, 8.535,-0.184,-9.262,-13.410};
static constexpr float Coupling[169]={ 6.331, 5.951, 4.853, 5.668, 9.597,13.966,21.427,13.966, 9.597, 5.668, 4.853, 5.951, 6.331, 5.645, 4.666, 5.023, 8.648,13.250,18.309,15.142, 9.130, 4.377, 2.899, 4.712, 6.081, 5.645, 4.333, 3.457, 5.621,10.432,15.718,18.485, 9.208, 2.179,-1.699, 0.142, 3.562, 4.969, 4.333, 2.604, 2.194, 5.158,10.530,15.330,11.107, 2.546,-5.915,-6.721,-2.935, 1.790, 3.453, 2.604, 0.248, 0.513, 3.309, 8.556, 9.372, 3.675,-7.395,-11.367,-8.648,-4.273, 0.026, 0.904, 0.248,-1.219,-0.279, 2.404, 4.462, 2.328,-5.560,-12.700,-12.308,-7.319,-2.714,-0.361,-0.648,-1.219,-1.948,-0.646, 0.771, 0.653,-3.904,-10.815,-14.526,-10.815,-3.904, 0.653, 0.771,-0.646,-1.948,-1.219,-0.648,-0.361,-2.714,-7.319,-12.308,-12.700,-5.560, 2.328, 4.462, 2.404,-0.279,-1.219, 0.248, 0.904, 0.026,-4.273,-8.648,-11.367,-7.395, 3.675, 9.372, 8.556, 3.309, 0.513, 0.248, 2.604, 3.453, 1.790,-2.935,-6.721,-5.915, 2.546,11.107,15.330,10.530, 5.158, 2.194, 2.604, 4.333, 4.969, 3.562, 0.142,-1.699, 2.179, 9.208,18.485,15.718,10.432, 5.621, 3.457, 4.333, 5.645, 6.081, 4.712, 2.899, 4.377, 9.130,15.142,18.309,13.250, 8.648, 5.023, 4.666, 5.645, 6.331, 5.951, 4.853, 5.668, 9.597,13.966,21.427,13.966, 9.597, 5.668, 4.853, 5.951, 6.331};

};

#endif
