#ifndef CALCW_H
#define CALCW_H

#include "traj.h"
#include "charges.h"
#include "vecManip.h"
#include "groFile.h"

#include <cmath>
#include <vector>
#include <cstdio>
#include <algorithm>

#define CM2PS 0.18836516 //cm/ps = 1/2pi*c
#define HART2CM 2.1947463e5 //convert hartree to wavenumber

class CalcW {
 public:
  CalcW(const int nchrom,const Charges &chg, const GroFile &gro);
  ~CalcW();
  void compute(const Traj &traj, const vector<int> &inds);
  void print(FILE *fFreq, FILE *fDip);

  float getCO(const int ii,const int jj) const {return CO[ii][jj];};
  float getFreq(const int ii) const {return freq[ii];};
  int getTS() const {return ts;};

 private:
  const int nchrom;
  const float *q;

  const string *type;
  const string *res;
  const int    *resnum;
  const int    *chain;

  int ts;

  float *freq;
  rvec  *CO;      //transition dipole moments

  //efield cutoffs
  const float cut,cut2,cutN,cutN2;

  vector<int> angleID;
  vector<int> angleCID;
  vector<float> phi;
  vector<float> psi;

  uint calcAngles(const int atomI, const int resI, const int chainI);
};

#endif
