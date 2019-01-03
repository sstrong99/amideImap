#ifndef CALCW_H
#define CALCW_H

#include "traj.h"
#include "charges.h"
#include "vecManip.h"

#include <cmath>
#include <vector>
#include <cstdio>

#define CM2PS 0.18836516 //cm/ps = 1/2pi*c
#define HART2CM 2.1947463e5 //convert hartree to wavenumber

class CalcW {
 public:
  CalcW(const int nchrom,const Charges &chg);
  ~CalcW();
  void compute(const Traj &traj, const vector<int> &inds);
  void print(FILE *fFreq, FILE *fDip);

  float getCO(const int ii,const int jj) const {return CO[ii][jj];};
  float getFreq(const int ii) const {return freq[ii];};
  int getTS() const {return ts;};

 private:
  const int nchrom;
  const float *q;
  const bool *exclude;
  int ts;

  float *freq;
  rvec  *CO;      //transition dipole moments

  //efield cutoffs
  const float cut,cut2,cutN,cutN2;
};

#endif
