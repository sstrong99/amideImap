#ifndef CALCW_H
#define CALCW_H

#include "traj.h"
#include "charges.h"
#include "vecManip.h"

#include <cmath>
#include <vector>

#define CM2PS 0.18836516 //cm/ps = 1/2pi*c
#define HART2CM 2.1947463e5 //convert hartree to wavenumber

class CalcW {
 public:
  CalcW(const int nchrom,const Charges &chg);
  ~CalcW();
  void compute(const Traj &traj, rvec *m, const vector<int> &inds);

  const float* getW() const { return freq; };

 private:
  const int nchrom;
  const float *q;
  const bool *exclude;
  
  float *En;      //scalar Efield at each N
  float *Ec;      //scalar Efield at each C
  float *freq;

  //efield cutoffs
  const float cut,cut2,cutN,cutN2;
};

#endif
