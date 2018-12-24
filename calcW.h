#ifndef CALCW_H
#define CALCW_H

#include "h2xind.h"
#include "traj.h"
#include "vecManip.h"

#include <cmath>
#include <vector>

#define CM2PS 0.18836516 //cm/ps = 1/2pi*c
#define HART2CM 2.1947463e5 //convert hartree to wavenumber

class CalcW {
public:
  CalcW(const int model,const int natoms,float avef=0.0);
  ~CalcW();
  void compute(Traj &traj,rvec *m,vector<int> inds);

  void getW(float* out) const ;
  const float* getW() const { return wMat; };
  inline float getW(int nn,int mm) const
  { return wMat[mm+nn*nH]; };

  float maxInterFreq() const;
  float avgIntraFreq() const;

private:
  void calcE(const Traj &traj);
  void mapE2W(rvec *m);

  inline void setNN(float *M, const float val, const int nn, const int mm)
  { M[mm+nn*nH]=val; };

  float *En;      //scalar Efield at each N
  float *Ec;      //scalar Efield at each C
  float *wMat;

  //map
  const float cut2;

};


#endif
