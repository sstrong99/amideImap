#ifndef VECMANIP_H
#define VECMANIP_H

#include <xdrfile_xtc.h>
#include <xdrfile.h>
#include <cstdlib>
#include <cmath>

inline void copyRvec(const rvec &M, rvec &out)
{ for (int ii=0; ii<DIM; ii++) out[ii]=M[ii]; };

inline void setRvec(rvec &v, const float val)
{ for (int ii=0; ii<DIM; ii++) v[ii]=val; };

inline void addRvec(const rvec &v1, const rvec &v2, rvec &v, const float &scl) {
//sgn=-1 for subtraction
  for (int ii=0; ii<DIM; ii++) v[ii]=v1[ii]+scl*v2[ii];
};

inline void addRvec(const rvec &v1, rvec &v, const float &scl) {
//sgn=-1 for subtraction
  for (int ii=0; ii<DIM; ii++) v[ii]+=scl*v1[ii];
};

inline void multRvec(rvec &v, float scalar)
{ for (int ii=0; ii<DIM; ii++) v[ii]*=scalar; };

inline float dot(const rvec &v1, const rvec &v2) {
  float out=0.0;
  for (int ii=0; ii<DIM; ii++)
    out+=v1[ii]*v2[ii];
  return out;
};

inline void cross(const rvec &v1, const rvec &v2, rvec &out) {
  out[0]=v1[1]*v2[2] - v1[2]*v2[1];
  out[1]=v1[2]*v2[0] - v1[0]*v2[2];
  out[2]=v1[0]*v2[1] - v1[1]*v2[0];
};

inline float norm2vec(const rvec &v) {
  float d2=0;
  for (int ii=0; ii<DIM; ii++)
    d2+=v[ii]*v[ii];
  return d2;
};

inline void pbc(rvec &v, const rvec &box) {
  for (int ii=0; ii<DIM; ii++)
    v[ii]-=box[ii]*round(v[ii]/box[ii]);
};

inline void pbcOther(rvec &v, const rvec &vTest, const rvec &box) {
  for (int ii=0; ii<DIM; ii++)
    v[ii]-=box[ii]*round(vTest[ii]/box[ii]);
};
#endif
