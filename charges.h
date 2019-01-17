#ifndef CHARGES_H
#define CHARGES_H

#include "input.h"
#include "groFile.h"
#include "itpFile.h"
#include <vector>

class Charges
{
public:
  Charges(const Input &input,const GroFile &gro);
  ~Charges();

  const float*      getCharges() const { return charges; }
  const vector<int> getCGst()    const { return cgSt; }
  const int         getNcg()     const { return nCG; }

private:
  int natom;
  float *charges;

  const string *type;
  const string *res;
  const int    *resnum;

  vector<int>  cgSt;  //vector of starting inds of each charge group
  int          nCG;
};
#endif
