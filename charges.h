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

  const float* getCharges() const { return charges; }

private:
  int natom;
  float *charges;

  const string *type;
  const string *res;
  const int    *resnum;

  vector<int>  cgSt;
};
#endif
