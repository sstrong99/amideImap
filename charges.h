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

  const float* getCharges() const { return charges; };
  const bool* getExclude() const { return exclude; };

 private:
  int natom;
  float *charges;
  bool  *exclude;

};
#endif
