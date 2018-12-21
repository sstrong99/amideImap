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

 private:
  int natom;
  float *charges;

};
#endif
