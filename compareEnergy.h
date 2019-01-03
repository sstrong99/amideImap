#ifndef COMPARE_ENERGY_H
#define COMPARE_ENERGY_H

#include "compare.h"

#include <sstream>
#include <vector>

class CompareEnergy : public Compare {
public:
  CompareEnergy(const string &filename, const int nchrom);
  ~CompareEnergy() {};

protected:
  void readline();
  void calcDiff(const CalcW &calcW);
  void print(FILE *fh);

private:
  vector<float> energy;
  vector<float> diff;
};

#endif
