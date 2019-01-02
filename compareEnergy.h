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

private:
  vector<float> energy;
};

#endif
