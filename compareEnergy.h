#ifndef COMPARE_ENERGY_H
#define COMPARE_ENERGY_H

#include "compare.h"

#include <sstream>
#include <vector>

class CompareEnergy : public Compare {
public:
  CompareEnergy(const string &reffilename, const string &outfilename,
		const int nchrom);
  ~CompareEnergy() {};

protected:
  void readline();
  void calcDiff(const CalcFreq &calcFreq);
  void print();

private:
  vector<float> energy;
  vector<float> diff;
};

#endif
