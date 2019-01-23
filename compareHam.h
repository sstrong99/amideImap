#ifndef COMPARE_HAM_H
#define COMPARE_HAM_H

#include "compare.h"

#include <sstream>
#include <vector>

class CompareHam : public Compare {
public:
  CompareHam(const string &reffilename, const string &outfilename,
		const int nchrom);
  ~CompareHam() {};

protected:
  void readline();
  void calcDiff(const CalcFreq &calcFreq);
  void print();

private:
  const int nham;

  vector<float> ham;
  vector<float> diff;
};

#endif
