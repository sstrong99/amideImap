#ifndef COMPARE_DIPOLE_H
#define COMPARE_DIPOLE_H

#include "compare.h"

#include <sstream>
#include <vector>

class CompareDipole : public Compare {
public:
  CompareDipole(const string &filename, const string &outfilename,
		const int nchrom);
  ~CompareDipole() {};

protected:
  void readline();
  void calcDiff(const CalcFreq &calcFreq);
  void print();

private:
  vector<rvec> dipole;
  vector<rvec> diff;
};

#endif
