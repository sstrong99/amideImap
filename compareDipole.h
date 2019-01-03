#ifndef COMPARE_DIPOLE_H
#define COMPARE_DIPOLE_H

#include "compare.h"

#include <sstream>
#include <vector>

class CompareDipole : public Compare {
public:
  CompareDipole(const string &filename, const int nchrom);
  ~CompareDipole() {};

protected:
  void readline();
  void calcDiff(const CalcW &calcW);
  void print(FILE *fh);

private:
  vector<rvec> dipole;
  vector<rvec> diff;
};

#endif
