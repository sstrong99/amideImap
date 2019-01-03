#ifndef COMPARE_H
#define COMPARE_H

#include "calcW.h"

#include <string>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdio>

using namespace std;
class Compare {
public:
  Compare(const string &filename,const int nchrom);
  ~Compare() {};

  void compare(const CalcW &calcW,FILE *fh);

protected:
  string filename;
  ifstream file;
  int nchrom;

  string line;
  int stat;
  int ts;

  virtual void readline() = 0;
  virtual void calcDiff(const CalcW &calcW)  = 0;
  virtual void print(FILE *fh) = 0;

};

#endif
