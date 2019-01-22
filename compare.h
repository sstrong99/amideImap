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
  Compare(const string &reffilename,const string &outfilename,const int nchrom);
  ~Compare();

  void compare(const CalcW &calcW);

protected:
  string reffilename,outfilename;
  ifstream reffile;
  FILE *outfile;
  int nchrom;

  bool doNothing;

  string line;
  int stat;
  int ts;

  virtual void readline() = 0;
  virtual void calcDiff(const CalcW &calcW)  = 0;
  virtual void print() = 0;

};

#endif
