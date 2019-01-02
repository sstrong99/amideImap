#ifndef COMPARE_H
#define COMPARE_H

#include <string>
#include <fstream>
#include <string>
#include <cstdlib>

using namespace std;
class Compare {
public:
  Compare(const string &filename,const int nchrom);
  ~Compare() {};

  void next();

protected:
  string filename;
  ifstream file;
  int nchrom;

  string line;
  int stat;
  int ts;

  virtual void readline() = 0;
};

#endif
