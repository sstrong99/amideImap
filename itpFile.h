#ifndef ITPFILE_H
#define ITPFILE_H

#include <sstream>
#include <fstream>
#include <string>

using namespace std;
class ItpFile {
 public:
  ItpFile(const string &filename);
  ~ItpFile();

  inline string getType(const int ii) const { return type[ii]; }
  inline string getRes(const int ii) const { return res[ii]; }
  inline float getQ(const int ii) const { return charge[ii]; }
  
  int findType(const string &s) const;
  
 private:
  string *res;   //TODO: might not need
  string *type;
  float *charge;

  int nTypes;
  
};

#endif
