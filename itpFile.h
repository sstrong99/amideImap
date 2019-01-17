#ifndef ITPFILE_H
#define ITPFILE_H

#include <sstream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;
class ItpFile {
 public:
  ItpFile(const string &filename);
  ~ItpFile() {};

  //inline string getType(const int ii) const { return type[ii]; }
  //inline string getRes(const int ii) const { return res[ii]; }
  inline float getQ(const int ii) const { return charge[ii]; }
  inline int getCG(const int ii) const { return cgnr[ii]; }

  int findType(const int whichnum, const string &whichtype,
	       const string &whichres) const;

 private:
  vector<string> res;
  vector<int>    resnum;
  vector<int>    cgnr;   //charge group number
  vector<string> type;
  vector<float>  charge;

  int nTypes;
  bool solvFlag;

};

#endif
