#ifndef GROFILE_H
#define GROFILE_H

#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <regex>

using namespace std;
class GroFile
{
 public:
  GroFile(const string &filename);
  ~GroFile();

  inline string getType(const int ii) const { return type[ii]; }
  inline string getRes(const int ii) const { return res[ii]; }
  inline int    getNatom() const {return natom;}
  vector<int> findCarboxyl(const string &whichRes,const int whichNum) const;
  
 private:
  string *type;
  string *res;
  int    *resnum;
  int     natom;

  static string extractAndTrim(const string &s, const int a,const int b);
};


#endif
