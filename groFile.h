#ifndef GROFILE_H
#define GROFILE_H

#include "input.h"

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

  const string* getType()   const { return type;   }
  const string* getRes()    const { return res;    }
  const int*    getResNum() const { return resnum; }
  const int*    getChain()  const { return chain;  }

  inline int    getNatom() const {return natom;}

  vector<int> getChromList(const Input &input, const int nchain) const;


 private:
  string *type;
  string *res;
  int    *resnum;
  int    *chain;
  int     natom;

  static string extractAndTrim(const string &s, const int a,const int b);

  int findCarboxyl(const string &whichRes,const int whichNum,
			   const int whichChain) const;
};


#endif
