#ifndef INPUT_H
#define INPUT_H

#include <sstream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;
class Input {
public:
  Input(const string &inputfile);
  ~Input() {};

  int    getNitp() const { return itpFiles.size(); }
  string getITPfile(const int ii) const;

  int    getNchrom() const { return resNames.size(); }
  string getResNames(const int ii) const;
  int    getResNums(const int ii)  const;

  string getTrajFile() const {return trajFile;};
  string getOutPostfix() const {return outPostfix;};
  string getGroFile() const {return groFile;};

  string getEnergyFile() const {return "Energy"+outPostfix+".txt";};
  string getDipoleFile() const {return "Dipole"+outPostfix+".txt";};

private:
  string trajFile;
  string outPostfix;
  string groFile;
  vector<string> itpFiles;
  vector<string> resNames;
  vector<int>    resNums;

  template<class T>
    T getVect(const uint ii,const vector<T> vect,const string &name) const {
    if (ii >= vect.size()) {
      string tmp = "ERROR: invalid "+name+" index\n";
      printf("%s",tmp.c_str());
      exit(EXIT_FAILURE);
    } else
      return vect[ii];
  };
};
#endif
