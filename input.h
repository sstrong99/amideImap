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

  int    getNitp()                 const { return itpFiles.size(); }
  string getITPfile(const int ii)  const;

  int    getNchrom()               const { return resNames.size(); }
  string getResNames(const int ii) const;
  int    getResNums(const int ii)  const;

  string getTrajFile()   const { return trajFile;   }
  string getGroFile()    const { return groFile;    }
  string getErefFile()   const { return eRefFile;   }
  string getDrefFile()   const { return dRefFile;   }
  string getEdiffFile()  const { return eDiffFile;  }
  string getDdiffFile()  const { return dDiffFile;  }

  string getEnergyFile() const {return eFile;}
  string getDipoleFile() const {return dFile;}

  void setTrajFile(const string &newTrajFile) { trajFile=newTrajFile; }

private:
  string trajFile,eFile,dFile,groFile,eRefFile,dRefFile,eDiffFile,dDiffFile;
  vector<string> itpFiles;
  vector<string> resNames;
  vector<int>    resNums;

  static void checkEmpty(const string &var,const string &name);
  static void checkDefault(const string &,const string &,const string &);
  static void eolErr(const string &name);


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
