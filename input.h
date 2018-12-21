#ifndef INPUT_H
#define INPUT_H

#include <sstream>
#include <fstream>
#include <string>

using namespace std;
class Input {
public:
  Input(const string &inputfile);
  ~Input();

  int    getNitp() const { return nITP; }
  string getITPfile(const int ii) const;
  
  string getTrajFile() const {return trajFile;};
  string getOutPostfix() const {return outPostfix;};
  string getGroFile() const {return groFile;};

private:
  string trajFile;
  string outPostfix;
  string groFile;
  string *itpFiles;

  int nITP;
};
#endif
