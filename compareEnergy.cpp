#include "compareEnergy.h"

CompareEnergy::CompareEnergy(const string &filename,const int nchrom) : Compare(filename,nchrom) {}

void CompareEnergy::readline() {
  string word;
  stringstream linestream(line);
  linestream >> word;
  ts=stoi(word);
  int skip=nchrom;
  for (int ii=0; ii<nchrom; ii++) {
    linestream >> word;
    energy.push_back(stof(word));
    for (int jj=0; jj<skip; jj++)
      linestream >> word; //skip couplings
    skip--;
  }
}
