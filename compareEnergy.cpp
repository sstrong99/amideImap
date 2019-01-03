#include "compareEnergy.h"

CompareEnergy::CompareEnergy(const string &filename,const int nchrom) : Compare(filename,nchrom), energy(nchrom), diff(nchrom) {}

void CompareEnergy::readline() {
  string word;
  stringstream linestream(line);
  linestream >> word;
  ts=stoi(word);
  int skip=nchrom;
  for (int ii=0; ii<nchrom; ii++) {
    linestream >> word;
    energy[ii] = stof(word);
    for (int jj=0; jj<skip; jj++)
      linestream >> word; //skip couplings
    skip--;
  }
}

void CompareEnergy::calcDiff(const CalcW &calcW) {
  if (ts != calcW.getTS()) {
    printf("ERROR: comparison with %s is not at the same timestep.\n",filename.c_str());
    exit(EXIT_FAILURE);
  }

  for (int ii=0; ii<nchrom; ii++)
    diff[ii]=calcW.getFreq(ii)-energy[ii];
}

void CompareEnergy::print(FILE *fh) {
  fprintf(fh,"%d",ts);
  for (int ii=0; ii<nchrom; ii++)
    fprintf(fh," %f",diff[ii]);
  fprintf(fh,"\n");
}
