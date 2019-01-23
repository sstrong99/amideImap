#include "compareEnergy.h"

CompareEnergy::CompareEnergy(const string &reffilename,
			     const string &outfilename, const int nchrom) :
  Compare(reffilename,outfilename,nchrom), energy(nchrom), diff(nchrom) {}

void CompareEnergy::readline() {
  string word;
  stringstream linestream(line);
  linestream >> word;
  ts=stoi(word);
  int skip=nchrom-1;
  for (int ii=0; ii<nchrom; ii++) {
    linestream >> word;
    energy[ii] = stof(word);
    for (int jj=0; jj<skip; jj++)
      linestream >> word; //skip couplings
    skip--;
  }
}

void CompareEnergy::calcDiff(const CalcFreq &calcFreq) {
  if (ts != calcFreq.getTS()) {
    printf("ERROR: comparison with %s is not at the same timestep.\n",reffilename.c_str());
    exit(EXIT_FAILURE);
  }

  for (int ii=0; ii<nchrom; ii++)
    diff[ii]=calcFreq.getHam(ii)-energy[ii];
}

void CompareEnergy::print() {
  fprintf(outfile,"%d",ts);
  for (int ii=0; ii<nchrom; ii++)
    fprintf(outfile," %f",diff[ii]);
  fprintf(outfile,"\n");
}
