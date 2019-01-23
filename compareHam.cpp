#include "compareHam.h"

CompareHam::CompareHam(const string &reffilename,
		       const string &outfilename, const int nchrom) :
  Compare(reffilename,outfilename,nchrom), nham(nchrom*(nchrom+1)/2),
  ham(nham), diff(nham) {}

void CompareHam::readline() {
  string word;
  stringstream linestream(line);
  linestream >> word;
  ts=stoi(word);
  for (int ii=0; ii<nham; ii++) {
    linestream >> word;
    ham[ii] = stof(word);
  }
}

void CompareHam::calcDiff(const CalcFreq &calcFreq) {
  if (ts != calcFreq.getTS()) {
    printf("ERROR: comparison with %s is not at the same timestep.\n",reffilename.c_str());
    exit(EXIT_FAILURE);
  }

  for (int ii=0; ii<nham; ii++)
    diff[ii]=calcFreq.getHam(ii)-ham[ii];
}

void CompareHam::print() {
  fprintf(outfile,"%d",ts);
  for (int ii=0; ii<nham; ii++)
    fprintf(outfile," %f",diff[ii]);
  fprintf(outfile,"\n");
}
