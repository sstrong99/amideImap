#include "compareDipole.h"

CompareDipole::CompareDipole(const string &reffilename,
			     const string &outfilename,const int nchrom) :
  Compare(reffilename,outfilename,nchrom), dipole(nchrom), diff(nchrom) {}

void CompareDipole::readline() {
  string word;
  stringstream linestream(line);
  linestream >> word;
  ts=stoi(word);
  for (int jj=0; jj<DIM; jj++) {
    for (int ii=0; ii<nchrom; ii++) {
      linestream >> word;
      dipole[ii][jj]=stof(word);
    }
  }
}

void CompareDipole::calcDiff(const CalcFreq &calcFreq) {
  if (ts != calcFreq.getTS()) {
    printf("ERROR: comparison with %s is not at the same timestep.\n",reffilename.c_str());
    exit(EXIT_FAILURE);
  }

  for (int ii=0; ii<nchrom; ii++)
    for (int jj=0; jj<DIM; jj++)
      diff[ii][jj]=calcFreq.getDip(ii,jj)/dipole[ii][jj];
}

void CompareDipole::print() {
  fprintf(outfile,"%d",ts);
  for (int ii=0; ii<nchrom; ii++)
    for (int jj=0; jj<DIM; jj++)
      fprintf(outfile," %f",diff[ii][jj]);
  fprintf(outfile,"\n");
}
