#include "compareDipole.h"

CompareDipole::CompareDipole(const string &filename,const int nchrom) : Compare(filename,nchrom), dipole(nchrom), diff(nchrom) {}

void CompareDipole::readline() {
  string word;
  stringstream linestream(line);
  linestream >> word;
  ts=stoi(word);
  for (int jj=0; jj<3; jj++) {
    for (int ii=0; ii<nchrom; ii++) {
      linestream >> word;
      dipole[ii][jj]=stof(word);
    }
  }
}

void CompareDipole::calcDiff(const CalcW &calcW) {
  if (ts != calcW.getTS()) {
    printf("ERROR: comparison with %s is not at the same timestep.\n",filename.c_str());
    exit(EXIT_FAILURE);
  }

  for (int ii=0; ii<nchrom; ii++)
    for (int jj=0; jj<3; jj++)
      diff[ii][jj]=calcW.getCO(ii,jj)/dipole[ii][jj];
}

void CompareDipole::print(FILE *fh) {
  fprintf(fh,"%d",ts);
  for (int ii=0; ii<nchrom; ii++)
    for (int jj=0; jj<3; jj++)
      fprintf(fh," %f",diff[ii][jj]);
  fprintf(fh,"\n");
}
