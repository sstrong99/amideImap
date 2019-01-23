#include "compare.h"

Compare::Compare(const string &reffilename,const string &outfilename,
		 const int nchrom) : reffilename(reffilename),
				     outfilename(outfilename),nchrom(nchrom)
{
  if (reffilename.empty() && outfilename.empty())
    doNothing=true;
  else {
    doNothing=false;

    reffile.open(reffilename);
    if (!reffile.good()) {
      printf("ERROR: Compare file %s cannot be read.\n",reffilename.c_str());
      exit(EXIT_FAILURE);
    }

    outfile=fopen(outfilename.c_str(),"w");
    if (outfile==NULL) {
      printf("ERROR: Output file %s is bad.\n",outfilename.c_str());
      exit(EXIT_FAILURE);
    }
  }
}

Compare::~Compare() {
  if (!doNothing) {
    reffile.close();
    fclose(outfile);
  }
}

void Compare::compare(const CalcFreq &calcFreq) {
  if (doNothing) return;

  if (!getline(reffile, line)) {
    printf("ERROR: no more lines in %s.\n",reffilename.c_str());
    exit(EXIT_FAILURE);
  }

  this->readline();
  this->calcDiff(calcFreq);
  this->print();
}
