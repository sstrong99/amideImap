#include "compare.h"

Compare::Compare(const string &filename,const int nchrom) : filename(filename),
							    file(filename),
							    nchrom(nchrom)
{
  if (!file.good()) {
    printf("ERROR: Compare file %s cannot be read.\n",filename.c_str());
    exit(EXIT_FAILURE);
  }
}

void Compare::next() {
  if (!getline(file, line)) {
    printf("ERROR: no more lines in %s.\n",filename.c_str());
    exit(EXIT_FAILURE);
  }

  readline();
}
