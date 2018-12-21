#include "itpFile.h"

ItpFile::ItpFile(const string &filename) {
  ifstream file(filename);
  if (!file.good()) {
    printf("ERROR: Itp file %s cannot be read.\n",filename.c_str());
    exit(EXIT_FAILURE);
  }

  string   line;
  int ii=0;
  while(getline(file, line)) {
    //skip lines until find [ atoms ]
    if (line.compare("[ atoms ]"))
      continue;

    //skip comment
    if (line.compare(0,1,";")==0)
      continue;

    //empty line indicates end of atoms section
    if (line.empty())
      break;

    string entry;
    stringstream   linestream(line);

    //skip first 3 columns
    linestream >> entry;
    linestream >> entry;
    linestream >> entry;
    
    //get residue
    linestream >> res[ii];

    //get type
    linestream >> type[ii];

    //skip line
    linestream >> entry;

    //get charge
    linestream >> entry;
    charge[ii] = stof(entry);

    ii++;
  }

  nTypes=ii;

  //TODO: could optimize this to remove redundant types

  //TODO: check if type is only thing necessary for idenetifying charge
}

ItpFile::~ItpFile() {
  delete[] res;
  delete[] type;
  delete[] charge;
}

int ItpFile::findType(const string &s) const {
  for (int ii=0; ii<nTypes; ii++) 
    if (s.compare(type[ii]) == 0) {
      return ii;
    }

  return -1;
}
