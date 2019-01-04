#include "itpFile.h"

ItpFile::ItpFile(const string &filename) {
  printf("Reading file %s...\n",filename.c_str());

  ifstream file(filename);
  if (!file.good()) {
    printf("ERROR: Itp file %s cannot be read.\n",filename.c_str());
    exit(EXIT_FAILURE);
  }

  string   line;
  int ii=0;
  bool readFlag=false;
  while(getline(file, line)) {
    //skip comment
    if (line.compare(0,1,";")==0 || line.compare(0,1,"#")==0)
      continue;

    //skip lines until find [ atoms ]
    if (!readFlag && line.compare("[ atoms ]")==0) {
      readFlag=true;
      continue;
    }

    if (readFlag) {
      //empty line indicates end of atoms section
      if (line.empty()) {
	readFlag=false;
	break;
      }

      string entry;
      stringstream   linestream(line);

      //skip first 3 columns
      linestream >> entry;
      linestream >> entry;
      linestream >> entry;

      //get residue
      linestream >> entry;
      //res.push_back(entry);
      //don't need residue, b/c each type has specific charge

      //get type
      linestream >> entry;
      type.push_back(entry);

      //skip line
      linestream >> entry;

      //get charge
      linestream >> entry;
      charge.push_back( stof(entry) );

      ii++;
    }
  }

  nTypes=ii;

  //TODO: could optimize this to remove redundant types
}

int ItpFile::findType(const string &s) const {
  for (int ii=0; ii<nTypes; ii++)
    if (s.compare(type[ii]) == 0) {
      return ii;
    }

  return -1;
}
