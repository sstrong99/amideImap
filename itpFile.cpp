#include "itpFile.h"

ItpFile::ItpFile(const string &filename) {
  ifstream file(filename);
  if (!file.good()) {
    printf("ERROR: Itp file %s cannot be read.\n",filename.c_str());
    exit(EXIT_FAILURE);
  }

  string   line;
  uint ii=0;
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

      //skip 2 cols
      linestream >> entry;
      linestream >> entry;

      //get res num
      linestream >> entry;
      resnum.push_back(stoi(entry));

      //get residue
      linestream >> entry;
      //if (entry.compare("SOL")==0)
      //	entry="HOH"; //rename SOL to HOH for compatibility with gro file
      res.push_back(entry);

      //get type
      linestream >> entry;
      type.push_back(entry);

      //get charge group number
      linestream >> entry;
      cgnr.push_back(stoi(entry));

      //get charge
      linestream >> entry;
      charge.push_back( stof(entry) );

      ii++;
    }
  }

  nTypes=ii;

  //check if this itp file has solvent or
  bool lastFlag;
  for (ii=0; ii<res.size(); ii++) {
    if (res[ii].compare("SOL")==0)
      lastFlag=true;
    else
      lastFlag=false;

    if (ii==0)
      solvFlag=lastFlag;
    else if (lastFlag!=solvFlag) {
      printf("ERROR: solvent itp files cannot have other molecules.\n");
      exit(EXIT_FAILURE);
    }
  }
}

int ItpFile::findType(const int whichnum, const string &whichtype,
		      const string &whichres) const {
  if (solvFlag) {    //if solvent, just compare resname and type
    for (int ii=0; ii<nTypes; ii++) {
      if (whichres.compare(res[ii])   == 0 &&
	  whichtype.compare(type[ii]) == 0   )
	return ii;
    }
  } else {           //otherwise compare exact resnum
    for (int ii=0; ii<nTypes; ii++) {
      if (whichnum==resnum[ii]             &&
	  whichres.compare(res[ii])   == 0 &&
	  whichtype.compare(type[ii]) == 0   )
	return ii;
    }
  }

  return -1;

  //TODO: check that multpile itp files don't define a charge?
}
