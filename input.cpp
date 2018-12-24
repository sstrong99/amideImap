#include "input.h"

Input::Input(const string &inputfile) : outPostfix("")
{
  ifstream file(inputfile);
  if (!file.good()) {
    printf("ERROR: Input file %s cannot be read.\n",inputfile.c_str());
    exit(EXIT_FAILURE);
  }

  string   line;
  string   key;
  string   tmpstr;

  while(getline(file, line))
  {
    //check if line is empty
    if (line.empty())
      continue;

    stringstream   linestream(line);

    linestream >> key;
    if (key.compare("trajFile")==0)
      linestream >> trajFile;
    else if (key.compare("outPostfix")==0) {
      linestream >> outPostfix;
      outPostfix.insert(0,"_"); }
    else if (key.compare("groFile")==0) 
      linestream >> groFile;
    else if (key.compare("itpFile")==0) {
      linestream >> tmpstr;
      itpFiles.push_back(tmpstr);
    }
    else if (key.compare("residue")==0) {
      linestream >> tmpstr;
      resNames.push_back(tmpstr);
      linestream >> tmpstr;
      resNums.push_back(stoi(tmpstr));
    }
    else if (key.compare(0,1,"#")==0)
      continue; //skip comment
    else {
      printf("ERROR: unrecognized keyword %s in %s\n",
	     key.c_str(),inputfile.c_str());
      exit(EXIT_FAILURE);
    }
  }
}

string Input::getITPfile(const int ii) const {
  return getVect(ii, itpFiles, "itp file");
}

string Input::getResNames(const int ii) const {
  return getVect(ii, resNames, "residue name");
}

int Input::getResNums(const int ii) const {
  return getVect(ii, resNums, "residue number");
}
