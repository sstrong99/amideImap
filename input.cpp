#include "input.h"

Input::Input(const string &inputfile)
{
  ifstream file(inputfile);
  if (!file.good()) {
    printf("ERROR: Input file %s cannot be read.\n",inputfile.c_str());
    exit(EXIT_FAILURE);
  }

  //default vaules
  string eFileDef="Energy.txt";
  string dFileDef="Dipole.txt";
  string eDiffDef="Energy_diff.txt";
  string dDiffDef="Dipole_diff.txt";
  trajFile="";
  groFile="";
  eFile=eFileDef;
  dFile=dFileDef;
  eRefFile="";
  dRefFile="";
  eDiffFile=eDiffDef;
  dDiffFile=dDiffDef;

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
    if (key.compare("trajFile")==0) {
      checkEmpty(trajFile,key);
      linestream >> trajFile;
      if (linestream.fail()) eolErr(key);
    } else if (key.compare("energyFile")==0) {
      checkDefault(eFile,eFileDef,key);
      linestream >> eFile;
      if (linestream.fail()) eolErr(key);
    } else if (key.compare("dipoleFile")==0) {
      checkDefault(dFile,dFileDef,key);
      linestream >> dFile;
      if (linestream.fail()) eolErr(key);
    } else if (key.compare("groFile")==0) {
      checkEmpty(groFile,key);
      linestream >> groFile;
      if (linestream.fail()) eolErr(key);
    } else if (key.compare("itpFile")==0) {
      linestream >> tmpstr;
      if (linestream.fail()) eolErr(key);
      itpFiles.push_back(tmpstr);
    }  else if (key.compare("residue")==0) {
      linestream >> tmpstr;
      if (linestream.fail()) eolErr(key);
      resNames.push_back(tmpstr);
      linestream >> tmpstr;
      if (linestream.fail()) eolErr(key);
      resNums.push_back(stoi(tmpstr));
    } else if (key.compare("energyRef")==0) {
      checkEmpty(eRefFile,key);
      linestream >> eRefFile;
      if (linestream.fail()) eolErr(key);
    } else if (key.compare("dipoleRef")==0) {
      checkEmpty(dRefFile,key);
      linestream >> dRefFile;
      if (linestream.fail()) eolErr(key);
    } else if (key.compare("energyDiff")==0) {
      checkDefault(eDiffFile,eDiffDef,key);
      linestream >> eDiffFile;
      if (linestream.fail()) eolErr(key);
    } else if (key.compare("dipoleDiff")==0) {
      checkDefault(dDiffFile,dDiffDef,key);
      linestream >> dDiffFile;
      if (linestream.fail()) eolErr(key);
    } else if (key.compare(0,1,"#")==0) { //skip comment
      continue;
    } else {
      printf("ERROR: unrecognized keyword %s in %s\n",
	     key.c_str(),inputfile.c_str());
      exit(EXIT_FAILURE);
    }
  }

  //check for bad input
  if (trajFile.empty()) {
    printf("ERROR: trajFile keyword not specified\n");
    exit(EXIT_FAILURE);
  } else if (groFile.empty()) {
    printf("ERROR: groFile keyword not specified\n");
    exit(EXIT_FAILURE);
  } else if (itpFiles.size()==0) {
    printf("ERROR: no itpFiles specified\n");
    exit(EXIT_FAILURE);
  } else if (resNames.size()==0) {
    printf("ERROR: no residues specified\n");
    exit(EXIT_FAILURE);
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

void Input::checkEmpty(const string &var,const string &name) {
  if (!var.empty()) {
    printf("ERROR: %s already set to %s\n",name.c_str(), var.c_str());
    exit(EXIT_FAILURE);
  }
}

void Input::checkDefault(const string &var,const string &def,
			 const string & name) {
  if (var.compare(def)) {
    printf("ERROR: %s already set to %s\n",name.c_str(),def.c_str());
    exit(EXIT_FAILURE);
  }
}

void Input::eolErr(const string &name) {
  printf("ERROR: value was not supplied for keyword %s\n",name.c_str());
  exit(EXIT_FAILURE);
}
