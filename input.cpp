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

  nITP=0;

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
      string file;
      //first count number of files
      while(linestream >> file)
	nITP++;
      if (nITP == 0) {
	printf("ERROR: no itp files found\n");
	exit(EXIT_FAILURE);
      }
      
      //now get files
      itpFiles = new string[nITP];
      stringstream   tmpstream(line);
      tmpstream >> key;
      for (int ii=0; ii<nITP; ii++)
	tmpstream >> itpFiles[ii];
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

Input::~Input() {
  delete[] itpFiles;
}

string Input::getITPfile(const int ii) const {
  if (ii >= nITP || ii < 0) {
    printf("ERROR: invalid ITP file index\n");
    exit(EXIT_FAILURE);
  } else
    return itpFiles[ii];
}
  
