#include "groFile.h"

#define RES_ST 0
#define RES_SP 7
#define TYP_ST 8
#define TYP_SP 14

GroFile::GroFile(const string &filename) {
  ifstream file(filename);
  if (!file.good()) {
    printf("ERROR: Gro file %s cannot be read.\n",filename.c_str());
    exit(EXIT_FAILURE);
  }

  //skip first line
  string line;
  string tmp;
  getline(file,line);

  //get number of atoms
  getline(file,line);
  natom=stoi(line);
  type = new string[natom];
  res  = new string[natom];

  int ii=0;
  while(getline(file, line)) {
    tmp = extractAndTrim(line,RES_ST,RES_SP);
    size_t jj;
    for(jj = 0; jj < tmp.size(); jj++)
      if (std::isalpha(tmp[jj]))
	break;
    res[ii]  = tmp.substr(jj);
    resnum[ii]= stoi(tmp.substr(0,jj-1));
    
    type[ii] = extractAndTrim(line,TYP_ST,TYP_SP);
    
    ii++;
  }
}

GroFile::~GroFile() {
  delete[] type;
  delete[] res;
}

string GroFile::extractAndTrim(const string &s, const int a, const int b) {
  string out = s.substr(a,b);

  //trim leading whitespace
  size_t start = out.find_first_not_of(" \n\r\t\f\v");
  return (start == string::npos) ? "" : out.substr(start);
}

//return index of C atom, i think O atom is always i+1 from C atom
vector<int> GroFile::findCarboxyl(const string &whichRes,const int whichNum) const {
  vector<int> list;
  for (int ii=0; ii<natom; ii++) 
    if ( whichNum==resnum[ii] && \
	 whichRes.compare(res[ii])==0 && \
	 type[ii].compare("C")==0 )
      list.push_back(ii);

  return list;
}
