#include "groFile.h"

#define RES_ST 0
#define RES_L  9
#define TYP_ST 9
#define TYP_L  6

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
  type   = new string[natom];
  res    = new string[natom];
  resnum = new int[natom];
  chain  = new int[natom];

  int ii=0;
  int chainid=0;
  int lastres=-1;
  int resdiff;
  while(getline(file, line)) {
    tmp = extractAndTrim(line,RES_ST,RES_L);
    size_t jj;
    for(jj = 0; jj < tmp.size(); jj++)
      if (std::isalpha(tmp[jj]))
	break;
    res[ii]  = tmp.substr(jj);
    resnum[ii]= stoi(tmp.substr(0,jj));

    type[ii] = extractAndTrim(line,TYP_ST,TYP_L);

    //number chains from 0-n
    //-1 = not part of protein chain
    if (res[ii].compare("HOH")==0  ||
	res[ii].compare("DPPC")==0 ||
	res[ii].compare("POT")==0  ||
	res[ii].compare("CL")==0 )
    {
      chain[ii]=-1;
    } else {
      //if first protein atom
      if (lastres==-1) {
	lastres=resnum[ii];
	chain[ii]=chainid;
      } else {
	resdiff=resnum[ii]-lastres;
	if (resdiff==0 || resdiff==1)
	  chain[ii]=chainid;
	else
	  chain[ii]=++chainid;
	lastres=resnum[ii];
      }
    }

    ii++;
  }
}

GroFile::~GroFile() {
  delete[] type;
  delete[] res;
  delete[] resnum;
  delete[] chain;
}

string GroFile::extractAndTrim(const string &s, const int a, const int b) {
  string out = s.substr(a,b);

  //remove trailing and leading whitespace
  return regex_replace(out, regex("^ +| +$|( ) +"), "$1");

  //trim leading whitespace
  //size_t start = out.find_first_not_of(" \n\r\t\f\v");
  //return (start == string::npos) ? "" : out.substr(start);
}

//return index of C atom, i think O atom is always i+1 from C atom
vector<int> GroFile::findCarboxyl(const string &whichRes,const int whichNum) const {
  vector<int> list;
  for (int ii=0; ii<natom; ii++)
    if ( whichNum==resnum[ii] &&
	 whichRes.compare(res[ii])==0 &&
	 type[ii].compare("C")==0 )
      list.push_back(ii);

  if (list.size()==0) {
    string err="ERROR: no "+whichRes+to_string(whichNum)+" residues found...\n";
    printf("%s",err.c_str());
    exit(EXIT_FAILURE);
  }

  return list;
}
