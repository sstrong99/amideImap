#include "groFile.h"

#define RES_ST 0
#define RES_L  9
#define TYP_ST 9
#define TYP_L  6

GroFile::GroFile(const string &filename) {
  printf("Reading file %s...\n",filename.c_str());

  ifstream file(filename);
  if (!file.good()) {
    printf("ERROR: Gro file %s cannot be read.\n",filename.c_str());
    exit(EXIT_FAILURE);
  }

  //skip first line
  string line;
  getline(file,line);

  //get number of atoms
  getline(file,line);
  natom=stoi(line);
  type   = new string[natom];
  res    = new string[natom];
  resnum = new int[natom];
  chain  = new int[natom];
  backbone = new bool[natom];

  resnumAll = new int[natom];

  string tmp,tmpres,tmptype,lastresname;
  int resdiff,lastres;
  int chainid=0;
  int nAtomIn=0;
  int thisResSt=0;
  nres=0;
  for (int ii=0; ii<natom; ii++) {
    if (!getline(file,line)) {
      printf("ERROR: number of atoms in gro file is incorrect.\n");
      exit(EXIT_FAILURE);
    }

    tmp = extractAndTrim(line,RES_ST,RES_L);
    size_t jj;
    for(jj = 0; jj < tmp.size(); jj++)
      if (std::isalpha(tmp[jj]))
	break;
    tmpres     = tmp.substr(jj);
    res[ii]    = tmpres;
    resnum[ii] = stoi(tmp.substr(0,jj));
    tmptype    = extractAndTrim(line,TYP_ST,TYP_L);
    type[ii]   = tmptype;

    if (tmptype.compare("C") ==0  ||
	tmptype.compare("N") ==0  ||
	tmptype.compare("CA")==0  ||
	tmptype.compare("O") ==0  ||
	tmptype.compare("H") ==0    ) //H is excluded too: see Carr JCP 140 2014
      backbone[ii]=true;
    else
      backbone[ii]=false;

    //init vars on first atom
    if (ii==0) {
      lastres=resnum[ii];
      lastresname=tmpres;
    }

    //number distinct residues/molecules
    if (tmpres.compare(lastresname)!=0 && lastres!=resnum[ii]) {
      resnumAll[ii]=nres++;
      lastresname=tmpres;  //only update if different
      atomsInRes.push_back(nAtomIn);
      resSt.push_back(thisResSt);
      thisResSt+=nAtomIn;
      nAtomIn=0;
    } else {
      resnumAll[ii]=nres;
      nAtomIn++;
    }

    //number protein chains from 0-n
    //-1 = not part of protein chain
    if (tmpres.compare("HOH")==0  ||
	tmpres.compare("DPPC")==0 ||
	tmpres.compare("POT")==0  ||
	tmpres.compare("CL")==0 )
    {
      chain[ii]=-1;
    } else {
      resdiff=resnum[ii]-lastres;
      if (resdiff==0 || resdiff==1)
	chain[ii]=chainid;
      else
	chain[ii]=++chainid;
      lastres=resnum[ii];
    }
  }
}

GroFile::~GroFile() {
  delete[] type;
  delete[] res;
  delete[] resnum;
  delete[] chain;
  delete[] resnumAll;
  delete[] backbone;
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
int GroFile::findCarboxyl(const string &whichRes,const int whichNum,
			  const int whichChain) const {
  for (int ii=0; ii<natom; ii++)
    if ( whichChain==chain[ii]        &&
	 whichNum==resnum[ii]         &&
	 whichRes.compare(res[ii])==0 &&
	 type[ii].compare("C")==0 )
      return ii;

  //if reach this point didn't find residue
  string err="ERROR: no "+whichRes+to_string(whichNum)+" residues found...\n";
  printf("%s",err.c_str());
  exit(EXIT_FAILURE);
}

vector<int> GroFile::getChromList(const Input &input, const int nchain) const {
  int nchrom=input.getNchrom();
  vector<int> indC(nchain*nchrom);
  for (int jj=0; jj<nchain; jj++)
    for (int ii=0; ii<nchrom; ii++)
      indC[ii+nchrom*jj] = findCarboxyl(input.getResNames(ii),
					input.getResNums(ii) , jj);

  return indC;
}
