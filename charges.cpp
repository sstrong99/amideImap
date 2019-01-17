#include "charges.h"

Charges::Charges(const Input &input,const GroFile &gro) :
  type(gro.getType()), res(gro.getRes()), resnum(gro.getResNum())
{
  printf("Reading itp files and finding Charges...\n");

  //read all itp files
  int nITP=input.getNitp();
  std::vector<ItpFile>  itps;
  for (int ii=0; ii<nITP; ii++)
    itps.push_back(ItpFile(input.getITPfile(ii)));


  //find charges for each atom in gro
  natom = gro.getNatom();
  charges = new float[natom];

  string tmptype,tmpres;
  int ii,jj,tmpresnum,cgThis;
  int thisI=-1;
  int cgLast=-1;
  int itpLast=-1;
  for (ii=0; ii<natom; ii++) {
    //loop through itp files to find charge
    tmptype   = type[ii];
    tmpres    = res[ii];
    tmpresnum = resnum[ii];
    for (jj=0; jj<nITP; jj++) {
      thisI=itps[jj].findType(tmpresnum, tmptype, tmpres);
      if (thisI >= 0) {
	charges[ii]=itps[jj].getQ(thisI);
	break;
      }
    }

    //test that jj==nITP for failure
    if (jj==nITP) {
      printf("ERROR: no charge found for %d%s: %s\n",
	     tmpresnum, tmpres.c_str(), tmptype.c_str());
      exit(EXIT_FAILURE);
    }

    cgThis=itps[jj].getCG(thisI);
    if (cgLast!=cgThis || itpLast!=jj) {
      cgSt.push_back(ii);
      cgLast=cgThis;
      itpLast=jj;
    }
  }

  nCG=cgSt.size();
  //add end of last cg
  cgSt.push_back(ii);
}

Charges::~Charges() {
  delete[] charges;
}
