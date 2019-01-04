#include "charges.h"

Charges::Charges(const Input &input,const GroFile &gro) :
  type(gro.getType()), res(gro.getRes())
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
  int jj,thisI;
  for (int ii=0; ii<natom; ii++) {
    tmpres=res[ii];

    //loop through itp files to find charge
    tmptype = type[ii];
    for (jj=0; jj<nITP; jj++) {
      thisI=itps[jj].findType(tmptype);
      if (thisI >= 0) {
	charges[ii]=itps[jj].getQ(thisI);
	break;
      }
    }

    //test that jj==nITP for failure
    if (jj==nITP) {
      printf("ERROR: no charge found for %s atom type\n",tmptype.c_str());
      exit(EXIT_FAILURE);
    }
  }
}

Charges::~Charges() {
  delete[] charges;
}
