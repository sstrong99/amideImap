#include "charges.h"

Charges::Charges(const Input &input,const GroFile &gro) {
  int nITP=input.getNitp();
  std::vector<ItpFile>  itps; 
  for (int ii=0; ii<nITP; ii++)
    itps.push_back(ItpFile(input.getITPfile(ii)));

  string tmptype,tmpres;
  natom = gro.getNatom();
  charges = new float[natom];
  exclude = new bool[natom];
  int jj,thisI;
  for (int ii=0; ii<natom; ii++) {
    //set exclude flag to exclude all peptide charges
    tmpres=gro.getRes(ii);
    if (tmpres.compare("HOH")==0 ||
	tmpres.compare("POT")==0 ||
	tmpres.compare("CL" )==0   ) {
      
      exclude[ii]=false;
      
      //loop through itp files to find charge
      tmptype = gro.getType(ii);
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
    } else {
      exclude[ii]=true;
      charges[ii]=0.0;
    }
  }
}

Charges::~Charges() {
  delete[] charges;
  delete[] exclude;
}
