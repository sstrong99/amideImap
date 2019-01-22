#include "timer.h"
#include "input.h"
#include "groFile.h"
#include "charges.h"
#include "traj.h"
#include "calcW.h"
#include "compareEnergy.h"
#include "compareDipole.h"

#include <cstdio>
#include <string>

//using namespace std;
int main(int argc, const char *argv[])
{
  string inputfile="in.amideI"; //default input file
  if (argc == 2)
    inputfile=argv[1];
  Input input(inputfile);

  Timer time_entire;

  //read files
  GroFile gro(input.getGroFile());
  Charges charges(input,gro);
  Traj traj((input.getTrajFile()).c_str());

  printf("Computing Frequencies...\n");

  //get indicies of C=O using gro.findCarboxyl(res,num)
  vector<int> indC = gro.getChromList(input,4);

  //init calcW
  int nchrom = indC.size();
  //use residues as cutoff groups
  //CalcW calcW(nchrom,charges,gro,gro.getResSt(),gro.getNres());

  //use charge groups for cutoff
  CalcW calcW(nchrom,charges,gro,charges.getCGst(),
	      charges.getCG(),charges.getNcg());

  //loop through timesteps
  CompareEnergy cmpE(input.getErefFile(),input.getEdiffFile(),nchrom);
  CompareDipole cmpD(input.getDrefFile(),input.getDdiffFile(),nchrom);
  FILE *fFreq=fopen(input.getEnergyFile().c_str(),"w");
  FILE *fDip=fopen(input.getDipoleFile().c_str(),"w");
  while(traj.next()==0) {
    calcW.compute(traj,indC);
    calcW.print(fFreq,fDip);

    cmpE.compare(calcW);
    cmpD.compare(calcW);
  }
  fclose(fFreq);
  fclose(fDip);

  string time = time_entire.getDHMS();
  printf("Completed in %s\n",time.c_str());

  return EXIT_SUCCESS;
};
