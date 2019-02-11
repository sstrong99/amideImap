#include "timer.h"
#include "input.h"
#include "groFile.h"
#include "charges.h"
#include "traj.h"
#include "calcFreq.h"
#include "compareHam.h"
#include "compareDipole.h"

#include <cstdio>
#include <string>

//using namespace std;
int main(int argc, const char *argv[])
{
  string inputfile="in.amideI"; //default input file
  if (argc == 2 || argc == 3)
    inputfile=argv[1];
  Input input(inputfile);

  //if supplied trajfile, override input file
  if (argc == 3) {
    string trajfile=argv[2];
    input.setTrajFile(trajfile);
  }

  Timer time_entire;

  //read files
  GroFile gro(input.getGroFile());
  Charges charges(input,gro);
  Traj traj((input.getTrajFile()).c_str());

  printf("Computing Frequencies...\n");

  //get indicies of C=O using gro.findCarboxyl(res,num)
  vector<int> indC = gro.getChromList(input,4);

  //init calcFreq
  int nchrom = indC.size();
  //use residues as cutoff groups
  //CalcFreq calcFreq(nchrom,charges,gro,gro.getResSt(),gro.getNres());

  //use charge groups for cutoff
  CalcFreq calcFreq(nchrom,charges,gro,charges.getCGst(),
	      charges.getCG(),charges.getNcg());

  //loop through timesteps
  CompareHam    cmpE(input.getErefFile(),input.getEdiffFile(),nchrom);
  CompareDipole cmpD(input.getDrefFile(),input.getDdiffFile(),nchrom);
  FILE *fFreq=fopen(input.getEnergyFile().c_str(),"w");
  FILE *fDip=fopen(input.getDipoleFile().c_str(),"w");
  while(traj.next()==0) {
    calcFreq.compute(traj,indC);
    calcFreq.print(fFreq,fDip);

    cmpE.compare(calcFreq);
    cmpD.compare(calcFreq);
  }
  fclose(fFreq);
  fclose(fDip);

  string time = time_entire.getDHMS();
  printf("Completed in %s\n",time.c_str());

  return EXIT_SUCCESS;
};
