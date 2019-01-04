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
  CalcW calcW(nchrom,charges,gro);

  //loop through timesteps
  CompareEnergy cmpE("Energy_carr.txt",nchrom);
  CompareDipole cmpD("Dipole_carr.txt",nchrom);
  FILE *fFreq=fopen(input.getEnergyFile().c_str(),"w");
  FILE *fDip=fopen(input.getDipoleFile().c_str(),"w");
  FILE *fEdiff=fopen("Energy_diff.txt","w");
  FILE *fDdiff=fopen("Dipole_diff.txt","w");
  while(traj.next()==0) {
    calcW.compute(traj,indC);
    calcW.print(fFreq,fDip);

    cmpE.compare(calcW,fEdiff);
    cmpD.compare(calcW,fDdiff);
  }
  fclose(fFreq);
  fclose(fDip);
  fclose(fEdiff);
  fclose(fDdiff);

  string time = time_entire.getDHMS();
  printf("Completed in %s\n",time.c_str());

  return EXIT_SUCCESS;
};
