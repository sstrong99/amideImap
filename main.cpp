#include "timer.h"
#include "input.h"
#include "groFile.h"
#include "charges.h"
#include "traj.h"
#include "calcW.h"

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

  //get indicies of C=O using gro.findCarboxyl(res,num)
  vector<int> indC;
  vector<int> tmp;
  for (int ii=0; ii<input.getNres(); ii++) {
    tmp=gro.findCarboxyl(input.getResNames(ii),input.getResNums(ii));
    indC.insert( indC.end(), tmp.begin(), tmp.end() );
  };

  //init calcW
  int nchrom = indC.size();
  CalcW calcW(nchrom,charges);

  //loop through timesteps
  FILE *fFreq=fopen(input.getEnergyFile().c_str(),"w");
  FILE *fDip=fopen(input.getDipoleFile().c_str(),"w");
  while(traj.next()==0) {
    calcW.compute(traj,indC);
    calcW.print(fFreq,fDip);
  }
  fclose(fFreq);
  fclose(fDip);

  //TODO: add output function to to calcW to print dipoles, energies
  
  string time = time_entire.getDHMS();
  printf("Completed in %s\n",time.c_str());

  return EXIT_SUCCESS;
};
