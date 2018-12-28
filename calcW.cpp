#include "calcW.h"

#define NDIST 0.2  //Distance between C and N atom in nm

CalcW::CalcW(const int nchrom, const Charges &chg) :
  nchrom(nchrom), q(chg.getCharges()), exclude(chg.getExclude()),
  cut(2.0), cut2(cut*cut), cutN(cut+NDIST), cutN2(cutN*cutN)
{
  //TODO: convert tip3p to tip4p?
  En     = new float[nchrom];
  Ec     = new float[nchrom];
  freq   = new float[nchrom];
}

CalcW::~CalcW() {
  delete[] En;
  delete[] Ec;
  delete[] freq;
}

void CalcW::compute(const Traj &traj, rvec *m, const vector<int> &inds) {
  const rvec  *x=traj.getCoords();
  rvec box;
  traj.getBox(box);
  int ii,jj;
  //check that box is larger than 2*cutoff
  for (ii=0; ii<DIM; ii++)
    if (box[ii]<2*cut) {
      printf("ERROR: Box is smaller than twice the cutoff\n");
      exit(EXIT_FAILURE);
    }

  rvec tmpCO, tmpvec, Cpos, tmpEn, tmpEc;
  float d2,d;
  //loop through labelled CO and calculate CO vec
  for (ii=0; ii<nchrom; ii++) {
    setRvec(Cpos,x[inds[ii]]);
    addRvec(Cpos,x[inds[ii]+1],tmpCO,-1);
    pbc(tmpCO,box);
    d2=norm2vec(tmpCO);
    d=sqrt(d2);
    multRvec(tmpCO,1.0/d);
    setRvec(m[ii],tmpCO);
    
    //loop through other atoms
    setRvec(tmpEn,0.0);
    setRvec(tmpEc,0.0);
    for (jj=0; jj<traj.getNatoms(); jj++) {
      if (exclude[jj])
	continue;

      //get distance to C atom
      addRvec(Cpos,x[jj],tmpvec,-1);
      pbc(tmpvec,box);
      d2=norm2vec(tmpvec);
      if (d2 < cutN2) { //1st check if either C or N could be within cutoff
	if (d2 < cut2) { //now check if C is within cutoff
	  d=sqrt(d2);
	  multRvec(tmpvec, q[jj]/(d*d*d) );
	  addRvec(tmpvec,tmpEc,+1);
	}

	//now test if N is within cutoff
	addRvec(x[inds[ii]+2],x[jj],tmpvec,-1); //N to jj
	pbc(tmpvec,box);
	d2=norm2vec(tmpvec);
	if (d2 < cut2) {
	  d=sqrt(d2);
	  multRvec(tmpvec, q[jj]/(d*d*d) );
	  addRvec(tmpvec,tmpEn,+1);
	}
      }
    }
    Ec[ii]=dot(tmpEc,tmpCO);
    En[ii]=dot(tmpEn,tmpCO);
  }

  for (ii=0; ii<nchrom; ii++)
    freq[ii]=1618 + 7729*Ec[ii] - 3576*En[ii];
}
