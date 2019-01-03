#include "calcW.h"

#define NDIST 0.2  //Distance between C and N atom in nm

CalcW::CalcW(const int nchrom, const Charges &chg, const GroFile &gro) :
  nchrom(nchrom), q(chg.getCharges()), type(gro.getType()),
  res(gro.getRes()), resnum(gro.getResNum()), chain(gro.getChain()),
  cut(2.0*A0INV), cut2(cut*cut), cutN(cut+NDIST*A0INV), cutN2(cutN*cutN)
{
  freq   = new float[nchrom];
  CO     = new rvec[nchrom];
}

CalcW::~CalcW() {
  delete[] freq;
  delete[] CO;
}

void CalcW::compute(const Traj &traj, const vector<int> &inds) {
  const rvec  *x=traj.getCoords();
  rvec box;
  traj.getBox(box);
  ts=traj.getNT()-1;
  int ii,jj;
  //check that box is larger than 2*cutoff
  for (ii=0; ii<DIM; ii++)
    if (box[ii]<2*cut) {
      printf("ERROR: Box is smaller than twice the cutoff\n");
      exit(EXIT_FAILURE);
    }

  rvec tmpCO, tmpvec, Cpos, tmpEn, tmpEc;
  float d2,d;
  int natoms = traj.getNatoms();
  int thisChain,thisRes,nnL,nnR;
  //loop through labelled CO and calculate CO vec
  for (ii=0; ii<nchrom; ii++) {
    setRvec(Cpos,x[inds[ii]]);
    addRvec(Cpos,x[inds[ii]+1],tmpCO,-1);
    pbc(tmpCO,box);
    setRvec(CO[ii],tmpCO);
    d2=norm2vec(tmpCO);
    d=sqrt(d2);
    multRvec(tmpCO,1.0/d);

    thisChain=chain[ii];

    //include NN peptide shifts
    thisRes=resnum[ii];
    nnL=thisRes-1;
    nnR=thisRes+1;
    calcAngles(ii,thisRes,thisChain);

    //loop through other atoms
    setRvec(tmpEn,0.0);
    setRvec(tmpEc,0.0);
    for (jj=0; jj<natoms; jj++) {
      //exclude NN peptides
      //TODO: should NN side chains be excluded?
      if (thisChain==chain[jj] &&
	  (resnum[jj]==nnL || resnum[jj]==nnR) )
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
    freq[ii] = 1618 + 7729*dot(tmpEc,tmpCO) - 3576*dot(tmpEn,tmpCO);
  }
}

void CalcW::print(FILE *fFreq, FILE *fDip) {
  fprintf(fFreq,"%d",ts);
  fprintf(fDip,"%d",ts);
  for (int ii=0; ii<nchrom; ii++) {
    fprintf(fFreq," %f",freq[ii]);
    for (int kk=0; kk<3; kk++) {
      fprintf(fDip," %f",CO[ii][kk]);
    }
  }
  fprintf(fFreq,"\n");
  fprintf(fDip,"\n");
}

uint CalcW::calcAngles(const int atomI, const int resI, const int chainI) {
  //check if already calculated this angle
  uint ind;
  for (ind=0; ind<angleID.size(); ind++)
    if (angleID[ind]==resI && angleCID[ind]==chainI)
      break;
  if (ind != angleID.size())
    return ind;

  angleID.push_back(resI);
  angleCID.push_back(chainI);
  calc1angle(
}
