#include "calcW.h"


CalcW::CalcW(const int nchrom, const Charges &chg, const GroFile &gro) :
  nchrom(nchrom), natoms(gro.getNatom()), type(gro.getType()),res(gro.getRes()),
  resnum(gro.getResNum()), chain(gro.getChain()),q(chg.getCharges())
{
  freq   = new float[nchrom];
  CO     = new rvec[nchrom];
}

CalcW::~CalcW() {
  delete[] freq;
  delete[] CO;
}

//needed due to a flaw in c++11
constexpr float CalcW::NtermShift[169];
constexpr float CalcW::CtermShift[169];

void CalcW::compute(const Traj &traj, const vector<int> &inds) {
  const rvec  *x=traj.getCoords();
  traj.getBox(box);
  ts=traj.getNT()-1;
  int ii,jj;
  //check that box is larger than 2*cutoff
  for (ii=0; ii<DIM; ii++)
    if (box[ii]<2*cut) {
      printf("ERROR: Box is smaller than twice the cutoff\n");
      exit(EXIT_FAILURE);
    }

  if (traj.getNatoms() != natoms) {
    printf("ERROR: gro file does not have the same number of atoms as traj.\n");
    exit(EXIT_FAILURE);
  }

  rvec tmpCO, tmpvec, Cpos, tmpEn, tmpEc;
  float d2,d;
  int thisAtom,thisChain,thisRes,nnL,nnR;
  //loop through labelled CO and calculate CO vec
  for (ii=0; ii<nchrom; ii++) {
    thisAtom=inds[ii];
    setRvec(Cpos,x[thisAtom]);
    addRvec(Cpos,x[thisAtom+1],tmpCO,-1);
    pbc(tmpCO,box);
    setRvec(CO[ii],tmpCO);
    normalize(tmpCO);



    //include NN peptide shifts
    float nnfs = 0.0;
    thisRes=resnum[thisAtom];
    thisChain=chain[thisAtom];
    nnL=thisRes-1;
    nnR=thisRes+1;
    uint a1,a2,Cnext;
    a1=calcAngles(thisAtom,thisRes,thisChain,x);
    //calculate angles for next resiude
    Cnext=search(thisAtom,thisRes+1,"C");
    a2=calcAngles(Cnext,thisRes+1,thisChain,x);

    nnfs += interp2(phi[a1],psi[a1],NtermShift);
    nnfs += interp2(phi[a2],psi[a2],CtermShift);

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
    freq[ii] = map(dot(tmpEc,tmpCO), dot(tmpEn,tmpCO)) + nnfs;
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

uint CalcW::calcAngles(const int atomI, const int resI, const int chainI,
		       const rvec *x) {
  //check if already calculated this angle
  uint ind;
  for (ind=0; ind<angleID.size(); ind++)
    if (angleID[ind]==resI && angleCID[ind]==chainI)
      break;
  if (ind != angleID.size())
    return ind;

  angleID.push_back(resI);
  angleCID.push_back(chainI);

  rvec x1,x2,x3,x4;

  //calculate psi: N-Ca-C-N
  setRvec(x1,x[search(atomI,resI,"N")  ]);
  setRvec(x2,x[search(atomI,resI,"CA") ]);
  setRvec(x3,x[atomI                   ]);
  setRvec(x4,x[search(atomI,resI+1,"N")]);
  psi.push_back(calcDihedral(x1,x2,x3,x4));

  //calculate phi: C-N-Ca-C
  setRvec(x1,x[search(atomI,resI-1,"C")]);
  setRvec(x2,x[search(atomI,resI,"N")  ]);
  setRvec(x3,x[search(atomI,resI,"CA") ]);
  setRvec(x4,x[atomI                   ]);
  phi.push_back(calcDihedral(x1,x2,x3,x4));

  return angleID.size()-1;
}

uint CalcW::search(const int st, const int resI, const string &whichtype) {
  if (resI>resnum[st]) {  //search forward
    for (int ii=st+1; ii<natoms; ii++)
      if (type[ii].compare(whichtype)==0 && chain[ii]==chain[st])
	return ii;
  } else {             //search backward
    for (int ii=st-1; ii>=0; ii--)
      if (type[ii].compare(whichtype)==0 && chain[ii]==chain[st])
	return ii;
  }

  //didn't find type
  printf("ERROR: failed to find dihedral. Might be at end of chain.\n");
  exit(EXIT_FAILURE);
}

float CalcW::calcDihedral(const rvec &x1, const rvec &x2,
			  const rvec &x3, const rvec &x4) {
  rvec b1,b2,b3;
  addRvec(x2,x1,b1,-1);
  addRvec(x3,x2,b2,-1);
  addRvec(x4,x3,b3,-1);

  //just for generality, even though in KcsA case protein is in center of box
  pbc(b1,box);
  pbc(b2,box);
  pbc(b3,box);

  normalize(b1);
  normalize(b2);
  normalize(b3);

  rvec n1,n2;
  cross(b1,b2,n1);
  cross(b2,b3,n2);

  //see  https://math.stackexchange.com/a/47084
  rvec m1;
  cross(n1,b2,m1);

  float d1,d2;
  d1=dot(n1,n2);
  d2=dot(m1,n2);

  //this returns answer in -180 to 180
  float ans = atan2(d2,d1)*180.0/PI;
  return ans;

  //this returns angle between 0 and 180;
  //return acosd(dot(n1,n2));
}

void CalcW::normalize(rvec &v) {
  float norm=sqrt(norm2vec(v));
  multRvec(v,1.0/norm);
}

float CalcW::interp2(const float &x, const float &y, const float z[]) {
  int nx=(int) (x+180.0)/dTheta;
  int ny=(int) (y+180.0)/dTheta;

  //account for case when x or y == +180
  if (nx==nTheta-1) nx--;
  if (ny==nTheta-1) ny--;

  if (nx<0 || nx>=nTheta || ny<0 || ny>=nTheta) {
    printf("ERROR: Ramachandran angles out of bounds.\n");
    exit(EXIT_FAILURE);
  }

  //corner points surrounding query point
  float xl,xh,yl,yh;
  xl=nx*dTheta-180.0;
  xh=xl+dTheta;
  yl=ny*dTheta-180.0;
  yh=yl+dTheta;

  //shift values at corners
  float z1,z2,z3,z4,u,t;
  z1=z[ny    *nTheta + nx  ];
  z2=z[(ny+1)*nTheta + nx  ];
  z3=z[(ny+1)*nTheta + nx+1];
  z4=z[ny    *nTheta + nx+1];

  u=(x-xl)/(xh-xl);
  t=(y-yl)/(yh-yl);

  return (1-t)*(1-u)*z1 + t*(1-u)*z2 + t*u*z3 + (1-t)*u*z4;
}
