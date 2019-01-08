#include "calcW.h"


CalcW::CalcW(const int nchrom, const Charges &chg, const GroFile &gro) :
  nchrom(nchrom), natoms(gro.getNatom()), type(gro.getType()),
  backbone(gro.getBackbone()), chain(gro.getChain()), nres(gro.getNres()),
  resnumAll(gro.getResNumAll()), atomsInRes(gro.getAtomsInRes()),
  resSt(gro.getResSt()), q(chg.getCharges())
{
  freq   = new float[nchrom];
  dip     = new rvec[nchrom];
}

CalcW::~CalcW() {
  delete[] freq;
  delete[] dip;
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

  //get centers of geometry
  rvec *cog = new rvec[nres];
  calcCOG(x,cog);

  rvec vecCO, vecCN, tmpvec, Cpos, Npos, tmpEn, tmpEc;
  rvec tmpcog;
  float d2,d;
  int thisAtom,thisRes,nnL,nnR;
  int kk;
  bool excludeBackbone;
  //loop through labelled CO and transition dipole
  for (ii=0; ii<nchrom; ii++) {
    thisAtom=inds[ii];
    copyRvec(x[thisAtom],Cpos);
    addRvec(x[thisAtom+1],Cpos,vecCO,-1); //points towards O (rO-rC)
    pbc(vecCO,box);
    normalize(vecCO);

    copyRvec(x[thisAtom+2],Npos);
    addRvec(Npos,Cpos,vecCN,-1); //points towards N (rN-rC)
    pbc(vecCN,box);
    normalize(vecCN);

    calcTD(vecCO,vecCN,tmpvec);
    copyRvec(tmpvec,dip[ii]);

    copyRvec(cog[resnumAll[ii]],tmpcog);

    //include NN peptide shifts
    float nnfs = 0.0;
    thisRes=resnumAll[thisAtom];
    nnL=thisRes-1;
    nnR=thisRes+1;
    uint a1,a2,Cnext;
    a1=calcAngles(thisAtom,thisRes,x);
    //calculate angles for next resiude
    Cnext=search(thisAtom,nnR,"C");
    a2=calcAngles(Cnext,nnR,x);

    nnfs += interp2(phi[a1],psi[a1],NtermShift);
    nnfs += interp2(phi[a2],psi[a2],CtermShift);

    //loop through other atoms
    setRvec(tmpEn,0.0);
    setRvec(tmpEc,0.0);
    excludeBackbone=false;
    for (jj=0; jj<nres; jj++) {
      //self and NN backbones are excluded, but side chains included
      //see Lin JCP 113 2009
      if ( jj==nnL || jj==nnR || jj==thisRes )
	excludeBackbone=true;

      //TODO: which atom is the cutoff with respect to?
      //get distance to C atom
      //addRvec(Cpos,x[jj],tmpvec,-1);
      //get distance to COG
      addRvec(tmpcog,cog[jj],tmpvec,-1);
      pbc(tmpvec,box);
      d2=norm2vec(tmpvec);
      if (d2 < cut2) {
	for (kk=resSt[jj]; kk<resSt[jj]+atomsInRes[jj]; kk++) {
	  if (q[kk]) {
	    if (backbone[kk] && excludeBackbone)
	      continue;

	    //calc Ec
	    addRvec(Cpos,x[kk],tmpvec,-1);
	    pbc(tmpvec,box);
	    d2=norm2vec(tmpvec);
	    d=sqrt(d2);
	    multRvec(tmpvec, q[kk]/(d*d*d) );
	    addRvec(tmpvec,tmpEc,+1);

	    //calc En
	    addRvec(Npos,x[kk],tmpvec,-1);
	    pbc(tmpvec,box);
	    d2=norm2vec(tmpvec);
	    d=sqrt(d2);
	    multRvec(tmpvec, q[kk]/(d*d*d) );
	    addRvec(tmpvec,tmpEn,+1);
	  }
	}
      }
    }
    freq[ii] = map(dot(tmpEc,vecCO), dot(tmpEn,vecCO)) + nnfs;
  }

  delete[] cog;
}

void CalcW::print(FILE *fFreq, FILE *fDip) {
  //frequencies
  fprintf(fFreq,"%d",ts);
  for (int ii=0; ii<nchrom; ii++)
    fprintf(fFreq," %f",freq[ii]);
  fprintf(fFreq,"\n");

  //dipole moments (order: all x components, all y comp, all z comp)
  fprintf(fDip,"%d",ts);
  for (int kk=0; kk<DIM; kk++)
    for (int ii=0; ii<nchrom; ii++)
      fprintf(fDip," %f",dip[ii][kk]);
  fprintf(fDip,"\n");
}

uint CalcW::calcAngles(const int atomI, const int resI, const rvec *x) {
  //check if already calculated this angle
  uint ind;
  for (ind=0; ind<angleID.size(); ind++)
    if (angleID[ind]==resI)
      break;
  if (ind != angleID.size())
    return ind;

  angleID.push_back(resI);

  rvec x1,x2,x3,x4;

  //calculate psi: N-Ca-C-N
  copyRvec(x[search(atomI,resI,"N")  ],x1);
  copyRvec(x[search(atomI,resI,"CA") ],x2);
  copyRvec(x[atomI                   ],x3);
  copyRvec(x[search(atomI,resI+1,"N")],x4);
  psi.push_back(calcDihedral(x1,x2,x3,x4));

  //calculate phi: C-N-Ca-C
  copyRvec(x[search(atomI,resI-1,"C")],x1);
  copyRvec(x[search(atomI,resI,"N")  ],x2);
  copyRvec(x[search(atomI,resI,"CA") ],x3);
  copyRvec(x[atomI                   ],x4);
  phi.push_back(calcDihedral(x1,x2,x3,x4));

  return angleID.size()-1;
}

uint CalcW::search(const int st, const int resI, const string &whichtype) {
  if (resI>resnumAll[st]) {  //search forward
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

inline void CalcW::normalize(rvec &v) {
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

  float ans=(1-t)*(1-u)*z1 + t*(1-u)*z2 + t*u*z3 + (1-t)*u*z4;
  return ans;
}

void CalcW::calcTD(const rvec &rCO, const rvec &rCN,rvec &out) {
  rvec n1,n2;
  cross(rCO,rCN,n1);  //n1: normal to OCN plane
  cross(n1,rCO,n2);   //n2: perpendicular to CO, in direction of N
  normalize(n2);

  copyRvec(n2,out);
  multRvec(out,sin(tdAngle));
  addRvec(rCO,out,-cos(tdAngle)); //reverse CO vector

  normalize(out);

  //test that angle is correct (10deg)
  //printf("angle = %f\n",acos(dot(out,rCO))*180.0/PI);

  multRvec(out,2.73); //magintude of td in D/A/amu^1/2
}

void CalcW::calcCOG(const rvec *x,rvec *cog) {
  int jj,thisSt,nthis;
  rvec tmpcog;
  for (int ii=0; ii<nres; ii++) {
    thisSt=resSt[ii];
    nthis=atomsInRes[ii];
    setRvec(tmpcog,0.0);
    for (jj=0; jj<nthis; jj++)
      addRvec(x[thisSt+jj],tmpcog,1.0);
    multRvec(tmpcog,1.0/nthis);
    copyRvec(tmpcog,cog[ii]);
  }
}
