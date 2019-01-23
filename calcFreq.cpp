#include "calcFreq.h"

CalcFreq::CalcFreq(const int nchrom, const Charges &chg, const GroFile &gro, const vector<int> &grpSt, const int *grpInd, const int nCutGrp) :
  nchrom(nchrom), nham(nchrom*(nchrom+1)/2), natoms(gro.getNatom()),
  type(gro.getType()), backbone(gro.getBackbone()), chain(gro.getChain()),
  resnumAll(gro.getResNumAll()), q(chg.getCharges()),
  grpSt(grpSt), grpInd(grpInd), nCutGrp(nCutGrp)
{
  ham     = new float[nham];
  dip     = new rvec[nchrom];
}

CalcFreq::~CalcFreq() {
  delete[] ham;
  delete[] dip;
}

//needed due to a flaw in c++11
constexpr float CalcFreq::NtermShift[169];
constexpr float CalcFreq::CtermShift[169];

void CalcFreq::compute(const Traj &traj, const vector<int> &inds) {
  const rvec  *x=traj.getCoords();
  traj.getBox(box);
  ts=traj.getNT()-1;
  int ii,jj,kk;
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
  rvec *cog = new rvec[nCutGrp];
  calcCOG(x,cog);

  //save angles so don't re-compute
  vector<int> angleID;
  vector<float> phi;
  vector<float> psi;

  uint          nnC[3];  //list of C atoms in two nn amides, and self
  rvec vecCO, vecCN, tmpvec, Cpos, Npos, tmpEn, tmpEc;
  float d2,d,d2n;
  uint thisC,thisRes,a1,a2,nextC,hii;
  rvec thisCcog,thisNcog;
  rvec *tdPos = new rvec[nchrom];
  //loop through labelled CO and transition dipole
  for (ii=0; ii<nchrom; ii++) {
    thisC=inds[ii];
    copyRvec(x[thisC],Cpos);
    addRvec(x[thisC+1],Cpos,vecCO,-1); //points towards O (rO-rC)
    pbc(vecCO,box);
    normalize(vecCO);

    copyRvec(x[thisC+2],Npos);
    addRvec(Npos,Cpos,vecCN,-1); //points towards N (rN-rC)
    pbc(vecCN,box);
    normalize(vecCN);

    calcTD(vecCO,vecCN,tmpvec);
    copyRvec(tmpvec,dip[ii]);
    calcTDpos(Cpos,vecCO,vecCN,tdPos[ii]);

    //include NN peptide shifts
    thisRes=resnumAll[thisC];
    a1=calcAngles(thisC,thisRes,x,angleID,phi,psi);
    //calculate angles for next resiude
    nextC=search(thisC,thisRes+1,"C");
    a2=calcAngles(nextC,thisRes+1,x,angleID,phi,psi);

    //for excluding previous amide bond
    nnC[0]=search(thisC,thisRes-1,"C");
    nnC[1]=thisC;
    nnC[2]=nextC;
    vector<uint> excludeGrp = getExcludes(nnC);

    //tested 1/21/19 that this gives same shifts as Jansen's code
    float nnfs = interp2(phi[a1],psi[a1],NtermShift);
    nnfs      += interp2(phi[a2],psi[a2],CtermShift);

    //loop through other atoms
    copyRvec(cog[grpInd[thisC]],thisCcog);
    copyRvec(cog[grpInd[thisC+2]],thisNcog);
    setRvec(tmpEn,0.0);
    setRvec(tmpEc,0.0);
    for (jj=0; jj<nCutGrp; jj++) {
      //self and NN backbones are excluded, but side chains included
      //see Lin JCP 113 2009
      if (std::find(excludeGrp.begin(), excludeGrp.end(), jj) !=
	  excludeGrp.end())
	continue;

      //get distance from COG to COG
      addRvec(thisCcog,cog[jj],tmpvec,-1);
      pbc(tmpvec,box);
      d2=norm2vec(tmpvec);
      if (d2 < cutExt2) {  //test if distance is within C or N of current amideI

	//get dist to N too
	addRvec(thisNcog,cog[jj],tmpvec,-1);
	pbc(tmpvec,box);
	d2n=norm2vec(tmpvec);

	//lots of duplicated code here, to avoid separate loop for each case
	if (d2 < cut2 && d2n < cut2) { //both C and N are in cutoff
	  for (kk=grpSt[jj]; kk<grpSt[jj+1]; kk++) {
	    if (q[kk]) {
	      //calc Ec
	      addRvec(Cpos,x[kk],tmpvec,-1);
	      pbc(tmpvec,box);
	      d2=norm2vec(tmpvec);
	      d=sqrt(d2);
	      multRvec(tmpvec, q[kk]/(d2*d) );
	      addRvec(tmpvec,tmpEc,+1);

	      //calc En
	      addRvec(Npos,x[kk],tmpvec,-1);
	      pbc(tmpvec,box);
	      d2=norm2vec(tmpvec);
	      d=sqrt(d2);
	      multRvec(tmpvec, q[kk]/(d2*d) );
	      addRvec(tmpvec,tmpEn,+1);
	    }
	  }
	} else if (d2 < cut2) {  //just C in cutoff
	  for (kk=grpSt[jj]; kk<grpSt[jj+1]; kk++) {
	    if (q[kk]) {
	      //calc Ec
	      addRvec(Cpos,x[kk],tmpvec,-1);
	      pbc(tmpvec,box);
	      d2=norm2vec(tmpvec);
	      d=sqrt(d2);
	      multRvec(tmpvec, q[kk]/(d2*d) );
	      addRvec(tmpvec,tmpEc,+1);
	    }
	  }
	} else if (d2n < cut2) { //just N in cutoff
	  for (kk=grpSt[jj]; kk<grpSt[jj+1]; kk++) {
	    if (q[kk]) {
	      //calc En
	      addRvec(Npos,x[kk],tmpvec,-1);
	      pbc(tmpvec,box);
	      d2=norm2vec(tmpvec);
	      d=sqrt(d2);
	      multRvec(tmpvec, q[kk]/(d2*d) );
	      addRvec(tmpvec,tmpEn,+1);
	    }
	  }
	}
      }
    }
    hii=ii*nchrom-(ii-1)*ii/2;
    ham[hii] = map(dot(tmpEc,vecCO), dot(tmpEn,vecCO)) + nnfs;
  }
  delete[] cog;

  //compute couplings
  float tmpcoup;
  for (ii=0; ii<nchrom-1; ii++) {
    for (jj=ii+1; jj<nchrom; jj++) {
      addRvec(tdPos[jj],tdPos[ii],tmpvec,-1);
      pbc(tmpvec,box);
      d2=norm2vec(tmpvec);
      d=sqrt(d2);

      tmpcoup=dot(dip[ii],dip[jj]) - 3*dot(dip[ii],tmpvec)*dot(dip[jj],tmpvec);

      hii=ii*nchrom - (ii-1)*ii/2 + jj - ii;
      ham[hii] = coupConst * tmpcoup / (d2*d2*d);
    }
    multRvec(dip[ii],tdMag);
  }
  multRvec(dip[ii],tdMag);
  delete[] tdPos;
}

void CalcFreq::print(FILE *fFreq, FILE *fDip) {
  //frequencies
  fprintf(fFreq,"%d",ts);
  for (int ii=0; ii<nham; ii++)
    fprintf(fFreq," %f",ham[ii]);
  fprintf(fFreq,"\n");

  //dipole moments (order: all x components, all y comp, all z comp)
  fprintf(fDip,"%d",ts);
  for (int kk=0; kk<DIM; kk++)
    for (int ii=0; ii<nchrom; ii++)
      fprintf(fDip," %f",dip[ii][kk]*tdMag);
  fprintf(fDip,"\n");
}

uint CalcFreq::calcAngles(const int atomI, const int resI, const rvec *x,vector<int> &angleID, vector<float> &phi, vector<float> &psi) {
  //check if already calculated this angle
  uint ind;
  for (ind=0; ind<angleID.size(); ind++)
    if (angleID[ind]==resI)
      break;
  if (ind != angleID.size())
    return ind;

  angleID.push_back(resI);

  int i1,i2,i3,i4;

  //calculate psi: N-Ca-C-N
  i1=search(atomI,resI,"N");
  i2=search(atomI,resI,"CA");
  i3=atomI;
  i4=search(atomI,resI+1,"N");
  psi.push_back(calcDihedral(x[i1],x[i2],x[i3],x[i4]));

  //calculate phi: C-N-Ca-C
  i1=search(atomI,resI-1,"C");
  i2=search(atomI,resI,"N");
  i3=search(atomI,resI,"CA");
  i4=atomI;
  phi.push_back(calcDihedral(x[i1],x[i2],x[i3],x[i4]));

  return angleID.size()-1;
}

uint CalcFreq::search(const int st, const int resI, const string &whichtype) {
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

float CalcFreq::calcDihedral(const rvec &x1, const rvec &x2,
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

inline void CalcFreq::normalize(rvec &v) {
  float norm=sqrt(norm2vec(v));
  multRvec(v,1.0/norm);
}

float CalcFreq::interp2(const float &x, const float &y, const float z[]) {
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

void CalcFreq::calcTD(const rvec &rCO, const rvec &rCN,rvec &out) {
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
}

void CalcFreq::calcTDpos(const rvec &Cpos, const rvec &vecCO, const rvec &vecCN,
			 rvec &tdPos) {
  copyRvec(Cpos,tdPos);
  addRvec(vecCO,tdPos,0.665);
  addRvec(vecCN,tdPos,0.258);
}


void CalcFreq::calcCOG(const rvec *x,rvec *cog) {
  int jj,nthis,thisSt;
  rvec tmpcog,tmpdiff,tmpvec,xref;
  for (int ii=0; ii<nCutGrp; ii++) {
    thisSt=grpSt[ii];
    nthis=grpSt[ii+1]-thisSt;
    copyRvec(x[thisSt],tmpcog);
    copyRvec(tmpcog,xref);
    for (jj=1; jj<nthis; jj++) {
      copyRvec(x[thisSt+jj],tmpvec);
      addRvec(tmpvec,xref,tmpdiff,-1.0);
      pbcOther(tmpvec,tmpdiff,box);   //get image closest to 1st point
      addRvec(tmpvec,tmpcog,1.0);
    }
    multRvec(tmpcog,1.0/nthis);
    copyRvec(tmpcog,cog[ii]);
  }
}

vector<uint> CalcFreq::getExcludes(const uint nnC[3]) {
  uint ii,jj,tmpExcl;
  vector<uint> excludeGrp;
  for (ii=0; ii<3; ii++) { //loop through neighboring Cs
    for (jj=0; jj<5; jj++) {  //loop through atoms in amide I bond to exclude
      tmpExcl=grpInd[nnC[ii]+jj];
      //if grp not already in list, add it
      if (std::find(excludeGrp.begin(), excludeGrp.end(), tmpExcl) ==
	  excludeGrp.end())
	excludeGrp.push_back(tmpExcl);
    }
  }

  return excludeGrp;
}
