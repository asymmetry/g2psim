/*
This function is provided to calculate quasi-elastic and resonance cross section for different kinds of particles
 Use EPC model for proton, neutron, pions
 Use QFS model for electron
 Do not support positron any more
 Use a subroutine to calculate photons
  The photon cross section is calculated from the pi0 decay so it could not be used with pi0 simultaneously
 
Usage: There are several forms of this getQElXS function
 getQElXS(PID,Z,N,Eb,theta,pf);
 getQElXS(PID,Z,N,Eb,theta,pf,EPS,EPSD,FP);
 getQElXS(PID,Z,N,Eb,theta,pf,Tb,Ta);
 getQElXS(PID,Z,N,Eb,theta,pf,EPS,EPSD,FP,Tb,Ta);
 
Meaning of parameters:
 PID: praticle ID, following the PDG definition
     =	2212	for p		;	2112	for n	;	211	for pi+	;
		-211	for pi-		;	111		for pi0	;	11	for e-	;
		22		for photon	;
 Z,N: proton and neutron number of the nucleus.
 Eb: incoming electron energy in GeV;
 theta: scattering angle for outgoing particle in radian;
 pf: outgoing particle momentum in GeV/c;
 
 EPS,EPSD,FP: nucleus parameters, which used in the QFS code. EPS - seperation energy in MeV, EPSD - delta seperation energy in MeV, FP - Fermi momentum in MeV/c; their values can be changed to fit experimetal results, they will be set with the recommended value if not provided (EPS=10, EPSD=-10, FP=220);
 Tb,Ta: target parameters, Tb - total radiative length before scattering, Ta - total radiative length after scattering, both in the unit of Radiation Length, if they are not provided, they will be set with 0, which means not taking target radiative length into account;
*/


extern "C" {
void epc_(int *, int *, int *, double *, double *, double *, double *);
void qfs_(int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
}

void photon(int tgt_Z,int tgt_N,double E1,double PTP,double THP,double *xs);

double getQElXS(int PID, int Z, int N, double Eb, double theta, double pf)
{
  double EPS,EPSD,FP;
  double Tb,Ta;
  double result;
  double pi;
  int A;

  pi=3.141592654;

  A=Z+N;

  theta=theta*180.0/pi;
  Eb=Eb*1000.0;
  pf=pf*1000.0;
  EPS=10;
  EPSD=-10;
  FP=220;
  if((Z+N)==1){EPS=0;EPSD=0;}
  Tb=0.0;
  Ta=0.0;
  
  if((PID==2212)||(PID==2112)||(PID==211)||(PID==-211)||(PID==111))
    epc_(&PID,&Z,&N,&Eb,&pf,&theta,&result);
  else if(PID==11)
    qfs_(&Z,&A,&Eb,&pf,&theta,&result,&EPS,&EPSD,&FP,&Tb,&Ta);
  else if(PID==22)
    photon(Z,N,Eb,pf,theta,&result);
  else
    result=-1.0;

  return result;
}

double getQElXS(int PID, int Z, int N, double Eb, double theta, double pf, double EPS, double EPSD, double FP)
{
  double Tb,Ta;
  double result;
  double pi;
  int A;

  pi=3.141592654;

  A=Z+N;

  theta=theta*180.0/pi;
  Eb=Eb*1000.0;
  pf=pf*1000.0;
  Tb=0.0;
  Ta=0.0;
  
  if((PID==2212)||(PID==2112)||(PID==211)||(PID==-211)||(PID==111))
    epc_(&PID,&Z,&N,&Eb,&pf,&theta,&result);
  else if(PID==11)
    qfs_(&Z,&A,&Eb,&pf,&theta,&result,&EPS,&EPSD,&FP,&Tb,&Ta);
  else if(PID==22)
    photon(Z,N,Eb,pf,theta,&result);
  else
    result=-1.0;

  return result;
}

double getQElXS(int PID, int Z, int N, double Eb, double theta, double pf, double Tb, double Ta)
{
  double EPS,EPSD,FP;
  double result;
  double pi;
  int A;

  pi=3.141592654;

  A=Z+N;

  theta=theta*180.0/pi;
  Eb=Eb*1000.0;
  pf=pf*1000.0;
  EPS=10;
  EPSD=-10;
  FP=220;
  if((Z+N)==1){EPS=0;EPSD=0;}
  
  if((PID==2212)||(PID==2112)||(PID==211)||(PID==-211)||(PID==111))
    epc_(&PID,&Z,&N,&Eb,&pf,&theta,&result);
  else if(PID==11)
    qfs_(&Z,&A,&Eb,&pf,&theta,&result,&EPS,&EPSD,&FP,&Tb,&Ta);
  else if(PID==22)
    photon(Z,N,Eb,pf,theta,&result);
  else
    result=-1.0;

  return result;
}

double getQElXS(int PID, int Z, int N, double Eb, double theta, double pf, double EPS, double EPSD, double FP, double Tb, double Ta)
{
  double result;
  double pi;
  int A;

  pi=3.141592654;

  A=Z+N;

  theta=theta*180.0/pi;
  Eb=Eb*1000.0;
  pf=pf*1000.0;
  
  if((PID==2212)||(PID==2112)||(PID==211)||(PID==-211)||(PID==111))
    epc_(&PID,&Z,&N,&Eb,&pf,&theta,&result);
  else if(PID==11)
    qfs_(&Z,&A,&Eb,&pf,&theta,&result,&EPS,&EPSD,&FP,&Tb,&Ta);
  else if(PID==22)
    photon(Z,N,Eb,pf,theta,&result);
  else
    result=-1.0;

  return result;
}
