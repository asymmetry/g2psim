#ifndef G2P_PHYS_H
#define G2P_PHYS_H

#include <vector>

#include "G2PAppsBase.hh"

class G2PPhysBase;

class G2PPhys : public G2PAppsBase
{
public:
    G2PPhys(const char *name);
    ~G2PPhys();

    void SetPars(double* array, int n) { fPars = array; nPars = n; }

    EStatus Init();
    void Clear() { }

    double GetXS(double Eb, double Ef, double theta);

    int GetSetting() { return iSetting; }

    static G2PPhys* GetInstance() { return pG2PPhys; }

protected:
    G2PPhys();

    int iSetting;

    int iZ, iA; // Define Target
    int iPID;

    double* fPars;
    int nPars;

    G2PPhysBase* pModel;

private:
    static G2PPhys* pG2PPhys;

    ClassDef(G2PPhys, 1)
};

#endif

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
 
 EPS,EPSD,FP: nucleus parameters, which used in the QFS code. EPS - seperation energy in MeV, 
 EPSD - delta seperation energy in MeV, FP - Fermi momentum in MeV/c; their values can be changed 
 to fit experimetal results, they will be set with the recommended value if not provided (EPS=10, EPSD=-10, FP=220);
 Tb,Ta: target parameters, Tb - total radiative length before scattering, Ta - total radiative 
 length after scattering, both in the unit of Radiation Length, if they are not provided, they will 
 be set with 0, which means not taking target radiative length into account;
*/
