// -*- C++ -*-

/* class G2PPhys
 * It will calculate cross sections at reaction point.
 *
 * Meaning of parameters:
 * fPID: incident particle ID, following the PDG definition:
 *       2212 for p        ;   2112 for n     ;   211 for pi+   ;
 *       -211 for pi-      ;   111  for pi0   ;   11  for e-    ;
 *       22   for photon   ;
 * fZ, fA: proton and mass number of the nucleus.
 *
 * Radiative correction parameters:
 * Tb: total radiative length before scattering in radiation length;
 * Ta: total radiative length after scattering in radiation length;
 * If they are not provided, they will be set to 0, which means not taking target radiative length into account;
 *
 * QFS model parameters:
 * EPS: separation energy in MeV;
 * EPSD: delta separation energy in MeV;
 * FP - Fermi momentum in MeV/c;
 * If they are not provided, they will be set to the recommended value (EPS=10, EPSD=-10, FP=220);
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Apr 2013, C. Gu, Add P. Bosted's model.
//   May 2013, C. Gu, Add L. Cardman's C12 elastic model.
//

#ifndef G2P_PHYS_H
#define G2P_PHYS_H

#include "G2PProcBase.hh"

class G2PPhysBase;

class G2PPhys : public G2PProcBase {
public:
    G2PPhys(const char *name);
    virtual ~G2PPhys();

    virtual int Begin();
    virtual int Process();
    virtual void Clear();

    // Gets

    // Sets
    void SetPars(double* array, int n);

protected:
    G2PPhys(); // Only for ROOT I/O

    double CalXS(const double* V5lab, const double* V5tr, double& scatangle);

    virtual int Configure(EMode mode = kTWOWAY);
    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    int fSetting;

    int fPID;

    int fZ, fA; // Define Target
    double fTargetMass;

    double* fPars;
    int fNPars;

    double fHRSAngle;
    double fHRSMomentum;

    double fBeamEnergy;

    double fXSreact;
    double fTHreact;
    double fXSrec;
    double fTHrec;

    G2PPhysBase* pModel;

private:
    static G2PPhys* pG2PPhys;

    ClassDef(G2PPhys, 1)
};

// inline functions

inline void G2PPhys::SetPars(double* array, int n) {
    fPars = array;
    fNPars = n;
}

#endif
