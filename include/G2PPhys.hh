// -*- C++ -*-

/* class G2PPhys
 * This file defines a class G2PPhys.
 * It is the interface class of G2PPhys package. It provides physics models.
 * G2PProcBase classes will call GetXS() to calculate cross sections.
 *
 * Meaning of parameters:
 * iPID: incident particle ID, following the PDG definition:
 *       2212 for p        ;   2112 for n     ;   211 for pi+   ;
 *       -211 for pi-      ;   111  for pi0   ;   11  for e-    ;
 *       22   for photon   ;
 * iZ, iA: proton and mass number of the nucleus.
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

#include <vector>

#include "G2PAppBase.hh"

class G2PPhysBase;

class G2PPhys : public G2PAppBase {
public:
    G2PPhys(const char *name);
    ~G2PPhys();

    void SetPars(double* array, int n) {
        fPars = array;
        nPars = n;
    }

    int Init();
    int Begin();

    double GetXS(double Eb, double Ef, double theta);

    static G2PPhys* GetInstance() {
        return pG2PPhys;
    }

protected:
    G2PPhys();

    int iSetting;

    int iZ, iA; // Define Target
    double fTargetMass;

    int iPID;

    double* fPars;
    int nPars;

    G2PPhysBase* pModel;

private:
    static G2PPhys* pG2PPhys;

    ClassDef(G2PPhys, 1)
};

#endif
