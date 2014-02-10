// -*- C++ -*-

/* class G2PFwdProc
 * It simulates the movement of the scatted particles in the spectrometers.
 * G2PDrift, G2PHRS and G2PSieve are used in this class.
 * Input variables: fV5tg_tr, fV5react_lab (register in gG2PVars).
 */

// History:
//   Apr 2013, C. Gu, First public version.
//   Oct 2013, J. Liu, Add Energy loss and Multiple scattering.
//

#ifndef G2P_FWDPROC_H
#define G2P_FWDPROC_H

#include "G2PProcBase.hh"

class G2PDrift;
class G2PHRS;
class G2PSieve;
class G2PMaterial;

class G2PFwdProc : public G2PProcBase {
public:
    G2PFwdProc();
    virtual ~G2PFwdProc();

    virtual int Init();
    virtual int Begin();
    virtual int Process();
    virtual void Clear(Option_t* /*option*/ = "");

protected:
    void RunType10(double* V5react_tr, double& z_tr, double* V5troj, double& dlentot, double& elosstot); // production target
    void RunType20(double thickness, double* V5react_tr, double& z_tr, double* V5troj, double& dlentot, double& elosstot); // carbon target, without LHe
    void RunType21(double thickness, double* V5react_tr, double& z_tr, double* V5troj, double& dlentot, double& elosstot); // carbon target, with LHe

    void ApplyVDCRes(double* V5fp);

    virtual int Configure(EMode mode = kTWOWAY);
    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    int fRunType;
    double fHRSAngle;
    double fHRSMomentum;

    bool fSieveOn;
    int fHoleID;

    double fV5sieve_tr[5];
    double fV5tpproj_tr[5];

    double fV5fp_tr[5];
    double fV5fp_rot[5];

    G2PDrift* pDrift;
    G2PHRS* pHRS;
    G2PSieve* pSieve;

private:
    static G2PFwdProc* pG2PFwdProc;

    ClassDef(G2PFwdProc, 1)
};

#endif
