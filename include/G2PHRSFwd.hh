// -*- C++ -*-

/* class G2PHRSFwd
 * It simulates the movement of the scatted particles in the spectrometers.
 * G2PDrift, G2PHRS and G2PSieve are used in this class.
 * Input variables: fV5tg_tr, fV5react_lab (register in gG2PVars).
 */

// History:
//   Apr 2013, C. Gu, First public version.
//   Oct 2013, J. Liu, Add Energy loss and Multiple scattering.
//

#ifndef G2P_HRSFWD_H
#define G2P_HRSFWD_H

#include "G2PProcBase.hh"

class G2PDrift;
class G2PSieve;
class HRSTransBase;

class G2PHRSFwd : public G2PProcBase {
public:
    G2PHRSFwd(const char* name);
    virtual ~G2PHRSFwd();

    virtual int Init();
    virtual int Begin();
    virtual int Process();
    virtual void Clear(Option_t* /*option*/ = "");

    // Gets

    // Sets
    void SetSieve(const char* opt);

protected:
    G2PHRSFwd(); // Only for ROOT I/O

    void RunType10(double* V5react_tr, double& z_tr, double* V5troj_tr, double& dlentot, double& elosstot); // production target
    void RunType20(double thickness, double* V5react_tr, double& z_tr, double* V5troj_tr, double& dlentot, double& elosstot); // carbon target, without LHe
    void RunType21(double thickness, double* V5react_tr, double& z_tr, double* V5troj_tr, double& dlentot, double& elosstot); // carbon target, with LHe
    void RunType31(double* V5react_tr, double& z_tr, double* V5troj_tr, double& dlentot, double& elosstot); // pure LHe

    bool Forward(const double* V5tp_tr, double* V5fp_tr);

    void ApplyVDCRes(double* V5fp_tr);

    virtual int Configure(EMode mode = kTWOWAY);
    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    int fRunType;
    double fHRSAngle;
    double fHRSMomentum;

    int fSetting;

    bool fSieveOn;
    int fHoleID;

    double fV5react_lab[5];
    double fV5react_tr[5];

    double fV5sieve_tr[5];
    double fV5tpproj_tr[5];

    double fV5fp_tr[5];
    double fV5fp_rot[5];

    G2PDrift* pDrift;
    G2PSieve* pSieve;
    HRSTransBase* pModel;

private:
    static G2PHRSFwd* pG2PHRSFwd;

    ClassDef(G2PHRSFwd, 1)
};

#endif
