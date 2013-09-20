// -*- C++ -*-

/* class G2PFwdProc
 * It simulates the movement of the scatted particles in the spectrometers.
 * G2PDrift, G2PHRS and G2PSieve are used in this class.
 * Input variables: fV5tg_tr, fV5react_lab (register in gG2PVars).
 */

// History:
//   Apr 2013, C. Gu, First public version.
//

#ifndef G2P_FWDPROC_H
#define G2P_FWDPROC_H

#include "G2PProcBase.hh"

class G2PDrift;
class G2PHRS;
class G2PSieve;

class G2PFwdProc : public G2PProcBase {
public:
    G2PFwdProc();
    virtual ~G2PFwdProc();

    virtual int Init();
    virtual int Begin();
    virtual int Process();
    virtual void Clear();

protected:
    void ApplyVDCRes(double* V5fp);

    virtual int Configure(EMode mode = kTWOWAY);
    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    double fHRSAngle;
    double fHRSMomentum;

    bool fSieveOn;

    double fV5sieve_tr[5];
    double fV5tpproj_tr[5];

    double fV5fp_tr[5];
    double fV5fp_rot[5];

    G2PDrift* pDrift;
    G2PHRS* pHRS;
    G2PSieve *pSieve;

private:
    static G2PFwdProc* pG2PFwdProc;

    ClassDef(G2PFwdProc, 1)
};

#endif
