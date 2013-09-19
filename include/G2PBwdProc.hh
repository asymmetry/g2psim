// -*- C++ -*-

/* class G2PBwdProc
 * It simulates the reconstruction of g2p kinematics.
 * G2PDrift, G2PHRS and G2PSieve are used in this class.
 * Input variables: fV5bpm_lab, fV5fp_tr (register in gG2PVars).
 */

// History:
//   Apr 2013, C. Gu, First public version.
//

#ifndef G2P_BWDPROC_H
#define G2P_BWDPROC_H

#include "G2PProcBase.hh"

class G2PDrift;
class G2PHRS;
class G2PSieve;

class G2PBwdProc : public G2PProcBase {
public:
    G2PBwdProc();
    virtual ~G2PBwdProc();

    virtual int Init();
    virtual int Begin();
    virtual int Process();
    virtual void Clear();

protected:
    double GetEffBPM(double xbpm_tr, const double* V5fp);

    virtual int Configure(EMode mode = kTWOWAY);
    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    double fBeamEnergy;
    double fHRSAngle;
    double fHRSMomentum;
    double fFieldRatio;

    double fV5rectg_tr[5];
    double fV5recsieve_tr[5];

    double fV5rec_tr[5];
    double fV5rec_lab[5];

    G2PDrift* pDrift;
    G2PHRS* pHRS;
    G2PSieve *pSieve;

private:
    static G2PBwdProc* pG2PBwdProc;

    ClassDef(G2PBwdProc, 1)
};

#endif
