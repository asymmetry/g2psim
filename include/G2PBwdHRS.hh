// -*- C++ -*-

/* class G2PBwdHRS
 * It simulates the reconstruction of g2p kinematics.
 * G2PDrift, G2PHRS and G2PGeoSieve are used in this class.
 * Input variables: fV5bpm_lab, fV5fp_tr (register in gG2PVars).
 */

// History:
//   Apr 2013, C. Gu, First public version.
//   Apr 2014, C. Gu, New effective bpm fitting.
//

#ifndef G2P_BWDHRS_H
#define G2P_BWDHRS_H

#include "G2PProcBase.hh"

class G2PSieve;
class HRSTransBase;

class G2PBwdHRS : public G2PProcBase
{
public:
    G2PBwdHRS(const char *name);
    virtual ~G2PBwdHRS();

    virtual int Begin();
    virtual int Process();
    virtual void Clear(Option_t *opt = "");

    // Gets

    // Sets
    void SetParsX(const double *pars);
    void SetParsY(const double *pars);
    void SetRecZ(double z);

protected:
    G2PBwdHRS(); // Only for ROOT I/O

    double GetEffBPM(int axis);

    bool Backward(const double *V5fp_tr, double *V5tp_tr);

    virtual int Configure(EMode mode = kTWOWAY);
    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    double fFieldRatio;

    int fSetting;

    double fFitPars[2][3];
    double frecz_lab;

    double fV5bpm_lab[5];
    double fV5bpm_tr[5];

    double fV5fp_tr[5];

    double fV5tpsnake_tr[5];
    double fV5sieveproj_tr[5];

    double fV5tprec_tr[5];
    double fV5tprec_lab[5];

    G2PSieve *pSieve;
    HRSTransBase *pModel;

private:
    static G2PBwdHRS *pG2PBwdHRS;

    ClassDef(G2PBwdHRS, 1)
};

#endif
