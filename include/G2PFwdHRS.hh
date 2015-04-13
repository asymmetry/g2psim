// -*- C++ -*-

/* class G2PFwdHRS
 * It simulates the movement of the scatted particles in the spectrometers.
 * G2PDrift, G2PHRS and G2PGeoSieve are used in this class.
 * Input variables: fV5tg_tr, fV5react_lab (register in gG2PVars).
 */

// History:
//   Apr 2013, C. Gu, First public version.
//   Oct 2013, J. Liu, Add Energy loss and Multiple scattering.
//   Jan 2015, C. Gu, Rewrite with geometry classes.
//

#ifndef G2P_FWDHRS_H
#define G2P_FWDHRS_H

#include "G2PProcBase.hh"

class G2PSieve;
class HRSTransBase;

class G2PFwdHRS : public G2PProcBase
{
public:
    G2PFwdHRS(const char *name);
    virtual ~G2PFwdHRS();

    virtual int Begin();
    virtual int Process();
    virtual void Clear(Option_t *opt = "");

    // Gets

    // Sets
    void SetSieve(const char *opt);
    void SetVDCRes(double x, double t, double y, double p);

protected:
    G2PFwdHRS(); // Only for ROOT I/O

    bool Forward(const double *V5tp_tr, double *V5fp_tr);

    void ApplyVDCRes(double *V5fp_tr);

    virtual int Configure(EMode mode = kTWOWAY);
    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    int fSetting;

    bool fSieveOn;
    int fHoleID;

    int fEndPlane;

    double fELoss;
    double fTa;

    double fVDCRes[4];

    double fV5beam_lab[5];

    double fV5react_lab[5];
    double fV5react_tr[5];
    double freactz_tr;

    double fV5sieve_tr[5];
    double fV5tpproj_tr[5];

    double fV5fp_tr[5];
    double fV5fp_rot[5];
    double fPlanePosX[30];
    double fPlanePosY[30];
    double fDumpFront[2];
    double fDumpBack[2];
    G2PSieve *pSieve;
    HRSTransBase *pModel;

private:
    static G2PFwdHRS *pG2PFwdHRS;

    ClassDef(G2PFwdHRS, 1)
};

#endif
