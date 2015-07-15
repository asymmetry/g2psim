// -*- C++ -*-

/* class G2PRec
 * The real class to do the reconstruction.
 */

// History:
//   Apr 2013, C. Gu, First public version.
//   Feb 2014, C. Gu, Modified for G2PRec.
//   Jan 2015, C. Gu, Move this file back to the g2psim package.
//

#ifndef G2P_REC_H
#define G2P_REC_H

#include "G2PProcBase.hh"

class G2PSieve;

class G2PRec : public G2PProcBase
{
public:
    G2PRec();
    virtual ~G2PRec();

    virtual int Begin();
    virtual int Process();
    virtual void Clear(Option_t *opt = "");

    // Gets

    // Sets

protected:
    void GetEffBPM(const double *V5tp_tr, const double *V5bpm_tr, double *V5bpmeff_tr);
    void Correct(const double *V5bpm_tr, const double *V5tp_tr, double *V5corr_tr);

    virtual int Configure(EMode mode = kTWOWAY);
    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    double fE0;
    double fFieldRatio;

    double frecz_lab;

    double fTgtXCorrT;
    double fTgtXCorrP;
    double fTgtXCorrD;

    double fTgtYCorrT;
    double fTgtYCorrP;
    double fTgtYCorrD;

    double fConsCorrT;
    double fConsCorrP;
    double fConsCorrD;

    double fFitPars[2][3];

    double fV5bpm_bpm[5];
    double fV5bpm_tr[5];
    double fbpmz_tr;

    double fV5tpmat_tr[5];
    double fV5tpcorr_tr[5];
    double fV5sieveproj_tr[5];

    double fV5tprec_tr[5];
    double fV5tprec_lab[5];

    G2PSieve *pSieve;

private:
    static G2PRec *pG2PRec;

    ClassDef(G2PRec, 1)
};

#endif
