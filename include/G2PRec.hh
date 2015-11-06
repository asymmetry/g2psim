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
    virtual int Process(const double *V5bpm_bpm, const double *V5tpmat_tr, double *V5rec_tr, double *V5rec_lab);
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

    double fFitPars[2][3];

    double fCorT[3];
    double fCorP[3];
    double fCorD[3];

    double fV5bpm_bpm[5];

    double fV5tpmat_tr[5];
    double fV5tpcorr_tr[5];
    double fV5sieveproj_tr[5];

    double fV5tprec_tr[5];
    double fV5tprec_lab[5];

    double fV5rec_tr[5];
    double fV5rec_lab[5];

    G2PSieve *pSieve;

private:
    static G2PRec *pG2PRec;

    ClassDef(G2PRec, 1)
};

#endif
