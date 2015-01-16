// -*- C++ -*-

/* class G2PBPM
 * Calculate beam readout at BPM and target using kinematics from event generator.
 * Transport functions defined in G2PBPMTrans is used in this class.
 * Orbits are defined in G2PBPMTrans.
 *
 * Variables ending with "_bpm" are defined in a special coordinates.
 * BPM2Lab() will transform it to lab coordinates.
 * In output, these variables are labeled as "b_".
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Apr 2013, C. Gu, Add Pengjia's fitting result.
//   Jul 2013, C. Gu, Treat optics (no field) case specially.
//   Sep 2013, C. Gu, Rewrite it as a G2PProcBase class.
//   Jan 2015, C. Gu, Use new Drift() function in G2PProcBase.
//

#ifndef G2P_BPM_H
#define G2P_BPM_H

#include "G2PProcBase.hh"

class G2PBPM : public G2PProcBase
{
public:
    G2PBPM();
    virtual ~G2PBPM();

    typedef void (G2PBPM::*pfGetBPM_)(const double *, double *, double *);

    virtual int Begin();
    virtual int Process();
    virtual void Clear(Option_t *opt = "");

    // Gets

    // Sets
    void SetBPMRes(double a, double b);

protected:
    void SetBPMPos();

    void BPM2Lab(const double *V5_bpm, double *V5_lab);

    void GetBPM(const double *V5beam_lab, double *V5bpm_bpm, double *V4bpmab_bpm);
    void GetBPM0(const double *V5beam_lab, double *V5bpm_bpm, double *V4bpmab_bpm);
    void GetBPM1(const double *V5beam_lab, double *V5bpm_bpm, double *V4bpmab_bpm);
    void GetBPM4(const double *V5beam_lab, double *V5bpm_bpm, double *V4bpmab_bpm);
    void GetBPM5(const double *V5beam_lab, double *V5bpm_bpm, double *V4bpmab_bpm);
    void GetBPM7(const double *V5beam_lab, double *V5bpm_bpm, double *V4bpmab_bpm);
    void GetBPM9(const double *V5beam_lab, double *V5bpm_bpm, double *V4bpmab_bpm);

    void GetBPMAB(const double *V5beam_lab, double *V4bpmab_bpm);

    virtual int Configure(EMode mode = kTWOWAY);
    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    double fE0;
    double fFieldRatio;

    double fBPMAX, fBPMAY;
    double fBPMBX, fBPMBY;
    double fBPMAZ, fBPMBZ;
    double fBPMARes, fBPMBRes;

    double fV5beam_lab[5];
    double fV5bpm_bpm[5];
    double fV5bpm_lab[5];
    double fV5bpm_tr[5];
    double fbpmz_tr;

    double fV4bpmab_bpm[4];

    pfGetBPM_ pfGetBPM;

private:
    static G2PBPM *pG2PBPM;

    ClassDef(G2PBPM, 1)
};

#endif
