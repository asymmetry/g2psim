// -*- C++ -*-

/* class G2PBPM
 * Calculate beam readout at BPM and target using kinematics from event generator.
 * Transport functions defined in G2PBPMTrans is used in this class.
 * Orbits are defined in G2PBPMTrans.
 *
 * Variables ending with "_bpm" are defined in a special coordinates.
 * TransBPM2Lab() will transform it to lab coordinates.
 * In output, these variables are labeled as "b_".
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Apr 2013, C. Gu, Add Pengjia's fitting result.
//   Jul 2013, C. Gu, Treat optics (no field) case specially.
//   Sep 2013, C. Gu, Rewrite it as a G2PProcBase class.
//

#ifndef G2P_BPM_H
#define G2P_BPM_H

#include "G2PProcBase.hh"

class G2PDrift;

class G2PBPM : public G2PProcBase {
public:
    G2PBPM();
    virtual ~G2PBPM();

    typedef void (G2PBPM::*pfGetBPM_)(const double*, double*, double*);

    virtual int Init();
    virtual int Begin();
    virtual int Process();
    virtual void Clear();

    // Gets

    // Sets
    void SetBPMRes(double a, double b);

protected:
    void GetBPM(const double* V5beam_lab, double* V5bpm_bpm, double* V4);
    void TransBPM2Lab(const double* V5_bpm, double* V5_lab);

    void GetBPM0(const double* V5beam_lab, double* V5bpm_bpm, double* V4);
    void GetBPM1(const double* V5beam_lab, double* V5bpm_bpm, double* V4);
    void GetBPM4(const double* V5beam_lab, double* V5bpm_bpm, double* V4);
    void GetBPM5(const double* V5beam_lab, double* V5bpm_bpm, double* V4);
    void GetBPM7(const double* V5beam_lab, double* V5bpm_bpm, double* V4);
    void GetBPM9(const double* V5beam_lab, double* V5bpm_bpm, double* V4);
    void GetBPMO(const double* V5beam_lab, double* V5bpm_bpm, double* V4);

    void GetBPMAB(const double* V5beam_lab, float* xout);

    void SetBPMPos();

    virtual int Configure(EMode mode = kTWOWAY);
    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    double fBeamEnergy;
    double fFieldRatio;

    double fBPMAX, fBPMAY;
    double fBPMBX, fBPMBY;
    double fBPMAZ, fBPMBZ;
    double fBPMARes, fBPMBRes;

    double fV5beam_lab[5];
    double fV5bpm_bpm[5];
    double fV5bpm_lab[5];
    double fV2bpma_bpm[2];
    double fV2bpmb_bpm[2];

    G2PDrift* pDrift;

    pfGetBPM_ pfGetBPM;

private:
    static G2PBPM* pG2PBPM;

    ClassDef(G2PBPM, 1)
};

#endif
