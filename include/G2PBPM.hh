// -*- C++ -*-

/* class G2PBPM
 * This file defines a class G2PBPM.
 * It calculates the beam position at BPM and target using kinematics from event generator.
 * G2PProcBase classes will call GetBPMValue() to get BPM readouts.
 * G2PDrift is used in this class.
 * Transport functions defined in G2PBPMTrans is used in this class.
 *
 * Variables ending with "_bpm" are defined in a special coordinates, TransBPM2Lab() will transform it to lab coordinates.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Apr 2013, C. Gu, Add Pengjia's fitting result.
//   Jul 2013, C. Gu, Treat optics (no field) case specially.
//

#ifndef G2P_BPM_H
#define G2P_BPM_H

#include "G2PAppBase.hh"

class G2PDrift;

class G2PBPM : public G2PAppBase {
public:
    G2PBPM();
    ~G2PBPM();

    typedef void (G2PBPM::*pfGetBPMValue_)(const double*, double*, double*);

    void SetBPMRes(double a, double b) {
        fBPMARes = a;
        fBPMBRes = b;
    }

    int Init();
    int Begin();

    void GetBPMValue(const double* V5beam_lab, double* V5bpm_bpm, double* V2bpma_bpm, double* V2bpmb_bpm);

    void TransBPM2Lab(const double* V5_bpm, double* V5_lab);

    static G2PBPM* GetInstance() {
        return pG2PBPM;
    }

protected:
    void SetBPM();

    void GetBPMValue0(const double* V5beam_lab, double* V5bpm_bpm, double* V4);
    void GetBPMValue1(const double* V5beam_lab, double* V5bpm_bpm, double* V4);
    void GetBPMValue4(const double* V5beam_lab, double* V5bpm_bpm, double* V4);
    void GetBPMValue5(const double* V5beam_lab, double* V5bpm_bpm, double* V4);
    void GetBPMValue7(const double* V5beam_lab, double* V5bpm_bpm, double* V4);
    void GetBPMValue9(const double* V5beam_lab, double* V5bpm_bpm, double* V4);
    void GetBPMValueO(const double* V5beam_lab, double* V5bpm_bpm, double* V4);

    void GetBPMAB(const double* V5beam_lab, float* xout);

    double fBeamEnergy;
    double fFieldRatio;

    double fBPMAX, fBPMAY;
    double fBPMBX, fBPMBY;
    double fBPMAZ, fBPMBZ;
    double fBPMARes, fBPMBRes;

    G2PDrift* pDrift;

    pfGetBPMValue_ pfGetBPMValue;

private:
    static G2PBPM* pG2PBPM;

    ClassDef(G2PBPM, 1)
};

#endif
