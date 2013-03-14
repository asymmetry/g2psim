// This file defines a class G2PBPM.
// This class is a tool class.
// G2PProcBase classes will call GetBPMValue() to get bpm readouts.
// The BPM values are in a special coords, TransBPM2Lab() will transform it to
//+lab coords
//
// History:
//   Mar 2013, C. Gu, First public version.
//

#ifndef G2P_BPM_H
#define G2P_BPM_H

#include "G2PAppBase.hh"

class G2PDrift;

class G2PBPM : public G2PAppBase
{
public:
    G2PBPM();
    ~G2PBPM();

    typedef void (G2PBPM::*pfGetBPMValue_)(const double*, double*);

    void SetBPMRes(double a, double b) { fBPMARes = a; fBPMBRes = b; }

    int Init();
    int Begin();

    void GetBPMValue(const double* V5beam_lab, double* V5bpm_bpm);
    
    void TransBPM2Lab(const double* V5_bpm, double* V5_lab);

    static G2PBPM* GetInstance() { return pG2PBPM; }

protected:
    void SetBPM();

    void GetBPMValue0(const double* V5beam_lab, double* V5bpm_bpm);
    void GetBPMValue1(const double* V5beam_lab, double* V5bpm_bpm);
    void GetBPMValue4(const double* V5beam_lab, double* V5bpm_bpm);
    void GetBPMValue5(const double* V5beam_lab, double* V5bpm_bpm);
    void GetBPMValue7(const double* V5beam_lab, double* V5bpm_bpm);
    void GetBPMValue9(const double* V5beam_lab, double* V5bpm_bpm);

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
