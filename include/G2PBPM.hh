#ifndef G2P_BPM_H
#define G2P_BPM_H

#include "G2PAppsBase.hh"

class G2PDrift;

class G2PBPM : public G2PAppsBase
{
public:
    G2PBPM();
    ~G2PBPM();

    void SetBPMRes(double value1, double value2) { fBPMARes = value1; fBPMBRes = value2; }

    EStatus Init();
    void Clear() { }

    void GetBPMValue(const double* V5beam_lab, double* V5bpm_lab);
    void TransBPM2Lab(const double*V5_bpm, double* V5_lab);

    double GetBPMARes() { return fBPMARes; }
    double GetBPMBRes() { return fBPMBRes; }

    static G2PBPM* GetInstance() { return pG2PBPM; }

    int RegisterModel();

protected:
    int SetBPM();

    double fBeamEnergy;

    double fBPMZ, fBPMAZ, fBPMBZ;
    double fBPMAX, fBPMAY;
    double fBPMBX, fBPMBY;
    double fBPMARes, fBPMBRes;

    G2PDrift* pDrift;

private:
    static G2PBPM* pG2PBPM;

    ClassDef(G2PBPM, 1)
};

#endif
