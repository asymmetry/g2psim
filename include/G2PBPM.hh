#ifndef G2P_BPM_H
#define G2P_BPM_H

#include <cstring>
#include <map>

#include "TROOT.h"
#include "TObject.h"

class G2PBPM : public TObject
{
public:
    G2PBPM();
    ~G2PBPM();

    void SetBPMRes(double value1, double value2) { fBPMARes = value1; fBPMBRes = value2; }

    double GetBPMARes() { return fBPMARes; }
    double GetBPMBRes() { return fBPMBRes; }

    void GetBPMValue(const double* V5beam_lab, double* V5bpm_lab);

private:
    std::map<std::string, int> mSetting;
    int iSetting;

    double fBPMARes, fBPMBRes;

    ClassDef(G2PBPM, 1);
};

#endif
