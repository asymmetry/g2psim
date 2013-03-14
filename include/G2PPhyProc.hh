#ifndef G2P_PHYPROC_H
#define G2P_PHYPROC_H

#include "G2PProcBase.hh"

class G2PPhys;

class G2PPhyProc : public G2PProcBase
{
public:
    G2PPhyProc();
    ~G2PPhyProc();

    int Init();
    int Begin();
    int Process();
    void Clear();

protected:
    int DefineVariables(EMode mode = kDefine);

    void MakePrefix();

    double CalXS(const double* V5lab, const double* V5tr, double& scatangle);

    double fBeamEnergy;
    double fHRSAngle;
    double fHRSMomentum;

    double fV5beam_lab[5];
    double fV5bpm_lab[5];

    double fV5react_tr[5];
    double fV5rec_tr[5];

    double fThinit, fThrec;
    double fXSinit, fXSrec;

    G2PPhys* pPhys;

private:
    static G2PPhyProc* pG2PPhyProc;

    ClassDef(G2PPhyProc, 1)
};

#endif
