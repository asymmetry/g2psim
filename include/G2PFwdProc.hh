#ifndef G2P_FWDPROC_H
#define G2P_FWDPROC_H

#include "G2PProcBase.hh"
#include "G2PSieve.hh"

class G2PDrift;
class G2PHRSTrans;
class G2PDBRec;

class G2PFwdProc : public G2PProcBase, public G2PSieve
{
public:
    G2PFwdProc();
    ~G2PFwdProc();

    int Init();
    int Begin();
    int Process();
    void Clear();

protected:
    int DefineVariables(EMode mode = kDefine);

    void MakePrefix();

    void ApplyVDCRes(double* V5fp);

    double fHRSAngle;
    double fHRSMomentum;

    double fV5tg_tr[5];
    double fV5react_lab[5];

    double fV5sieve_tr[5];
    double fV5projtg_tr[5];

    double fV5fp_tr[5];
    double fV5fp_rot[5];

    G2PDrift* pDrift;
    G2PHRSTrans* pHRS;
    G2PDBRec* pDBRec;

private:
    static G2PFwdProc* pG2PFwdProc;

    ClassDef(G2PFwdProc, 1)
};

#endif
