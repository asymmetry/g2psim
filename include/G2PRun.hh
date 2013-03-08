#ifndef G2P_RUN_H
#define G2P_RUN_H

#include "G2PRunBase.hh"
#include "G2PSieve.hh"

class TTree;

class G2PRun : public G2PRunBase, public G2PSieve
{
public:
    G2PRun();
    ~G2PRun();

    typedef int (G2PRun::*pfRun_)();

    void SetUseEffBPM(bool b) { bUseEffBPM = b; }

    EStatus Init();
    int Begin() { return 0; }
    int End() { return 0; }
    void Clear();

    int Run() { return (this->*pfRun)(); }

    int DefineVariables(TTree *t);

protected:
    int RunSim();
    int RunData();

    void SmearVDC(double* V5_fp);
    double GetEffBPM(double xbpm_tr, const double* V5fp);
    double CalXS(const double* V5lab, const double* V5tr, double& scatangle);
    //double DriftPath();

    bool bUseEffBPM;

    bool bIsGood;

    double fV5beam_lab[5];
    double fV5react_tr[5];
    double fV5tg_tr[5];
    double fV5tg_lab[5];

    double fV5fpdata_tr[5];
    double fV5fpdata_rot[5];

    double fV5bpm_lab[5];
    double fV5bpm_tr[5];

    double fV5fp_tr[5];
    double fV5fp_rot[5];

    double fV5rec_tr[5];
    double fV5rec_lab[5];
    double fV5recdb_tr[5];
    double fV5recdb_lab[5];

    double fV5sieve_tr[5];
    double fV5tgproj_tr[5];
    double fV5rectg_tr[5];
    double fV5recsieve_tr[5];

    double fEffXbeam, fRealXbeam;

    double fThetainit, fThetarec;
    double fXSinit, fXSrec;

    TTree* pTree;

private:
    pfRun_ pfRun;

    ClassDef(G2PRun, 1)
};

#endif
