// -*- C++ -*-

/* class G2PBwdProc
 * This file defines a class G2PBwdProc.
 * It simulates the reconstruction of g2p kinematics.
 * G2PBPM, G2PDBRec, G2PDrift, G2PHRSTrans are used in this class.
 * Input variables: fV5bpm_bpm, fV5projtg_tr, fV5fp_tr (register in G2PRun).
 */

// History:
//   Apr 2013, C. Gu, First public version.
//

#ifndef G2P_BWDPROC_H
#define G2P_BWDPROC_H

#include "G2PProcBase.hh"
#include "G2PSieve.hh"

class G2PBPM;
class G2PDrift;
class G2PHRSTrans;
class G2PDBRec;

class G2PBwdProc : public G2PProcBase, public G2PSieve {
public:
    G2PBwdProc();
    ~G2PBwdProc();

    int Init();
    int Begin();
    int Process();
    void Clear();

protected:
    int DefineVariables(EMode mode = kDefine);

    void MakePrefix();

    double GetEffBPM(double xbpm_tr, const double* V5fp);

    double fBeamEnergy;
    double fHRSAngle;
    double fHRSMomentum;
    double fFieldRatio;

    double fV5bpm_bpm[5];
    double fV5projtg_tr[5];

    double fV5fp_tr[5];

    double fV5rectg_tr[5];
    double fV5recsiv_tr[5];

    double fV5rec_tr[5];
    double fV5rec_lab[5];

    G2PBPM* pBPM;
    G2PDrift* pDrift;
    G2PHRSTrans* pHRS;
    G2PDBRec* pDBRec;

private:
    static G2PBwdProc* pG2PBwdProc;

    ClassDef(G2PBwdProc, 1)
};

#endif
