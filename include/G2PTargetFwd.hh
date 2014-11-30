// -*- C++ -*-

/* class G2PTargetFwd
 * It simulates the movement of the scatted particles only in the target field without any cuts and energy loss.
 * Input variables: fV5tp_tr, fV5react_lab (register in gG2PVars).
 */

// History:
//   May 2014, C. Gu, First public version.
//

#ifndef G2P_TARGETFWD_H
#define G2P_TARGETFWD_H

#include "G2PProcBase.hh"

class G2PDrift;
class G2PGeoSieve;

class G2PTargetFwd : public G2PProcBase {
public:
    G2PTargetFwd();
    virtual ~G2PTargetFwd();

    virtual int Init();
    virtual int Begin();
    virtual int Process();
    virtual void Clear(Option_t* /*option*/ = "");

    // Gets

    // Sets
    void SetSieve(const char* opt);

protected:
    virtual int Configure(EMode mode = kTWOWAY);
    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    double fHRSAngle;
    double fHRSMomentum;

    bool fSieveOn;
    int fHoleID;

    double fV5react_lab[5];
    double fV5react_tr[5];

    double fV5sieve_tr[5];
    double fV5tpproj_tr[5];

    G2PDrift* pDrift;
    G2PGeoSieve* pSieve;

private:
    static G2PTargetFwd* pG2PTargetFwd;

    ClassDef(G2PTargetFwd, 1)
};

#endif
