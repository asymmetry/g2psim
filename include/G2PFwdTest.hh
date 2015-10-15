// -*- C++ -*-

/* class G2PFwdTest
 * It simulates the movement of the scatted particles only in the target field without any cuts and energy loss.
 * Input variables: fV5tp_tr, fV5react_lab (register in gG2PVars).
 */

// History:
//   May 2014, C. Gu, First public version.
//

#ifndef G2P_FWDTEST_H
#define G2P_FWDTEST_H

#include "G2PProcBase.hh"

class G2PSieve;

class G2PFwdTest : public G2PProcBase
{
public:
    G2PFwdTest();
    virtual ~G2PFwdTest();

    virtual int Begin();
    virtual int Process();
    virtual void Clear(Option_t *opt = "");

    // Gets

    // Sets
    void SetSieve(const char *opt);

protected:
    virtual int Configure(EMode mode = kTWOWAY);
    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    bool fSieveOn;
    int fHoleID;

    double fV5react_tr[5];
    double freactz_tr;

    double fV5sieve_tr[5];
    double fV5tpproj_tr[5];

    G2PSieve *pSieve;

private:
    static G2PFwdTest *pG2PFwdTest;

    ClassDef(G2PFwdTest, 1)
};

#endif
