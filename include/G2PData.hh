// -*- C++ -*-

/* class G2PData
 * It loads real data to the simulation.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Sep 2013, C. Gu, Rewrite G2PDataGun as a standalone G2PProcBase class.
//

#ifndef G2P_DATA_H
#define G2P_DATA_H

#include <vector>

#include "G2PProcBase.hh"

using namespace std;

class G2PData : public G2PProcBase {
public:
    G2PData(const char* filename);
    virtual ~G2PData();

    virtual int Begin();
    virtual int Process();
    virtual void Clear();

protected:
    G2PData(); // Only for ROOT I/O

    int LoadData();

    virtual int Configure(EMode mode = kTWOWAY);
    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    const char* fDataFile;

    typedef struct {
        int ind;
        double xb, tb, yb, pb, zb, xf, tf, yf, pf;
    } sData;

    vector<sData> fData;

    double fHRSAngle;

    double fV5bpm_lab[5];
    double fV5fp_tr[5];
    double fV5fp_rot[5];

private:
    static G2PData* pG2PData;

    ClassDef(G2PData, 1)
};

#endif
