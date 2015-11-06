// -*- C++ -*-

/* class G2PData
 * It loads real data to the simulation.
 * The file should be generated with tree2ascii program.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Sep 2013, C. Gu, Rewrite G2PDataGun as a standalone G2PProcBase class.
//

#ifndef G2P_DATA_H
#define G2P_DATA_H

#include <queue>

#include "G2PProcBase.hh"

using namespace std;

class G2PData : public G2PProcBase
{
public:
    G2PData(const char *filename);
    virtual ~G2PData();

    virtual int Begin();
    virtual int Process();
    virtual void Clear(Option_t *opt = "");

protected:
    G2PData(); // Only for ROOT I/O

    int LoadData();

    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    const char *fDataFile;

    typedef struct {
        int ind;
        double xf, tf, yf, pf, eb, xb, tb, yb, pb;
    } sData;

    queue<sData> fData;

    int fHoleID;
    double fE;

    double fV5bpm_bpm[5];
    double fV5bpm_lab[5];
    double fV5bpm_tr[5];
    double fbpmz_tr;

    double fV5fp_det[5];
    double fV5fp_tr[5];
    double fV5fp_rot[5];

private:
    static G2PData *pG2PData;

    ClassDef(G2PData, 1)
};

#endif
