// -*- C++ -*-

/* class G2POptics
 * Special class to treat optics data.
 */

// History:
//   Jan 2014, C. Gu, First public version.
//

#ifndef G2P_OPTICS_H
#define G2P_OPTICS_H

#include <queue>
#include <vector>

#include "G2PProcBase.hh"

using namespace std;

class G2PSieve;

class G2POptics : public G2PProcBase
{
public:
    G2POptics(const char *filename);
    virtual ~G2POptics();

    virtual int Begin();
    virtual int Process();
    virtual void Clear(Option_t *opt = "");

    // Gets

    // Sets
    void SetHRSMomentum(int n, double *value);
    void SetBPMZ(double value);
    void SetFoilZ(int n, double *value);
    void SetTargetMass(int n, double *value);
    void SetEnergyLoss(int n, double *value);

protected:
    G2POptics(); // Only for ROOT I/O

    void CalPos(double *ang, double *pos);
    double Distance(double *V2a, double *V2b);

    int LoadData();

    virtual int Configure(EMode mode = kTWOWAY);
    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    const char *fDataFile;

    typedef struct {
        int ind;
        double xf, tf, yf, pf, eb, xb, tb, yb, pb;
    } sData;

    queue<sData> fData;

    double fE0;
    double fm;
    double fM0;
    double fELoss;
    double fBPMZ;
    double fFoilZ;

    int fNFoil;
    vector<double> fHRSMomentumV;
    vector<double> fFoilZV;
    vector<double> fM0V;
    vector<double> fELossV;

    int fHoleID;

    double fV5bpm_bpm[5];
    double fV5bpm_lab[5];
    double fV5bpm_tr[5];
    double fbpmz_tr;

    double fV5react_tr[5];

    double fV5sieve_tr[5];
    double fV5tpproj_tr[5];

    double fV5fp_det[5];
    double fV5fp_tr[5];
    double fV5fp_rot[5];

    G2PSieve *pSieve;

private:
    static G2POptics *pG2POptics;

    ClassDef(G2POptics, 1)
};

#endif
