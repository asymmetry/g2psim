#ifndef G2P_DATAGUN_H
#define G2P_DATAGUN_H

#include <vector>

#include "G2PGunBase.hh"
#include "G2PSieve.hh"

using namespace std;

class G2PDataGun : public G2PGunBase, public G2PSieve
{
public:
    G2PDataGun(const char* filename);
    ~G2PDataGun();

    typedef bool (G2PDataGun::*pfGun_)(double*, double*, double*);

    void SetOpticsData(bool b) { bIsOptics = b; }

    EStatus Init();
    bool Shoot(double* V51, double* V52, double* V53 = NULL) { return (this->*pfGun)(V51, V52, V53); }

    bool UseData() { return true; }

protected:
    G2PDataGun(); // Only for ROOT I/O

    bool ShootNormal(double* V5bpm_lab, double* V5react_tr, double* V5fp_tr);
    bool ShootOptics(double* V5bpm_lab, double* V5react_tr, double* V5fp_tr);

    bool LoadData();

    bool bIsOptics;

    const char* pDataFileName;

    double fTargetMass;
    double fEnergyLoss;

    typedef struct {
        int ind;
        double xb, tb, yb, pb, zb, xf, tf, yf, pf;
    } sData;
    
    vector<sData> fData;

private:
    pfGun_ pfGun;

    ClassDef(G2PDataGun, 1)
};

#endif
