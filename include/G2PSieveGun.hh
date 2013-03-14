#ifndef G2P_SIEVEGUN_H
#define G2P_SIEVEGUN_H

#include <vector>

#include "G2PGunBase.hh"
#include "G2PSieve.hh"

using namespace std;

class G2PSieveGun : public G2PGunBase, public G2PSieve
{
public:
    G2PSieveGun();
    ~G2PSieveGun();

    typedef int (G2PSieveGun::*pfGun_)(double*, double*, double*);

    void SetUseFast(bool fast) { bUseFast = fast; }

    int Begin();
    int Shoot(double* V51, double* V52, double* V53 = NULL) { return (this->*pfGun)(V51, V52, V53); }

    bool UseData() { return false; }

protected:
    int ShootNormal(double* V5beam_lab, double* V5react_tr, double* reserved);
    int ShootFast(double* V5beam_lab, double* V5react_tr, double* reserved);

    bool bUseFast;

    double fTargetMass;
    double fEnergyLoss;

    double fThreshold;

private:
    pfGun_ pfGun;

    ClassDef(G2PSieveGun, 1)
};

#endif
