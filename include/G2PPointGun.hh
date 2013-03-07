#ifndef G2P_POINTGUN_H
#define G2P_POINTGUN_H

#include "G2PGunBase.hh"

class G2PPointGun : public G2PGunBase
{
public:
    G2PPointGun();
    ~G2PPointGun();
   
    bool Shoot(double* V51, double* V52, double* V53 = NULL);

    bool UseData() { return false; }

private:
    ClassDef(G2PPointGun, 1)
};

#endif
