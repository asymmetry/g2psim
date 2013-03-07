#ifndef G2P_FLATGUN_H
#define G2P_FLATGUN_H

#include "G2PGunBase.hh"

class G2PFlatGun : public G2PGunBase
{
public:
    G2PFlatGun();
    ~G2PFlatGun();
   
    bool Shoot(double* V51, double* V52, double* V53 = NULL);

    bool UseData() { return false; }

private:
    ClassDef(G2PFlatGun, 1)
};

#endif
