// -*- C++ -*-

/* class G2PFlatGun
 * It generates events in flat distribution.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//

#ifndef G2P_FLATGUN_H
#define G2P_FLATGUN_H

#include "G2PGun.hh"

class G2PFlatGun : public G2PGun {
public:
    G2PFlatGun();
    virtual ~G2PFlatGun();

protected:
    virtual int Shoot(double* V51, double* V52);

private:
    ClassDef(G2PFlatGun, 1)
};

#endif
