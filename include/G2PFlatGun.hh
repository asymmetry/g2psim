// -*- C++ -*-

/* class G2PFlatGun
 * It generates events in flat distribution.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Jun 2014, J. Liu, Generate flat distribution in lab coordinates.
//

#ifndef G2P_FLATGUN_H
#define G2P_FLATGUN_H

#include "G2PGun.hh"

class G2PFlatGun : public G2PGun
{
public:
    G2PFlatGun(const char *coords);
    virtual ~G2PFlatGun();

    typedef int (G2PFlatGun::*pfShoot_)(double *, double *);

    virtual int Shoot(double *V5beam_lab, double *V5react_tr);

protected:
    G2PFlatGun(); // Only for ROOT I/O

    virtual int ShootTCS(double *V51, double *V52);
    virtual int ShootHCS(double *V51, double *V52);

    pfShoot_ pfShoot;

private:
    ClassDef(G2PFlatGun, 1)
};

#endif
