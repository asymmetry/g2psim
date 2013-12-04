// -*- C++ -*-

/* class G2PSieveGun
 * It generates special events which can pass through the holes on sieve slit.
 * Only work for no target field situation.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//

#ifndef G2P_SIEVEGUN_H
#define G2P_SIEVEGUN_H

#include <vector>

#include "G2PGun.hh"

class G2PSieve;

class G2PSieveGun : public G2PGun {
public:
    G2PSieveGun();
    virtual ~G2PSieveGun();

    virtual int Init();

protected:
    virtual int Shoot(double* V5beam_lab, double* V5react_tr);

    G2PSieve* pSieve;

private:
    ClassDef(G2PSieveGun, 1)
};

#endif
