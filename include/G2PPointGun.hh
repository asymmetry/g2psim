// -*- C++ -*-

/* class G2PPointGun
 * This file defines a class G2PPointGun.
 * It generates events around one particular 4-D point using Gaussian distribution.
 * G2PProcBase classes will call Shoot() to get reaction point kinematics.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//

#ifndef G2P_POINTGUN_H
#define G2P_POINTGUN_H

#include "G2PGun.hh"

class G2PPointGun : public G2PGun {
public:
    G2PPointGun();
    ~G2PPointGun();

    int Shoot(double* V51, double* V52, double* V53 = NULL);

    bool UseData() {
        return false;
    }

private:
    ClassDef(G2PPointGun, 1)
};

#endif
