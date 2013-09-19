// -*- C++ -*-

/* class G2PRand
 * A wrapper class of TRandom2 to allow only one instance at a time.
 */

// History:
//   Jan 2013, C. Gu, First public version.
//   Feb 2013, C. Gu, Change it to a class.
//   Sep 2013, C. Gu, Derived directly from TRandom2.
//

#ifndef G2P_RAND_H
#define G2P_RAND_H

#include "TRandom2.h"

class G2PRand : public TRandom2 {
public:
    static G2PRand* GetInstance(); // Only allow one instance at a time

private:
    G2PRand(); // Only for ROOT I/O
    virtual ~G2PRand();

    static G2PRand* pG2PRand;

    ClassDef(G2PRand, 1)
};

#endif
