// -*- C++ -*-

/* class G2PRand
 * A wrapper class of TRandom2 to allow only one instance at a time.
 */

// History:
//   Jan 2013, C. Gu, First public version.
//   Feb 2013, C. Gu, Change it to a class.
//   Sep 2013, C. Gu, Derived directly from TRandom2.
//

#include <cstdlib>
#include <cstdio>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TRandom2.h"

#include "G2PRand.hh"

using namespace std;

G2PRand* G2PRand::pG2PRand = NULL;

G2PRand::G2PRand() {
    // Nothing to do
}

G2PRand::~G2PRand() {
    // Nothing to do
}

G2PRand* G2PRand::GetInstance() {
    if (!pG2PRand) {
        static G2PRand instance;
        pG2PRand = &instance;
    }
    return pG2PRand;
}

ClassImp(G2PRand)
