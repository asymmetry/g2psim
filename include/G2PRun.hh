// -*- C++ -*-

/* class G2PRun
 * This file defines a class G2PRun.
 * It describes the procedure of a simulation.
 * G2PProcBase classes and their input variables should be registered here.
 */

// History:
//   Jan 2013, C. Gu, First public version.
//   May 2013, C. Gu, Add G2PProcBase classes, G2PRun uses these classes to describe the simulation procedure.
//

#ifndef G2P_RUN_H
#define G2P_RUN_H

#include "G2PRunBase.hh"

class G2PRun : public G2PRunBase {
public:
    G2PRun();
    ~G2PRun();

    int Init();

private:

    ClassDef(G2PRun, 1)
};

#endif
