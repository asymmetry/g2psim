// -*- C++ -*-

/* class G2PTarget
 * Fake class to initialize and configure the geometries used in the simulation.
 * The geometries are stored in a G2PAppList.
 */

// History:
//   Dec 2014, C. Gu, Add this class for g2p geometries.
//

#ifndef G2P_RUNTYPE_H
#define G2P_RUNTYPE_H

#include "G2PAppBase.hh"

class G2PAppList;

class G2PTarget : public G2PAppBase
{
public:
    G2PTarget();
    virtual ~G2PTarget();

    virtual int Begin();

protected:
    virtual int Configure(EMode mode = kTWOWAY);
    virtual void MakePrefix();

    int fFieldType;
    int fTargetType;

    double fPF;

    G2PAppList *fMats;
    G2PAppList *fGeos;

private:
    ClassDef(G2PTarget, 1)
};

#endif
