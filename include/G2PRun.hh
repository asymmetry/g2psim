#ifndef G2P_RUN_H
#define G2P_RUN_H

#include "G2PRunBase.hh"

class G2PRun : public G2PRunBase
{
public:
    G2PRun();
    ~G2PRun();

    int Init();

private:

    ClassDef(G2PRun, 1)
};

#endif
