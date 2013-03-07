#ifndef G2P_MAPFIELD_H
#define G2P_MAPFIELD_H

#include "G2PFieldBase.hh"

class G2PMapField : public G2PFieldBase
{
public:
    G2PMapField(const char* name);
    ~G2PMapField();

    EStatus Init();

protected:
    G2PMapField(); // Only for ROOT I/O

    bool ReadMap();

private:
    ClassDef(G2PMapField, 1)
};

#endif
