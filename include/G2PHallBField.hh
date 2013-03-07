#ifndef G2P_HALLBFIELD_H
#define G2P_HALLBFIELD_H

#include "G2PFieldBase.hh"

class G2PHallBField : public G2PFieldBase
{
public:
    G2PHallBField();
    ~G2PHallBField();

    EStatus Init();

protected:
    bool ReadMap();
    bool CreateMap();

private:
    ClassDef(G2PHallBField, 1)
};

#endif
