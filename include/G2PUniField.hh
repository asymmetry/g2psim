#ifndef G2P_UNIFIELD_H
#define G2P_UNIFIELD_H

#include "G2PFieldBase.hh"

class G2PUniField : public G2PFieldBase
{
public:
    G2PUniField();
    ~G2PUniField();

    int Begin();

protected:
    int CreateMap();

private:
    ClassDef(G2PUniField, 1)
};

#endif
