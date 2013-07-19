// -*- C++ -*-

/* class G2PUniField
 * This file defines a class G2PUniField.
 * It is derived from G2PField class.
 * The interpolation function is defined in G2PField class, this class only sets field map.
 * A uniform field along positive z direction is set in this class.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//

#ifndef G2P_UNIFIELD_H
#define G2P_UNIFIELD_H

#include "G2PField.hh"

class G2PUniField : public G2PField {
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
