// -*- C++ -*-

/* class G2PHallBField
 * This file defines a class G2PHallBField.
 * It is derived from G2PField class.
 * The interpolation function is defined in G2PField class, this class only sets field map.
 * If field map "hallbfield.map" exists, it will read this map.
 * If field map does not exist, it will use sda_ptf_() in hallbfield.f to calculate the field map and store it in "hallbfield.map".
 */

// History:
//   Mar 2013, C. Gu, First public version.
//

#ifndef G2P_HALLBFIELD_H
#define G2P_HALLBFIELD_H

#include "G2PField.hh"

class G2PHallBField : public G2PField {
public:
    G2PHallBField();
    ~G2PHallBField();

    int Begin();

protected:
    int ReadMap();
    int CreateMap();

private:
    ClassDef(G2PHallBField, 1)
};

#endif
