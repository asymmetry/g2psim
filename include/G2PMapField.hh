// -*- C++ -*-

/* class G2PMapField
 * This file defines a class G2PMapField.
 * It is derived from G2PField class.
 * The interpolation function is defined in G2PField class, this class only sets field map.
 * The file name of the field map is set when initializing this class.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//

#ifndef G2P_MAPFIELD_H
#define G2P_MAPFIELD_H

#include "G2PField.hh"

class G2PMapField : public G2PField {
public:
    G2PMapField(const char* name);
    ~G2PMapField();

    int Begin();

protected:
    G2PMapField(); // Only for ROOT I/O

    int ReadMap();

private:
    ClassDef(G2PMapField, 1)
};

#endif
