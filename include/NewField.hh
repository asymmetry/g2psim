// -*- C++ -*-

/* class NewField
 * Generate field map for New Oxford magnets for SoLID NH3 target.
 * Calculate field strength of a particular point from the field map using Lagrange polynomial interpolation, default is 2nd order.
 *
 * Field map may have an angle respect to the lab coordinates.
 * Use SetEulerAngle() to set this angle and the program will rotate the field map to correct direction.
 */

// History:
//   Sep 2013, C. Gu, Put New Oxford field map into NewField.
//

#ifndef New_FIELD_H
#define New_FIELD_H

#define DEBUGWITHROOT

#include <vector>

#include "G2PField.hh"

using namespace std;

class NewField : public G2PField
{
public:
    NewField();
    virtual ~NewField();

protected:
    virtual int CreateMap();

private:
    static NewField *pNewField;

    ClassDef(NewField, 1)
};

#endif
