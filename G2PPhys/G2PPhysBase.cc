// -*- C++ -*-

/* class G2PPhysBase
 * Abstract base class of G2PPhys classes.
 * It provides interface functions.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//

#include <cstdlib>
#include <vector>
#include <cmath>

#include "G2PPhysBase.hh"

using namespace std;

static const double kU = 931.49406121;

G2PPhysBase::G2PPhysBase() :
iZ(1), iA(1), fTargetMass(0.0), iPID(11)
{
    fPars.clear();
}

G2PPhysBase::~G2PPhysBase()
{
    // Nothing to do
}

void G2PPhysBase::SetTarget(int Z, int A)
{
    iZ = Z;
    iA = A;

    if (fabs(fTargetMass) < 1.0e-8) SetTargetMass();
}

void G2PPhysBase::SetPars(double* array, int n)
{
    fPars.clear();

    for (int i = 0; i < n; i++) fPars.push_back(array[i]);
}

void G2PPhysBase::SetTargetMass()
{
    double atomicmass = 0.0;

    if (iZ == 1) atomicmass = 1.00794;
    if (iZ == 2) atomicmass = 4.002602;
    if (iZ == 6) atomicmass = 12.0107;
    if (iZ == 7) atomicmass = 14.0067;
    if (iZ == 8) atomicmass = 15.9994;
    if (iZ == 26) atomicmass = 55.845;
    if (iZ == 29) atomicmass = 63.546;
    if (iZ == 74) atomicmass = 183.84;
    if (fabs(atomicmass) < 1.0e-8) atomicmass = iA; // only an estimation

    fTargetMass = atomicmass*kU;
}
