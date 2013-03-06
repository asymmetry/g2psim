#include <vector>

#include "G2PPhysBase.hh"

using namespace std;

G2PPhysBase::G2PPhysBase() :
    iZ(1), iA(1), iPID(11)
{
    fPars.clear();
}

G2PPhysBase::~G2PPhysBase()
{
    // Nothing to do
}

void G2PPhysBase::SetPars(double* array, int n)
{
    for (int i = 0; i<n; i++) fPars.push_back(array[i]);
}
