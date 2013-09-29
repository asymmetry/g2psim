// -*- C++ -*-

/* class G2PMaterial
 * Calculate energy loss and multi-scattering for defined material.
 * Most of the algorithm is from SAMC package.
 * Thanks to the author of SAMC.
 */

// History:
//   Sep 2013, C. Gu, First public version.
//

#ifndef G2P_MATERIAL_H
#define	G2P_MATERIAL_H

#include "G2PAppBase.hh"

class G2PMaterial : public G2PAppBase {
public:
    G2PMaterial(const char* name);
    virtual ~G2PMaterial();

    virtual int Begin();

    virtual double EnergyLoss(double E);
    virtual double MultiScattering(double E);

protected:
    G2PMaterial(); // Only for ROOT I/O

    double Ionization(double E);
    double Bremsstrahlung(double E);
    double b();

    virtual int Configure(EMode mode = kTWOWAY);
    virtual void MakePrefix();

    const char* fName;
    int fZ, fA;
    double fMass;
    double fLength;
    double fDensity; // density in g/cm^3
    double fX0; // radlen in g/cm^2
    double fThickness; // thickness in g/cm^2
    double fThicknessR; // thickness in radlen
    double fBT;

private:
    ClassDef(G2PMaterial, 1)
};

#endif

