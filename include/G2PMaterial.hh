// -*- C++ -*-

/* class G2PMaterial
 * Calculate energy loss and multi-scattering for defined material.
 * Most of the algorithm is from SAMC package.
 * Thanks to the author of SAMC.
 */

// History:
//   Sep 2013, C. Gu, First public version.
//   Oct 2013, J. Liu, Check formulas.
//

#ifndef G2P_MATERIAL_H
#define	G2P_MATERIAL_H

#include "G2PAppBase.hh"

class G2PMaterial : public G2PAppBase {
public:
    G2PMaterial(const char* name, double z, double a, double x0, double density);
    virtual ~G2PMaterial();

    virtual double EnergyLoss(double E, double l);
    virtual double MultiScattering(double E, double l);

protected:
    G2PMaterial(); // Only for ROOT I/O

    double Ionization(double E, double l);
    double Bremsstrahlung(double E, double l);
    double b();

    virtual int Configure(EMode mode = kTWOWAY);
    virtual void MakePrefix();

    const char* fName;
    int fZ, fA;
    double fMass;
    double fDensity; // density in g/cm^3
    double fX0; // radiation length in g/cm^2

private:
    static G2PAppList* pG2PMaterial;

    ClassDef(G2PMaterial, 1)
};

#endif

