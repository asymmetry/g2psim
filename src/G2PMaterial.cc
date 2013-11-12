// -*- C++ -*-

/* class G2PMaterial
 * Calculate energy loss and multi-scattering for defined material.
 * Most of the algorithm is from SAMC package.
 * Thanks to the author of SAMC.
 */

// History:
//   Sep 2013, C. Gu, First public version.
//   Oct 2013, J. Liu, Definiton of material
//

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PMaterial.hh"
#include "G2PRand.hh"
#include "G2PVarDef.hh"

using namespace std;

static const double kELECTRONMASS = 0.510998918; //MeV

G2PAppList* G2PMaterial::pG2PMaterial = NULL;

G2PMaterial::G2PMaterial() :
fName(NULL)
{
    // Only for ROOT I/O
}

G2PMaterial::G2PMaterial(const char* name, double z, double a, double x0, double density) :
fName(name), fZ(z), fA(a), fMass(0), fDensity(density), fX0(x0)
{
    // Constructor

    pG2PMaterial->Add(this);
}

G2PMaterial::~G2PMaterial()
{
    // Destructor

    pG2PMaterial->Remove(this);
}

double G2PMaterial::EnergyLoss(double E, double l)
{
    double EMeV = E * 1000; // MeV

    double result = 0;
    result += Ionization(EMeV, l);
    result += Bremsstrahlung(EMeV, l);

    return result / 1000.0; // GeV
}

double G2PMaterial::MultiScattering(double E, double l)
{
    // only for electron

    double EMeV = E * 1000;

    double thicknessr = l * fDensity / fX0;

    double lPsq = EMeV * EMeV - kELECTRONMASS*kELECTRONMASS;
    double bcp = lPsq / EMeV;
    double ltheta0 = 13.6 / bcp * sqrt(thicknessr)*(1 + 0.038 * log(thicknessr));
    if (thicknessr != 0) {
        //return pRand->Gaus(0, ltheta0 / 2.3548); // rad // sigma=width/(2*sqrt(2))
        return pRand->Gaus(0, ltheta0 / 1.3548); // rad
    } else
        return 0;
}

double G2PMaterial::Ionization(double E, double l)
{
    // Particle Data Group Booklet Equ (27.9)

    double thickness = l*fDensity;

    double lK = 0.307075; // cm^2/g for A=1 g/mol
    double lbetasq = 1 - kELECTRONMASS * kELECTRONMASS / (E * E);
    double lxi = lK / 2 * fZ / fA * thickness / lbetasq; // fThickness: g/cm^2
    double lhbarwsq = 28.816 * 28.816 * fDensity * fZ / fA * 1e-12; // MeV // fDensity is density of absorber
    double j = 0.200;
    double Delta_p = lxi * (log(2 * kELECTRONMASS * lxi / lhbarwsq) + j);
    double lw = 4 * lxi;
    double result = 0;
    if (fZ != 0 && fA != 0 && thickness != 0 && fDensity != 0)
        result = pRand->Landau(Delta_p, lw);
    if (result > (E - kELECTRONMASS))
        result = E - kELECTRONMASS;
    if (result < 0)
        result = 0;

    return result; // GeV!
}

double G2PMaterial::Bremsstrahlung(double E, double l)
{
    // Bremsstrahlung Energy Loss for external and internal(equivalent radiator)
    // Xiaodong Jiang, PhD.thesis Equ (5.15)
    // http://filburt.mit.edu/oops/Html/Pub/theses/xjiang.ps
    // *0.999 to avoid lose all energy

    double bt = l * fDensity / fX0 * b();

    double result = 0;
    if (bt != 0)
        result = E * pow(pRand->Uniform()*0.999, 1. / bt);
    if (result > (E - kELECTRONMASS)) // GeV!
        result = E - kELECTRONMASS;
    if (result < 0)
        result = 0;

    return result; // GeV!
}

double G2PMaterial::b()
{
    // Phys.Rev.D 12,1884 A45

    if (fZ != 0) {
        double eta = log(1440 * pow(fZ, -2 / 3.)) / log(183 * pow(fZ, -1 / 3.));
        return 4. / 3. * (1 + 1. / 9. * (fZ + 1) / (fZ + eta) / log(183 * pow(fZ, -1 / 3.)));
    }

    return 0;
}

int G2PMaterial::Configure(EMode mode)
{
    if (mode == kREAD || mode == kTWOWAY) {
        if (fIsInit) return 0;

        else fIsInit = true;
    }

    ConfDef confs[] = {
        {"z", "Z", kINT, &fZ},
        {"a", "A", kINT, &fA},
        {"mass", "Mass", kDOUBLE, &fMass},
        {"density", "Density", kDOUBLE, &fDensity},
        {"radlen", "Radiation Length", kDOUBLE, &fX0},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

void G2PMaterial::MakePrefix()
{
    const char* base = "material";

    G2PAppBase::MakePrefix(Form("%s.%s", base, fName));
}

ClassImp(G2PMaterial)
