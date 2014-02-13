// -*- C++ -*-

/* class G2PPhysEl
 * Class to calculate elastic cross section.
 * Unit is ub/sr.
 * 
 * Elastic cross section models.
 * * All: Form factors from K. C. Stansfield et al., Phys. Rev. C, 3(1971)1448
 * * 1H : Form factors from J. Arrington, Phys. Rev. C, 69(2004)022201
 * * 4He: Charge and magnetization densities from De Jager, At. Data Nucl. Data Tables, 14(1974)
 * * 12C: Charge distribution from L. S. Cardman et al., Phys. Lett. B, 91(1970)203 
 * * 14N: Charge and magnetization densities from De Jager, At. Data Nucl. Data Tables, 14(1974)
 */

// History:
//   Mar 2013, C. Gu, First public version, with C12 models.
//   Apr 2013, C. Gu, Add L. Cardman's C12 charge densities.
//   Nov 2013, C. Gu, Add He charge and magnetization densities from D. Jager, original coded by M. Friedman.
//

#include <cstdio>
#include <cmath>
#include <vector>

#include "gauss_legendre.h"
#include "G2PPhysBase.hh"

#include "G2PPhysEl.hh"

#define QUADRATURE_ORDER 128

using namespace std;

extern "C" {
void carbon_(double* Ei, double* ang, double* xs);
}

static const double kPi = 3.1415926535897932384626;
static const double kDEG = kPi / 180.0;
static const double kMEV = 1.0e-3;
static const double kAlpha = 1 / 137.035999679;
static const double kC = 299792458;
static const double kQe = 1.602176565e-19;
static const double kHbar = 1.054571726e-34;
static const double kFm = kHbar*kC / 1e-15 / 1e9 / kQe;

static const double kIntegralMin = 0.0;
static const double kIntegralMax = 10.0;
static double kRho0GE_He4, kRho0GE_N14, kRho0GM_N14;

static double rhoGE_He4(double x, void* data)
{
    return 4 * kPi * x * x * (1 + 0.445 * x * x / 1.016) / (1 + exp((x - 1.008) / 0.327));
}

static double funcGE_He4(double x, void* data)
{
    double qfm = *(double*) data;
    return kRho0GE_He4 * rhoGE_He4(x, NULL) * sin(qfm * x) / (qfm * x);
}

static double GE_He4(double q2)
{
    // calculate GE from density function

    double q = sqrt(q2);
    double qfm = q / kFm;

    double params[1] = {qfm};
    return fabs(gauss_legendre(QUADRATURE_ORDER, funcGE_He4, params, kIntegralMin, kIntegralMax));
}

static double Carbon(double Ei, double theta)
{
    double e = Ei * 1000;
    double t = theta / kDEG;
    double XS;
    carbon_(&e, &t, &XS);

    XS *= 1.0e4; // fm^2 (1e-30m^2) to microbarn (1e-34m^2)
    return XS;
}

static double rhoGE_N14(double x, void* data)
{
    return 4 * kPi * x * x * (1 + 1.234 * x * x / 3.0976) * exp(-x * x / 3.0976);
}

static double funcGE_N14(double x, void* data)
{
    double qfm = *(double*) data;
    return kRho0GE_N14 * rhoGE_N14(x, NULL) * sin(qfm * x) / (qfm * x);
}

static double rhoGM_N14(double x, void* data)
{
    // double Z = 7, A = 14, a0 = 1.61; 
    // double alpha0 = (Z-2.)/3.;
    // double rp     = 0.810; // proton magnetic radius fm (De Jager 1974)
    // double ap     = (2/3.)*rp*rp;	
    // double a      = (A-1)/A*a0*a0+ap*ap;
    // double alpha  = alpha0*a0*a0/(a*a+3/2.*alpha0*(a*a-a0*a0));

    return 4 * kPi * x * x * (1 + 0.25193 * x * x / 6.751) * exp(-x * x / 6.751);
}

static double funcGM_N14(double x, void* data)
{
    double qfm = *(double*) data;
    return kRho0GM_N14 * rhoGM_N14(x, NULL) * sin(qfm * x) / (qfm * x);
}

static double GE_N14(double q2)
{
    // calculate GE from density function

    double q = sqrt(q2);
    double qfm = q / kFm;

    double params[1] = {qfm};
    return fabs(gauss_legendre(QUADRATURE_ORDER, funcGE_N14, params, kIntegralMin, kIntegralMax));
}

static double GM_N14(double q2)
{
    // calculate GM from density function
    double q = sqrt(q2);
    double qfm = q / kFm;

    double params[1] = {qfm};
    return fabs(gauss_legendre(QUADRATURE_ORDER, funcGM_N14, params, kIntegralMin, kIntegralMax));
}

static double FF_All(int Z, double Q2)
{
    // K. C. Stansfield et al., PRC 3(1971)1448
    // Only valid for Q2 in between the first diffraction minimum
    // For C12, only valid for Q2 in between the second diffraction minimum

    double XA = 1.64; // unit in Fm
    double FF;

    double XALPHA = (Z - 2.0) / 3.0;
    double Q2FM = Q2 / (kFm * kFm); // fm^-2 // Corrected on Feb 17, 2011

    if (Z == 6) {
        if (Q2FM < 3.2)
            XA = 1.64; // fm
        else if (Q2FM > 3.5) // second diffraction minimum
            XA = 1.68; // fm
        else
            XA = 0; // fm
    }

    FF = 1.0 - XALPHA / (2.0 * (2. + 3. * XALPHA)) * Q2FM * XA*XA;
    FF = FF * exp(-(Q2FM * XA * XA) / 4.0);

    if (FF < 1.0e-6) FF = 1.0e-6;

    return FF;
}

G2PPhysEl::G2PPhysEl() : iSetting(1)
{
    // Constructor

    double temp = gauss_legendre(QUADRATURE_ORDER, rhoGE_He4, NULL, kIntegralMin, kIntegralMax);
    kRho0GE_He4 = 1. / temp;
    temp = gauss_legendre(QUADRATURE_ORDER, rhoGE_N14, NULL, kIntegralMin, kIntegralMax);
    kRho0GE_N14 = 1. / temp;
    temp = gauss_legendre(QUADRATURE_ORDER, rhoGM_N14, NULL, kIntegralMin, kIntegralMax);
    kRho0GM_N14 = 1. / temp;
}

G2PPhysEl::~G2PPhysEl()
{
    // Nothing to do
}

void G2PPhysEl::SetPars(double* array, int n)
{
    G2PPhysBase::SetPars(array, n);

    switch (n) {
    case 0:
        break;
    case 1:
        iSetting = int(fPars[0]);
        break;
    default:
        printf("Error: G2PPhysEl::SetPars(): Invalid number of pars.\n");
        break;
    }
}

double G2PPhysEl::GetXS(double Ei, double Ef, double theta)
{
    if (iSetting == 2) return GetXS_All(Ei, theta);

    if ((fZ == 1)&&(fA == 1)) return GetXS_H1(Ei, theta);
    if ((fZ == 2)&&(fA == 4)) return GetXS_He4(Ei, theta);
    if ((fZ == 6)&&(fA == 12)) return Carbon(Ei, theta);
    if ((fZ == 7)&&(fA == 14)) return GetXS_N14(Ei, theta);

    return GetXS_All(Ei, theta);
}

double G2PPhysEl::GetXS_H1(double Ei, double theta)
{
    // reference: Xiaohui thesis, eq. 1.39

    int Z = 1;
    double M = 0.93889; // GeV
    double eP = Ei / (1 + Ei / M * (1 - cos(theta))); // scattered electron energy
    double Q2 = 4 * Ei * eP * sin(theta / 2) * sin(theta / 2); // units GeV2
    double tau = Q2 / 4.0 / M / M;

    double x1 = Q2;
    double x2 = x1*Q2;
    double x3 = x2*Q2;
    double x4 = x3*Q2;
    double x5 = x4*Q2;
    double x6 = x5*Q2;
    double GE = 1.0 / (1.0 + 2.94 * x1 + 3.04 * x2 - 2.255 * x3 + 2.002 * x4 - 0.5338 * x5 + 0.04875 * x6);
    double GM = 2.793 / (1 + 3.0 * x1 + 1.39 * x2 + 0.122 * x3 - 0.00834 * x4 + 4.25e-4 * x5 - 7.79e-6 * x6);

    double sigma = Z * Z * kAlpha * kAlpha / Q2 * (eP / Ei)*(eP / Ei)*(2 * tau * GM * GM + 1 / tan(theta / 2) / tan(theta / 2) / (1 + tau)*(GE * GE + tau * GM * GM));

    return sigma * kFm * kFm * 1e4; // microbarn
}

double G2PPhysEl::GetXS_He4(double Ei, double theta)
{
    // reference: Xiaohui thesis, eq. 1.39

    int Z = 2;
    double M = 3.7284; // GeV
    double eP = Ei / (1 + Ei / M * (1 - cos(theta))); // scattered electron energy
    double Q2 = 4 * Ei * eP * sin(theta / 2) * sin(theta / 2); // units GeV2
    double tau = Q2 / 4.0 / M / M;
    double GE = GE_He4(Q2);
    double GM = 0.0;
    double sigma = Z * Z * kAlpha * kAlpha / Q2 * (eP / Ei)*(eP / Ei)*(2 * tau * GM * GM + 1 / tan(theta / 2) / tan(theta / 2) / (1 + tau)*(GE * GE + tau * GM * GM));

    return sigma * kFm * kFm * 1e4; // microbarn
}

double G2PPhysEl::GetXS_N14(double Ei, double theta)
{
    // reference: Xiaohui thesis, eq. 1.39

    int Z = 7;
    double M = 13.04378; // GeV
    double eP = Ei / (1 + Ei / M * (1 - cos(theta))); // scattered electron energy
    double Q2 = 4 * Ei * eP * sin(theta / 2) * sin(theta / 2); // units GeV2
    double tau = Q2 / 4.0 / M / M;
    double GE = GE_N14(Q2);
    double GM = GM_N14(Q2);
    double sigma = Z * Z * kAlpha * kAlpha / Q2 * (eP / Ei)*(eP / Ei)*(2 * tau * GM * GM + 1 / tan(theta / 2) / tan(theta / 2) / (1 + tau)*(GE * GE + tau * GM * GM));

    return sigma * kFm * kFm * 1e4; // microbarn
}

double G2PPhysEl::GetXS_All(double Ei, double theta)
{
    double Recoil = 1.0 / (1.0 + Ei / fTargetMass * (1. - cos(theta)));
    double Ef = Recoil*Ei;
    double Q2 = 2. * fTargetMass * (Ei - Ef);

    // Calculated cross sections are often written in units of hbar2c2/GeV2
    double hbc2 = 0.38938; // turn GeV to mbar

    double C = cos(theta / 2.0);
    double S = fZ * kAlpha * C / (2. * Ei * (1. - C * C));
    double Mott = S * S * hbc2 * 1000.0; // microbarn

    double FF = FF_All(fZ, Q2);

    return Recoil * Mott * FF*FF; // microbarn
}
