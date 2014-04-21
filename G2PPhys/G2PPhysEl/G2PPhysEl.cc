// -*- C++ -*-

/* class G2PPhysEl
 * Class to calculate elastic cross section.
 * Unit is ub/sr.
 *
 * Elastic cross section models.
 * * All: Form factors from K. C. Stansfield et al., Phys. Rev. C, 3(1971)1448
 * * H1 : Form factors from S. Venkat et al., Phys. Rev. C, 83(2011)015203 (global fit, with TPE correction)
 *                          J. Arrington et al., Phys. Rev. C 76(2007)035201 (low Q2, with/without TPE correction)
 * * He4: Charge and magnetization densities from De Jager, At. Data Nucl. Data Tables, 14(1974)
 * * C12: Charge distribution from L. S. Cardman et al., Phys. Lett. B, 91(1970)203
 * * N14: Charge and magnetization densities from De Jager, At. Data Nucl. Data Tables, 14(1974)
 *
 * How to set parameters:
 * H1:
 * If set 1 parameters with SetPars(pars,1), pars[0]=2 means to use 2007 low Q2 fit without TPE correction,
 *   pars[0]=3 means to use 2007 low Q2 fit with TPE correction,
 *   default is to use 2011 global fit with TPE correction;
 * C12:
 * If set 1 parameters with SetPars(pars,1), pars[0]=2 means to use Stansfield's form factors,
 *   default is to use Cardman's fit;
 * Other uses will be considered as invalid.
 */

// History:
//   Mar 2013, C. Gu, First public version, with C12 models.
//   Apr 2013, C. Gu, Add L. Cardman's C12 charge densities.
//   Nov 2013, C. Gu, Add He4 and N14 charge and magnetization densities from D. Jager, original coded by M. Friedman.
//   Apr 2014, C. Gu, Add H1 form factors from J. Arrington.(Thanks to M. Cummings)
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
    // Three-parameter Fermi model

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
    // Harmonic-oscillator model

    return 4 * kPi * x * x * (1 + 1.234 * x * x / 3.0976) * exp(-x * x / 3.0976);
}

static double funcGE_N14(double x, void* data)
{
    double qfm = *(double*) data;
    return kRho0GE_N14 * rhoGE_N14(x, NULL) * sin(qfm * x) / (qfm * x);
}

static double rhoGM_N14(double x, void* data)
{
    // Harmonic-oscillator model

    // double Z = 7, A = 14, a0 = 1.6120;
    // double alpha0 = (Z-2.)/3. = 1.667;
    // double rp     = 0.810; // proton magnetic radius fm (De Jager 1974)
    // double ap     = sqrt((2/3.)*rp*rp) = 0.661;
    // double a      = sqrt((A-1)/A*a0*a0+ap*ap) = 1.688;
    // double alpha  = alpha0*a0*a0/(a*a+3/2.*alpha0*(a*a-a0*a0)) = 1.244;

    return 4 * kPi * x * x * (1 + 1.244582 * x * x / 2.8503) * exp(-x * x / 2.8503);
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

G2PPhysEl::G2PPhysEl() : fSetting(1)
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
        fSetting = int(fPars[0]);
        break;
    default:
        printf("Error: G2PPhysEl::SetPars(): Invalid number of pars.\n");
        break;
    }
}

double G2PPhysEl::GetXS(double Ei, double Ef, double theta)
{
    if ((fZ == 1)&&(fA == 1)) return GetXS_H1(Ei, theta);
    if ((fZ == 2)&&(fA == 4)) return GetXS_He4(Ei, theta);
    if ((fZ == 6)&&(fA == 12)) {
        if (fSetting == 2) return GetXS_All(Ei, theta);
        else return Carbon(Ei, theta);
    }
    if ((fZ == 7)&&(fA == 14)) return GetXS_N14(Ei, theta);

    return GetXS_All(Ei, theta);
}

double G2PPhysEl::GetXS_H1(double Ei, double theta)
{
    // reference: Xiaohui thesis, eq. 1.39 for sigma

    int Z = 1;
    double M = 0.938783; // GeV
    double eP = Ei / (1 + Ei / M * (1 - cos(theta))); // scattered electron energy
    double Q2 = 4 * Ei * eP * sin(theta / 2) * sin(theta / 2); // units GeV2
    double tau = Q2 / 4.0 / M / M;

    double GE = 0.0, GM = 0.0;
    // J. Arrington, Phys. Rev. C, 69(2004)022201
    // Use Rosenbluth form factors in TABLE I because it yields the correct cross sections
    //double x1 = Q2;
    //double x2 = x1*Q2;
    //double x3 = x2*Q2;
    //double x4 = x3*Q2;
    //double x5 = x4*Q2;
    //double x6 = x5*Q2;
    //GE = 1.0 / (1.0 + 3.226 * x1 + 1.508 * x2 - 0.3773 * x3 + 0.611 * x4 - 0.1853 * x5 + 0.01596 * x6);
    //GM = 2.79292 / (1 + 3.19 * x1 + 1.355 * x2 + 0.151 * x3 - 0.0114 * x4 + 5.33e-4 * x5 - 9.00e-6 * x6);

    // J. Arrington, W. Melnitchouk, and J. A. Tjon, Phys. Rev. C, 76(2007)035205
    //double x1 = tau;
    //double x2 = x1*tau;
    //double x3 = x2*tau;
    //double x4 = x3*tau;
    //double x5 = x4*tau;
    //GE = (1.0 + 3.439 * x1 - 1.602 * x2 + 0.068 * x3) / (1 + 15.055 * x1 + 48.061 * x2 + 99.304 * x3 + 0.012 * x4 + 8.650 * x5);
    //GM = 2.79292 * (1.0 - 1.465 * x1 + 1.260 * x2 + 0.262 * x3) / (1 + 9.627 * x1 + 11.179 * x4 + 13.245 * x5);

    switch (fSetting) {
    case 2:
        // J. Arrington and I. Sick, Phys. Rev. C 76(2007)035201
        // Without Two-Photon-Exchange correction
        GE = 1 / (1 + 3.366 * Q2 / (1 - 0.189 * Q2 / (1 - 1.263 * Q2 / (1 + 1.351 * Q2 / (1 - 0.301 * Q2)))));
        GM = 2.79292 / (1 + 3.205 * Q2 / (1 - 0.318 * Q2 / (1 - 1.228 * Q2 / (1 + 5.619 * Q2 / (1 - 1.116 * Q2)))));
        break;
    case 3:
        // J. Arrington and I. Sick, Phys. Rev. C 76(2007)035201
        // With Two-Photon-Exchange correction
        GE = 1 / (1 + 3.478 * Q2 / (1 - 0.140 * Q2 / (1 - 1.311 * Q2 / (1 + 1.128 * Q2 / (1 - 0.233 * Q2)))));
        GM = 2.79292 / (1 + 3.224 * Q2 / (1 - 0.313 * Q2 / (1 - 0.868 * Q2 / (1 + 4.278 * Q2 / (1 - 1.102 * Q2)))));
        break;
    case 1:
    default:
        // S. Venkat, J. Arrington, G. A. Miller and X. Zhan, Phys. Rev. C, 83(2011)015203
        // With Two-Photon-Exchange correction
        double x1 = tau;
        double x2 = x1*tau;
        double x3 = x2*tau;
        double x4 = x3*tau;
        double x5 = x4*tau;
        GE = (1.0 + 2.90966 * x1 - 1.11542229 * x2 + 3.866171e-2 * x3) / (1 + 14.5187212 * x1 + 40.88333 * x2 + 99.999998 * x3 + 4.579e-5 * x4 + 10.3580447 * x5);
        GM = 2.792782 * (1.0 - 1.43573 * x1 + 1.19052066 * x2 + 2.5455841e-1 * x3) / (1 + 9.70703681 * x1 + 3.7357e-4 * x2 + 6.0e-8 * x3 + 9.9527277 * x4 + 12.7977739 * x5);
    }

    double sigma = Z * Z * kAlpha * kAlpha / Q2 * (eP / Ei)*(eP / Ei)*(2 * tau * GM * GM + 1 / tan(theta / 2) / tan(theta / 2) / (1 + tau)*(GE * GE + tau * GM * GM));

    return sigma * kFm * kFm * 1e4; // microbarn
}

double G2PPhysEl::GetXS_He4(double Ei, double theta)
{
    // reference: Xiaohui thesis, eq. 1.39 for sigma

    int Z = 2;
    double M = 3.7284; // GeV
    double eP = Ei / (1 + Ei / M * (1 - cos(theta))); // scattered electron energy
    double Q2 = 4 * Ei * eP * sin(theta / 2) * sin(theta / 2); // units GeV2
    double tau = Q2 / 4.0 / M / M;

    // De Jager, At. Data Nucl. Data Tables, 14(1974)
    double GE = GE_He4(Q2);
    double GM = 0.0;

    double sigma = Z * Z * kAlpha * kAlpha / Q2 * (eP / Ei)*(eP / Ei)*(2 * tau * GM * GM + 1 / tan(theta / 2) / tan(theta / 2) / (1 + tau)*(GE * GE + tau * GM * GM));

    return sigma * kFm * kFm * 1e4; // microbarn
}

double G2PPhysEl::GetXS_N14(double Ei, double theta)
{
    // reference: Xiaohui thesis, eq. 1.39 for sigma

    int Z = 7;
    double M = 13.04378; // GeV
    double eP = Ei / (1 + Ei / M * (1 - cos(theta))); // scattered electron energy
    double Q2 = 4 * Ei * eP * sin(theta / 2) * sin(theta / 2); // units GeV2
    double tau = Q2 / 4.0 / M / M;

    // De Jager, At. Data Nucl. Data Tables, 14(1974)
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
