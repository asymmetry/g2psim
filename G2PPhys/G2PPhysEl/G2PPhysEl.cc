#include <cstdio>
#include <cmath>
#include <vector>

#include "G2PPhysBase.hh"

#include "G2PPhysEl.hh"

using namespace std;

static const double kDEG = 3.14159265358979323846/180.0;
static const double kMEV = 1.0e-3;

static double FF_C12(double Q2)
{
    double XA;
    double Fc12;
    /*
    C     K. Slifer 10/04/02
    C     12C Form Factor
    C     See K. C. Stansfield et al. PRC 3, 1448 (1971)
    */
    double HBARC = 0.197327053; // GEV-FM
    double ZT = 6.0;            // ATOMIC NUMBER|CHARGE

    double XALPHA = (ZT-2.0)/3.0;
    double Q2FM = Q2/(HBARC*HBARC); // fm^-2 // Corrected on Feb 17, 2011

    if (Q2FM<3.2)
        XA = 1.64;              // fm
    else if(Q2FM>3.5)
        XA = 1.68;              // fm
    else
        XA = 0;                 // fm

    Fc12 = 1.0-XALPHA/(2.0*(2.+3.*XALPHA))*Q2FM*XA*XA;
    Fc12 = Fc12*exp(-(Q2FM*XA*XA)/4.0);

    if ((Q2FM>3.2)&&(Q2FM<3.5)) // DIFFRACTION MINIMUM
        Fc12 = 1.0/100000.0;
    return Fc12;
}

static double FF_All(double Q2, double Z)
{
    //input: Q2 in GeV2, Z is atomic number

	double XA = 1.64; //unit in Fm
	double FF;

	double HBARC = 0.197327053; // GEV-FM
	double XALPHA = (Z-2.0)/3.0;
	double Q2FM = Q2/(HBARC*HBARC); // fm^-2 // Corrected on Feb 17, 2011

	FF = 1.0-XALPHA/(2.0*(2.+3.*XALPHA))*Q2FM*XA*XA;
	FF = FF*exp(-(Q2FM*XA*XA)/4.0); 
	if (FF<1.0e-6) FF=1.0e-6;

	return FF;
}

G2PPhysEl::G2PPhysEl()
{
    // Nothing to do
}

G2PPhysEl::~G2PPhysEl()
{
    // Nothing to do
}

void G2PPhysEl::SetPars(double* array, int n)
{
    G2PPhysBase::SetPars(array, n);

    // switch (n) {
    // case 0:
    //     break;
    // case 2:
    //     fTb = fPars[0]; fTa = fPars[1];
    //     break;
    // case 3:
    //     fEPS = fPars[0]; fEPSD = fPars[1]; fFP = fPars[2];
    //     break;
    // case 5:
    //     fEPS = fPars[0]; fEPSD = fPars[1]; fFP = fPars[2];
    //     fTb = fPars[3]; fTa = fPars[4];
    //     break;
    // default:
    //     printf("Error: G2PPhysEl::SetPars(): Invalid number of pars.\n");
    //     break;
    // }
}

double G2PPhysEl::GetXS(double Ei, double Ef, double theta)
{
    if ((iZ==6)&&(iA==12)) return GetXS_C12(Ei, theta);

    return 0.0;
}

double G2PPhysEl::GetXS_C12(double Ei, double theta)
{
    double Recoil, Ef, Q2;
    double Mtg = 11.188;

    /*
    C     Ei        - INCIDENT ENERGY (GEV)
    C     theta     - SCATTERING ANGLE (RADIANS)
    C     Mtg       - TARGET MASS (GEV)
    C     Recoil    - RECOIL FACTOR
    C     Ef        - SCATTERED ENERGY (GEV)
    C     Q2        - MOMENTUM TRANSFER SQUARED (GEV**2)
    C     MOTT      - MOTT CROSS SECTION
    C     FF_C12    - NUCLEAR FORM FACTOR
    */

    Recoil = 1.0/(1.0+Ei/Mtg*(1.-cos(theta)));
    Ef = Recoil*Ei;
    Q2 = 2.*Mtg*(Ei-Ef);

    double Z = 6;
    double alpha = 1.0/137.035989561;
    // Calculated cross sections are often written in units of hbar2c2/GeV2 (approximately 0.3894 mb).
    double hbc2 = 0.38938;  // turn GeV to mbar

    double C = cos(theta/2.0);
    double S = Z*alpha*C/(2.*Ei*(1.-C*C));
    double Mott = S*S*hbc2*1000.0; // microbarn

    double FF = FF_C12(Q2);

    return Recoil*Mott*FF*FF; // microbarn
}

double G2PPhysEl::GetXS_All(double Ei, double theta)
{
    double Recoil, Ef, Q2;
    double Mtg = fTargetMass;

    Recoil = 1.0/(1.0+Ei/Mtg*(1.-cos(theta)));
    Ef = Recoil*Ei;
    Q2 = 2.*Mtg*(Ei-Ef);

    double Z = iZ;
    double alpha = 1.0/137.035989561;
    // Calculated cross sections are often written in units of hbar2c2/GeV2 (approximately 0.3894 mb).
    double hbc2 = 0.38938;  // turn GeV to mbar

    double C = cos(theta/2.0);
    double S = Z*alpha*C/(2.*Ei*(1.-C*C));
    double Mott = S*S*hbc2*1000.0; // microbarn

    double FF = FF_All(Q2, iZ);

    return Recoil*Mott*FF*FF; // microbarn
}
