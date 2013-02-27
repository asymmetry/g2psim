#include <cstdio>
#include <cstdlib>

#include "TROOT.h"
#include "TObject.h"

#include "G2PDrift.hh"
#include "HRSTransTCSNHCS.hh"

#include "G2PBPM.hh"

ClassImp(G2PBPM);

G2PBPM::G2PBPM()
{
    // Nothing to do
}

G2PBPM::~G2PBPM()
{
    // Nothing to do
}

void G2PBPM::GetBPMValue(const double* V5beam_lab, double* V5bpm_lab)
{
    V5bpm_lab[0] = V5beam_lab[0];
    V5bpm_lab[1] = V5beam_lab[1];
    V5bpm_lab[2] = V5beam_lab[2];
    V5bpm_lab[3] = V5beam_lab[3];
    V5bpm_lab[4] = V5beam_lab[4];
}

// void G2PGun::GetBPMValue(const double* V5beam_lab, double* V5bpm_lab)
// {
//     double x[3] = { V5beam_lab[0], V5beam_lab[2], V5beam_lab[4] };
//     double p[3] = { fBeamEnergy*sin(V5beam_lab[1])*cos(V5beam_lab[3]),
//                     fBeamEnergy*sin(V5beam_lab[1])*sin(V5beam_lab[3]),
//                     fBeamEnergy*cos(V5beam_lab[1]) };

//     if (bUseField) {
//         double z_lab = x[2];
//         G2PDrift::Drift(x, p, fBPMZ_lab, 10.0, x, p);

//         double theta_tr = G2PRand::Gaus(atan(p[0]/p[2]), fBPMAngRes);
//         double phi_tr = G2PRand::Gaus(atan(p[1]/p[2]), fBPMAngRes);
//         p[0] = p[2]*tan(theta_tr);
//         p[1] = p[2]*tan(phi_tr);
//         double normF = fBeamEnergy/sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
//         p[0] *= normF;
//         p[1] *= normF;
//         p[2] *= normF;
//         x[0] = G2PRand::Gaus(x[0], fBPMPosRes);
//         x[1] = G2PRand::Gaus(x[1], fBPMPosRes);
//         G2PDrift::Drift(x, p, z_lab, 10.0, x, p);
//         V5bpm_lab[0] = x[0];
//         V5bpm_lab[1] = acos(p[2]/fBeamEnergy);
//         V5bpm_lab[2] = x[1];
//         V5bpm_lab[3] = atan(p[1]/p[0]);
//         V5bpm_lab[4] = x[2];
//     }
//     else {
//         double theta_tr = atan(p[0]/p[2]);
//         double phi_tr = atan(p[1]/p[2]);
//         double z_lab = x[2];
//         HRSTransTCSNHCS::Project(x[0], x[1], x[2], fBPMZ_lab-x[2], theta_tr, phi_tr);
//         theta_tr = G2PRand::Gaus(theta_tr, fBPMAngRes);
//         phi_tr = G2PRand::Gaus(phi_tr, fBPMAngRes);
//         x[0] = G2PRand::Gaus(x[0], fBPMPosRes);
//         x[1] = G2PRand::Gaus(x[1], fBPMPosRes);
//         HRSTransTCSNHCS::Project(x[0], x[1], x[2], z_lab-x[2], theta_tr, phi_tr);
//         p[0] = p[2]*tan(theta_tr);
//         p[1] = p[2]*tan(phi_tr);
//         double normF = fBeamEnergy/sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
//         p[0] *= normF;
//         p[1] *= normF;
//         p[2] *= normF;
//         V5bpm_lab[0] = x[0];
//         V5bpm_lab[1] = acos(p[2]/fBeamEnergy);
//         V5bpm_lab[2] = x[1];
//         V5bpm_lab[3] = atan(p[1]/p[0]);
//         V5bpm_lab[4] = x[2];
//     }

// #ifdef GUN_DEBUG
//     printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5bpm_lab[0], V5bpm_lab[1], V5bpm_lab[2], V5bpm_lab[3], V5bpm_lab[4]);
// #endif
// }
