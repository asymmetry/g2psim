#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "TROOT.h"
#include "TObject.h"

#include "G2PDrift.hh"
#include "HRSTransTCSNHCS.hh"
#include "G2PRand.hh"

#include "G2PBPM.hh"

ClassImp(G2PBPM);

G2PBPM::G2PBPM()
    :bIsInit(false), iSetting(0), bUseField(false), fBeamEnergy(2.254)
{
    // Nothing to do
}

G2PBPM::~G2PBPM()
{
    // Nothing to do
}

void G2PBPM::Init()
{
    if (G2PDrift::HasField()) {
        bUseField = true;
    }
    SetBPM();
    bIsInit = true;
}

void G2PBPM::GetBPMValue(const double* V5beam_lab, double* V5bpm_lab)
{
    double x[3] = { V5beam_lab[0], V5beam_lab[2], V5beam_lab[4] };
    double p[3] = { fBeamEnergy*sin(V5beam_lab[1])*cos(V5beam_lab[3]),
                    fBeamEnergy*sin(V5beam_lab[1])*sin(V5beam_lab[3]),
                    fBeamEnergy*cos(V5beam_lab[1]) };

    if (bUseField) {
        G2PDrift::Drift(x, p, fBPMAZ, 10.0, x, p);
        double res = fBPMBRes/(fBPMBZ-fBPMAZ);
        // here theta phi are special coords defined by Pengjia
        // theta is dy/dz in hall coords
        // phi is dx/dz in hall coords
        double theta = G2PRand::Gaus(atan(p[1]/p[2]), res);
        double phi = G2PRand::Gaus(atan(p[0]/p[2]), res);
        p[1] = p[2]*tan(theta);
        p[0] = p[2]*tan(phi);
        double normF = fBeamEnergy/sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
        p[0] *= normF;
        p[1] *= normF;
        p[2] *= normF;
        x[0] = G2PRand::Gaus(x[0], fBPMARes);
        x[1] = G2PRand::Gaus(x[1], fBPMARes);
        G2PDrift::Drift(x, p, 0.0, 10.0, x, p);
        V5bpm_lab[0] = x[0];
        V5bpm_lab[1] = atan(p[1]/p[2]);
        V5bpm_lab[2] = x[1];
        V5bpm_lab[3] = atan(p[0]/p[2]);
        V5bpm_lab[4] = x[2];
    }
    else {
        double theta_tr = atan(p[0]/p[2]);
        double phi_tr = atan(p[1]/p[2]);
        double z_lab = x[2];
        HRSTransTCSNHCS::Project(x[0], x[1], x[2], fBPMAZ-x[2], theta_tr, phi_tr);
        double res = fBPMBRes/(fBPMBZ-fBPMAZ);
        theta_tr = G2PRand::Gaus(theta_tr, res);
        phi_tr = G2PRand::Gaus(phi_tr, res);
        x[0] = G2PRand::Gaus(x[0], fBPMARes);
        x[1] = G2PRand::Gaus(x[1], fBPMARes);
        HRSTransTCSNHCS::Project(x[0], x[1], x[2], z_lab-x[2], theta_tr, phi_tr);
        p[0] = p[2]*tan(theta_tr);
        p[1] = p[2]*tan(phi_tr);
        double normF = fBeamEnergy/sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
        p[0] *= normF;
        p[1] *= normF;
        p[2] *= normF;
        V5bpm_lab[0] = x[0];
        V5bpm_lab[1] = atan(p[1]/p[2]);
        V5bpm_lab[2] = x[1];
        V5bpm_lab[3] = atan(p[0]/p[2]);
        V5bpm_lab[4] = x[2];
    }
}

void G2PBPM::SetBPM()
{
    fBPMAZ = -956e-3;
    fBPMBZ = -690e-3;
    fBPMZ = (fBPMAZ+fBPMBZ)/2.0;
}
