#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <map>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"

#include "G2PAppBase.hh"
#include "G2PGlobals.hh"
#include "G2PPhys.hh"
#include "G2PProcBase.hh"
#include "G2PRunBase.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PPhyProc.hh"

static const double kDEG = 3.14159265358979323846/180.0;

G2PPhyProc::G2PPhyProc() :
    fBeamEnergy(0.0), fHRSAngle(0.0), fHRSMomentum(0.0), 
    pPhys(NULL)
{
    Clear();
}

G2PPhyProc::~G2PPhyProc()
{
    // Nothing to do
}

int G2PPhyProc::Init()
{
    static const char* const here = "Init()";

    if (G2PProcBase::Init()!=0) return fStatus;

    pPhys = G2PPhys::GetInstance();
    if (!pPhys) {
        Error(here, "Cannot initialize, no G2PPhys found.");
        return (fStatus = kINITERROR);
    }

    fApps->Add(pPhys);

    return (fStatus = kOK);
}

int G2PPhyProc::Begin()
{
    //static const char* const here = "Begin()";

    if (G2PProcBase::Begin()!=0) return fStatus;

    fBeamEnergy = gG2PRun->GetBeamEnergy();
    fHRSAngle = gG2PRun->GetHRSAngle();
    fHRSMomentum = gG2PRun->GetHRSMomentum();

    mName["fV5beam_lab"] = fV5beam_lab; mLength["fV5beam_lab"] = 5;
    mName["fV5bpm_lab"] = fV5bpm_lab; mLength["fV5bpm_lab"] = 5;
    mName["fV5react_tr"] = fV5react_tr; mLength["fV5react_tr"] = 5;
    mName["fV5rec_tr"] = fV5rec_tr; mLength["fV5rec_tr"] = 5;

    return (fStatus = kOK);
}

int G2PPhyProc::Process()
{
    static const char* const here = "Process()";

    //double V5[5];

    fXSinit = CalXS(fV5beam_lab, fV5react_tr, fThinit);
    fXSrec = CalXS(fV5bpm_lab, fV5rec_tr, fThrec);

    if (fDebug>1) {
        Info(here, "phys_init : %10.3e %10.3e", fThinit/kDEG, fXSinit);
        Info(here, "phys_rec  : %10.3e %10.3e", fThrec/kDEG, fXSrec);
    }

    return 0;
}

void G2PPhyProc::Clear()
{
    memset(fV5beam_lab, 0, sizeof(fV5beam_lab));
    memset(fV5bpm_lab, 0, sizeof(fV5bpm_lab));
    memset(fV5react_tr, 0, sizeof(fV5react_tr));
    memset(fV5rec_tr, 0, sizeof(fV5rec_tr));

    fThinit = 0.0; fXSinit = 0.0;
    fThrec = 0.0; fXSrec = 0.0;
}

int G2PPhyProc::DefineVariables(EMode mode)
{
    if (mode==kDefine&&bIsSetup) return 0;
    bIsSetup = (mode==kDefine);

    VarDef vars[] = {
        { "react.angle", "Real scattering angle", kDouble, &fThinit },
        { "react.xs", "Cross section with real kins", kDouble, &fXSinit },
        { "rec.angle", "Rec scattering angle", kDouble, &fThrec },
        { "rec.xs", "Cross section with rec kins", kDouble, &fXSrec },
        { 0 }
    };

    return DefineVarsFromList(vars, mode);
}

void G2PPhyProc::MakePrefix()
{
    const char* basename = "phy";

    G2PAppBase::MakePrefix(basename);
}

double G2PPhyProc::CalXS(const double* V5lab, const double* V5tr, double& scatangle)
{
    double Eb[3] = { sin(V5lab[1])*cos(V5lab[3]),
                     sin(V5lab[1])*sin(V5lab[3]),
                     cos(V5lab[1]) };

    double theta, phi;
    TCS2HCS(V5tr[1], V5tr[3], fHRSAngle, theta, phi);

    double Ef[3] = { sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) };

    scatangle = acos(Eb[0]*Ef[0]+Eb[1]*Ef[1]+Eb[2]*Ef[2]);

    double Ebval = fBeamEnergy;
    double Efval = (1+V5tr[4])*fHRSMomentum;

    return pPhys->GetXS(Ebval, Efval, scatangle);
}

ClassImp(G2PPhyProc)
