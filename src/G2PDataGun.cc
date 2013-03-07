#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"

#include "G2PDrift.hh"
#include "G2PFieldBase.hh"
#include "G2PGlobals.hh"
#include "G2PGunBase.hh"
#include "G2PRand.hh"
#include "G2PRunBase.hh"
#include "G2PSieve.hh"

#include "G2PDataGun.hh"

using namespace std;

G2PDataGun::G2PDataGun(const char* filename) :
    bIsOptics(false), pDataFileName(filename),
    fTargetMass(0.0), fEnergyLoss(0.0), pfGun(NULL)
{
    fData.clear();
}

G2PDataGun::~G2PDataGun()
{
    // Nothing to do
}

G2PAppsBase::EStatus G2PDataGun::Init()
{
    static const char* const here = "Init()";

    if (G2PGunBase::Init()) return fStatus;

    if (pDrift->GetField()) bIsOptics = false;

    if (bIsOptics) {
        pfGun = &G2PDataGun::ShootOptics;
        fTargetMass = gG2PRun->GetTargetMass();
        fEnergyLoss = gG2PRun->GetEnergyLoss();
    }
    else pfGun = &G2PDataGun::ShootNormal;

    if (!LoadData()) {
        Error(here, "Cannot initialize, cannot read data.");
        return (fStatus = kINITERROR);
    }

    SetSieve(fHRSAngle);

    return (fStatus = kOK);
}

bool G2PDataGun::ShootNormal(double* V5bpm_lab, double* V5react_tr, double* V5fp_tr)
{
    static const char* const here = "Shoot()";

    bool noerror = true;

    sData tempdata;

    if (fData.empty()) return false;
    tempdata = fData.back();
    V5bpm_lab[0] = tempdata.xb;
    V5bpm_lab[1] = tempdata.tb;
    V5bpm_lab[2] = tempdata.yb;
    V5bpm_lab[3] = tempdata.pb;
    V5bpm_lab[4] = tempdata.zb;
    V5fp_tr[0] = tempdata.xf;
    V5fp_tr[1] = tempdata.tf;
    V5fp_tr[2] = tempdata.yf;
    V5fp_tr[3] = tempdata.pf;
    V5fp_tr[4] = 0.0;

    if (fDebug>1) {
        Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e", V5bpm_lab[0], V5bpm_lab[1], V5bpm_lab[2], V5bpm_lab[3], V5bpm_lab[4]);
        Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e", V5fp_tr[0], V5fp_tr[1], V5fp_tr[2], V5fp_tr[3], V5fp_tr[4]);
    }

    fData.pop_back();
    return noerror;
}

bool G2PDataGun::ShootOptics(double* V5bpm_lab, double* V5react_tr, double* V5fp_tr)
{
    static const char* const here = "Shoot()";

    sData tempdata;

    if (fData.empty()) return false;
    tempdata = fData.back();
    int index = tempdata.ind;
    V5bpm_lab[0] = tempdata.xb;
    V5bpm_lab[1] = tempdata.tb;
    V5bpm_lab[2] = tempdata.yb;
    V5bpm_lab[3] = tempdata.pb;
    V5bpm_lab[4] = tempdata.zb;
    V5fp_tr[0] = tempdata.xf;
    V5fp_tr[1] = tempdata.tf;
    V5fp_tr[2] = tempdata.yf;
    V5fp_tr[3] = tempdata.pf;
    V5fp_tr[4] = 0.0;

    int col = index/(fSieve.nRow);
    int row = index%(fSieve.nRow);

    double Xreact_tr, Yreact_tr, Zreact_tr;
    HCS2TCS(V5bpm_lab[0], V5bpm_lab[2], V5bpm_lab[4], fHRSAngle, Xreact_tr, Yreact_tr, Zreact_tr);

    double V3sieve_tr[3];
    double V3pd_tr[3];

    V3sieve_tr[0] = fSieve.fXOffset+fSieve.fX[row];
    V3sieve_tr[1] = fSieve.fYOffset+fSieve.fY[col];
    V3sieve_tr[2] = fSieve.fZ;

    V3pd_tr[0] = V3sieve_tr[0]-Xreact_tr;
    V3pd_tr[1] = V3sieve_tr[1]-Yreact_tr;
    V3pd_tr[2] = V3sieve_tr[2]-Zreact_tr;

    double Thetareact_tr = atan(V3pd_tr[0]/V3pd_tr[2]);
    double Phireact_tr = atan(V3pd_tr[1]/V3pd_tr[2]);

    // Calculate delta based on angle
    double V3pd_lab[3];
    TCS2HCS(V3pd_tr[0], V3pd_tr[1], V3pd_tr[2], fHRSAngle, V3pd_lab[0], V3pd_lab[1], V3pd_lab[2]);

    double cosscatangle = V3pd_lab[2]/(sqrt(V3pd_lab[0]*V3pd_lab[0]+V3pd_lab[1]*V3pd_lab[1]+V3pd_lab[2]*V3pd_lab[2]));

    double scatmom = (fTargetMass*fBeamEnergy)/(fTargetMass+fBeamEnergy-fBeamEnergy*cosscatangle);

    double Delta = scatmom/fHRSMomentum-1-fEnergyLoss/fHRSMomentum;

    V5react_tr[0] = Xreact_tr;
    V5react_tr[1] = Thetareact_tr;
    V5react_tr[2] = Yreact_tr;
    V5react_tr[3] = Phireact_tr;
    V5react_tr[4] = Delta;

    if (fDebug>1) {
        Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e", V5bpm_lab[0], V5bpm_lab[1], V5bpm_lab[2], V5bpm_lab[3], V5bpm_lab[4]);
        Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e", V5fp_tr[0], V5fp_tr[1], V5fp_tr[2], V5fp_tr[3], V5fp_tr[4]);
    }

    return true;
}

bool G2PDataGun::LoadData()
{
    FILE *fp;

    if ((fp=fopen(pDataFileName, "r"))==NULL) return false;

    sData temp;
    fscanf(fp, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf", &temp.ind, &temp.xb, &temp.tb, &temp.yb, &temp.pb, &temp.zb, &temp.xf, &temp.tf, &temp.yf, &temp.pf);
    while (!feof(fp)) {
        fData.push_back(temp);
        fscanf(fp, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf", &temp.ind, &temp.xb, &temp.tb, &temp.yb, &temp.pb, &temp.zb, &temp.xf, &temp.tf, &temp.yf, &temp.pf);
    }

    fclose(fp);

    if (!fData.empty()) return true;
    else return false;
}

ClassImp(G2PDataGun)
