// -*- C++ -*-

/* class G2PRun
 * This file defines a class G2PRun.
 * It describes the procedure of a simulation.
 * G2PProcBase classes and their input variables should be registered here.
 */

// History:
//   Jan 2013, C. Gu, First public version.
//   May 2013, C. Gu, Add G2PProcBase classes, G2PRun uses these classes to describe the simulation procedure.
//

#include <cstdlib>
#include <cmath>
#include <vector>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"

#include "G2PBwdProc.hh"
#include "G2PFwdProc.hh"
#include "G2PGlobals.hh"
#include "G2PGunProc.hh"
#include "G2PRunBase.hh"
#include "G2PPhyProc.hh"
#include "G2PProcBase.hh"

#include "G2PRun.hh"

using namespace std;

static const double kDEG = 3.14159265358979323846 / 180.0;

G2PRun::G2PRun() {
    // Nothing to do
}

G2PRun::~G2PRun() {
    // Nothing to do
}

int G2PRun::Init() {
    //static const char* const here = "Init()";

    // Notice: the order is important
    fProcs->Add(new G2PGunProc());
    fProcs->Add(new G2PFwdProc());
    fProcs->Add(new G2PBwdProc());
    fProcs->Add(new G2PPhyProc());

    vector<const char*> gunreqs;
    gunreqs.clear();
    fProcReqs.push_back(gunreqs);

    vector<const char*> fwdreqs;
    fwdreqs.clear();
    fwdreqs.push_back("fV5tg_tr");
    fwdreqs.push_back("fV5react_lab");
    fProcReqs.push_back(fwdreqs);

    vector<const char*> bwdreqs;
    bwdreqs.clear();
    bwdreqs.push_back("fV5bpm_bpm");
    bwdreqs.push_back("fV5projtg_tr");
    bwdreqs.push_back("fV5fp_tr");
    fProcReqs.push_back(bwdreqs);

    vector<const char*> phyreqs;
    phyreqs.clear();
    phyreqs.push_back("fV5beam_lab");
    phyreqs.push_back("fV5bpm_lab");
    phyreqs.push_back("fV5react_tr");
    phyreqs.push_back("fV5rec_tr");
    fProcReqs.push_back(phyreqs);

    if (G2PRunBase::Init() != 0) return fStatus;

    return (fStatus = kOK);
}

// double G2PRun::DriftPath()
// {
//     if (bUseField) {
//         double x[3] = { fV5bpm_lab[0], fV5bpm_lab[2], fV5bpm_lab[4] };
//         double p[3] = { fBeamEnergy*sin(fV5bpm_lab[1])*cos(fV5bpm_lab[3]),
//                         fBeamEnergy*sin(fV5bpm_lab[1])*sin(fV5bpm_lab[3]),
//                         fBeamEnergy*cos(fV5bpm_lab[1]) };

//         double xrec[3];
//         HRSTransTCSNHCS::X_TCS2HCS(fV5rec_tr[0], fV5rec_tr[2], 0.0, fHRSAngle, xrec[0], xrec[1], xrec[2]);
//         double theta,phi;
//         HRSTransTCSNHCS::P_TCS2HCS(fV5rec_tr[1], fV5rec_tr[3], fHRSAngle, theta, phi);
//         double prec[3] = { fHRSMomentum*(1+fV5rec_tr[4])*sin(theta)*cos(phi),
//                            fHRSMomentum*(1+fV5rec_tr[4])*sin(theta)*sin(phi),
//                            fHRSMomentum*(1+fV5rec_tr[4])*cos(theta) };
//         pDrift->Drift(xrec, prec, 0, 10.0, xrec, prec);

//         double z = 0.0;
//         double zmin = 0.0;
//         double distmin = sqrt((x[0]-xrec[0])*(x[0]-xrec[0])+(x[1]-xrec[1])*(x[1]-xrec[1]));
//         for (int i = 1; i<150; i++) {
//             z = 0.1e-3*i;
//             pDrift->Drift(x, p, z, 10.0, x, p);
//             pDrift->Drift(xrec, prec, z, 10.0, xrec, prec);
//             //printf("%e\t%e\t%e\t%e\t%e\n", z, x[0], x[1], xrec[0], xrec[1]);
//             double distance = sqrt((x[0]-xrec[0])*(x[0]-xrec[0])+(x[1]-xrec[1])*(x[1]-xrec[1]));
//             if (distance<distmin) {
//                 zmin = z;
//                 distmin = distance;
//             }
//         }

//         // for (int i = 20; i<200; i++) {
//         //     z = 1e-3*i;
//         //     pDrift->Drift(xrec, prec, z, 10.0, xrec, prec);
//         //     printf("%e\t%e\t%e\t%e\t%e\n", z, 1000.0, 1000.0, xrec[0], xrec[1]);
//         // }

//         for (int i = 1; i<150; i++) {
//             z = -0.1e-3*i;
//             pDrift->Drift(x, p, z, 10.0, x, p);
//             pDrift->Drift(xrec, prec, z, 10.0, xrec, prec);
//             //printf("%e\t%e\t%e\t%e\t%e\n", z, x[0], x[1], xrec[0], xrec[1]);
//             double distance = sqrt((x[0]-xrec[0])*(x[0]-xrec[0])+(x[1]-xrec[1])*(x[1]-xrec[1]));
//             if (distance<distmin) {
//                 zmin = z;
//                 distmin = distance;
//             }
//         }

//         // for (int i = -20; i>-200; i--) {
//         //     z = 1e-3*i;
//         //     pDrift->Drift(x, p, z, 10.0, x, p);
//         //     printf("%e\t%e\t%e\t%e\t%e\n", z, x[0], x[1], 1000.0, 1000.0);
//         // }

//         pDrift->Drift(x, p, zmin, 10.0, x, p);
//         pDrift->Drift(xrec, prec, zmin, 10.0, xrec, prec);

//         double x_tr,y_tr,z_tr;
//         HRSTransTCSNHCS::X_HCS2TCS(xrec[0], xrec[1], xrec[2], fHRSAngle, x_tr, y_tr, z_tr);
//         return (z_tr);
//     }
//     else return 0;
// }

ClassImp(G2PRun)
