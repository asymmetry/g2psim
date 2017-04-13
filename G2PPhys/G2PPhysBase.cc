// -*- C++ -*-

/* class G2PPhysBase
 * Abstract base class of G2PPhys classes.
 * It provides interface functions.
 *
 * Elastic cross section models (G2PPhysEl):
 * * All: Form factors from K. C. Stansfield et al., Phys. Rev. C, 3(1971)1448
 * * H1 : Form factors from S. Venkat et al., Phys. Rev. C, 83(2011)015203 (global fit, with TPE correction)
 *                          J. Arrington et al., Phys. Rev. C 76(2007)035201 (low Q2, with/without TPE correction)
 * * He4: Charge and magnetization densities from De Jager, At. Data Nucl. Data Tables, 14(1974)
 * * C12: Charge and magnetization densities from De Jager, At. Data Nucl. Data Tables, 14(1974)
 *        Charge distribution from L. S. Cardman et al., Phys. Lett. B, 91(1970)203
 * * N14: Charge and magnetization densities from De Jager, At. Data Nucl. Data Tables, 14(1974)
 *
 * Inelastic cross section models:
 * * G2PPhysEPC: J. W. Lightbody et al, Computers in Physics, 2(1988)57
 * * G2PPhysPB: P. E. Bosted et al, Phys. Rev. C, 78(2008)015202 and arXiv:1203.2262
 * * G2PPhysQFS: J. W. Lightbody et al, Computers in Physics, 2(1988)57
 * * G2PPhysWISER: D. E. Wiser, Ph.D. Thesis
 *
 * Radiative correction added for P. Bosted model and QFS.
 *
 * Meaning of parameters:
 * fPID: incident particle ID, following the PDG definition:
 *       2212 for p        ;   2112 for n     ;   211 for pi+   ;
 *       -211 for pi-      ;   111  for pi0   ;   11  for e-    ;
 *       22   for photon   ;
 * fZ, fA: proton and mass number of the nucleus.
 *
 * Please also read headers of QFS, PBosted, EPC and WISER models. They contains very important usage information!
 *
 * Unit of cross section is ub/MeV-sr.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Nov 2016, C. Gu, Rewrite parameter parser.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>

#include "G2PPhysBase.hh"

using namespace std;

static const double kU = 0.93149406121; // MeV

G2PPhysBase::G2PPhysBase() : fZ(1), fA(1), fTargetMass(0.0), fPID(11)
{
    // Nothing to do
}

G2PPhysBase::~G2PPhysBase()
{
    // Nothing to do
}

void G2PPhysBase::SetTarget(int Z, int A)
{
    fZ = Z;
    fA = A;

    if (fabs(fTargetMass) < 1.0e-8) SetTargetMass();
}

void G2PPhysBase::SetTargetMass(double value)
{
    fTargetMass = value;
}

void G2PPhysBase::SetParticle(int pid)
{
    fPID = pid;
}

void G2PPhysBase::SetTargetMass()
{
    double atomicmass = 0.0;

    if (fZ == 1) atomicmass = 1.00794;
    else if (fZ == 2) atomicmass = 4.002602;
    else if (fZ == 6) atomicmass = 12.0107;
    else if (fZ == 7) atomicmass = 14.0067;
    else if (fZ == 8) atomicmass = 15.9994;
    else if (fZ == 13) atomicmass = 26.982;
    else if (fZ == 26) atomicmass = 55.845;
    else if (fZ == 29) atomicmass = 63.546;
    else if (fZ == 74) atomicmass = 183.84;
    else if (fabs(atomicmass) < 1.0e-8) atomicmass = fA; // only an estimation

    fTargetMass = atomicmass * kU;
}
