// -*- C++ -*-

/* class G2PDrift
 * Use Nystrom-Runge-Kutta method to derive the trajectory of a charged particle in static magnetic field.
 * G2PProcBase classes will call Drift() to get the end point position and momentum of the trajectory.
 * Drift() has 2 prototypes, one for lab coordinates and one for HRS transportation coordinates.
 */

// History:
//   Feb 2013, C. Gu, First public version.
//   Feb 2013, C. Gu, Change algorithm to Nystrom-Runge-Kutta method.
//   Mar 2013, C. Gu, Add flexible step length and boundary check.
//   Oct 2013, J. Liu, Add drift function to stop at a cylinder boundary.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PField.hh"
#include "G2PGlobals.hh"

#include "G2PDrift.hh"

using namespace std;

static const double c = 2.99792458e8;
static const double e = 1.60217656535e-19;
static const double kDEG = 3.14159265358979323846 / 180.0;
static const double kGEV = 1.0e9 * e / c / c;
static const double kOneSixth = 1.0 / 6.0;

G2PDrift* G2PDrift::pG2PDrift = NULL;

G2PDrift::G2PDrift() :
fM0(0.51099892811e-3), fQ(-1 * e), fQSave(-1 * e), fStep(1.0e-3), fStepLimit(1.0e-6), fErrLoLimit(1.0e-7), fErrHiLimit(1.0e-6), fVelocity(0.0), fVelocity2(0.0), fGamma(0.0), fCof(0.0), pField(NULL), pfDriftHCS(NULL), pfDriftTCS(NULL)
{
    if (pG2PDrift) {
        Error("G2PDrift()", "Only one instance of G2PDrift allowed.");
        MakeZombie();
        return;
    }
    pG2PDrift = this;

    Clear();
}

G2PDrift::~G2PDrift()
{
    if (pG2PDrift == this) pG2PDrift = NULL;
}

int G2PDrift::Init()
{
    //static const char* const here = "Init()";

    if (G2PAppBase::Init() != 0) return fStatus;

    pField = static_cast<G2PField*> (gG2PApps->Find("G2PField"));

    return (fStatus = kOK);
}

int G2PDrift::Begin()
{
    //static const char* const here = "Begin()";

    if (G2PAppBase::Begin() != 0) return fStatus;

    if (pField) {
        pfDriftHCS = &G2PDrift::DriftHCS;
        pfDriftTCS = &G2PDrift::DriftTCS;
        pfDriftTCSL = &G2PDrift::DriftTCSL;
        pfDriftTCSV = &G2PDrift::DriftTCSV;
    } else {
        pfDriftHCS = &G2PDrift::DriftHCSNF;
        pfDriftTCS = &G2PDrift::DriftTCSNF;
    }

    return (fStatus = kOK);
}

void G2PDrift::Clear(Option_t* option)
{
    memset(fField, 0, sizeof (fField));
    memset(fIPoint, 0, sizeof (fIPoint));
    memset(fMPoint, 0, sizeof (fMPoint));
    memset(fEPoint, 0, sizeof (fEPoint));

    G2PAppBase::Clear(option);
}

double G2PDrift::Drift(const double* x, const double* p, double zlimit, double llimit, double *xout, double *pout)
{
    // Drift in lab coordinate

    static const char* const here = "Drift()";

    double xx[3] = {x[0], x[1], x[2]};
    double pp[3] = {p[0], p[1], p[2]};

    double result = (this->*pfDriftHCS)(x, p, zlimit, llimit, xout, pout);

    if (fDebug > 2) {
        Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e", xx[0], xx[1], xx[2], pp[0], pp[1], pp[2], xout[0], xout[1], xout[2], pout[0], pout[1], pout[2]);
    }

    return result;
}

double G2PDrift::Drift(const double* x, double p, double z_tr, double angle, double zlimit, double llimit, double* xout)
{
    // Drift in transport coordinate

    static const char* const here = "Drift()";

    double xx[5] = {x[0], x[1], x[2], x[3], x[4]};

    double result = (this->*pfDriftTCS)(x, p, z_tr, angle, zlimit, llimit, xout);

    if (fDebug > 2) {
        Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e", xx[0], xx[1], xx[2], xx[3], xx[4], xout[0], xout[1], xout[2], xout[3], xout[4]);
    }

    return result;
}

double G2PDrift::Drift(const double* x, double p, double z_tr, double angle, double zlimit, double llimit, double rlimit, double* xout, int cylinder_type)
{
    // J. Liu: added for cylinder boundary, like target
    // vertical cylinder    : cylinder_type = 1;
    // longitudinal cylinder: cylinder_type = 0;

    static const char* const here = "Drift()";

    double xx[5] = {x[0], x[1], x[2], x[3], x[4]};

    double result = 0;
    if (cylinder_type == 0) {
        result = (this->*pfDriftTCSL)(x, p, z_tr, angle, zlimit, llimit, rlimit, xout);
    }
    if (cylinder_type == 1) {
        result = (this->*pfDriftTCSV)(x, p, z_tr, angle, zlimit, llimit, rlimit, xout);
    }

    if (fDebug > 2) {
        Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e", xx[0], xx[1], xx[2], xx[3], xx[4], xout[0], xout[1], xout[2], xout[3], xout[4]);
    }

    return result;
}

void G2PDrift::SetLimit(double lo, double hi)
{
    fErrLoLimit = lo;
    fErrHiLimit = hi;

    fConfigIsSet.insert((unsigned long) &fErrLoLimit);
    fConfigIsSet.insert((unsigned long) &fErrHiLimit);
}

double G2PDrift::DriftHCS(const double* x, const double* p, double zlimit, double llimit, double *xout, double *pout)
{
    // Drift in lab coordinate

    double xi[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double xf[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    double M = sqrt(fM0 * fM0 + p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);

    double v[3] = {p[0] * c / M, p[1] * c / M, p[2] * c / M};
    fVelocity2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    fVelocity = sqrt(fVelocity2);
    fGamma = 1 / sqrt(1 - (fVelocity2 / c / c));

    xi[0] = x[0];
    xi[1] = x[1];
    xi[2] = x[2];
    xi[3] = v[0];
    xi[4] = v[1];
    xi[5] = v[2];

    // A trick to drift along the whole trajectory
    // 4 different situation:
    // z > zlimit, p > 0: true  ^ false = true  (backward)
    // z > zlimit, p < 0: true  ^ true  = false (forward)
    // z < zlimit, p > 0: false ^ false = false (forward)
    // z < zlimit, p < 0: false ^ true  = true  (backward)
    bool sign = (x[2] > zlimit) ? true : false;
    if ((p[2] < 0)^(sign)) {
        xi[3] *= -1.0;
        xi[4] *= -1.0;
        xi[5] *= -1.0;
        fQ = -1.0 * fQSave;
    }

    double dxdt[6], err[6];
    double dt = fStep / fVelocity;
    double dtlimit = fStepLimit / fVelocity;
    double error = (fErrHiLimit + fErrLoLimit) / 2.0;
    double l = 0.0;
    bool flag = true;

    for (int i = 0; i < 6; i++) xf[i] = xi[i];
    while ((l < llimit)&&((xf[2] < zlimit)^(sign))) {
        if ((error > fErrHiLimit)&&(dt > dtlimit)) {
            dt /= 2.0;
            flag = false;
        } else {
            for (int i = 0; i < 6; i++) xi[i] = xf[i];
            if (error < fErrLoLimit) dt *= 2.0;
            l += fVelocity*dt;
            flag = true;
        }
        ComputeRHS(xi, dxdt);
        NystromRK4(xi, dxdt, dt, xf, err);
        error = DistChord();
    }

    if (!flag) l += fVelocity * dt;
    if (l < llimit) {
        while (dt > dtlimit) {
            if ((xf[2] > zlimit)^(sign)) {
                l = l - fVelocity*dt;
                dt /= 2.0;
            } else {
                for (int i = 0; i < 6; i++) xi[i] = xf[i];
            }
            l += fVelocity*dt;
            ComputeRHS(xi, dxdt);
            NystromRK4(xi, dxdt, dt, xf, err);
        }
    }

    if ((p[2] < 0)^(sign)) {
        xf[3] *= -1.0;
        xf[4] *= -1.0;
        xf[5] *= -1.0;
        fQ = fQSave;
    }

    xout[0] = xf[0];
    xout[1] = xf[1];
    xout[2] = xf[2];
    pout[0] = xf[3] * M / c;
    pout[1] = xf[4] * M / c;
    pout[2] = xf[5] * M / c;

    return l;
}

double G2PDrift::DriftHCSNF(const double* x, const double* p, double zlimit, double llimit, double *xout, double *pout)
{
    xout[0] = x[0]+(zlimit - x[2]) * p[0] / p[2];
    xout[1] = x[1]+(zlimit - x[2]) * p[1] / p[2];
    xout[2] = zlimit;

    pout[0] = p[0];
    pout[1] = p[1];
    pout[2] = p[2];

    double dx = xout[0] - x[0];
    double dy = xout[1] - x[1];
    double dz = xout[2] - x[2];

    return sqrt(dx * dx + dy * dy + dz * dz);
}

double G2PDrift::DriftTCS(const double* x, double p, double z_tr, double angle, double zlimit, double llimit, double* xout)
{
    // Drift in transport coordinate

    double xi[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double xf[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    double sinang = sin(angle);
    double cosang = cos(angle);

    TCS2HCS(x[0], x[2], z_tr, angle, xi[0], xi[1], xi[2]);
    double theta, phi;
    TCS2HCS(x[1], x[3], angle, theta, phi);

    double pp = (1 + x[4]) * p;
    xi[3] = pp * sin(theta) * cos(phi);
    xi[4] = pp * sin(theta) * sin(phi);
    xi[5] = pp * cos(theta);
    double M = sqrt(fM0 * fM0 + pp * pp);
    xi[3] *= c / M;
    xi[4] *= c / M;
    xi[5] *= c / M;
    fVelocity2 = xi[3] * xi[3] + xi[4] * xi[4] + xi[5] * xi[5];
    fVelocity = sqrt(fVelocity2);
    fGamma = 1 / sqrt(1 - (fVelocity2 / c / c));

    bool sign = (z_tr > zlimit) ? true : false;
    if (sign) {
        xi[3] *= -1.0;
        xi[4] *= -1.0;
        xi[5] *= -1.0;
        fQ = -1.0 * fQSave;
    }

    double dxdt[6], err[6];
    double error = (fErrHiLimit + fErrLoLimit) / 2.0;
    double dt = fStep / fVelocity;
    double dtlimit = fStepLimit / fVelocity;
    double l = 0.0;
    bool flag = true;

    for (int i = 0; i < 6; i++) xf[i] = xi[i];
    double newz_tr = z_tr;
    while ((l < llimit)&&((newz_tr < zlimit)^(sign))) {
        if ((error > fErrHiLimit)&&(dt > dtlimit)) {
            dt /= 2.0;
            flag = false;
        } else {
            for (int i = 0; i < 6; i++) xi[i] = xf[i];
            if (error < fErrLoLimit) dt *= 2.0;
            l += fVelocity*dt;
            flag = true;
        }
        ComputeRHS(xi, dxdt);
        NystromRK4(xi, dxdt, dt, xf, err);
        error = DistChord();
        newz_tr = xf[0] * sinang + xf[2] * cosang;
    }

    if (!flag) l += fVelocity * dt;
    if (l < llimit) {
        while (dt > dtlimit) {
            if ((newz_tr > zlimit)^(sign)) {
                l = l - fVelocity*dt;
                dt /= 2.0;
            } else {
                for (int i = 0; i < 6; i++) xi[i] = xf[i];
            }
            l += fVelocity*dt;
            ComputeRHS(xi, dxdt);
            NystromRK4(xi, dxdt, dt, xf, err);
            newz_tr = xf[0] * sinang + xf[2] * cosang;
        }
    }

    if (sign) {
        xf[3] *= -1.0;
        xf[4] *= -1.0;
        xf[5] *= -1.0;
        fQ = fQSave;
    }

    theta = acos(xf[5] / fVelocity);
    phi = atan2(xf[4], xf[3]);

    HCS2TCS(theta, phi, angle, xout[1], xout[3]);
    double temp;
    HCS2TCS(xf[0], xf[1], xf[2], angle, xout[0], xout[2], temp);
    xout[4] = x[4];

    return l;
}

double G2PDrift::DriftTCSL(const double* x, double p, double z_tr, double angle, double zlimit, double llimit, double rlimit, double* xout)
{
    // J. Liu: added for longitudinal cylinder
    // zlimit is the limit in lab coordinate

    double xi[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double xf[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    TCS2HCS(x[0], x[2], z_tr, angle, xi[0], xi[1], xi[2]);
    double theta, phi;
    TCS2HCS(x[1], x[3], angle, theta, phi);

    double pp = (1 + x[4]) * p;
    xi[3] = pp * sin(theta) * cos(phi);
    xi[4] = pp * sin(theta) * sin(phi);
    xi[5] = pp * cos(theta);
    double M = sqrt(fM0 * fM0 + pp * pp);
    xi[3] *= c / M;
    xi[4] *= c / M;
    xi[5] *= c / M;
    fVelocity2 = xi[3] * xi[3] + xi[4] * xi[4] + xi[5] * xi[5];
    fVelocity = sqrt(fVelocity2);
    fGamma = 1 / sqrt(1 - (fVelocity2 / c / c));

    double dxdt[6], err[6];
    double error = (fErrHiLimit + fErrLoLimit) / 2.0;
    double dt = fStep / fVelocity;
    double dtlimit = fStepLimit / fVelocity;
    double l = 0.0;
    bool flag = true;

    double r_lab = sqrt(xi[0] * xi[0] + xi[1] * xi[1]);
    for (int i = 0; i < 6; i++) xf[i] = xi[i];
    double newz_lab = xi[2];
    double newr_lab = r_lab;
    while ((l < llimit)&&((newz_lab < zlimit))&&(newr_lab < rlimit)) {
        if ((error > fErrHiLimit)&&(dt > dtlimit)) {
            dt /= 2.0;
            flag = false;
        } else {
            for (int i = 0; i < 6; i++) xi[i] = xf[i];
            if (error < fErrLoLimit) dt *= 2.0;
            l += fVelocity*dt;
            flag = true;
        }
        ComputeRHS(xi, dxdt);
        NystromRK4(xi, dxdt, dt, xf, err);
        error = DistChord();
        newz_lab = xf[2];
        newr_lab = sqrt(xf[0] * xf[0] + xf[1] * xf[1]);
    }

    // updated before subtraction;
    if (newz_lab > zlimit) xout[6] = 0;
    if (newr_lab > rlimit) xout[6] = 1;
    if (!flag) l += fVelocity * dt;
    if (l < llimit) {
        while (dt > dtlimit) {
            if ((newz_lab > zlimit) || (newr_lab > rlimit)) {
                l = l - fVelocity*dt;
                dt /= 2.0;
            } else {
                for (int i = 0; i < 6; i++) xi[i] = xf[i];
            }
            l += fVelocity*dt;
            ComputeRHS(xi, dxdt);
            NystromRK4(xi, dxdt, dt, xf, err);
            newz_lab = xf[2];
            newr_lab = sqrt(xf[0] * xf[0] + xf[1] * xf[1]);
        }
    }

    theta = acos(xf[5] / fVelocity);
    phi = atan2(xf[4], xf[3]);

    HCS2TCS(theta, phi, angle, xout[1], xout[3]);
    HCS2TCS(xf[0], xf[1], xf[2], angle, xout[0], xout[2], xout[5]);
    xout[4] = x[4];

    return l;
}

double G2PDrift::DriftTCSV(const double* x, double p, double z_tr, double angle, double zlimit, double llimit, double rlimit, double* xout)
{
    // J. Liu: added for vertical cylinder
    // zlimit is the limit in lab coordinate

    double xi[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double xf[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    double sinang = sin(angle);
    double cosang = cos(angle);

    TCS2HCS(x[0], x[2], z_tr, angle, xi[0], xi[1], xi[2]);
    double theta, phi;
    TCS2HCS(x[1], x[3], angle, theta, phi);

    double pp = (1 + x[4]) * p;
    xi[3] = pp * sin(theta) * cos(phi);
    xi[4] = pp * sin(theta) * sin(phi);
    xi[5] = pp * cos(theta);
    double M = sqrt(fM0 * fM0 + pp * pp);
    xi[3] *= c / M;
    xi[4] *= c / M;
    xi[5] *= c / M;
    fVelocity2 = xi[3] * xi[3] + xi[4] * xi[4] + xi[5] * xi[5];
    fVelocity = sqrt(fVelocity2);
    fGamma = 1 / sqrt(1 - (fVelocity2 / c / c));

    double dxdt[6], err[6];
    double error = (fErrHiLimit + fErrLoLimit) / 2.0;
    double dt = fStep / fVelocity;
    double dtlimit = fStepLimit / fVelocity;
    double l = 0.0;
    bool flag = true;

    double r_tr = sqrt(xi[0] * xi[0] + xi[2] * xi[2]);
    for (int i = 0; i < 6; i++) xf[i] = xi[i];
    double newz_tr = z_tr;
    double newr_tr = r_tr;
    while ((l < llimit)&&((newz_tr < zlimit))&&(newr_tr < rlimit)) {
        if ((error > fErrHiLimit)&&(dt > dtlimit)) {
            dt /= 2.0;
            flag = false;
        } else {
            for (int i = 0; i < 6; i++) xi[i] = xf[i];
            if (error < fErrLoLimit) dt *= 2.0;
            l += fVelocity*dt;
            flag = true;
        }
        ComputeRHS(xi, dxdt);
        NystromRK4(xi, dxdt, dt, xf, err);
        error = DistChord();
        newz_tr = xf[0] * sinang + xf[2] * cosang;
        newr_tr = sqrt(xf[0] * xf[0] + xf[2] * xf[2]);
    }

    if (!flag) l += fVelocity * dt;
    if (l < llimit) {
        while (dt > dtlimit) {
            if ((newz_tr > zlimit) || (newr_tr > rlimit)) {
                l = l - fVelocity*dt;
                dt /= 2.0;
            } else {
                for (int i = 0; i < 6; i++) xi[i] = xf[i];
            }
            l += fVelocity*dt;
            ComputeRHS(xi, dxdt);
            NystromRK4(xi, dxdt, dt, xf, err);
            newz_tr = xf[0] * sinang + xf[2] * cosang;
            newr_tr = sqrt(xf[0] * xf[0] + xf[2] * xf[2]);
        }
    }

    theta = acos(xf[5] / fVelocity);
    phi = atan2(xf[4], xf[3]);

    HCS2TCS(theta, phi, angle, xout[1], xout[3]);
    HCS2TCS(xf[0], xf[1], xf[2], angle, xout[0], xout[2], xout[5]);
    xout[4] = x[4];

    return l;
}

double G2PDrift::DriftTCSNF(const double* x, double p, double z_tr, double angle, double zlimit, double llimit, double* xout)
{
    xout[0] = x[0]+(zlimit - z_tr) * tan(x[1]);
    xout[1] = x[1];
    xout[2] = x[2]+(zlimit - z_tr) * tan(x[3]);
    xout[3] = x[3];
    xout[4] = x[4];

    double dx = xout[0] - x[0];
    double dy = xout[2] - x[2];
    double dz = zlimit - z_tr;

    return sqrt(dx * dx + dy * dy + dz * dz);
}

void G2PDrift::NystromRK4(const double* x, const double* dxdt, double step, double* xo, double* err)
{
    // Using Nystrom-Runge-Kutta method to solve the motion equation of the 
    // electron in static magnetic field numerically
    // Using 4th order Nystrom-Runge-Kutta method

    static const char* const here = "NystromRK4()";

    double S = step;
    double S5 = step * 0.5;
    double S4 = step * 0.25;
    double S6 = step*kOneSixth;

    double R[3] = {x[0], x[1], x[2]};
    double A[3] = {dxdt[0], dxdt[1], dxdt[2]};

    fIPoint[0] = R[0];
    fIPoint[1] = R[1];
    fIPoint[2] = R[2];

    // Point 1
    double K1[3] = {dxdt[3], dxdt[4], dxdt[5]};

    // Point 2
    double p[3] = {R[0] + S5 * (A[0] + S4 * K1[0]),
        R[1] + S5 * (A[1] + S4 * K1[1]),
        R[2] + S5 * (A[2] + S4 * K1[2])};

    pField->GetField(p, fField);

    double A2[3] = {A[0] + S5 * K1[0], A[1] + S5 * K1[1], A[2] + S5 * K1[2]};
    double K2[3] = {(A2[1] * fField[2] - A2[2] * fField[1]) * fCof,
        (A2[2] * fField[0] - A2[0] * fField[2]) * fCof,
        (A2[0] * fField[1] - A2[1] * fField[0]) * fCof};

    fMPoint[0] = p[0];
    fMPoint[1] = p[1];
    fMPoint[2] = p[2];

    // Point 3
    double A3[3] = {A[0] + S5 * K2[0], A[1] + S5 * K2[1], A[2] + S5 * K2[2]};
    double K3[3] = {(A3[1] * fField[2] - A3[2] * fField[1]) * fCof,
        (A3[2] * fField[0] - A3[0] * fField[2]) * fCof,
        (A3[0] * fField[1] - A3[1] * fField[0]) * fCof};

    // Point 4
    p[0] = R[0] + S * (A[0] + S5 * K3[0]);
    p[1] = R[1] + S * (A[1] + S5 * K3[1]);
    p[2] = R[2] + S * (A[2] + S5 * K3[2]);

    pField->GetField(p, fField);

    double A4[3] = {A[0] + S * K3[0], A[1] + S * K3[1], A[2] + S * K3[2]};
    double K4[3] = {(A4[1] * fField[2] - A4[2] * fField[1]) * fCof,
        (A4[2] * fField[0] - A4[0] * fField[2]) * fCof,
        (A4[0] * fField[1] - A4[1] * fField[0]) * fCof};

    // New position
    xo[0] = R[0] + S * (A[0] + S6 * (K1[0] + K2[0] + K3[0]));
    xo[1] = R[1] + S * (A[1] + S6 * (K1[1] + K2[1] + K3[1]));
    xo[2] = R[2] + S * (A[2] + S6 * (K1[2] + K2[2] + K3[2]));

    fEPoint[0] = xo[0];
    fEPoint[1] = xo[1];
    fEPoint[2] = xo[2];

    // New direction
    xo[3] = A[0] + S6 * (K1[0] + K4[0] + 2. * (K2[0] + K3[0]));
    xo[4] = A[1] + S6 * (K1[1] + K4[1] + 2. * (K2[1] + K3[1]));
    xo[5] = A[2] + S6 * (K1[2] + K4[2] + 2. * (K2[2] + K3[2]));

    // Errors
    err[3] = S * fabs(K1[0] - K2[0] - K3[0] + K4[0]);
    err[4] = S * fabs(K1[1] - K2[1] - K3[1] + K4[1]);
    err[5] = S * fabs(K1[2] - K2[2] - K3[2] + K4[2]);
    err[0] = S * err[3];
    err[1] = S * err[4];
    err[2] = S * err[5];

    double normF = sqrt(fVelocity2 / (xo[3] * xo[3] + xo[4] * xo[4] + xo[5] * xo[5]));
    xo[3] *= normF;
    xo[4] *= normF;
    xo[5] *= normF;

    if (fDebug > 4) {
        Info(here, " err: %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e", err[0] / fabs(xo[0] - x[0]), err[1] / fabs(xo[1] - x[1]), err[2] / fabs(xo[2] - x[2]), err[3] / fabs(xo[3] - x[3]), err[4] / fabs(xo[4] - x[4]), err[5] / fabs(xo[5] - x[5]));
    }
}

double G2PDrift::DistChord()
{
    // Calculate the deviation of the trajectory

    double ax = fEPoint[0] - fIPoint[0];
    double ay = fEPoint[1] - fIPoint[1];
    double az = fEPoint[2] - fIPoint[2];
    double dx = fMPoint[0] - fIPoint[0];
    double dy = fMPoint[1] - fIPoint[1];
    double dz = fMPoint[2] - fIPoint[2];
    double d2 = ax * ax + ay * ay + az*az;

    if (d2 != 0.) {
        double s = (ax * dx + ay * dy + az * dz) / d2;
        dx -= (s * ax);
        dy -= (s * ay);
        dz -= (s * az);
    }

    return sqrt(dx * dx + dy * dy + dz * dz);
}

void G2PDrift::ComputeRHS(const double* x, double* dxdt)
{
    // Calculate the right hand side of the motion equation

    static const char* const here = "ComputeHRS()";

    pField->GetField(x, fField);

    double M = fM0 * fGamma*kGEV;
    dxdt[0] = x[3];
    dxdt[1] = x[4];
    dxdt[2] = x[5];
    fCof = fQ / M;
    dxdt[3] = (x[4] * fField[2] - fField[1] * x[5]) * fCof;
    dxdt[4] = (x[5] * fField[0] - fField[2] * x[3]) * fCof;
    dxdt[5] = (x[3] * fField[1] - fField[0] * x[4]) * fCof;

    if (fDebug > 4) {
        Info(here, "   x: %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e", x[0], x[1], x[2], x[3], x[4], x[5]);
        Info(here, "dxdt: %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e", dxdt[0], dxdt[1], dxdt[2], dxdt[3], dxdt[4], dxdt[5]);
    }
}

int G2PDrift::Configure(EMode mode)
{
    if (mode == kREAD || mode == kTWOWAY) {
        if (fIsInit) return 0;
        else fIsInit = true;
    }

    ConfDef confs[] = {
        {"run.particle.mass", "Particle Mass", kDOUBLE, &fM0},
        {"run.particle.charge", "Particle Charge", kDOUBLE, &fQ},
        {"step", "Particle Charge", kDOUBLE, &fStep},
        {"step.limit", "Particle Charge", kDOUBLE, &fStepLimit},
        {"error.low", "Lower limit", kDOUBLE, &fErrLoLimit},
        {"error.high", "Upper Limit", kDOUBLE, &fErrHiLimit},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

void G2PDrift::MakePrefix()
{
    const char* base = "drift";

    G2PAppBase::MakePrefix(base);
}

ClassImp(G2PDrift)
