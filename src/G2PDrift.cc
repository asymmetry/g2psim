#include <cmath>

#include "TROOT.h"
#include "TObject.h"

#include "G2PTargetField.hh"

#include "G2PDrift.hh"

using namespace std;

ClassImp(G2PDrift);

const double c = 2.99796458e+08;
const double kGEV = 1.78261428e-27;
const double e = 1.60217657e-19;
const double kOneSixth = 1./6.;

G2PDrift::G2PDrift()
    :fM0(0.000510998928), fQ(-1.0*e), fStep(1.0e-3), pField(NULL)
{
    // Nothing to do
}

G2PDrift::G2PDrift(G2PTargetField *field)
    :fM0(0.000510998928), fQ(-1.0*e), fStep(1.0e-3)
{
    pField = field;
}

G2PDrift::~G2PDrift()
{
    // Nothing to do
}

void G2PDrift::Drift(double* x, double* p, double zlimit, double llimit)
{
    double xi[6], xf[6];

    double M = sqrt(fM0*fM0+p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    double v[3];
    v[0] = p[0]*(c/M);
    v[1] = p[1]*(c/M);
    v[2] = p[2]*(c/M);
    double vv = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

    xi[0] = x[0];
    xi[1] = x[1];
    xi[2] = x[2];
    xi[3] = v[0];
    xi[4] = v[1];
    xi[5] = v[2];

    double dt = fStep/vv;
    double l = 0;

    double dxdt[6], err[6];
    bool sign = (x[2]>zlimit)?true:false;
    while ((l<llimit)&&((xf[2]>zlimit)^(sign))) {
        ComputeRHS(xi, dxdt);
        NystromRK4(xi, dxdt, dt, xf, err);
        for (int i = 0; i<6; i++) xi[i] = xf[i];
        l+=vv*dt;
    }

    x[0] = xf[0];
    x[1] = xf[1];
    x[2] = xf[2];
    p[0] = xf[3]*M/c;
    p[1] = xf[4]*M/c;
    p[2] = xf[5]*M/c;
}

void G2PDrift::NystromRK4(const double* x, const double* dxdt, double step, double* xo, double* err)
{
    double S = step;
    double S5 = step*0.5;
    double S4 = step*0.25;
    double S6 = step*kOneSixth;

    double R[3] = { x[0], x[1], x[2] };
    double A[3] = { dxdt[0], dxdt[1], dxdt[2] };

    fInPoint[0] = R[0];
    fInPoint[1] = R[1];
    fInPoint[2] = R[2];

    // Point 1
    double K1[3] = { x[3], x[4], x[5] };

    // Point 2
    double p[3] = { R[0]+S5*(A[0]+S4*K1[0]),
                    R[1]+S5*(A[1]+S4*K1[1]),
                    R[2]+S5*(A[2]+S4*K1[2])};
    pField->GetField(p, fField);

    double A2[3] = { A[0]+S5*K1[0], A[1]+S5*K1[1], A[2]+S5*K1[2] };
    double K2[3] = { (A2[1]*fField[2]-A2[2]*fField[1])*fCof,
                     (A2[2]*fField[0]-A2[0]*fField[2])*fCof,
                     (A2[0]*fField[1]-A2[1]*fField[0])*fCof};

    fMidPoint[0] = p[0];
    fMidPoint[1] = p[1];
    fMidPoint[2] = p[2];

    // Point 3
    double A3[3] = { A[0]+S5*K2[0], A[1]+S5*K2[1], A[2]+S5*K2[2] };
    double K3[3] = { (A3[1]*fField[2]-A3[2]*fField[1])*fCof,
                     (A3[2]*fField[0]-A3[0]*fField[2])*fCof,
                     (A3[0]*fField[1]-A3[1]*fField[0])*fCof};

    // Point 4
    p[0] = R[0]+S*(A[0]+S5*K3[0]);
    p[1] = R[1]+S*(A[1]+S5*K3[1]);
    p[2] = R[2]+S*(A[2]+S5*K3[2]);

    pField->GetField(p, fField);

    double A4[3] = { A[0]+S*K3[0], A[1]+S*K3[1], A[2]+S*K3[2] };
    double K4[3] = { (A4[1]*fField[2]-A4[2]*fField[1])*fCof,
                     (A4[2]*fField[0]-A4[0]*fField[2])*fCof,
                     (A4[0]*fField[1]-A4[1]*fField[0])*fCof};

    // New position
    xo[0] = R[0]+S*(A[0]+S6*(K1[0]+K2[0]+K3[0]));
    xo[1] = R[1]+S*(A[1]+S6*(K1[1]+K2[1]+K3[1]));
    xo[2] = R[2]+S*(A[2]+S6*(K1[2]+K2[2]+K3[2]));

    fFinPoint[0] = xo[0];
    fFinPoint[1] = xo[1];
    fFinPoint[2] = xo[2];

    // New direction
    xo[3] = A[0]+S6*(K1[0]+K4[0]+2.*(K2[0]+K3[0]));
    xo[4] = A[1]+S6*(K1[1]+K4[1]+2.*(K2[1]+K3[1]));
    xo[5] = A[2]+S6*(K1[2]+K4[2]+2.*(K2[2]+K3[2]));

    // Errors
    err[3] = S*fabs(K1[0]-K2[0]-K3[0]+K4[0]);
    err[4] = S*fabs(K1[1]-K2[1]-K3[1]+K4[1]);
    err[5] = S*fabs(K1[2]-K2[2]-K3[2]+K4[2]);
    err[0] = S*err[3];
    err[1] = S*err[4];
    err[2] = S*err[5];

    double normF = sqrt(fVelocity2/(xo[3]*xo[3]+xo[4]*xo[4]+xo[5]*xo[5]));
    xo[3]*=normF;
    xo[4]*=normF;
    xo[5]*=normF;
}

double G2PDrift::DistChord() const
{
    double ax = fFinPoint[0]-fInPoint[0];
    double ay = fFinPoint[1]-fInPoint[1];
    double az = fFinPoint[2]-fInPoint[2];
    double dx = fMidPoint[0]-fInPoint[0];
    double dy = fMidPoint[1]-fInPoint[1];
    double dz = fMidPoint[2]-fInPoint[2];
    double d2 = ax*ax+ay*ay+az*az;

    if (d2!=0.) {
        double s = (ax*dx+ay*dy+az*dz)/d2;
        dx-=(s*ax);
        dy-=(s*ay);
        dz-=(s*az);
    }

    return sqrt(dx*dx+dy*dy+dz*dz);
}

void G2PDrift::ComputeRHS(const double* x, double* dxdt)
{
    pField->GetField(x, fField);
    fVelocity2 = x[3]*x[3]+x[4]*x[4]+x[5]*x[5];
    double gamma = 1.0/sqrt(1.0-(fVelocity2)/c/c);
    double m = fM0*gamma;
    dxdt[0] = x[3];
    dxdt[1] = x[4];
    dxdt[2] = x[5];
    fCof = fQ/m;
    dxdt[3] = (x[4]*fField[2]-fField[1]*x[5])*fCof;
    dxdt[4] = (x[5]*fField[0]-fField[2]*x[3])*fCof;
    dxdt[5] = (x[3]*fField[1]-fField[0]*x[4])*fCof;
}
