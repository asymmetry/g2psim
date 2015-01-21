// -*- C++ -*-

/* class G2PMaterial
 * Calculate energy loss and multi-scattering for defined material.
 * Most of the algorithm is from SAMC package.
 * Thanks to the author of SAMC.
 */

// History:
//   Sep 2013, C. Gu, First public version.
//   Oct 2013, J. Liu, Check formulas
//   Nov 2014, J. Liu, Add ionization fluctuation, internal and external bremsstrahlung.
//

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TClass.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PMaterial.hh"
#include "G2PRand.hh"
#include "G2PVarDef.hh"

using namespace std;

static const double kMe = 0.510998918; // MeV

static double b_func(int Z)
{
    // Rev. Mod. Phys. 46(1974)815

    double Lrad, Lradp;

    switch (Z) {
    case 0:
        return 0;

    case 1:
    case 2:
        Lrad = (4.79 - 5.31) * (Z - 1) + 5.31;
        Lradp = (5.621 - 6.144) * (Z - 1) + 6.144;
        break;

    case 3:
        Lrad = (4.74 - 4.79) * (Z - 2) + 4.79;
        Lradp = (5.805 - 5.621) * (Z - 2) + 5.621;
        break;

    case 4:
        Lrad = (4.71 - 4.74) * (Z - 3) + 4.74;
        Lradp = (5.924 - 5.805) * (Z - 3) + 5.805;
        break;

    default:
        Lrad = log(184.15 * pow(Z, -1.0 / 3.0));
        Lradp = log(1194.0 * pow(Z, -2.0 / 3.0));
        break;
    }

    return (4.0 / 3.0) * (1.0 + (1.0 / 9.0) * (Z + 1) / (Lrad * Z + Lradp));
}

G2PMaterial::G2PMaterial()
{
    // Only for ROOT I/O
}

G2PMaterial::G2PMaterial(const char *name, double z, double a, double x0, double density, double ion, double c) : fName(name), fZ(z), fA(a), fMass(0), fDensity(density), fX0(x0), fIon(ion), fC(c)
{
    fB = b_func(fZ);
}

G2PMaterial::~G2PMaterial()
{
    // Nothing to do
}

double G2PMaterial::EnergyLoss(double E, double l)
{
    double EMeV = E * 1000; // MeV
    double lcm = l * 100; // cm

    assert(lcm >= 0);

    double result = 0;
    result += Ionization(EMeV, lcm);
    result += ExterBremsstrahlung(EMeV, lcm);

    if (result > (E - 2 * kMe))
        result = E - 2 * kMe;

    if ((result < 0) || (E < 2 * kMe))
        result = 0;

    if (result > 0)
        return result / 1000.0; // GeV

    return 0;
}

double G2PMaterial::MultiScattering(double E, double l)
{
    // Only for electron

    double EMeV = E * 1000; // MeV
    double thicknessr = (l * 100) * fDensity / fX0; // l: m -> cm

    assert(l >= 0);

    double lPsq = EMeV * EMeV - kMe * kMe;
    double bcp = lPsq / EMeV;
    double ltheta0 = 13.6 / bcp * sqrt(thicknessr) * (1 + 0.038 * log(thicknessr));

    if (thicknessr != 0)
        return pRand->Gaus(0, ltheta0); // rad

    return 0.0;
}

const char *G2PMaterial::GetName()
{
    return fName;
}

double G2PMaterial::GetRadLen()
{
    return fX0;
}

double G2PMaterial::GetDensity()
{
    return fDensity;
}

double G2PMaterial::Ionization(double E, double l)
{
    // Calculate ionization dE/dx by Seltzer-Berger formula

    double fcut = 6;
    double thickness = l * fDensity;
    double lK = 0.307075; // cm^2/g for A=1 g/mol
    double tau = E / kMe - 1;
    double gamma = tau + 1.0;
    double gamma2 = gamma * gamma;
    double bg2 = tau * (tau + 2);
    double beta2 = bg2 / gamma2;

    double eexc = fIon * 1e-6 / kMe;
    // double I = kIONI[fZ]* 1e-6;
    double eexc2 = eexc * eexc;
    double d = fcut / kMe;
    double dedx = log(2.0 * (tau + 2.0) / eexc2) - 1.0 - beta2 + log((tau - d) * d) + tau / (tau - d) + (0.5 * d * d + (2.0 * tau + 1.) * log(1. - d / tau)) / gamma2;

    // density correction
    double x = 0.5 * log10(bg2);

    double delta = 4.6052 * x - fC;
    dedx -= delta;

    // now you can compute the total ionization loss
    dedx *= lK / 2 * fZ / fA * thickness / beta2;

    if (dedx < 0.0)
        dedx = 0.0;

    double result = FluctIonization(E, dedx);

    if (result > (E - kMe))
        result = E - kMe;

    if (result < 0)
        result = 0;

    return result; // MeV
}

double G2PMaterial::FluctIonization(double E, double meanloss)
{
    // Calculate actual ionization energy loss from mean loss
    // Urban model

    double fcut = 6;
    double tau = E / kMe - 1;
    double gamma = tau + 1.0;
    double gamma2 = gamma * gamma;
    double bg2 = tau * (tau + 2);
    double beta2 = bg2 / gamma2;

    double eV = 1e-6;

    int nmax = 16;
    double fw = 4;
    double loss = 0;
    double e1 = 0;
    double e2 = 0;

    double f2fluct = 2. / fZ;

    if (fZ < 3)
        f2fluct = 0.;

    double f1fluct = 1. - f2fluct;
    double e2fluct = 10. * fZ * fZ * eV;
    double loge2fluct = log(e2fluct);
    double ipotfluct = fIon * eV;
    double logipotfluct = log(fIon * eV);
    double loge1fluct = (logipotfluct - f2fluct * loge2fluct) / f1fluct;
    double e1fluct = exp(loge1fluct);
    double e0 = 10. * eV;
    double rate = 0.55;
    double tmax = fcut;
    double esmall = 0.5 * sqrt(e0 * ipotfluct);

    if (tmax <= e0)
        return meanloss;

    double losstot = 0.;
    int nstep = 1;

    if (meanloss < 25. * ipotfluct) {
        if (pRand->Uniform() < 0.04 * meanloss / ipotfluct)
            nstep = 1;
        else {
            nstep = 2;
            meanloss *= 0.5;
        }
    }

    for (int istep = 0; istep < nstep; istep++) {
        loss = 0.;
        double a1 = 0., a2 = 0., a3 = 0.;

        if (tmax > ipotfluct) {
            double w2 = log(2. * kMe * beta2 * gamma2) - beta2;

            if (w2 > logipotfluct) {
                double C = meanloss * (1. - rate) / (w2 - logipotfluct);
                a1 = C * f1fluct * (w2 - loge1fluct) / e1fluct;

                if (w2 > loge2fluct)
                    a2 = C * f2fluct * (w2 - loge2fluct) / e2fluct;

                if (a1 < nmax) {
                    // small energy loss
                    double sa1 = sqrt(a1);

                    if (pRand->Uniform() < exp(-sa1)) {
                        e1 = esmall;
                        a1 = meanloss * (1. - rate) / e1;
                        a2 = 0.;
                        e2 = e2fluct;
                    } else {
                        a1 = sa1;
                        e1 = sa1 * e1fluct;
                        e2 = e2fluct;
                    }
                } else {
                    // not small energy loss
                    // correction to get better FWHM value
                    a1 /= fw;
                    e1 = fw * e1fluct;
                    e2 = e2fluct;
                }
            }
        }

        double w1 = tmax / e0;

        if (tmax > e0)
            a3 = rate * meanloss * (tmax - e0) / (e0 * tmax * log(w1));

        // "nearly" Gaussian fluctuation if (a1 > nmax && a2 > nmax && a3 > nmax)
        double emean = 0.;
        double sig2e = 0., sige = 0.;
        double p1 = 0., p2 = 0., p3 = 0.;

        // excitation of type 1
        if (a1 > nmax) {
            emean += a1 * e1;
            sig2e += a1 * e1 * e1;
        } else if (a1 > 0.) {
            p1 = double(pRand->Poisson(a1));
            loss += p1 * e1;

            if (p1 > 0.)
                loss += (1. - 2. * pRand->Uniform()) * e1;
        }

        // excitation of type 2
        if (a2 > nmax) {
            emean += a2 * e2;
            sig2e += a2 * e2 * e2;
        } else if (a2 > 0.) {
            p2 = double(pRand->Poisson(a2));
            loss += p2 * e2;

            if (p2 > 0.)
                loss += (1. - 2. * pRand->Uniform()) * e2;
        }

        if (emean > 0.) {
            sige = sqrt(sig2e);
            loss += max(0., pRand->Gaus(emean, sige));
        }

        // ionization
        double lossc = 0.;

        if (a3 > 0.) {
            emean = 0.;
            sig2e = 0.;
            sige = 0.;
            p3 = a3;
            double alfa = 1.;

            if (a3 > nmax) {
                alfa = w1 * (nmax + a3) / (w1 * nmax + a3);
                double alfa1 = alfa * log(alfa) / (alfa - 1.);
                double namean = a3 * w1 * (alfa - 1.) / ((w1 - 1.) * alfa);
                emean += namean * e0 * alfa1;
                sig2e += e0 * e0 * namean * (alfa - alfa1 * alfa1);
                p3 = a3 - namean;
            }

            double w2 = alfa * e0;
            double w = (tmax - w2) / tmax;
            int nb = pRand->Poisson(p3);

            if (nb > 0) {
                for (int k = 0; k < nb; k++)
                    lossc += w2 / (1. - w * pRand->Uniform());
            }
        }

        if (emean > 0.) {
            sige = sqrt(sig2e);
            lossc += max(0., pRand->Gaus(emean, sige));
        }

        loss += lossc;
        losstot += loss;
    }

    return losstot;
}

double G2PMaterial::ExterBremsstrahlung(double E, double l)
{
    // Taken from gener_cone MC, modified by R.M
    // Modified by Jie Liu, Nov 25 2014

    // L.Van Hoorebeke, University of Gent, e-mail: Luc.VanHoorebeke@UGent.be, luc@inwfsun1.rug.ac.be
    // This function generates external bremsstrahlung in the target
    // the distribution used contains multiple emission effects
    // k: electron momentum
    // fracrl = fraction of radiation lengths distance, traveled through material.

    double cut, Ekin, bt, prob, prob_sample;
    double sample, xtry, env, value, ref;
    double Euler = 0.5772157;

    // Initialization of lower limit of bremsstrahlung (1 keV)
    cut = 0.001;
    Ekin = E - kMe;

    bt = l * fDensity / fX0 * fB;

    // Calculation of probability to have bremsstrahlung effect above 1 keV
    prob = 1. - pow(cut / Ekin, bt) - bt / (bt + 1.) * (1. - pow(cut / Ekin, bt + 1.)) + 0.75 * bt / (2. + bt) * (1. - pow(cut / Ekin, bt + 2.));
    prob = prob / (1. - bt * Euler + bt * bt / 2.*(Euler * Euler + 3.1415926 * 3.1415926 / 6.)); /* Gamma function */

    prob_sample = pRand->Uniform();

    if (prob_sample > prob)
        return 0.;

    // Bremsstrahlung has taken place! Generate photon energy with sample and reject
    // using 1/x as envelope

    do {
        sample = pRand->Uniform();
        xtry = cut * pow(Ekin / cut, sample);

        env = 1. / xtry;
        value = 1. / xtry * (1. - xtry / Ekin + 0.75 * pow(xtry / Ekin, 2)) * pow(xtry / Ekin, bt);

        sample = pRand->Uniform();
        ref = value / env;
    } while (sample > ref);

    return xtry;
}

void G2PMaterial::MakePrefix()
{
    const char *base = "material";

    G2PAppBase::MakePrefix(Form("%s.%s", base, fName));
}

ClassImp(G2PMaterial)
