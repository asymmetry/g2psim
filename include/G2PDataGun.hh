// -*- C++ -*-

/* class G2PDataGun
 * This file defines a class G2PDataGun.
 * It reads kinematics from real data and generate events for simulation.
 * G2PProcBase classes will call Shoot() to get reaction point kinematics.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//

#ifndef G2P_DATAGUN_H
#define G2P_DATAGUN_H

#include <vector>

#include "G2PGun.hh"
#include "G2PSieve.hh"

using namespace std;

class G2PDataGun : public G2PGun, public G2PSieve {
public:
    G2PDataGun(const char* filename);
    ~G2PDataGun();

    typedef int (G2PDataGun::*pfGun_)(double*, double*, double*);

    void SetOpticsData(bool b) {
        bIsOptics = b;
    }

    int Begin();

    int Shoot(double* V51, double* V52, double* V53 = NULL) {
        return (this->*pfGun)(V51, V52, V53);
    }

    bool UseData() {
        return true;
    }

protected:
    G2PDataGun(); // Only for ROOT I/O

    int ShootNormal(double* V5bpm_lab, double* V5react_tr, double* V5fp_tr);
    int ShootOptics(double* V5bpm_lab, double* V5react_tr, double* V5fp_tr);

    int LoadData();

    bool bIsOptics;

    const char* pDataFileName;

    double fTargetMass;
    double fEnergyLoss;

    typedef struct {
        int ind;
        double xb, tb, yb, pb, zb, xf, tf, yf, pf;
    } sData;

    vector<sData> fData;

private:
    pfGun_ pfGun;

    ClassDef(G2PDataGun, 1)
};

#endif
