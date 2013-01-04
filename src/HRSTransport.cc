#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "HRSTransport.hh"
#include "Rand.hh"

const double deg = acos(0.0)/90.0;

//#define DEBUG_HRS_FORWARD
//#define DEBUG_HRS_BACKWARD

// Transport particles through HRS using SNAKE model
// Use iSetting to identify which SNAKE model to be used
// 10: 484816 no shim, 5.65 deg, Wrong Bx, by JJL (g2p test run)
// 11: 484816 with shim, 5.65 deg, 3 cm raster, by JJL 
// 12: 403216 with shim, 5.65 deg, SNAKE Model not ready yet 
// 13: 400016 with shim, 5.65 deg, SNAKE Model not ready yet 
// 19: 484816 with shim, 5.65 deg, Wrong Bx, 2 cm raster, by Min
// May add more HRS packages later

bool SNAKEForward(int pIsLeftArm, int iSetting, const double* pV5_tg, double* pV5_fp)
{
    // Definition of variables
    // pV5_tg = {x_tg, theta_tg, y_tg, phi_tg, delta@tg};
    // pV5_fp = {x_fg, theta_fp, y_fp, phi_fp, delta@tg};
    // delta does not change
    double pV5[5];
    
    pV5[0]=pV5_tg[0];
    pV5[1]=tan(pV5_tg[1]);
    pV5[2]=pV5_tg[2];
    pV5[3]=tan(pV5_tg[3]);
    pV5[4]=pV5_tg[4];

#ifdef DEBUG_HRS_FORWARD
	printf("FORWARDIN:  %8.4f %8.4f %8.4f %8.4f %8.4f", pV5[0], pV5[1], pV5[2], pV5[3], pV5[4]);
#endif

    bool bGoodParticle=false;

    switch (pIsLeftArm*100+iSetting) {
    case 10: // No shim, Wrong Bx, by JLL
        bGoodParticle=TransportRightHRS_g2_(pV5);
        break;
    case 11: // 484816 with shim, 3cm raster, by JLL
        bGoodParticle=TransportRightHRS_Shim_484816(pV5); 
        break;
    case 12: // 403216 with shim
        bGoodParticle=TransportRightHRS_Shim_403216(pV5);
        break;
    case 13: // 400016 with shim
        bGoodParticle=TransportRightHRS_Shim_400016(pV5);
        break;
    case 19: // 484816 with shim, Wrong Bx, 2cm raster, by Min
        bGoodParticle=TransportRightHRS_Shim_484816_WrongBx(pV5);
        break;

    case 110: // No shim, Wrong Bx, by JLL
        bGoodParticle=TransportLeftHRS_g2_(pV5);
        break;
    case 111: // 484816 with shim, 3cm raster, by JLL
        bGoodParticle=TransportLeftHRS_Shim_484816(pV5); 
        break;
    case 112: // 403216 with shim
        bGoodParticle=TransportLeftHRS_Shim_403216(pV5);
        break;
    case 113: // 400016 with shim
        bGoodParticle=TransportLeftHRS_Shim_400016(pV5);
        break;
    case 119: // 484816 with shim, Wrong Bx, 2cm raster, by Min
        bGoodParticle=TransportLeftHRS_Shim_484816_WrongBx(pV5);
        break;

    default:
        printf("iSetting = %d, which is not valid.\n", iSetting);
        exit(0);
    }

    pV5_fp[0]=pV5[0];
    pV5_fp[1]=atan(pV5[1]);
    pV5_fp[2]=pV5[2];
    pV5_fp[3]=atan(pV5[3]);
    pV5_fp[4]=pV5[4];
    
    return bGoodParticle;
}

bool SNAKEBackward(int pIsLeftArm, int iSetting, const double* pV5_fp, double* pV5_tg)
{
    // Definition of variables
    // pV5_fp = {x_fg, theta_fp, y_fp, phi_fp, x_tg};
    // pV5_tg = {x_tg, theta_tg, y_tg, phi_tg, delta@tg};
    // x_tg does not change
    
    double pV5[5];
    
    pV5[0]=pV5_fp[0];
    pV5[1]=tan(pV5_fp[1]);
    pV5[2]=pV5_fp[2];
    pV5[3]=tan(pV5_fp[3]);
    pV5[4]=pV5_fp[4];

#ifdef DEBUG_HRS_BACKWARD
	printf("IN:  %8.4f %8.4f %8.4f %8.4f %8.4f", pV5[0], pV5[1], pV5[2], pV5[3], pV5[4]);
#endif

    switch (pIsLeftArm*100+iSetting) {
    case 10: // No shim, Wrong Bx, by JLL
        ReconstructRightHRS_g2_(pV5);
        break;
    case 11: // 484816 with shim, 3cm raster, by JLL
        ReconstructRightHRS_Shim_484816(pV5); 
        break;
    case 12: // 403216 with shim
        ReconstructRightHRS_Shim_403216(pV5);
        break;
    case 13: // 400016 with shim
        ReconstructRightHRS_Shim_400016(pV5);
        break;
    case 19: // 484816 with shim, Wrong Bx, 2cm raster, by Min
        ReconstructRightHRS_Shim_484816_WrongBx(pV5);
        break;

    case 110: // No shim, Wrong Bx, by JLL
        ReconstructLeftHRS_g2_(pV5);
        break;
    case 111: // 484816 with shim, 3cm raster, by JLL
        ReconstructLeftHRS_Shim_484816(pV5); 
        break;
    case 112: // 403216 with shim
        ReconstructLeftHRS_Shim_403216(pV5);
        break;
    case 113: // 400016 with shim
        ReconstructLeftHRS_Shim_400016(pV5);
        break;
    case 119: // 484816 with shim, Wrong Bx, 2cm raster, by Min
        ReconstructLeftHRS_Shim_484816_WrongBx(pV5);
        break;

    default:
        printf("iSetting = %d, which is not valid.\n", iSetting);
        exit(0);
    }

    pV5_tg[0]=pV5[0];
    pV5_tg[1]=atan(pV5[1]);
    pV5_tg[2]=pV5[2];
    pV5_tg[3]=atan(pV5[3]);
    pV5_tg[4]=pV5[4];
    
    return true;
}

void DeltaCorrection(double &pDelta, double &pP_rec){;}
void XtgCorrection(double &pX,double pP_rec){;}
void ThetatgCorrection(double &pTheta,double pP_rec){;}
void YtgCorrection(double &pY,double pP_rec){;}
void PhitgCorrection(double &pPhi,double pP_rec){;}
