#include "HRSTransport.hh"
#include "GlobalDebuger.hh"

#include <iostream>
#include <iomanip>
#include <math.h>
#include "stdio.h"
#include "stdlib.h"

using namespace std;

#define DEBUG_HRS_TRANSPORT 2
#define DEBUG_HRS_RECONSTRUCTION 1
//input of forward:
//pVector_tg[0] = x_tg;
//pVector_tg[1] = theta_tg;
//pVector_tg[2] = y_tg;
//pVector_tg[3] = phi_tg;
//pVector_tg[4] = delta;
//output
//pVector_fp[0] = x_fp;
//pVector_fp[1] = theta_fp;
//pVector_fp[2] = y_fp;
//pVector_fp[3] = phi_fp;
//pVector_fp[4] = delta;		// = delta @ tg; no change
//input of backward:
//pVector_fp[0] = x_fp;
//pVector_fp[1] = theta_fp;
//pVector_fp[2] = y_fp;
//pVector_fp[3] = phi_fp;
//pVector_fp[4] = x_tg;
//output:
//pVector_rec_tg[0] = x_rec ;	// = x_tg; no change
//pVector_rec_tg[1] = theta_rec;
//pVector_rec_tg[2] = y_rec;
//pVector_rec_tg[3] = phi_rec;
//pVector_rec_tg[4] = delta;

void CMP_HRS_TRANSPORTATION(double *pV5_tg, double *pV5_g2p)
{
	bool bGoodParticleA=false;
	bool bGoodParticleB=false;
	bool bGoodParticleC=false;
	bool bGoodParticleD=false;
	bool bGoodParticleE=false;
	bool bGoodParticleF=false;
	bool bGoodParticleG=false;
	bool bGoodParticleH=false;
	bool bGoodParticleI=false;
	bool bGoodParticleJ=false;

	//Use the fortran codes
	double pV5A[5]={pV5_tg[0],pV5_tg[1],pV5_tg[2],pV5_tg[3],pV5_tg[4]};
	bGoodParticleA=TransportRightHRS6(pV5A);	
	double pV5B[5]={pV5_tg[0],pV5_tg[1],pV5_tg[2],pV5_tg[3],pV5_tg[4]};			
	bGoodParticleB=TransportLeftHRS6(pV5B);			
	double pV5C[5]={pV5_tg[0],pV5_tg[1],pV5_tg[2],pV5_tg[3],pV5_tg[4]};
	bGoodParticleC=TransportRightHRS6_LargeX0(pV5C);			
	double pV5D[5]={pV5_tg[0],pV5_tg[1],pV5_tg[2],pV5_tg[3],pV5_tg[4]};
	bGoodParticleD=TransportLeftHRS6_LargeX0(pV5D);
	double pV5E[5]={pV5_tg[0],pV5_tg[1],pV5_tg[2],pV5_tg[3],pV5_tg[4]};
	bGoodParticleE=TransportRightHRS(pV5E);			
	double pV5F[5]={pV5_tg[0],pV5_tg[1],pV5_tg[2],pV5_tg[3],pV5_tg[4]};
	bGoodParticleF=TransportLeftHRS(pV5F);

	//use thrown variables to test g2p package
	//double pV5G[5]={pV5_g2p[0],pV5_g2p[1],pV5_g2p[2],pV5_g2p[3],pV5_g2p[4]};
	//bGoodParticleG=TransportRightHRS_Shim_484816_WrongBx(pV5G);
	//double pV5H[5]={pV5_g2p[0],pV5_g2p[1],pV5_g2p[2],pV5_g2p[3],pV5_g2p[4]};
	//bGoodParticleH=TransportLeftHRS_Shim_484816_WrongBx(pV5H);

	double pV5G[5]={pV5_tg[0],pV5_tg[1],pV5_tg[2],pV5_tg[3],pV5_tg[4]};
	bGoodParticleG=TransportRightHRS_Shim_484816_WrongBx(pV5G);
	double pV5H[5]={pV5_tg[0],pV5_tg[1],pV5_tg[2],pV5_tg[3],pV5_tg[4]};
	bGoodParticleH=TransportLeftHRS_Shim_484816_WrongBx(pV5H);

	double pV5I[5]={pV5_tg[0],pV5_tg[1],pV5_tg[2],pV5_tg[3],pV5_tg[4]};
	bGoodParticleI=TransportRightHRS_12p5_Min_(pV5I);
	double pV5J[5]={pV5_tg[0],pV5_tg[1],pV5_tg[2],pV5_tg[3],pV5_tg[4]};
	bGoodParticleJ=TransportLeftHRS_12p5_Min_(pV5J);

	cout.setf(ios_base::fixed);
	int prec=cout.precision(5);

	cout<<endl<<"DEBUG THE TRANSPORTATION:"<<endl;
	cout<<"Input: x_tg="<<pV5_tg[0]<<"  theta_tg="<<pV5_tg[1]<<"  y_tg="<<pV5_tg[2]
	<<"  phi_tg="<<pV5_tg[3]<<"  delta_tg="<<pV5_tg[4]<<endl<<endl;
	cout<<"Input: x_g2p="<<pV5_g2p[0]<<"  theta_g2p="<<pV5_g2p[1]<<"  y_g2p="<<pV5_g2p[2]
	<<"  phi_g2p="<<pV5_g2p[3]<<"  delta_g2p="<<pV5_g2p[4]<<endl<<endl;
	cout<<"Fortran R6:      x_fp="<<pV5A[0]<<"  theta_fp="<<pV5A[1]<<"  y_fp="<<pV5A[2]
	<<"  phi_fp="<<pV5A[3]<<(bGoodParticleA?"  succeed":"   fail")<<endl;
	cout<<"Fortran L6:      x_fp="<<pV5B[0]<<"  theta_fp="<<pV5B[1]<<"  y_fp="<<pV5B[2]
	<<"  phi_fp="<<pV5B[3]<<(bGoodParticleB?"  succeed":"   fail")<<endl;
	cout<<"Fortran R6_v2:   x_fp="<<pV5C[0]<<"  theta_fp="<<pV5C[1]<<"  y_fp="<<pV5C[2]
	<<"  phi_fp="<<pV5C[3]<<(bGoodParticleC?"  succeed":"   fail")<<endl;
	cout<<"Fortran L6_v2:   x_fp="<<pV5D[0]<<"  theta_fp="<<pV5D[1]<<"  y_fp="<<pV5D[2]
	<<"  phi_fp="<<pV5D[3]<<(bGoodParticleD?"  succeed":"   fail")<<endl;
	cout<<"Fortran Std R:   x_fp="<<pV5E[0]<<"  theta_fp="<<pV5E[1]<<"  y_fp="<<pV5E[2]
	<<"  phi_fp="<<pV5E[3]<<(bGoodParticleE?"  succeed":"   fail")<<endl;
	cout<<"Fortran Std L:   x_fp="<<pV5F[0]<<"  theta_fp="<<pV5F[1]<<"  y_fp="<<pV5F[2]
	<<"  phi_fp="<<pV5F[3]<<(bGoodParticleF?"  succeed":"   fail")<<endl;
	cout<<"Fortran G2P_S R: x_fp="<<pV5G[0]<<"  theta_fp="<<pV5G[1]<<"  y_fp="<<pV5G[2]
	<<"  phi_fp="<<pV5G[3]<<(bGoodParticleG?"  succeed":"   fail")<<endl;
	cout<<"Fortran G2P_S L: x_fp="<<pV5H[0]<<"  theta_fp="<<pV5H[1]<<"  y_fp="<<pV5H[2]
	<<"  phi_fp="<<pV5H[3]<<(bGoodParticleH?"  succeed":"   fail")<<endl;
	cout<<"Fortran G2P12.5R x_fp="<<pV5I[0]<<"  theta_fp="<<pV5I[1]<<"  y_fp="<<pV5I[2]
	<<"  phi_fp="<<pV5I[3]<<(bGoodParticleI?"  succeed":"   fail")<<endl;
	cout<<"Fortran G2P12.5L x_fp="<<pV5J[0]<<"  theta_fp="<<pV5J[1]<<"  y_fp="<<pV5J[2]
	<<"  phi_fp="<<pV5J[3]<<(bGoodParticleJ?"  succeed":"   fail")<<endl;


	if (!bGoodParticleA && !bGoodParticleE) return;
	//reconstruction
	pV5A[4] = pV5_tg[0];
	ReconstructRightHRS6(pV5A);	
	pV5B[4] = pV5_tg[0];
	ReconstructLeftHRS6(pV5B);			
	pV5C[4] = pV5_tg[0];
	ReconstructRightHRS6_LargeX0(pV5C);		
	pV5D[4] = pV5_tg[0];
	ReconstructLeftHRS6_LargeX0(pV5D);
	pV5E[4] = pV5_tg[0];
	ReconstructRightHRS(pV5E);			
	pV5F[4] = pV5_tg[0];
	ReconstructLeftHRS(pV5F);

	pV5G[4] = pV5_g2p[0];
	ReconstructRightHRS_Shim_484816_WrongBx(pV5G);			
	pV5H[4] = pV5_g2p[0];
	ReconstructLeftHRS_Shim_484816_WrongBx(pV5H);
	pV5I[4] = pV5_tg[0];
	ReconstructRightHRS_12p5_Min_(pV5I);			
	pV5J[4] = pV5_tg[0];
	ReconstructLeftHRS_12p5_Min_(pV5J);

	cout<<endl<<"DEBUG THE RECONSTRUCTION"<<endl;
	cout<<"Fortran R6:      theta_rec="<<pV5A[1]<<"  y_rec="<<pV5A[2]
	<<"  phi_rec="<<pV5A[3]<<"  delta="<<pV5A[4]<<endl;
	cout<<"Fortran L6:      theta_rec="<<pV5B[1]<<"  y_rec="<<pV5B[2]
	<<"  phi_rec="<<pV5B[3]<<"  delta="<<pV5B[4]<<endl;
	cout<<"Fortran R6_v2:   theta_rec="<<pV5C[1]<<"  y_rec="<<pV5C[2]
	<<"  phi_rec="<<pV5C[3]<<"  delta="<<pV5C[4]<<endl;
	cout<<"Fortran L6_v2:   theta_rec="<<pV5D[1]<<"  y_rec="<<pV5D[2]
	<<"  phi_rec="<<pV5D[3]<<"  delta="<<pV5D[4]<<endl;
	cout<<"Fortran Std R:   theta_rec="<<pV5E[1]<<"  y_rec="<<pV5E[2]
	<<"  phi_rec="<<pV5E[3]<<"  delta="<<pV5E[4]<<endl;
	cout<<"Fortran Std L:   theta_rec="<<pV5F[1]<<"  y_rec="<<pV5F[2]
	<<"  phi_rec="<<pV5F[3]<<"  delta="<<pV5F[4]<<endl;
	cout<<"Fortran G2P_S R: theta_rec="<<pV5G[1]<<"  y_rec="<<pV5G[2]
	<<"  phi_rec="<<pV5G[3]<<"  delta="<<pV5G[4]<<endl;
	cout<<"Fortran G2P_S L: theta_rec="<<pV5H[1]<<"  y_rec="<<pV5H[2]
	<<"  phi_rec="<<pV5H[3]<<"  delta="<<pV5H[4]<<endl;
	cout<<"Fortran G2P12.5R theta_rec="<<pV5I[1]<<"  y_rec="<<pV5I[2]
	<<"  phi_rec="<<pV5I[3]<<"  delta="<<pV5I[4]<<endl;
	cout<<"Fortran G2P12.5L theta_rec="<<pV5J[1]<<"  y_rec="<<pV5J[2]
	<<"  phi_rec="<<pV5J[3]<<"  delta="<<pV5J[4]<<endl;

	cout.unsetf(ios_base::fixed);
	cout.precision(prec);
#ifdef DEBUG_HRS_TRANSPORT 
	STOP4DEBUG;
#endif
}


void CMP_HRS_TRANSPORTATION(float *pV5_tg,float *pV5_g2p)
{
	double dV5_tg[5], dV5_g2p[5];
	for(int i=0;i<5;i++) 
	{
		dV5_tg[i]=pV5_tg[i];
		dV5_g2p[i]=pV5_g2p[i];
	}
	CMP_HRS_TRANSPORTATION(dV5_tg,dV5_g2p);
	//reset the vector and return it back to the caller
	for(int i=0;i<5;i++) 
	{
		pV5_tg[i]=(float)dV5_tg[i];
		pV5_g2p[i]=(float)dV5_g2p[i];
	}
}


//flat random number generator between [0,1)
double fRand()
{
	return double(rand())/double(RAND_MAX);
}

//boxmuller gauss number generator
double fRandGaus(double m=0, double s=1.0)	/* normal random variate generator */
{				        /* mean m, standard deviation s */
	if(s==0.0) return m;

	double x1, x2, w, y1;
	static double y2;
	static int use_last = 0;

	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {;
		x1 = 2.0 * double(rand())/double(RAND_MAX) - 1.0;
		x2 = 2.0 * double(rand())/double(RAND_MAX) - 1.0;
		w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	return( m + y1 * s );
}

//return random number in [low,High) following a*x+c prob density
double fLinearRand(double a,double c,double low,double high)
{
	double m=double(RAND_MAX);
	double x,y;
	
	double ylow=a*low+c,yhigh=a*high+c;
	if(ylow<yhigh) 
	{
		double tmp=ylow;
		ylow=yhigh;
		yhigh=tmp;
	}

	do{
		x=low+(high-low)*((double)rand())/m;
		y=ylow+(yhigh-ylow)*((double)rand())/m;

		if((a*x+c)>=y) break;
	}while(true);
	return x;
}



void DeltaCorrection(double &pDelta, double &pP_rec){;}
void XtgCorrection(double &pX,double pP_rec){;}
void ThetatgCorrection(double &pTheta,double pP_rec){;}
void YtgCorrection(double &pY,double pP_rec){;}
void PhitgCorrection(double &pPhi,double pP_rec){;}


//transport particles through HRS, use iExperiment to identify which HRS packages
//will be use, iExperiment==[10,20) is for g2p, 20 for E97110 GDH iExperiment
//Here is a list of experiment candidates:
//experiment=10: g2p septum, NO shim, 5.65 deg, Created by John, use this for g2p test run
//experiment=11: g2p septum 484816+shim, 5.65 deg, 3 cm raster, by John, 
//experiment=12: g2p septum 403216+shim, 5.65 deg, SNAKE Model not ready yet 
//experiment=13: g2p septum 400016+shim, 5.65 deg, SNAKE Model not ready yet 
//experiment=19: g2p septum 403216+shim, 5.65 deg, 2 cm raster with Wrong Bx in septum, by Min
//experiment=20: GDH exp with large X0 version
//experiment=all-other-values: Standard HRS setting, no septum field, no target field
//will add more HRS packages later

bool SNAKEThruHRS(int pIsLeftArm, double pEndPlaneAngle, double pX0_tg_m, 
				  double pHRSMomentum, int iFieldRotation, int iExperiment,
				  double* pV5_tg, double* pV5_fp)
{
	const double deg = acos(0.0)/90.0;

	//pV5_tg is the vector to be used in forward transportation functions 
	//in unit of meter, rad

	//Note that the 'anlge' defined in snake is tan(angle)
	pV5_tg[1]=tan(pV5_tg[1]);
	pV5_tg[3]=tan(pV5_tg[3]);


#if defined DEBUG_HRS_RECONSTRUCTION && (DEBUG_HRS_RECONSTRUCTION>=2)
	char str[1024];
	sprintf(str,"IN:  %8.4f %8.4f %8.4f %8.4f %8.4f",
		pV5_tg[0],pV5_tg[1],pV5_tg[2],pV5_tg[3],pV5_tg[4]);
	cout<<str<<endl;
#endif

#if defined DEBUG_HRS_RECONSTRUCTION && (DEBUG_HRS_RECONSTRUCTION>=3)
	//convert to meter and rad
	double pV5[5]={pV5_tg[0],pV5_tg[1],pV5_tg[2],pV5_tg[3],pV5_tg[4]};
	CMP_HRS_TRANSPORTATION(pV5_tg,pV5);
#endif


	bool bGoodParticle=false;
	int iSnakeFlag=(iExperiment%10);

	if( iExperiment>=10 && iExperiment<20 )
	{	
		//use the transportation without target field, which is delta dependent, not momentum dependent	
		if(pIsLeftArm) 
		{
			if(cos(pEndPlaneAngle)>cos(12.0*deg)) 
			{
				if(iSnakeFlag==0) bGoodParticle=TransportLeftHRS_g2_(pV5_tg);	//from JLL, No shim, Wrong Bx
				else if(iSnakeFlag==1) bGoodParticle=TransportLeftHRS_Shim_484816(pV5_tg); //From John, 3cm raster
				else if(iSnakeFlag==2) bGoodParticle=TransportLeftHRS_Shim_403216(pV5_tg);
				else if(iSnakeFlag==3) bGoodParticle=TransportLeftHRS_Shim_400016(pV5_tg);
				else if(iSnakeFlag==9) bGoodParticle=TransportLeftHRS_Shim_484816_WrongBx(pV5_tg); //from Min, 2cm raster, wrong Bx in septum
			}
			else
			{
				if(iSnakeFlag==9) bGoodParticle=TransportLeftHRS_12p5_Min_(pV5_tg);  //not ready yet
				else bGoodParticle=TransportLeftHRS(pV5_tg);
			}
		}
		else 
		{
			if(cos(pEndPlaneAngle)>cos(12.0*deg)) 
			{
				if(iSnakeFlag==0) bGoodParticle=TransportRightHRS_g2_(pV5_tg);	//from JLL, No shim, Wrong Bx
				else if(iSnakeFlag==1) bGoodParticle=TransportRightHRS_Shim_484816(pV5_tg); //From John, 3cm raster
				else if(iSnakeFlag==2) bGoodParticle=TransportRightHRS_Shim_403216(pV5_tg);
				else if(iSnakeFlag==3) bGoodParticle=TransportRightHRS_Shim_400016(pV5_tg);
				else if(iSnakeFlag==9) bGoodParticle=TransportRightHRS_Shim_484816_WrongBx(pV5_tg); //from Min, 2cm raster, wrong Bx in septum
			}
			else
			{
				if(iSnakeFlag==9) bGoodParticle=TransportRightHRS_12p5_Min_(pV5_tg); //not ready yet
				else bGoodParticle=TransportRightHRS(pV5_tg);
			}
		}
	}
	else if( iExperiment == 20 )
	{
		if(pIsLeftArm) 
		{
			if(cos(pEndPlaneAngle)>cos(12.0*deg)) bGoodParticle=TransportLeftHRS6_LargeX0(pV5_tg);
			else bGoodParticle=TransportLeftHRS(pV5_tg);
		}
		else 
		{
			if(cos(pEndPlaneAngle)>cos(12.0*deg)) bGoodParticle=TransportRightHRS6_LargeX0(pV5_tg);
			else bGoodParticle=TransportRightHRS(pV5_tg);
		}
	}
	else
	{
		if(pIsLeftArm) 
		{
			if(cos(pEndPlaneAngle)>cos(12.0*deg)) bGoodParticle=TransportLeftHRS6(pV5_tg);
			else bGoodParticle=TransportLeftHRS(pV5_tg);
		}
		else 
		{
			if(cos(pEndPlaneAngle)>cos(12.0*deg)) bGoodParticle=TransportRightHRS6(pV5_tg);
			else bGoodParticle=TransportRightHRS(pV5_tg);
		}
	}

	//Note that after calling forward routine, pV5_tg stores variable of focal plane
	bool bApplyVDCSmearing=true;
#if defined DEBUG_HRS_RECONSTRUCTION && (DEBUG_HRS_RECONSTRUCTION>=4)
	bApplyVDCSmearing=false;
#endif
	if(bApplyVDCSmearing)
	{
		/* Resolutions */
		//double mHRSres = 2.9e-4;
		//double mBeamRes = 3.0e-5;
		double mWireChamberRes_x = 0.0013; //m;
		double mWireChamberRes_y = 0.0013; //m;
		double mWireChamberRes_theta = 0.0003; //rad;
		double mWireChamberRes_phi = 0.0003; //rad;
		pV5_tg[0] += mWireChamberRes_x * fRandGaus(0,1.0);
		pV5_tg[2] += mWireChamberRes_y * fRandGaus(0,1.0);
		pV5_tg[1] += mWireChamberRes_theta * fRandGaus(0,1.0);
		pV5_tg[3] += mWireChamberRes_phi * fRandGaus(0,1.0);
	}

	pV5_fp[0]=pV5_tg[0];
	pV5_fp[1]=pV5_tg[1];
	pV5_fp[2]=pV5_tg[2];
	pV5_fp[3]=pV5_tg[3];
	pV5_fp[4]=pV5_tg[4];


	//now prepare for backward
	double pV5_rec[5];
	for(int k=0;k<4;k++) pV5_rec[k]=pV5_fp[k];
	pV5_rec[4]=pX0_tg_m;


	//now to smear with VDC resolution and then go backward
	if(bGoodParticle)
	{
		/*********************************************************************/
		/*                       Wire chamber smearing                       */
		/*********************************************************************/

		//now do the smearing and recostruction only if it is a good particle
		

		/*************************************************************************/
		/*                        Electron Reconstruction                        */
		/*************************************************************************/


		// Transportation from the focal plane to the target plane 
		if( iExperiment>=1 && iExperiment<20)
		{
			//use the transportation without target field, which is delta dependent, not momentum dependent	
			if(pIsLeftArm)
			{
				if(cos(pEndPlaneAngle)>cos(12.0*deg))	
				{//with septa
					if(iSnakeFlag==0) ReconstructLeftHRS_g2_(pV5_rec);  //from JLL, No shim, Wrong Bx
					else if(iSnakeFlag==1) ReconstructLeftHRS_Shim_484816(pV5_rec);    //By John, 3cm raster, 
					else if(iSnakeFlag==2) ReconstructLeftHRS_Shim_403216(pV5_rec);    
					else if(iSnakeFlag==3) ReconstructLeftHRS_Shim_400016(pV5_rec);    
					else if(iSnakeFlag==9) ReconstructLeftHRS_Shim_484816_WrongBx(pV5_rec);   //using Min's 5.69 version, 2 cm raster and Wrong Bx
				}
				else 
				{//no septa
					//John does not provide 12.5 degree resonstruction for g2p
					if(iSnakeFlag==9) ReconstructLeftHRS_12p5_Min_(pV5_rec);  //using Min's tuned reconstruction, this package need to be updated
					else ReconstructLeftHRS(pV5_rec);
				}
			}
			else 
			{
				if(cos(pEndPlaneAngle)>cos(12.0*deg))	
				{
					if(iSnakeFlag==0) ReconstructRightHRS_g2_(pV5_rec);  //from JLL, No shim, Wrong Bx
					else if(iSnakeFlag==1) ReconstructRightHRS_Shim_484816(pV5_rec);    //By John, 3cm raster,
					else if(iSnakeFlag==2) ReconstructRightHRS_Shim_403216(pV5_rec);    
					else if(iSnakeFlag==3) ReconstructRightHRS_Shim_400016(pV5_rec);    
					else if(iSnakeFlag==9) ReconstructRightHRS_Shim_484816_WrongBx(pV5_rec);   //using Min's 5.69 version, 2 cm raster and Wrong Bx
				}
				else 
				{
					//John does not provide 12.5 degree resonstruction for g2p
					if(iSnakeFlag==9) ReconstructRightHRS_12p5_Min_(pV5_rec);  //using Min's tuned reconstruction, this package need to be updated
					else ReconstructRightHRS(pV5_rec);
				}
			}
		}
		else if( iExperiment == 20 )
		{
			if(pIsLeftArm)
			{
				if(pEndPlaneAngle/deg<12.0)	ReconstructLeftHRS6_LargeX0(pV5_rec);
				else ReconstructLeftHRS(pV5_rec);
			}
			else 
			{
				if(pEndPlaneAngle/deg<12.0)	ReconstructRightHRS6_LargeX0(pV5_rec);
				else ReconstructRightHRS(pV5_rec);
			}
		}
		else
		{
			if(pIsLeftArm)
			{
				if(pEndPlaneAngle/deg<12.0)	ReconstructLeftHRS6(pV5_rec);
				else ReconstructLeftHRS(pV5_rec);
			}
			else 
			{
				if(pEndPlaneAngle/deg<12.0)	ReconstructRightHRS6(pV5_rec);
				else ReconstructRightHRS(pV5_rec);
			}
		}

		//store the reconsructed values into pV5_tg
		for(int k=0;k<5;k++)  pV5_tg[k]=pV5_rec[k];

#ifdef ThisPartNotReadyYet
		//if there is target field, x,theta,y,phi and delta all need to be corrected, 
		//these corrections are momentum dependent
		//need to find these correction later
		if(iFieldRotation>1)
		{			
			double pP_rec=pHRSMomentum*(1.0+pV5_tg[4]);
			DeltaCorrection(pV5_tg[4],pP_rec);
			XtgCorrection(pV5_tg[0],pP_rec);
			ThetatgCorrection(pV5_tg[1],pP_rec);
			YtgCorrection(pV5_tg[2],pP_rec);
			PhitgCorrection(pV5_tg[3],pP_rec);
		}
#endif

#if defined DEBUG_HRS_RECONSTRUCTION && (DEBUG_HRS_RECONSTRUCTION>=2)
		sprintf(str,"OUT: %8.4f %8.4f %8.4f %8.4f %8.4f \n",
			pV5_tg[0],pV5_tg[1],pV5_tg[2],pV5_tg[3],pV5_tg[4]);
		cout<<str<<endl;
#endif
	}


	//Note that the 'anlge' defined in snake is tan(angle)
	pV5_tg[1]=atan(pV5_tg[1]);
	pV5_tg[3]=atan(pV5_tg[3]);
	pV5_fp[1]=atan(pV5_fp[1]);
	pV5_fp[3]=atan(pV5_fp[3]);

	return bGoodParticle;
}
