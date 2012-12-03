#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include <iostream>
#include <fstream>

#include "HRSTransform_TCSNHCS.hh"
#include "HRSRecUseDB.hh"

using namespace std;

void Ascii2Mac(const char *file)
{
	ifstream fin(file);

	char outfile[255];
	sprintf(outfile,"%s.mac",file);
	ofstream fout(outfile);

	const double deg=acos(0.0)/90;

	double pHRSAngle=5.65*deg;
	char buff[256];
	int pHoleIndex;
	double pVar[13];  
	double pXtg_lab,pYtg_lab,pThetatg_lab,pPhitg_lab, pZtg_lab;
	
	double pV5tg_tr[5]={0,0,0,0,0};
	double pMomentum_eff=0, HRSP0=2.251;
	HRSRecUseDB *pRecDBL=new HRSRecUseDB("L","db_L.vdc.dat");

	//0           1         2     3       4     5       6       7         8     9    10    11  12
	//HoleIndex, x_fp, theta_fp, y_fp, phi_fp, urb_x, urb_y, theta_tg, phi_tg, delta, p, beam, run

	fout<<"/mydet/particleNum 1"<<endl; 
	//fout<<"/tracking/particleOnly electron"<<endl;
	fout<<"/mydet/particle1/particlePDGCode 11"<<endl;
	fout<<"/mydet/particle1/detectorAngle 5.65 deg"<<endl;
	fout<<"/mydet/particle1/randomizeInTCS 1"<<endl;
	fout<<"/mydet/leftHRSMomentum "<<HRSP0<<" GeV"<<endl;
	fout<<"/mydet/rightHRSMomentum "<<HRSP0<<" GeV"<<endl<<endl;


	int n=0;
	while (fin.good())
	{
		//get one line
		fin>>pHoleIndex;
		for(int i=1;i<12;i++) fin>>pVar[i];
		fin.getline(buff,256); //eat the rest of this line
	
		//It turns out that these tg plane variables in the ascii file is from the initial  replay
		//I have to run the reconstruction here to get some good values
		//and also the table includes multiple delta data sets, I can not use the momentum!
		//The beam position in this table is not calibrated yet
		//only the focal plane variable can be used
		
		pRecDBL->CalcTargetCoords(&pVar[1],pV5tg_tr);

		//convert TCS to HCS
		pXtg_lab=0; pYtg_lab=0; 
		Transform::X_TCS2HCS(pV5tg_tr[0],pV5tg_tr[2],0.0,pHRSAngle,pXtg_lab,pYtg_lab,pZtg_lab);
		Transform::P_TCS2HCS(pV5tg_tr[1],pV5tg_tr[3],pHRSAngle,pThetatg_lab,pPhitg_lab);

		//the X,Y reconstruction is bad
		pXtg_lab=0; pYtg_lab=0; 		
		pZtg_lab=-0.888;

		pMomentum_eff=HRSP0*(1.0+pVar[9]);

		//now build G4 cmd
		fout<<"/mydet/position3V "<<pXtg_lab<<" "<<pYtg_lab<<" "<<pZtg_lab<<" m "<<endl;
		fout<<"/mydet/particle1/theta "<<pThetatg_lab<<" rad"<<endl;
		fout<<"/mydet/particle1/phi "<<pPhitg_lab<<"  rad"<<endl;
		fout<<"/mydet/particle1/momentum  "<<pMomentum_eff<<"  GeV"<<endl; 
		fout<<"/run/beamOn \n"<<endl;
		n++;
		if(!(n%1000)) cout<<n<<" events processed\n";
		//if(n>1000) break;
	}

	fin.close();
	fout.close();
}