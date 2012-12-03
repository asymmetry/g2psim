#include <iostream>
using namespace std;

//input:
//iExperiment=10: g2p normal, 11: g2p 484816_shim; 12 g2p 403216_shim; 13 g2p 400016_shim 
//iSmearX0=0,1,2 means X0=0, X0=X0_BPM, X0=gaussian_smeared_X0_BMP
//iSourceDistr=0,1,2 means pV5_tg[5] is in delta, flat, gaussian distribution
//iArm = 0,1,2 means random, left, right
void TestSNAKE(int iNEvent,int iExperiment=11, int iSmearX0=2, int iSourceDistr=2, int iArm=0);
void PlotSNAKE( int iExperiment, int iSmearX0, int iSourceDistr);

void Ascii2Mac(const char *file="Dp.full.f51");
int main(int argc, char** argv)
{	
  //Ascii2Mac("../HRSMC/RunHRSMC/g2psievedata/Dp.full.f51"); return 0;
	if(argc<4)
	{
		cout<<"Usage: "<<argv[0]<<" <iNEvent> <iExperiment(=10|11|12|13|19|20)> <iSmearX0(=0|1|2|3)> \n"
			<<"\t [iSourceDistr=1(0|1|2)] [iArm=0(0|1|2)]  [PlotFigureOnly=0(0|1)]\n"<<endl
			<<"\t 10<=iExperiment<20 is for g2p|gep experiment, iExperiment=20 for GDH experiment\n"
			<<"\t iExperiment=10,11,12,13 are for normal 484816 septa, 484816+shim,403216+shim,400016+shim\n"
			<<"\t iSmearX0=0,1,2 and 3 means X0_BPM=0, X0_BPM, 1mm_gaussian_smeared, 0-2.5mm_gaussian_smeared\n"
			<<"\t iSourceDistr=0,1,2 means delta,flat, gaussian distribution for x,theta,y,phi,delta\n"
			<<"\t iArm = 0,1,2 means random, left arm, right arm\n"
			<<"\t PlotFigureOnly=1 means do not create the ntuple, just use it to create figure\n"
			<<endl;
		exit(-1);
	}

	int iNEvent= atoi(argv[1]);
	int iExperiment = atoi(argv[2]);
	int iSmearX0 = atoi(argv[3]);
	int iSourceDistr = 1;
	if(argc>4) iSourceDistr = atoi(argv[4]);
	int iArm = 0; 
	if(argc>5) iArm = atoi(argv[5]);
	int iPlotFigureOnly = 0; 
	if(argc>6) iPlotFigureOnly = atoi(argv[6]);

	if(iPlotFigureOnly)
		PlotSNAKE(iExperiment,iSmearX0,iSourceDistr);
	else
		TestSNAKE(iNEvent,iExperiment,iSmearX0,iSourceDistr,iArm);

	return 0;
}
