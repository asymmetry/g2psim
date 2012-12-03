//script to plot variants for one single sieve hole by a single click
//left click is for 2x3.5 elipse, middle click is for 3x6 elipse

#include <iostream>
#include "math.h"
using namespace std;

#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TQObject.h>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TString.h"
#include "TCut.h"
#include "TCutG.h"
#include "TPaveText.h"
#include "TText.h"
#include "TPad.h"
#include "TF1.h"
#include "TStyle.h"
#include "TTree.h"

/////////////////////////////////////////////////////////////////////////////////
void DoPlot(TCut theCut="",int holeindex=33);
void DoPlotMapping(TCut theCut="",int holeindex=9999);

//find out the first and the last bin that its entries not less than a given fraction of the maximum
void GetValidBinIndex(TH1 *h1,double fraction2max, int &firstxbin, int &lastxbin);
void GetValidBinIndexX(TH2 *h2,double fraction2max, int &firstxbin, int &lastxbin);
void GetValidBinIndexY(TH2 *h2,double fraction2max, int &firstybin, int &lastybin);

//dependence on x0_tr,y0_tr,theta0_tr,phi0_tr,z0_tr,delta
void DrawDependence(int iDrawContent=1,const char *cut="TrackClass>5");

//draw elipse centered at the the location of the mouse point
void DrawElipse(double px,double py,double a=2.0,double b=3.5);

//Plot TH1 using given strings (title,target) and cuts, will use the given RMS_factor
//to determine the histo range
TH1* DrawH1_RMS(char *hname, char *title, char *target, TCut &cut, char* hrange, 
				const char *option="", int color=1, double RMS_factor=3.0);
TH2* DrawH2_RMS(char *hname, char *title, char *target, TCut &cut, char* hrange, 
				const char *option="", int color=0, int marker=0, double RMS_factor=3.0);

//Plot TH1 using given strings (title,target,histo binning) and cuts
TH1* DrawH1(char *hname, char *title, char *target, TCut &cut, char* hrange, 
			const char *option="", int color=4);

//find out the first and the last bin that its entries not less than a given fraction of the maximum
void GetValidBinIndex(TH1 *h1,double fraction2max, int &firstxbin, int &lastxbin);
/////////////////////////////////////////////////////////////////////////////////

char			gKey[255];
double			gX0,gY0;

//Declaration of leaves types tree
Int_t           Run;
Int_t           SkimLevel;
Int_t           BookTrees;
Double_t        Beam;
Double_t        TargetM;
Double_t        TargetAtomicNumber;
Double_t        LHRSAngle;
Double_t        RHRSAngle;
Double_t        TargetXOffset;
Double_t        TargetYOffset;
Double_t        TargetZOffset;
Double_t        PivotXOffset;
Double_t        PivotYOffset;
Double_t        PivotZOffset;
Double_t        LHRSMomentum;
Double_t        RHRSMomentum;
Int_t           UseHelmField;
Double_t        HelmXOffset;
Double_t        HelmYOffset;
Double_t        HelmZOffset;
Double_t        HelmRotAxis1;
Double_t        HelmRotAxis2;
Double_t        HelmRotAxis3;
Double_t        HelmRotAngle1;
Double_t        HelmRotAngle2;
Double_t        HelmRotAngle3;
Double_t        HelmCurrentRatio=1.0;
Int_t           UseSeptumField;
Double_t        SeptumXOffset;
Double_t        SeptumYOffset;
Double_t        SeptumZOffset;
Double_t        SeptumRotAxis1;
Double_t        SeptumRotAxis2;
Double_t        SeptumRotAxis3;
Double_t        SeptumRotAngle1;
Double_t        SeptumRotAngle2;
Double_t        SeptumRotAngle3;
Double_t        SeptumCurrentRatioL;
Double_t        SeptumCurrentRatioR;
Double_t        BigBiteAngle;
Double_t        BigBiteTiltAngle;
Double_t        Pivot2BigBiteFace;

bool			bIsCombinedTree=false;	//to ideneify if this is a combined ntuple
bool            bIsLeft=true;

TCanvas *pCan=0, *pRecCan=0;

bool ReadConfig()
{
	TTree *config = (TTree*)gDirectory->Get("config");

	// Set branch addresses.
	config->SetBranchAddress("Run",&Run);
	config->SetBranchAddress("SkimLevel",&SkimLevel);
	config->SetBranchAddress("BookTrees",&BookTrees);
	config->SetBranchAddress("Beam",&Beam);
	config->SetBranchAddress("TargetM",&TargetM);
	config->SetBranchAddress("TargetAtomicNumber",&TargetAtomicNumber);
	config->SetBranchAddress("LHRSAngle",&LHRSAngle);
	config->SetBranchAddress("RHRSAngle",&RHRSAngle);
	config->SetBranchAddress("TargetXOffset",&TargetXOffset);
	config->SetBranchAddress("TargetYOffset",&TargetYOffset);
	config->SetBranchAddress("TargetZOffset",&TargetZOffset);
	config->SetBranchAddress("PivotXOffset",&PivotXOffset);
	config->SetBranchAddress("PivotYOffset",&PivotYOffset);
	config->SetBranchAddress("PivotZOffset",&PivotZOffset);
	config->SetBranchAddress("LHRSMomentum",&LHRSMomentum);
	config->SetBranchAddress("RHRSMomentum",&RHRSMomentum);
	if(config->GetBranch("HelmCurrentRatio"))
	{
		config->SetBranchAddress("UseHelmField",&UseHelmField);
		config->SetBranchAddress("HelmXOffset",&HelmXOffset);
		config->SetBranchAddress("HelmYOffset",&HelmYOffset);
		config->SetBranchAddress("HelmZOffset",&HelmZOffset);
		config->SetBranchAddress("HelmRotAxis1",&HelmRotAxis1);
		config->SetBranchAddress("HelmRotAxis2",&HelmRotAxis2);
		config->SetBranchAddress("HelmRotAxis3",&HelmRotAxis3);
		config->SetBranchAddress("HelmRotAngle1",&HelmRotAngle1);
		config->SetBranchAddress("HelmRotAngle2",&HelmRotAngle2);
		config->SetBranchAddress("HelmRotAngle3",&HelmRotAngle3);
		config->SetBranchAddress("HelmCurrentRatio",&HelmCurrentRatio);
		config->SetBranchAddress("UseSeptumField",&UseSeptumField);
		config->SetBranchAddress("SeptumXOffset",&SeptumXOffset);
		config->SetBranchAddress("SeptumYOffset",&SeptumYOffset);
		config->SetBranchAddress("SeptumZOffset",&SeptumZOffset);
		config->SetBranchAddress("SeptumRotAxis1",&SeptumRotAxis1);
		config->SetBranchAddress("SeptumRotAxis2",&SeptumRotAxis2);
		config->SetBranchAddress("SeptumRotAxis3",&SeptumRotAxis3);
		config->SetBranchAddress("SeptumRotAngle1",&SeptumRotAngle1);
		config->SetBranchAddress("SeptumRotAngle2",&SeptumRotAngle2);
		config->SetBranchAddress("SeptumRotAngle3",&SeptumRotAngle3);
		config->SetBranchAddress("SeptumCurrentRatioL",&SeptumCurrentRatioL);
		config->SetBranchAddress("SeptumCurrentRatioR",&SeptumCurrentRatioR);
		config->SetBranchAddress("BigBiteAngle",&BigBiteAngle);
		config->SetBranchAddress("BigBiteTiltAngle",&BigBiteTiltAngle);
		config->SetBranchAddress("Pivot2BigBiteFace",&Pivot2BigBiteFace);
	}
	//     This is the loop skeleton
	//       To read only selected branches, Insert statements like:
	// config->SetBranchStatus("*",0);  // disable all branches
	// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

	Long64_t nentries = config->GetEntries();

	bIsCombinedTree=false;
	if(nentries==0) return false;
	Long64_t nbytes = 0;
	double tmpAtg,tmpE,tmpZOff;
	for (Long64_t i=0; i<nentries;i++) 
	{
		nbytes += config->GetEntry(i);
		//check if this is a combined tree with various beam energies
		if(i>0 && (fabs(tmpAtg-TargetAtomicNumber)>0.01 || fabs(tmpE-Beam)>0.01 ||
			fabs(tmpZOff-PivotZOffset)>30.0 ) )  
		{
			bIsCombinedTree=true;
			break;
		}
		else
		{
			tmpAtg=TargetAtomicNumber;
			tmpE=Beam;
			tmpZOff=PivotZOffset;
		}
	}

	TH1F *hXvb=new TH1F("hXvb","",100,-500,500);

	TTree *track0=(TTree*)gROOT->FindObject("track0");

	TCanvas *pCan=new TCanvas("pCanCon","Sieve Slit",630,800);
	pCan->cd();
	track0->Draw("Xvb>>hXvb","TrackClass>5");
	double MeanXvb=hXvb->GetMean();  
	if(MeanXvb>40) bIsLeft=true; 
	else if(MeanXvb<-40) bIsLeft=false;
	else 
	{
		bIsLeft=true;
		cout<<"***warning: Left and right arm data are mixed in this ntuple!!!"<<endl;
	}

	delete pCan;
	return true;
}

//get the index and return the center position
int GetHoleIndex(double x, double y, double &x0, double &y0)
{
	double xbinsL[8]={-24.0,-18.0,-12.0,-6.0,0.0,4.5,9.0,13.5};
	double xbinsR[8]={-13.5,-9.0,-4.5,0.0,6.0,12.0,18.0,24.0};
	double *xbins=(bIsLeft)?xbinsL:xbinsR;

	//vertical span=13.3mm x 7=93.1mm, from -43.6 to 49.6
	double ybins[8]={-43.6,-30.3,-17.0,-3.7,9.6,22.9,36.2,49.5};
	int h=-1,v=-1;
	for(int i=0;i<8;i++)
	{
		if(x>=xbins[i] && x<xbins[i+1]) {h=i;x0=(xbins[h]+xbins[h+1])/2;break;} 
	}
	for(int j=0;j<8;j++)
	{
		if(y>=ybins[j] && y<ybins[j+1]) {v=j;y0=(ybins[v]+ybins[v+1])/2;break;} 
	}

	int idx=(h<0 || v<0)?-1:10*h+v; 
	std::cout<<"x="<<x<<" y="<<y<<" ==> holeindex="<<idx<<std::endl;
	return idx;
}

void SetThisStyle()
{
	
	gStyle->SetPalette(1);

	//gStyle->SetFillColor(0);  //this line cause problem for 2-D histo, use the next  
	//Fill area color
	//gStyle->SetFrameFillColor(0);
	//gStyle->SetTitleFillColor(0);
	//gStyle->SetStatColor(0);
	//gStyle->SetCanvasColor(0);
	//gStyle->SetHistFillColor(0);  //I may change this for peticular histogram

	//pad margins	
	double pPadRightMargin=0.05;
	gStyle->SetPadRightMargin(pPadRightMargin);
	gStyle->SetPadTopMargin(0.16);
	gStyle->SetPadBottomMargin(0.13);	
	gStyle->SetPadLeftMargin(0.13);

	//Title box
	//gStyle->SetTitleH(0.12);
	gStyle->SetTitleH(gStyle->GetPadTopMargin()-0.01);
	gStyle->SetTitleX(gStyle->GetPadLeftMargin()+0.01);
	gStyle->SetTitleW(1.0-gStyle->GetTitleX()-0.01);
	//gStyle->SetTitleW(1.0-gStyle->GetPadRightMargin()-gStyle->GetTitleX());
	gStyle->SetTitleBorderSize(0);


	//for stat pavetext
	gStyle->SetOptStat(0);
	//	The type of information about fit parameters printed in the histogram
	//		statistics box can be selected via the parameter mode.
	//		The parameter mode can be = pcev  (default = 0111)
	//		p = 1;  print Probability
	//		c = 1;  print Chisquare/Number of degress of freedom
	//		e = 1;  print errors (if e=1, v must be 1)
	//		v = 1;  print name/values of parameters
	//Example: gStyle->SetOptFit(1011);
	//	print fit probability, parameter names/values and errors.
	//		When "v"=1 is specified, only the non-fixed parameters are shown.
	//		When "v"=2 all parameters are shown.
	//
	//Note: gStyle->SetOptFit(1) means "default value", so it is equivalent to
	//	  gStyle->SetOptFit(111)

	//Note that the input parameter is a int in octal number system,
	//in C++, an octal number is the octal number should stat with '0'

	//gStyle->SetOptFit(0010);	//only show name and value, no uncertainties, 
	//the vertical size of stat pavetext got doubled somehow
	//need to set it to a small value
	gStyle->SetOptFit(0);

	gStyle->SetStatX(1.0-gStyle->GetPadRightMargin());
	gStyle->SetStatY(1.0-gStyle->GetPadTopMargin());
	if(gStyle->GetOptFit()==8 || gStyle->GetOptFit()==9) gStyle->SetStatH(0.08);
	else gStyle->SetStatH(0.16);
	gStyle->SetStatW(0.25);
	//gStyle->SetStatStyle(4000);

	//pad grid
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);

	//axis label
	gStyle->SetNdivisions(505,"X"); 
	gStyle->SetNdivisions(505,"Y");
	gStyle->SetLabelSize(0.07,"X");
	gStyle->SetLabelSize(0.07,"Y");
	//gStyle->SetLabelOffset(0.001,"X");
	//gStyle->SetLabelOffset(0.001,"Y");

	//axis title
	gStyle->SetTitleOffset(0.8,"X"); 
	gStyle->SetTitleOffset(0.8,"Y"); 
	gStyle->SetTitleSize(0.07,"X");
	gStyle->SetTitleSize(0.07,"Y");
	//gStyle->SetTitleColor(1,"X");
	//gStyle->SetTitleColor(1,"Y");
}

TF1* FitGaus(TH1* h1, double range_in_sigma=1.0)
{
	//cout<<"h1->GetEntries()="<<h1->GetEntries()<<endl;
	if(h1->GetEntries()<50) return NULL;

	double xmin=h1->GetMean()-1.0*h1->GetRMS();
	double xmax=h1->GetMean()+1.0*h1->GetRMS();
	h1->Fit("gaus","RQ","",xmin,xmax);
	TF1 *f=(TF1*)h1->GetListOfFunctions()->FindObject("gaus");
	double mean=f->GetParameter(1);
	double sigma=f->GetParameter(2);
	xmin=mean-range_in_sigma*sigma;
	xmax=mean+range_in_sigma*sigma;

	h1->Fit("gaus","RQ","",xmin,xmax);
	f=(TF1*)h1->GetListOfFunctions()->FindObject("gaus");
	mean=f->GetParameter(1);
	sigma=f->GetParameter(2);

	if(gStyle->GetOptFit()==0)
	{
		char str[100];
		TText *text=0;

		double xx=gStyle->GetPadLeftMargin()+0.03; 
		TPaveText *pt = new TPaveText(xx,0.20,xx+0.45,0.45,"brNDC");
		pt->SetBorderSize(0);
		pt->SetFillColor(0);
		sprintf(str,"Mean = %.3G",mean);
		text=pt->AddText(str);
		text->SetTextColor(2);
		sprintf(str,"Sigma = %.3G",sigma);
		text=pt->AddText(str);
		text->SetTextColor(2);
		pt->Draw("same");
	}

	return f;
}

TF1* FitGaus(char* hName, double range_in_sigma=1.0)
{
	TH1* h1=(TH1 *)gROOT->FindObject(hName);
	return FitGaus(h1,range_in_sigma);
}

//////////////////////////////////////////////////////////////////
//find out the first and the last bin that its entries not less than a given fraction of the maximum
void GetValidBinIndex(TH1 *h1,double fraction2max, int &firstxbin, int &lastxbin)
{
	firstxbin=0;
	lastxbin=h1->GetNbinsX(); 
	double maxentries=h1->GetMaximum(); 
	for(int ib=0;ib<lastxbin;ib++)
	{
		if(h1->GetBinContent(ib)>=fraction2max*maxentries) 
		{
			firstxbin=ib;
			break;
		}
	}	
	for(int ib=lastxbin;ib>firstxbin;ib--)
	{
		if(h1->GetBinContent(ib)>=fraction2max*maxentries) 
		{
			lastxbin=ib;
			break;
		}
	}
	return;
}


void GetValidBinIndexX(TH2 *h2,double fraction2max, int &firstxbin, int &lastxbin)
{
	cout<<"GetValidBinIndexX() h2Name="<<h2->GetName()<<endl;
	TH1D* h1= (TH1D*) gROOT->FindObject("tmp_projx");
	if(h1) delete h1;
	h2->ProjectionX("tmp_projx");
	h1= (TH1D*) gROOT->FindObject("tmp_projx");
	//cout<<"GetValidBinIndexX() h1Name="<<h1->GetName()<<endl;
	GetValidBinIndex(h1,fraction2max,firstxbin,lastxbin);
	delete h1;
	return;
}

void GetValidBinIndexY(TH2 *h2,double fraction2max, int &firstybin, int &lastybin)
{
	TH1* h1=(TH1*) (h2->ProjectionX());
	GetValidBinIndex(h1,fraction2max,firstybin,lastybin);
	delete h1;
	return;
}

///////////////////////////////////////////////////////////////////////////////////
//dependence on x0_tr,y0_tr,theta0_tr,phi0_tr,z0_tr,delta
void DrawDependence(int iDrawContent,const char *cut)
{
	//int iDrawContent=1;  //can be 1(sigma only), 2(mean+sigma),3(2-D histo+mean+sigma)
	int OldOptFitStyle=gStyle->GetOptFit();
	double OldStatY=gStyle->GetStatY();
	
	bool DoFitPol1=true;
	if(DoFitPol1)
	{
		gStyle->SetOptFit(0010);
		if(gStyle->GetOptFit()==8 || gStyle->GetOptFit()==9) gStyle->SetStatH(0.08);
		else gStyle->SetStatH(0.16);
		gStyle->SetStatY(0.4);
	}
	//TCut theCut="TrackClass>5 && abs(Xfp_tr)<800 && abs(Yvb)<100";
	system("mkdir -p graph");

	TCut theCut=cut;

	const int nitemx=6;
	const int nitemy=6;
	char *ylist[nitemy]={"X0-X_rec","Y0-Y_rec","(Theta0-Theta_rec)*1000",
		"(Phi0-Phi_rec)*1000","Z0-Z_rec","(Delta-Delta_rec)*10000"};
	char *ynamelist[nitemy]={"dX","dY","dTheta","dPhi","dZ","dDelta"};
	char *yunit[nitemy]={"(mm)","(mm)","(mrad)","(mrad)","mm","(x 10^{-4})"};

	//char *xlist[nitemx]={"X0_tr","Y0_tr","Theta0_tr*1000","Phi0_tr*1000","Z_rec","Delta*10000"};
	//char *xnamelist[nitemx]={"X0_tr","Y0_tr","Theta0_tr","Phi0_tr","Z_rec","Delta"};
	char *xlist[nitemx]={"Z0_tr","Y0_tr","Theta0_tr*1000","Phi0_tr*1000","Z_rec",
		"(Ytg_rec_tr-Y_proj2tg_tr)*10-877"
		//"Ytg_rec_tr*cos(Phitg_rec_tr)/sin(0.099+Phitg_rec_tr)-X0/tan(0.099+Phitg_rec_tr)"
		};
	char *xnamelist[nitemx]={"Z0_tr","Y0_tr","Theta0_tr","Phi0_tr","Z_rec","My_Z_rec"};
	char *xunit[nitemx]={"(mm)","(mm)","(mrad)","(mrad)","(mm)","mm"/*"(x 10^{-4})"*/};


	int ncol=iDrawContent;
	int nrow=nitemy;
	if(iDrawContent<3 || nrow>5)
	{
		if((nrow%2)==0) 
		{
			ncol*=2;
			nrow=nrow/2;
		}
		else if((nrow%3)==0) 
		{
			ncol*=3;
			nrow=nrow/3;
		}
	}
	TCanvas *pCanR=new TCanvas("pCanR","Resolution",320*ncol,240*nrow);
	pCanR->Divide(ncol,nrow,0.001,0.001);

	char hname[255],hdest[500],htitle[255];
	
	TH1 *h1mean,*h1sigma;
	TH2 *h2=0;
	TTree *track0 = (TTree*) gROOT->FindObject("track0");

	for (int j=0;j<nitemy;j++)
	{
		for (int i=0;i<nitemx;i++)
		{
			pCanR->cd(i*iDrawContent+1);
			sprintf(hname,"h2%sVs%s",ynamelist[j],xnamelist[i]);
			sprintf(htitle,"%s Vs %s; %s; %s",ynamelist[j],xnamelist[i],xlist[i],ylist[j]);
			sprintf(hdest,"%s:%s >> %s",ylist[j],xlist[i],hname);
			h2=(TH2*) gROOT->FindObject(hname);
			if(h2) delete h2;

			//cout<<"Draw 2-D histo "<<hname<<"(\""<<ylist[i]<<":"<<xlist[i]<<"\")"
			//	<<" and title=\""<<htitle<<"\"\n"; 

			track0->Draw(hdest,theCut,"");	
			h2=(TH2*) gROOT->FindObject(hname);
			h2->SetTitle(htitle);

			//search for start bin index and end bin index
			int firstxbin,lastxbin; 			
			GetValidBinIndexX(h2,0.15,firstxbin,lastxbin); 
			
			h2->FitSlicesY(0,firstxbin,lastxbin);
			
			h1mean  = (TH1*)gDirectory->Get(Form("%s_1",hname));
			h1mean->SetTitle(Form("%s Offset Vs %s; %s %s; %s Offset %s",
				ynamelist[j],xnamelist[i],xnamelist[i],xunit[i],ynamelist[j],yunit[j]));
			h1mean->SetMarkerStyle(20); 

			h1sigma = (TH1*)gDirectory->Get(Form("%s_2",hname));
			h1sigma->SetMarkerStyle(20); 
			h1sigma->SetTitle(Form("%s Resolution Vs %s; %s %s; %s Resolution %s",
				ynamelist[j],xnamelist[i],xnamelist[i],xunit[i],ynamelist[j],yunit[j]));

			
			if(DoFitPol1) h1sigma->Fit("pol1","Q");

			if(iDrawContent==1) 
			{
				pCanR->cd(i*1+1); h1sigma->Draw();	
			}
			else if(iDrawContent==2)
			{
				pCanR->cd(i*2+1); h1mean->Draw();
				pCanR->cd(i*2+2); h1sigma->Draw();	
			}
			else if(iDrawContent==3) 
			{
				pCanR->cd(i*3+1); gPad->SetRightMargin(0.15);h2->Draw("contz");
				pCanR->cd(i*3+2); h1mean->Draw();
				pCanR->cd(i*3+3); h1sigma->Draw();	
			}	

		}
		pCanR->cd();

		pCanR->SaveAs(Form("graph/Resolution_E%.3f_Helm%.0fdeg_R%.1f_%s.png",
			Beam,HelmRotAngle1*180/3.141593,HelmCurrentRatio,ynamelist[j]));	

	}
	
	if(DoFitPol1)
	{
		gStyle->SetOptFit(OldOptFitStyle);
		gStyle->SetStatY(OldStatY);
	}
}

////////////////////////////////////////////////////////////////

TH2* DrawH2_RMS(char *hname, char *title, char *target, TCut &cut, char* hrange, 
				const char *option, int color, int marker, double RMS_factor)
{
	char hdest[512];
	TTree *track0=(TTree*)gROOT->FindObject("track0");
	double pMeanX,pRMSX,pMeanY,pRMSY;

	TH2 *h2 = (TH2*) (gROOT->FindObject(hname));
	if(h2)  {delete h2;}
	sprintf(hdest,"%s >> %s",target,hname);
	track0->Draw(hdest,cut,"");
	h2 = (TH2*) (gROOT->FindObject(hname));
	/*
	TH1 *h1X=(TH1*) h2->ProjectionX();
	TH1 *h1Y=(TH1*) h2->ProjectionY();
	pMeanX=h1X->GetMean();
	pMeanY=h1X->GetMean();
	pRMSX=h1X->GetRMS(); 
	pRMSY=h1Y->GetRMS();
	delete h1X;
	delete h1Y;
	*/
	
	pMeanX=h2->GetMean(1);
	pMeanY=h2->GetMean(2);
	pRMSX=h2->GetRMS(1); 
	pRMSY=h2->GetRMS(2);

	sprintf(hrange,"(60,%.3f,%.3f,60,%.3f,%.3f)",pMeanX-RMS_factor*pRMSX,
		pMeanX+RMS_factor*pRMSX,pMeanY-RMS_factor*pRMSY, pMeanY+2.5*RMS_factor*pRMSY);
	//cout<<"DrawH2_RMS() histo range: "<<hrange<<endl;

	sprintf(hdest,"%s >> %s %s",target,hname,hrange);
	h2 = (TH2*) (gROOT->FindObject(hname));
	if(h2)  {delete h2;}
	track0->Draw(hdest,cut,option);
	h2 = (TH2*) (gROOT->FindObject(hname));
	h2->SetTitle(title);
	if(marker>0) h2->SetMarkerStyle(marker); 
	if(color>0) h2->SetMarkerColor(color);
	return h2;
}


TH2* DrawH2(char *hname, char *title, char *target, TCut &cut, char* hrange, 
				const char *option, int color, int marker)
{
	char hdest[512];
	TTree *track0=(TTree*)gROOT->FindObject("track0");
	TH2 * h2=0;
	//cout<<"DrawH2_RMS() histo range: "<<hrange<<endl;

	sprintf(hdest,"%s >> %s %s",target,hname,hrange);
	h2 = (TH2*) (gROOT->FindObject(hname));
	if(h2)  {delete h2;}
	track0->Draw(hdest,cut,option);
	h2 = (TH2*) (gROOT->FindObject(hname));
	h2->SetTitle(title);
	if(marker>0) h2->SetMarkerStyle(marker); 
	if(color>0) h2->SetMarkerColor(color);
	return h2;
}


//Plot TH1 using given strings (title,target) and cuts, will use the given RMS_factor
//to determine the histo range
TH1* DrawH1_RMS(char *hname, char *title, char *target, TCut &cut, char* hrange, 
				const char *option, int color, double RMS_factor)
{
	char hdest[512];
	TH1 *h1=0;
	TTree *track0=(TTree*)gROOT->FindObject("track0");

	//if the option contains 'same' or 'SAME', do not search the range
	if(!(strstr(option,"same") || strstr(option,"SAME")))  
	{
		double pMean,pRMS;
		h1 = (TH1*) (gROOT->FindObject(hname));
		if(h1)  {delete h1;}
		sprintf(hdest,"%s >> %s",target,hname);
		track0->Draw(hdest,cut,"");
		h1 = (TH1*) (gROOT->FindObject(hname));
		pMean=h1->GetMean();
		pRMS=h1->GetRMS();
		sprintf(hrange,"(60,%.3f,%.3f)",pMean-RMS_factor*pRMS,pMean+RMS_factor*pRMS);
	}

	sprintf(hdest,"%s >> %s %s",target,hname,hrange);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,cut,option);
	h1 = (TH1*) (gROOT->FindObject(hname));
	h1->SetTitle(title);
	if(color>0) h1->SetLineColor(color); 
	h1->GetYaxis()->SetRangeUser(h1->GetMinimum(),1.2*h1->GetMaximum());  
	return h1;
}

//Plot TH1 using given strings (title,target,histo binning) and cuts
TH1* DrawH1(char *hname, char *title, char *target, TCut &cut, char* hrange, 
			const char *option, int color)
{
	char hdest[512];
	TTree *track0=(TTree*)gROOT->FindObject("track0");

	sprintf(hdest,"%s >> %s %s",target,hname,hrange);
	TH1* h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,cut,option);
	h1 = (TH1*) (gROOT->FindObject(hname));
	h1->SetTitle(title);
	h1->SetLineColor(color);
	if(!(strstr(option,"same") && strstr(option,"SAME")))   
		h1->GetYaxis()->SetRangeUser(h1->GetMinimum(),1.2*h1->GetMaximum());  
	return h1;
}



///////////////////////////////////////////////////////////////////
//Draw the reconstruction results
void PlotRec(char* extraCut="(1>0)")
{

	TCut all="TrackClass>5 && abs(Xfp_tr)<800 && abs(Yvb)<100";
	TCut theCut=all;

	if(gROOT->FindObject("CUTGVBXY")) theCut+="CUTGVBXY";	
	if(gROOT->FindObject("CUTG")) theCut+="CUTG";
	theCut+=extraCut;


	DoPlot(theCut,0);return;
	/*
	///////////////////////////////////////////////////////////////////	

	TTree *track0=(TTree*)gROOT->FindObject("track0");

	if(!pRecCan) pRecCan=new TCanvas("pTgCan","Sieve Slit",650,30,1350,860);
	pRecCan->Clear();
	pRecCan->Divide(5,5,0.001,0.001);

	TH1 *h1, *h12; h1=0; h12=0;
	TH2 *h2, *h22; h2=0; h22=0;
	char hdest[255],hname[255],target[255],htitle[255];
	double pMean=0,pRMS=0;
	char pHistoRange[255];

	////////////////////////////////////
	int holeindex=0;

	pRecCan->cd(1);
	gPad->SetRightMargin(0.15); 

	sprintf(pHistoRange,"(12,%.1f,%.1f,12,%.1f,%.1f)",gX0-4.5,gX0+4.5,gY0-13.3,gY0+13.3);

	sprintf(hname,"h2VH_vb_%02d",holeindex);
	sprintf(hdest,"-Xvb_tr:-Yvb_tr >> %s %s",hname,pHistoRange);
	h2 = (TH2*) (gROOT->FindObject(hname));
	if(h2)  {delete h2;}
	track0->Draw(hdest,theCut,"contz");
	h2 = (TH2*) (gROOT->FindObject(hname));
	h2->SetTitle(Form("Sieve Hole %02d: Hits at VB; -Y_{vb}^{tr} (mm) ; -X_{vb}^{tr} (mm) ",holeindex));

	TText *text=0;
	double xx=1.0-gPad->GetRightMargin(); 
	double yy=1.0-gPad->GetTopMargin()-0.01; 
	TPaveText *pt = new TPaveText(xx-0.5,yy-0.3,xx,yy,"brNDC");
	pt->SetBorderSize(0);
	pt->SetFillColor(0);
	//text=pt->AddText(Form("Sieve Hole %02d",holeindex));
	//text->SetTextColor(2);
	if(UseHelmField)
	{
		text=pt->AddText(Form("TargetFieldRatio=%0.1f",HelmCurrentRatio));
	}
	else
	{
		text=pt->AddText("NO Target Field");
	}
	text->SetTextColor(2);
	text=pt->AddText(Form("Beam = %.3f",Beam));
	text->SetTextColor(2);
	if(!bIsCombinedTree)
	{  
		text=pt->AddText(Form("Target_A = %.0f",TargetAtomicNumber));
	} 
	else
	{  
		text=pt->AddText("Multi-Targets");	
	}
	text->SetTextColor(2);
	pt->Draw("same");


	////////////////////////////////////
	pRecCan->cd(2);

	sprintf(hname,"h1X_vb_%02d_tmp",holeindex);
	sprintf(hdest,"-Xvb_tr >> %s",hname);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	pMean=h1->GetMean();
	delete h1;

	sprintf(pHistoRange,"(50,%.3f,%.3f)",pMean-5,pMean+5);

	sprintf(hname,"h1X_vb_%02d",holeindex);
	sprintf(hdest,"-Xvb_tr >> %s %s",hname,pHistoRange);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	h1->SetTitle("-X_{vb}^{tr} ; -X_{vb}^{tr} (mm) ");

	sprintf(hname,"h1X_proj2sl_%02d",holeindex);
	sprintf(hdest,"-X_proj2sl_tr >> %s",hname);
	h12 = (TH1*) (gROOT->FindObject(hname));
	if(h12)  {delete h12;}
	track0->Draw(hdest,theCut,"same");
	h12 = (TH1*) (gROOT->FindObject(hname));
	h12->SetLineColor(4);
	h1->SetTitle("-X_{vb}^{tr}(black), -X_{img2sl}^{tr}(blue) (mm); -X_{vb}^{tr} or -X_{img2sl}^{tr} (mm) ");


	////////////////////////////////////
	pRecCan->cd(3);
	sprintf(hname,"h1Y_vb_%02d_tmp",holeindex);
	sprintf(hdest,"-Yvb_tr >> %s",hname);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	pMean=h1->GetMean();
	delete h1;

	sprintf(pHistoRange,"(50,%.3f,%.3f)",pMean-5,pMean+5);

	sprintf(hname,"h1Y_vb_%02d",holeindex);
	sprintf(hdest,"-Yvb_tr >> %s %s",hname,pHistoRange);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	h1->SetTitle("-Y_{vb}^{tr} ; -Y_{vb}^{tr} (mm) ");

	sprintf(hname,"h1Y_proj2sl_%02d",holeindex);
	sprintf(hdest,"-Y_proj2sl_tr >> %s %s",hname,pHistoRange);
	h12 = (TH1*) (gROOT->FindObject(hname));
	if(h12)  {delete h12;}
	track0->Draw(hdest,theCut,"same");
	h12 = (TH1*) (gROOT->FindObject(hname));
	h12->SetLineColor(4);
	h1->SetTitle("-Y_{vb}^{tr}(black), -Y_{img2sl}^{tr}(blue) (mm); -Y_{vb}^{tr} or -Y_{img2sl}^{tr} (mm) ");


	////////////////////////////////////
	pRecCan->cd(4);

	gPad->SetRightMargin(0.15); 
	sprintf(hname,"h2VH_fp_%02d",holeindex);
	sprintf(hdest,"-Xfp_tr:-Yfp_tr >> %s",hname);
	h2 = (TH2*) (gROOT->FindObject(hname));
	if(h2)  {delete h2;}
	track0->Draw(hdest,theCut,"contz");
	h2 = (TH2*) (gROOT->FindObject(hname));
	h2->SetTitle("Focus Plane: -X_{fp}^{tr} Vs -Y_{fp}^{tr} (mm); -Yfp (mm) ; -Xfp (mm) ");

	////////////////////////////////////
	pRecCan->cd(5);

	gPad->SetRightMargin(0.15); 
	sprintf(hname,"h2TP_fp_%02d",holeindex);
	sprintf(hdest,"Thetafp_tr*1000:Phifp_tr*1000 >> %s",hname);
	h2 = (TH2*) (gROOT->FindObject(hname));
	if(h2)  {delete h2;}
	track0->Draw(hdest,theCut,"contz");
	h2 = (TH2*) (gROOT->FindObject(hname));
	h2->SetTitle("Focus Plane: #theta_{fp}^{tr} Vs #phi_{fp}^{tr} (mrad); #phi_{fp}^{tr} (mrad) ; #theta_{fp}^{tr} (mrad) ");

	/////////////////////////////////////
	pRecCan->cd(6);

	if(gROOT->FindObject("h1Pvb")) delete gROOT->FindObject("h1Pvb");
	track0->Draw("Pvb >> h1Pvb",theCut);
	h12 = (TH1*) (gROOT->FindObject("h1Pvb"));
	double MeanPvb=h12->GetMean();

	sprintf(hname,"h1Pvb_%02d",holeindex);
	sprintf(hdest,"Pvb >> %s",hname);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	if(!bIsCombinedTree)
	{
		sprintf(hdest,"Pvb >> %s",hname);
		track0->Draw(Form("Pvb >> %s(50,%.3f,%.3f)",hname,
			MeanPvb-0.005,MeanPvb+0.005),theCut);
	}
	else
	{
		if(Beam<1.3)
			track0->Draw(Form("Pvb >> %s(75,%.3f,%.3f)",hname,
			MeanPvb-0.01,MeanPvb+0.005),theCut);
		else if (Beam<1.8)
			track0->Draw(Form("Pvb >> %s(125,%.3f,%.3f)",hname,
			MeanPvb-0.02,MeanPvb+0.005),theCut);
		else
			track0->Draw(Form("Pvb >> %s(175,%.3f,%.3f)",hname,
			MeanPvb-0.03,MeanPvb+0.005),theCut);
	}
	track0->Draw(hdest,theCut);
	h1 = (TH1*) (gROOT->FindObject(hname));
	h1->SetTitle("P_{0}(black), P_{vb}(red) and P_{rec}(blue); P (GeV)");
	h1->SetLineColor(2);

	track0->Draw("P0",theCut,"same");

	TH1 *h1Prec=0; 
	h1Prec = (TH1*) (gROOT->FindObject("h1Prec"));
	if(h1Prec) delete h1Prec;
	track0->Draw("P_rec >> h1Prec",theCut,"same");
	h1Prec = (TH1*) (gROOT->FindObject("h1Prec"));
	h1Prec->SetLineColor(4);


	////////////////////////////////////
	pRecCan->cd(7);

	sprintf(hname,"h1X0_%02d_tmp",holeindex);
	sprintf(hdest,"X0 >> %s",hname);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	pMean=h1->GetMean();
	delete h1;

	sprintf(pHistoRange,"(50,%.3f,%.3f)",pMean-5,pMean+5);

	sprintf(hname,"h1X0_%02d",holeindex);
	sprintf(hdest,"X0 >> %s %s",hname,pHistoRange);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	h1->SetTitle("X_{0} ; X_{0} (mm) ");

	sprintf(hname,"h1X_rec_%02d",holeindex);
	sprintf(hdest,"X_rec >> %s %s",hname,pHistoRange);
	h12 = (TH1*) (gROOT->FindObject(hname));
	if(h12)  {delete h12;}
	track0->Draw(hdest,theCut,"same");
	h12 = (TH1*) (gROOT->FindObject(hname));
	h12->SetLineColor(4);
	h1->SetTitle("X_{0}(black), X_{rec}(blue) (mm); X_{0} or X_{rec} (mm) ");

	////////////////////////////////////
	pRecCan->cd(8);

	sprintf(hname,"h1Y0_%02d_tmp",holeindex);
	sprintf(hdest,"Y0 >> %s",hname);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	pMean=h1->GetMean();
	delete h1;

	sprintf(pHistoRange,"(50,%.3f,%.3f)",pMean-5,pMean+5);

	sprintf(hname,"h1Y0_%02d",holeindex);
	sprintf(hdest,"Y0 >> %s %s",hname,pHistoRange);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	h1->SetTitle("Y_{0} ; Y_{0} (mm) ");

	sprintf(hname,"h1Y_rec_%02d",holeindex);
	sprintf(hdest,"Y_rec >> %s %s",hname,pHistoRange);
	h12 = (TH1*) (gROOT->FindObject(hname));
	if(h12)  {delete h12;}
	track0->Draw(hdest,theCut,"same");
	h12 = (TH1*) (gROOT->FindObject(hname));
	h12->SetLineColor(4);
	h1->SetTitle("Y_{0}(black), Y_{rec}(blue) (mm); Y_{0} or Y_{rec} (mm) ");

	////////////////////////////////////
	pRecCan->cd(9);

	sprintf(hname,"h1Theta0_%02d_tmp",holeindex);
	sprintf(hdest,"Theta0*1000 >> %s",hname);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	pMean=h1->GetMean();
	delete h1;

	sprintf(pHistoRange,"(60,%.3f,%.3f)",pMean-15,pMean+15);

	sprintf(hname,"h1Theta0_%02d",holeindex);
	sprintf(hdest,"Theta0*1000 >> %s %s",hname,pHistoRange);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	h1->SetTitle("#theta_{0} ; #theta_{0} (mrad) ");

	sprintf(hname,"h1Theta_rec_%02d",holeindex);
	sprintf(hdest,"Theta_rec*1000 >> %s %s",hname,pHistoRange);
	h12 = (TH1*) (gROOT->FindObject(hname));
	if(h12)  {delete h12;}
	track0->Draw(hdest,theCut,"same");
	h12 = (TH1*) (gROOT->FindObject(hname));
	h12->SetLineColor(4);
	h1->SetTitle("#theta_{0}(black), #theta_{rec}(blue) (mrad); #theta_{0} or #theta_{rec} (mrad) ");

	////////////////////////////////////
	pRecCan->cd(10);

	sprintf(hname,"h1Phi0_%02d_tmp",holeindex);
	sprintf(hdest,"Phi0*1000 >> %s",hname);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	pMean=h1->GetMean();
	delete h1;

	sprintf(pHistoRange,"(50,%.3f,%.3f)",pMean-100,pMean+100);

	sprintf(hname,"h1Phi0_%02d",holeindex);
	sprintf(hdest,"Phi0*1000 >> %s %s",hname,pHistoRange);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	h1->SetTitle("#phi_{0} ; #phi_{0} (mrad) ");

	sprintf(hname,"h1Phi_rec_%02d",holeindex);
	sprintf(hdest,"Phi_rec*1000 >> %s %s",hname,pHistoRange);
	h12 = (TH1*) (gROOT->FindObject(hname));
	if(h12)  {delete h12;}
	track0->Draw(hdest,theCut,"same");
	h12 = (TH1*) (gROOT->FindObject(hname));
	h12->SetLineColor(4);
	h1->SetTitle("#phi_{0}(black), #phi_{rec}(blue) (mrad); #phi_{0} or #phi_{rec} (mrad) ");


	////////////////////////////////////
	pRecCan->cd(11);

	sprintf(hname,"h1Z_rec_%02d_tmp",holeindex);
	sprintf(hdest,"Z_rec >> %s",hname);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	pMean=h1->GetMean();
	pRMS=h1->GetRMS();
	delete h1;

	if(pRMS<1.0) sprintf(pHistoRange,"(60,%.3f,%.3f)",pMean-3,pMean+3);
	else sprintf(pHistoRange,"(60,%.3f,%.3f)",pMean-3*pRMS,pMean+3*pRMS);
	//sprintf(pHistoRange,"(60,%.3f,%.3f)",pMean-3*pRMS,pMean+3*pRMS);
	//cout<<"Vertex Z histo range: "<<pHistoRange<<endl;

	sprintf(hname,"h1Z0_%02d",holeindex);
	sprintf(hdest,"Z0 >> %s %s",hname,pHistoRange);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	h1->SetTitle("Z_{0} ; Z_{0} (mm) ");

	sprintf(hname,"h1Z_rec_%02d",holeindex);
	sprintf(hdest,"Z_rec >> %s %s",hname,pHistoRange);
	h12 = (TH1*) (gROOT->FindObject(hname));
	if(h12)  {delete h12;}
	track0->Draw(hdest,theCut,"same");
	h12 = (TH1*) (gROOT->FindObject(hname));
	h12->SetLineColor(4);
	h1->SetTitle("Z_{0}(black), Z_{rec}(blue) (mm); Z_{0} or Z_{rec} (mm) ");

	////////////////////////////////////
	pRecCan->cd(12);

	sprintf(hname,"h1X0_tr_%02d_tmp",holeindex);
	sprintf(hdest,"X0_tr >> %s",hname);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	pMean=h1->GetMean();
	delete h1;

	sprintf(pHistoRange,"(60,%.3f,%.3f)",pMean-6,pMean+6);

	sprintf(hname,"h1X0_tr_%02d",holeindex);
	sprintf(target,"X0_tr");
	sprintf(htitle,"X_{0}^{tr} (mm); X_{0}^{tr} (mm) ");
	h1=DrawH1(hname,htitle,target,theCut,pHistoRange,"",1);

	sprintf(hname,"h1X_rec_tr_%02d",holeindex);
	sprintf(target,"X_rec_tr");
	sprintf(htitle,"X_{rec}^{tr} (mm); X_{rec}^{tr} (mm) ");
	h12=DrawH1(hname,htitle,target,theCut,pHistoRange,"same",4);

	h1->SetTitle("X_{0}^{tr}(black), X_{rec}^{tr}(blue) (mm); X_{0}^{tr} or X_{rec}^{tr} (mm) ");


	////////////////////////////////////
	pRecCan->cd(13);

	sprintf(hname,"h1Y0_tr_%02d_tmp",holeindex);
	sprintf(hdest,"Y0_tr >> %s",hname);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	pMean=h1->GetMean();
	delete h1;

	sprintf(pHistoRange,"(60,%.3f,%.3f)",pMean-6,pMean+6);

	sprintf(hname,"h1Y0_tr_%02d",holeindex);
	sprintf(target,"Y0_tr");
	sprintf(htitle,"Y_{0}^{tr} (mm); Y_{0}^{tr} (mm) ");
	h1=DrawH1(hname,htitle,target,theCut,pHistoRange,"",1);

	sprintf(hname,"h1Y_rec_tr_%02d",holeindex);
	sprintf(target,"Y_rec_tr");
	sprintf(htitle,"Y_{rec}^{tr} (mm); Y_{rec}^{tr} (mm) ");
	h12=DrawH1(hname,htitle,target,theCut,pHistoRange,"same",4);

	h1->SetTitle("Y_{0}^{tr}(black), Y_{rec}^{tr}(blue) ; Y_{0}^{tr} or Y_{rec}^{tr} (mm) ");


	////////////////////////////////////
	pRecCan->cd(14);

	sprintf(hname,"h1Theta0_tr_%02d_tmp",holeindex);
	sprintf(hdest,"Theta0_tr*1000 >> %s",hname);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	pMean=h1->GetMean();
	delete h1;

	sprintf(pHistoRange,"(60,%.3f,%.3f)",pMean-15,pMean+15);

	sprintf(hname,"h1Theta0_tr_%02d",holeindex);
	sprintf(target,"Theta0_tr*1000");
	sprintf(htitle,"#theta_{0}^{tr} (mrad); #theta_{0}^{tr} (mrad) ");
	h1=DrawH1(hname,htitle,target,theCut,pHistoRange,"",1);

	sprintf(hname,"h1Theta_rec_tr_%02d",holeindex);
	sprintf(target,"Theta_rec_tr*1000");
	sprintf(htitle,"#theta_{rec}^{tr} (mrad); #theta_{rec}^{tr} (mrad) ");
	h12=DrawH1(hname,htitle,target,theCut,pHistoRange,"same",4);

	h1->SetTitle("#theta_{0}^{tr}(black), #theta_{rec}^{tr}(blue); #theta_{0}^{tr} or #theta_{rec}^{tr} (mrad) ");

	////////////////////////////////////
	pRecCan->cd(15);

	sprintf(hname,"h1Phi0_tr_%02d_tmp",holeindex);
	sprintf(hdest,"Phi0_tr >> %s",hname);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	pMean=h1->GetMean();
	delete h1;

	sprintf(pHistoRange,"(60,%.3f,%.3f)",pMean-15,pMean+15);

	sprintf(hname,"h1Phi0_tr_%02d",holeindex);
	sprintf(target,"Phi0_tr*1000");
	sprintf(htitle,"#phi_{0}^{tr} (mrad); #phi_{0}^{tr} (mrad) ");
	h1=DrawH1(hname,htitle,target,theCut,pHistoRange,"",1);

	sprintf(hname,"h1Phi_rec_tr_%02d",holeindex);
	sprintf(target,"Phi_rec_tr*1000");
	sprintf(htitle,"#phi_{rec}^{tr} (mrad); #phi_{rec}^{tr} (mrad) ");
	h12=DrawH1(hname,htitle,target,theCut,pHistoRange,"same",4);

	h1->SetTitle("#phi_{0}^{tr}(black), #phi_{rec}^{tr}(blue); #phi_{0}^{tr} or #phi_{rec}^{tr} (mrad) ");

	////////////////////////////////////
	pRecCan->cd(16);

	sprintf(hname,"h1DD_%02d",holeindex);
	sprintf(target,"(Delta-Delta_rec)*10000");
	sprintf(htitle,"(Delta-Delta_rec) ( 10^{-4}); Delta-Delta_rec ( 10^{-4})");
	sprintf(pHistoRange,"(100,-10,10)");

	h1=DrawH1(hname,htitle,target,theCut,pHistoRange,"",1);
	FitGaus(h1);

	////////////////////////////////////
	pRecCan->cd(17);
	sprintf(hname,"h1DX_snake_%02d",holeindex);
	sprintf(target,"X_proj2tg_tr - X_rec2tg_tr ");
	sprintf(htitle,"X_{proj2tg}^{tr} - X_{snake2tg}^{tr} (mm); X_{proj2tg}^{tr} - X_{snake2tg}^{tr} (mm) ");
	sprintf(pHistoRange,"(100,-5,5)");

	h1=DrawH1(hname,htitle,target,theCut,pHistoRange,"",1);
	//for delta distribution, gaus fit willnot work
	pRMS=h1->GetRMS();
	if(pRMS>0.01) FitGaus(h1);

	////////////////////////////////////
	pRecCan->cd(18);

	sprintf(hname,"h1DY_snake_%02d",holeindex);
	sprintf(target,"Y_proj2tg_tr - Y_rec2tg_tr");
	sprintf(htitle,"Y_{proj2tg}^{tr} - Y_{snake2tg}^{tr} (mm); Y_{proj2tg}^{tr} - Y_{snake2tg}^{tr} (mm) ");
	sprintf(pHistoRange,"(100,-5,5)");

	h1=DrawH1(hname,htitle,target,theCut,pHistoRange,"",1);
	FitGaus(h1);

	////////////////////////////////////
	pRecCan->cd(19);

	sprintf(hname,"h1DTheta0_snake_%02d",holeindex);
	sprintf(target,"(Thetavb_tr-Theta_rec2tg_tr)*1000");
	sprintf(htitle,"#theta_{proj2tg}^{tr} - #theta_{snake2tg}^{tr} (mrad); #theta_{proj2tg}^{tr} - #theta_{snake2tg}^{tr} (mrad)");
	sprintf(pHistoRange,"(100,-10,10)");

	h1=DrawH1(hname,htitle,target,theCut,pHistoRange,"",1);
	FitGaus(h1);

	////////////////////////////////////
	pRecCan->cd(20);

	sprintf(hname,"h1DPhi0_snake_%02d",holeindex);
	sprintf(target,"(Phivb_tr-Phi_rec2tg_tr)*1000");
	sprintf(htitle,"#phi_{proj2tg}^{tr} - #phi_{snake2tg}^{tr} (mrad); #phi_{proj2tg}^{tr} - #phi_{snake2tg}^{tr} (mrad)");
	sprintf(pHistoRange,"(100,-10,10)");

	h1=DrawH1(hname,htitle,target,theCut,pHistoRange,"",1);
	FitGaus(h1);

	////////////////////////////////////
	pRecCan->cd(21);

	sprintf(hname,"h1DZ0_%02d",holeindex);
	sprintf(target,"Z0 - Z_rec");
	sprintf(htitle,"Z_{0} - Z_{rec} (mm); Z_{0} - Z_{rec0} (mm) ");

	h1=DrawH1_RMS(hname,htitle,target,theCut,pHistoRange,"",1,3.0);
	FitGaus(h1);

	////////////////////////////////////
	pRecCan->cd(22);

	sprintf(hname,"h1DX0_%02d",holeindex);
	sprintf(target,"X0 - X_rec");
	sprintf(htitle,"X_{0} - X_{rec} (mm); X_{0} - X_{rec0} (mm) ");

	h1=DrawH1_RMS(hname,htitle,target,theCut,pHistoRange,"",1,3.0);
	FitGaus(h1);

	////////////////////////////////////
	pRecCan->cd(23);

	sprintf(hname,"h1DY0_%02d",holeindex);
	sprintf(target,"Y0 - Y_rec");
	sprintf(htitle,"Y_{0} - Y_{rec} (mm); Y_{0} - Y_{rec0} (mm) ");

	h1=DrawH1_RMS(hname,htitle,target,theCut,pHistoRange,"",1,3.0);
	FitGaus(h1);

	////////////////////////////////////
	pRecCan->cd(24);

	sprintf(hname,"h1DTheta0_%02d",holeindex);
	sprintf(target,"(Theta0-Theta_rec)*1000");
	sprintf(htitle,"#theta_{0} - #theta_{rec} (mrad); #theta_{0} - #theta_{rec} (mrad)");

	h1=DrawH1_RMS(hname,htitle,target,theCut,pHistoRange,"",1,3.0);
	FitGaus(h1);

	////////////////////////////////////
	pRecCan->cd(25);

	sprintf(hname,"h1DPhi0_%02d",holeindex);
	sprintf(htitle,"#phi_{0} - #phi_{rec} (mrad); #phi_{0} - #phi_{rec} (mrad)");
	sprintf(target,"(Phi0-Phi_rec)*1000");

	h1=DrawH1_RMS(hname,htitle,target,theCut,pHistoRange,"",1,3.0);
	FitGaus(h1);


	////////////////////////////////////
	pRecCan->cd();
	system("mkdir -p graph");
	cout<<"gKey="<<gKey<<endl;
	if(!bIsCombinedTree)
	{  
		pRecCan->SaveAs(Form("graph/SieveHole%02d_E%.3f_Helm%.0fdeg_R%.1f_A%.0f_%s.png",holeindex,
			Beam,HelmRotAngle1*180/3.141593,HelmCurrentRatio,TargetAtomicNumber,gKey));	
	} 
	else
	{ 
		pRecCan->SaveAs(Form("graph/SieveHole%02d_E%.3f_Helm%.0fdeg_R%.1f_combined_%s.png",holeindex,
			Beam,HelmRotAngle1*180/3.141593,HelmCurrentRatio,gKey));
	}
	*/
}



///////////////////////////////////////////////////////////////////

void VarInHole(int holeindex=33)
{
	if(!gROOT->FindObject("CUTGVBXY")) return;


	TCut all="TrackClass>5 && abs(Xfp_tr)<800 && abs(Yvb)<100";
	TCut theCut="CUTGVBXY";
	theCut += all;
	if(gROOT->FindObject("CUTG")) theCut += "CUTG";

	//DoPlot(theCut,holeindex);	
	DoPlotMapping( theCut, holeindex);
}

//plot for mapping purpose,
//just plot x,theta,y,phi for fp plane, vb plane, target plane
//and also z0,delta,thtea0,phi0
void DoPlotMapping(TCut theCut,int holeindex)
{
	TTree *track0=(TTree*)gROOT->FindObject("track0");
	TCanvas *pMapCan=0;
	if(!pMapCan) pMapCan=new TCanvas("pMapCan","mapping",650,30,1350,860);
	pMapCan->Clear();
	pMapCan->Divide(5,5,0.001,0.001);

	TH1 *h1, *h12; h1=0; h12=0;
	TH2 *h2, *h22; h2=0; h22=0;
	char hdest[255],hname[255],target[255],htitle[255];
	double pMean=0,pRMS=0;
	char pHistoRange[255];

	const char *strKey[5]={"X","Y","Theta","Phi","Delta"};
	const char *strKeyTitle[5]={"X","Y","#theta","#phi","#delta"};
	const char *strUnit[5]={"(mm)","(mm)","(mrad)","(mrad)","(x 10^{-4})"};

	bool bPlotTargetPlane=true;
	////////////////////////////////////

	pMapCan->cd(1);
	gPad->SetRightMargin(0.15); 

	if(holeindex>0 && holeindex<78)
	{
		sprintf(pHistoRange,"(12,%.1f,%.1f,18,%.1f,%.1f)",gX0-4.5,gX0+4.5,gY0-13.3,gY0+13.3 + 13.3);
		sprintf(hname,"h2VH_vb_%02d",holeindex);
		sprintf(target,"-Xvb_tr:-Yvb_tr");
		sprintf(htitle,"Sieve Hole %02d: Hits at Sieve; -Y_{vb}^{tr} (mm) ; -X_{vb}^{tr} (mm)",holeindex);
	}
	else
	{
		sprintf(hname,"h1VH_vb_%02d",holeindex);
		sprintf(target,"-Xvb_tr:-Yvb_tr");
		sprintf(htitle,"%s; -Y_{vb}^{tr} (mm) ; -X_{vb}^{tr} (mm) ",gKey);
	}
		
	h2=DrawH2_RMS(hname,htitle,target,theCut,pHistoRange,"contz",1,4);
	

	TText *text=0;
	double xx=1.0-gPad->GetRightMargin()-0.01; 
	double yy=1.0-gPad->GetTopMargin()-0.01; 
	TPaveText *pt = new TPaveText(gPad->GetLeftMargin()+0.01,yy-0.3,xx,yy,"brNDC");
	pt->SetBorderSize(0);
	pt->SetFillColor(0);
	pt->SetTextColor(2);
 
	//if(holeindex>0 && holeindex<78) text=pt->AddText("Sieve Hole %02d",holeindex);
	if(UseHelmField)
	{
		text=pt->AddText(Form("TargetFieldRatio=%0.1f",HelmCurrentRatio));
	}
	else
	{
		text=pt->AddText("NO Target Field");
	}
	text=pt->AddText(Form("Beam = %.3f",Beam));
	if(!bIsCombinedTree)
	{  
		text=pt->AddText(Form("Target_A = %.0f",TargetAtomicNumber));
	} 
	else
	{  
		text=pt->AddText("Multi-Targets");	
	}
	pt->Draw("same");


	////////////////////////////////////
	//pMapCan->cd(2);
	
	for(int i=0; i<4;i++)
	{
		pMapCan->cd(i+2);

		sprintf(htitle,"%s_{sl}^{tr}(black) %s_{proj2sl}^{tr}(blue) %s",
			strKeyTitle[i],strKeyTitle[i],strUnit[i]);
		
		//get range
		sprintf(hname,"h1%s_proj2sl_tr_%02d",strKey[i],holeindex);
		if(i<2) sprintf(target,"%s_proj2sl_tr",strKey[i]);
		else sprintf(target,"%s_rec2tg_tr*1000",strKey[i]);

		h12=DrawH1_RMS(hname,htitle,target,theCut,pHistoRange,"",4,3.0);

		sprintf(hname,"h1%s_vb_tr_%02d",strKey[i],holeindex);
		if(i<2) sprintf(target,"%svb_tr",strKey[i]);
		else sprintf(target,"%svb_tr*1000",strKey[i]);

		h1=DrawH1(hname,htitle,target,theCut,pHistoRange,"",1);
		
		pRMS=h1->GetRMS();
		if(!(i<2 && pRMS<0.01)) FitGaus(h1);
		
		h12->Draw("same");
	} 

	/////////////////////////////////////
	pMapCan->cd(6);
	
	sprintf(hname,"h1Pvb_%02d",holeindex);
	sprintf(hdest,"Pvb >> %s %s",hname,pHistoRange);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	pMean=h1->GetMean();
	pRMS=h1->GetRMS();
	delete h1;
	sprintf(pHistoRange,"(50,%.3f,%.3f)",pMean-3*pRMS,pMean+2*pRMS);

	sprintf(hname,"h1P0_%02d",holeindex);
	sprintf(hdest,"P0 >> %s %s",hname,pHistoRange);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	//h1->SetTitle("P_{0} ; P_{0} (GeV/c) ");
	h1->SetTitle("P_{0}(black), P_{vb}(red) ; P_{0} or P_{vb}(GeV/c)");

	sprintf(hname,"h1Pvb_%02d",holeindex);
	sprintf(hdest,"Pvb >> %s %s",hname,pHistoRange);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"same");
	h1 = (TH1*) (gROOT->FindObject(hname));
	h1->SetTitle("P_{vb} ; P_{vb} (GeV/c) ");
	h1->SetLineColor(2);


	////////////////////////////////////
	//pMapCan->cd(7);

	for(int i=0; i<4;i++)
	{
		pMapCan->cd(i+7);

		sprintf(htitle,"%s_{0}^{tr}(black), %s_{rec}^{tr}(blue) %s",
			strKeyTitle[i],strKeyTitle[i],strUnit[i]);

		//this part to get the range
		sprintf(hname,"h1%s_rec_tr",strKey[i]);
		if(i<2) sprintf(target,"%s_rec_tr",strKey[i]);
		else sprintf(target,"%s_rec_tr*1000",strKey[i]);
		h1=DrawH1_RMS(hname,htitle,target,theCut,pHistoRange,"",4,3.0);

		sprintf(hname,"h1%s0_tr",strKey[i]);
		if(i<2) sprintf(target,"%s0_tr",strKey[i]);
		else sprintf(target,"%s0_tr*1000",strKey[i]);
		h12=DrawH1(hname,htitle,target,theCut,pHistoRange,"same",1);

		if(h12->GetMaximum() > 1.2*h1->GetMaximum())
		{
			h12->Draw();
			h1->Draw("same");
		}
	} 

	////////////////////////////////////
	pMapCan->cd(11);

	sprintf(hname,"h1Z0_%02d",holeindex);
	sprintf(target,"Z0");
	//sprintf(htitle,"Z_{0} (mm); Z_{0} (mm) ");
	sprintf(htitle,"Z_{0}(black), Z_{rec}(blue) (mm); Z_{0} or Z_{rec} (mm) ");
	h1=DrawH1_RMS(hname,htitle,target,theCut,pHistoRange,"",1,3.0);

	sprintf(hname,"h1Z_rec_%02d",holeindex);
	sprintf(target,"Z_rec");
	h12=DrawH1(hname,htitle,target,theCut,pHistoRange,"same",4);
	
	if(h12->GetMaximum() > 1.2*h1->GetMaximum())
	{
		h12->Draw();
		h1->Draw("same");
	}

	////////////////////////////////////
	//pMapCan->cd(12);

	for(int i=0; i<4;i++)
	{
		pMapCan->cd(i+12);

		sprintf(htitle,"%s_{0}(black), %s_{rec}(blue) %s",
			strKeyTitle[i],strKeyTitle[i],strUnit[i]);

		//this part to get the range
		sprintf(hname,"h1%s_rec",strKey[i]);
		if(i<2) sprintf(target,"%s_rec",strKey[i]);
		else sprintf(target,"%s_rec*1000",strKey[i]);
		h1=DrawH1_RMS(hname,htitle,target,theCut,pHistoRange,"",4,3.0);

		sprintf(hname,"h1%s0",strKey[i]);
		if(i<2) sprintf(target,"%s0",strKey[i]);
		else sprintf(target,"%s0*1000",strKey[i]);
		h12=DrawH1(hname,htitle,target,theCut,pHistoRange,"same",1);
		if(h12->GetMaximum() > 1.2*h1->GetMaximum())
		{
			h12->Draw();
			h1->Draw("same");
		}
	} 


	////////////////////////////////////
	pMapCan->cd(16);

	sprintf(hname,"h1Delta_%02d",holeindex);
	sprintf(target,"Delta*10000");
	sprintf(htitle,"Delta (10^{-4}); Delta(10^{-4})");
	sprintf(pHistoRange,"(100,-500,500)");

	h1=DrawH1(hname,htitle,target,theCut,pHistoRange,"",1);

	////////////////////////////////////
	//pMapCan->cd(17);
	
	for(int i=0; i<4;i++)
	{
		pMapCan->cd(i+17);

		sprintf(htitle,"%s_{fp}^{tr} %s",strKeyTitle[i],strUnit[i]);
		sprintf(hname,"h1%s_fp_tr_%02d",strKey[i],holeindex);
		sprintf(target,"%sfp_tr",strKey[i]);
		
		h1=DrawH1_RMS(hname,htitle,target,theCut,pHistoRange,"",1,3.0);
		
	} 
	
	////////////////////////////////////
	pMapCan->cd(21);

	sprintf(hname,"h1DD_%02d",holeindex);
	sprintf(target,"(Delta-Delta_rec)*10000");
	sprintf(htitle,"(Delta-Delta_rec)  (10^{-4}); Delta-Delta_rec (10^{-4})");
	sprintf(pHistoRange,"(100,-10,10)");

	h1=DrawH1(hname,htitle,target,theCut,pHistoRange,"",1);
	FitGaus(h1);


	////////////////////////////////////
	//pRecCan->cd(22);

	for(int i=0; i<4;i++)
	{
		pMapCan->cd(i+22);

		sprintf(htitle,"%s_{0}^{tr} - %s_{rec}^{tr} %s",
			strKeyTitle[i],strKeyTitle[i],strUnit[i]);

		//this part to get the range
		sprintf(hname,"h1D%s_tr_%02d",strKey[i],holeindex);
		if(i<2) sprintf(target,"%s0_tr-%s_rec_tr",strKey[i],strKey[i]);
		else sprintf(target,"(%s0_tr-%s_rec_tr)*1000",strKey[i],strKey[i]);
		h1=DrawH1_RMS(hname,htitle,target,theCut,pHistoRange,"",4,3.0);
		pRMS=h1->GetRMS();
		if(!(i<2 && pRMS<0.01)) FitGaus(h1);
	} 

	////////////////////////////////////

	pMapCan->cd();
	system("mkdir -p graph");
	//cout<<"gKey="<<gKey<<endl;
	if(!bIsCombinedTree)
	{  
		pMapCan->SaveAs(Form("graph/SieveHole%02d_E%.3f_Helm%.0fdeg_R%.1f_A%.0f_%s.png",holeindex,
			Beam,HelmRotAngle1*180/3.141593,HelmCurrentRatio,TargetAtomicNumber,gKey));	
	} 
	else
	{ 
		pMapCan->SaveAs(Form("graph/SieveHole%02d_E%.3f_Helm%.0fdeg_R%.1f_combined_%s.png",holeindex,
			Beam,HelmRotAngle1*180/3.141593,HelmCurrentRatio,gKey));
	}
}

///////////////////////////////////////////////////////////////////	
void DoPlot(TCut theCut,int holeindex)
{
	TTree *track0=(TTree*)gROOT->FindObject("track0");
	TCanvas *pRecCan=0;
	if(!pRecCan) pRecCan=new TCanvas("pTgCan","Sieve Slit",650,30,1350,860);
	pRecCan->Clear();
	pRecCan->Divide(5,5,0.001,0.001);

	TH1 *h1, *h12; h1=0; h12=0;
	TH2 *h2, *h22; h2=0; h22=0;
	char hdest[255],hname[255],target[255],htitle[255];
	double pMean=0,pRMS=0;
	char pHistoRange[255];

	const char *strKey[5]={"X","Y","Theta","Phi","Delta"};
	const char *strKeyTitle[5]={"X","Y","#theta","#phi","#delta"};
	const char *strUnit[5]={"(mm)","(mm)","(mrad)","(mrad)","(x 10^{-4})"};

	bool bPlotTargetPlane=true;
	////////////////////////////////////

	pRecCan->cd(1);
	gPad->SetRightMargin(0.15); 

	if(holeindex>0 && holeindex<78)
	{
		sprintf(pHistoRange,"(12,%.1f,%.1f,18,%.1f,%.1f)",gX0-4.5,gX0+4.5,gY0-13.3,gY0+13.3 + 13.3);
		sprintf(hname,"h2VH_vb_%02d",holeindex);
		sprintf(target,"-Xvb_tr:-Yvb_tr");
		sprintf(htitle,"Sieve Hole %02d: Hits at Sieve; -Y_{vb}^{tr} (mm) ; -X_{vb}^{tr} (mm)",holeindex);
	}
	else
	{
		sprintf(hname,"h1VH_vb_%02d",holeindex);
		sprintf(target,"-Xvb_tr:-Yvb_tr");
		sprintf(htitle,"%s; -Y_{vb}^{tr} (mm) ; -X_{vb}^{tr} (mm) ",gKey);
	}
		
	h2=DrawH2_RMS(hname,htitle,target,theCut,pHistoRange,"contz",1,4);
	

	TText *text=0;
	double xx=1.0-gPad->GetRightMargin()-0.01; 
	double yy=1.0-gPad->GetTopMargin()-0.01; 
	TPaveText *pt = new TPaveText(gPad->GetLeftMargin()+0.01,yy-0.3,xx,yy,"brNDC");
	pt->SetBorderSize(0);
	pt->SetFillColor(0);
	pt->SetTextColor(2);
 
	//if(holeindex>0 && holeindex<78) text=pt->AddText("Sieve Hole %02d",holeindex);
	if(UseHelmField)
	{
		text=pt->AddText(Form("TargetFieldRatio=%0.1f",HelmCurrentRatio));
	}
	else
	{
		text=pt->AddText("NO Target Field");
	}
	text=pt->AddText(Form("Beam = %.3f",Beam));
	if(!bIsCombinedTree)
	{  
		text=pt->AddText(Form("Target_A = %.0f",TargetAtomicNumber));
	} 
	else
	{  
		text=pt->AddText("Multi-Targets");	
	}
	pt->Draw("same");


	////////////////////////////////////
	if(!bPlotTargetPlane)
	{
	pRecCan->cd(2);

	sprintf(hname,"h1X_vb_%02d_tmp",holeindex);
	sprintf(hdest,"-Xvb_tr >> %s",hname);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	pMean=h1->GetMean();
	delete h1;

	sprintf(pHistoRange,"(50,%.3f,%.3f)",pMean-5,pMean+5);

	sprintf(hname,"h1X_vb_%02d",holeindex);
	sprintf(hdest,"-Xvb_tr >> %s %s",hname,pHistoRange);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	h1->SetTitle("-X_{vb}^{tr} ; -X_{vb}^{tr} (mm) ");

	sprintf(hname,"h1X_proj2sl_%02d",holeindex);
	sprintf(hdest,"-X_proj2sl_tr >> %s",hname);
	h12 = (TH1*) (gROOT->FindObject(hname));
	if(h12)  {delete h12;}
	track0->Draw(hdest,theCut,"same");
	h12 = (TH1*) (gROOT->FindObject(hname));
	h12->SetLineColor(4);
	h1->SetTitle("-X_{vb}^{tr}(black), -X_{img2sl}^{tr}(blue) (mm); -X_{vb}^{tr} or -X_{img2sl}^{tr} (mm) ");


	////////////////////////////////////
	pRecCan->cd(3);
	sprintf(hname,"h1Y_vb_%02d_tmp",holeindex);
	sprintf(hdest,"-Yvb_tr >> %s",hname);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	pMean=h1->GetMean();
	delete h1;

	sprintf(pHistoRange,"(50,%.3f,%.3f)",pMean-5,pMean+5);

	sprintf(hname,"h1Y_vb_%02d",holeindex);
	sprintf(hdest,"-Yvb_tr >> %s %s",hname,pHistoRange);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	h1->SetTitle("-Y_{vb}^{tr} ; -Y_{vb}^{tr} (mm) ");

	sprintf(hname,"h1Y_proj2sl_%02d",holeindex);
	sprintf(hdest,"-Y_proj2sl_tr >> %s %s",hname,pHistoRange);
	h12 = (TH1*) (gROOT->FindObject(hname));
	if(h12)  {delete h12;}
	track0->Draw(hdest,theCut,"same");
	h12 = (TH1*) (gROOT->FindObject(hname));
	h12->SetLineColor(4);
	h1->SetTitle("-Y_{vb}^{tr}(black), -Y_{img2sl}^{tr}(blue) (mm); -Y_{vb}^{tr} or -Y_{img2sl}^{tr} (mm) ");


	////////////////////////////////////
	pRecCan->cd(4);

	gPad->SetRightMargin(0.15); 
	sprintf(hname,"h2VH_fp_%02d",holeindex);
	sprintf(hdest,"-Xfp_tr:-Yfp_tr >> %s",hname);
	h2 = (TH2*) (gROOT->FindObject(hname));
	if(h2)  {delete h2;}
	track0->Draw(hdest,theCut,"contz");
	h2 = (TH2*) (gROOT->FindObject(hname));
	h2->SetTitle("Focus Plane: -X_{fp}^{tr} Vs -Y_{fp}^{tr} (mm); -Yfp (mm) ; -Xfp (mm) ");

	////////////////////////////////////
	pRecCan->cd(5);

	gPad->SetRightMargin(0.15); 
	sprintf(hname,"h2TP_fp_%02d",holeindex);
	sprintf(hdest,"Thetafp_tr*1000:Phifp_tr*1000 >> %s",hname);
	h2 = (TH2*) (gROOT->FindObject(hname));
	if(h2)  {delete h2;}
	track0->Draw(hdest,theCut,"contz");
	h2 = (TH2*) (gROOT->FindObject(hname));
	h2->SetTitle("Focus Plane: #theta_{fp}^{tr} Vs #phi_{fp}^{tr} (mrad); #phi_{fp}^{tr} (mrad) ; #theta_{fp}^{tr} (mrad) ");
}
	/////////////////////////////////////
	pRecCan->cd(6);
	
	sprintf(hname,"h1P_rec_%02d_tmp",holeindex);
	sprintf(hdest,"P_rec >> %s",hname);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	pMean=h1->GetMean();
	pRMS=h1->GetRMS();
	delete h1;
	sprintf(pHistoRange,"(50,%.3f,%.3f)",pMean-3*pRMS,pMean+2*pRMS);

	sprintf(hname,"h1P0_%02d",holeindex);
	sprintf(hdest,"P0 >> %s %s",hname,pHistoRange);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	//h1->SetTitle("P_{0} ; P_{0} (GeV/c) ");
	h1->SetTitle("P_{0}(black), P_{vb}(red), P_{rec}(blue) ; P_{0}, P_{vb}, or P_{rec} (GeV/c)");

	
	sprintf(hname,"h1Pvb_%02d",holeindex);
	sprintf(hdest,"Pvb >> %s %s",hname,pHistoRange);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"same");
	h1 = (TH1*) (gROOT->FindObject(hname));
	h1->SetTitle("P_{vb} ; P_{vb} (GeV/c) ");
	h1->SetLineColor(2);

	
	sprintf(hname,"h1P_rec_%02d",holeindex);
	sprintf(hdest,"P_rec >> %s %s",hname,pHistoRange);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"same");
	h1 = (TH1*) (gROOT->FindObject(hname));
	h1->SetTitle("P_{rec} ; P_{rec} (GeV/c) ");
	h1->SetLineColor(4);


	////////////////////////////////////
	//pRecCan->cd(7);

	for(int i=0; i<4;i++)
	{
		pRecCan->cd(i+7);

		sprintf(htitle,"%s_{0}^{tr}(black), %s_{rec}^{tr}(blue) %s",
			strKeyTitle[i],strKeyTitle[i],strUnit[i]);

		//this part to get the range
		sprintf(hname,"h1%s_rec_tr",strKey[i]);
		if(i<2) sprintf(target,"%s_rec_tr",strKey[i]);
		else sprintf(target,"%s_rec_tr*1000",strKey[i]);
		h1=DrawH1_RMS(hname,htitle,target,theCut,pHistoRange,"",4,3.0);

		sprintf(hname,"h1%s0_tr",strKey[i]);
		if(i<2) sprintf(target,"%s0_tr",strKey[i]);
		else sprintf(target,"%s0_tr*1000",strKey[i]);
		h12=DrawH1(hname,htitle,target,theCut,pHistoRange,"same",1);
		if(h12->GetMaximum() > 1.2*h1->GetMaximum())
		{
			h12->Draw();
			h1->Draw("same");
		}
	} 

	////////////////////////////////////
	pRecCan->cd(11);

	sprintf(hname,"h1Z_rec_%02d_tmp",holeindex);
	sprintf(hdest,"Z_rec >> %s",hname);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	track0->Draw(hdest,theCut,"");
	h1 = (TH1*) (gROOT->FindObject(hname));
	pMean=h1->GetMean();
	pRMS=h1->GetRMS();
	delete h1;

	if(pRMS<1.0) sprintf(pHistoRange,"(60,%.3f,%.3f)",pMean-15,pMean+15);
	else sprintf(pHistoRange,"(60,%.3f,%.3f)",pMean-3*pRMS,pMean+3*pRMS);
	//sprintf(pHistoRange,"(60,%.3f,%.3f)",pMean-3*pRMS,pMean+3*pRMS);
	//cout<<"Vertex Z histo range: "<<pHistoRange<<endl;

	
	sprintf(hname,"h1Z0_%02d",holeindex);
	sprintf(target,"Z0");
	//sprintf(htitle,"Z_{0} (mm); Z_{0} (mm) ");
	sprintf(htitle,"Z_{0}(black), Z_{rec}(blue) (mm); Z_{0} or Z_{rec} (mm) ");
	h1=DrawH1(hname,htitle,target,theCut,pHistoRange,"",1);

	sprintf(hname,"h1Z_rec_%02d",holeindex);
	sprintf(target,"Z_rec");
	h12=DrawH1(hname,htitle,target,theCut,pHistoRange,"same",4);
	
	if(h12->GetMaximum() > 1.2*h1->GetMaximum())
	{
		h12->Draw();
		h1->Draw("same");
	}

	////////////////////////////////////
	//pRecCan->cd(12);

	for(int i=0; i<4;i++)
	{
		pRecCan->cd(i+12);

		sprintf(htitle,"%s_{0}(black), %s_{rec}(blue) %s",
			strKeyTitle[i],strKeyTitle[i],strUnit[i]);

		//this part to get the range
		sprintf(hname,"h1%s_rec",strKey[i]);
		if(i<2) sprintf(target,"%s_rec",strKey[i]);
		else sprintf(target,"%s_rec*1000",strKey[i]);
		h1=DrawH1_RMS(hname,htitle,target,theCut,pHistoRange,"",4,3.0);

		sprintf(hname,"h1%s0",strKey[i]);
		if(i<2) sprintf(target,"%s0",strKey[i]);
		else sprintf(target,"%s0*1000",strKey[i]);
		h12=DrawH1(hname,htitle,target,theCut,pHistoRange,"same",1);
		if(h12->GetMaximum() > 1.2*h1->GetMaximum())
		{
			h12->Draw();
			h1->Draw("same");
		}
	} 


	////////////////////////////////////
	pRecCan->cd(16);

	sprintf(hname,"h1DD_%02d",holeindex);
	sprintf(target,"(Delta-Delta_rec)*10000");
	sprintf(htitle,"(Delta-Delta_rec)  (10^{-4}); Delta-Delta_rec (10^{-4})");
	sprintf(pHistoRange,"(100,-10,10)");

	h1=DrawH1(hname,htitle,target,theCut,pHistoRange,"",1);
	FitGaus(h1);


	////////////////////////////////////
	//pRecCan->cd(17);
	
	for(int i=0; i<4;i++)
	{
		pRecCan->cd(i+17);

		sprintf(htitle,"%s_{proj2tg}^{tr} - %s_{snakerec}^{tr} %s",
			strKeyTitle[i],strKeyTitle[i],strUnit[i]);

		//this part to get the range
		sprintf(hname,"h1D%s_snake_tr_%02d",strKey[i],holeindex);
		if(i<2)
		{
		sprintf(target,"%s_proj2tg_tr-%s_rec2tg_tr",strKey[i],strKey[i]);
		}
		else
		{
		sprintf(target,"(%svb_tr-%s_rec2tg_tr)*1000",strKey[i],strKey[i]);
		}
		h1=DrawH1_RMS(hname,htitle,target,theCut,pHistoRange,"",1,3.0);
		pRMS=h1->GetRMS();
		if(!(i<2 && pRMS<0.01)) FitGaus(h1);
	} 
	

	////////////////////////////////////
	pRecCan->cd(21);

	sprintf(hname,"h1DZ0_%02d",holeindex);
	sprintf(target,"Z0 - Z_rec");
	sprintf(htitle,"Z_{0} - Z_{rec} (mm); Z_{0} - Z_{rec0} (mm) ");

	h1=DrawH1_RMS(hname,htitle,target,theCut,pHistoRange,"",1,3.0);
	FitGaus(h1);

	////////////////////////////////////
	//pRecCan->cd(22);

	for(int i=0; i<4;i++)
	{
		pRecCan->cd(i+22);

		sprintf(htitle,"%s_{0} - %s_{rec} %s",
			strKeyTitle[i],strKeyTitle[i],strUnit[i]);

		//this part to get the range
		sprintf(hname,"h1D%s_lab_%02d",strKey[i],holeindex);
		if(i<2) sprintf(target,"%s0-%s_rec",strKey[i],strKey[i]);
		else sprintf(target,"(%s0-%s_rec)*1000",strKey[i],strKey[i]);
		h1=DrawH1_RMS(hname,htitle,target,theCut,pHistoRange,"",4,3.0);
		pRMS=h1->GetRMS();
		if(!(i<2 && pRMS<0.01)) FitGaus(h1);
	} 


	////////////////////////////////////
	//Just plot the target plate
	////////////////////////////////////
	if (bPlotTargetPlane)
	{
		//pRecCan->cd(2);
		for(int i=0; i<4;i++)
		{
			pRecCan->cd(i+2);

			sprintf(htitle,"%s_{tg}^{tr} - %s_{drift2tg}^{tr} %s",
				strKeyTitle[i],strKeyTitle[i],strUnit[i]);

			//this part to get the range
			sprintf(hname,"h1D%s_tgplane_tr_%02d",strKey[i],holeindex);
			if(i<2)
			{
				sprintf(target,"%stg_tr-%stg_rec_tr",strKey[i],strKey[i]);
			}
			else
			{
				sprintf(target,"(%stg_tr-%stg_rec_tr)*1000",strKey[i],strKey[i]);
			}
			h1=DrawH1_RMS(hname,htitle,target,theCut,pHistoRange,"",4,3.0);
			pRMS=h1->GetRMS();
			if(!(i<2 && pRMS<0.01)) FitGaus(h1);
		} 
	}
	////////////////////////////////////

	pRecCan->cd();
	system("mkdir -p graph");
	//cout<<"gKey="<<gKey<<endl;
	if(!bIsCombinedTree)
	{  
		pRecCan->SaveAs(Form("graph/SieveHole%02d_E%.3f_Helm%.0fdeg_R%.1f_A%.0f_%s.png",holeindex,
			Beam,HelmRotAngle1*180/3.141593,HelmCurrentRatio,TargetAtomicNumber,gKey));	
	} 
	else
	{ 
		pRecCan->SaveAs(Form("graph/SieveHole%02d_E%.3f_Helm%.0fdeg_R%.1f_combined_%s.png",holeindex,
			Beam,HelmRotAngle1*180/3.141593,HelmCurrentRatio,gKey));
	}
}


///////////////////////////////////////////////////////////////////
//draw elipse centered at the the location of the mouse point
void DrawElipse(double px,double py,double a,double b)
{
	double x = gPad->AbsPixeltoX(px);
	double y = gPad->AbsPixeltoY(py);

	//redefine CUTGVBXY
	if(gROOT->FindObject("CUTGVBXY")) delete gROOT->FindObject("CUTGVBXY");
	const int kNpt=9;
	TCutG *cutg = new TCutG("CUTGVBXY",kNpt);
	cutg->SetVarX("-Yvb_tr");
	cutg->SetVarY("-Xvb_tr");
	cutg->SetTitle("SelectedHoleCut");
	cutg->SetFillColor(0);
	cutg->SetMarkerStyle(20);
	cutg->SetLineWidth(2);

	for(int i=0;i<=kNpt;i++)
	{
		double phi=i*4*acos(0.0)/kNpt;       
		cutg->SetPoint(i,x+a*cos(phi),y+b*sin(phi));
	}

	cutg->Draw("same");

	gPad->Update();
}


/////////////////////////////////////////////////////////////////
void ReDrawMainCan(TCanvas *pCan)
{
	TCut theCut="TrackClass>5";
	
	if(gROOT->FindObject("CUTGVBXY")) theCut+="CUTGVBXY";	
	if(gROOT->FindObject("CUTG")) theCut+="CUTG";

	TTree* track0=(TTree*)gROOT->FindObject("track0");

	pCan->cd(3);
	gPad->SetRightMargin(0.15);
	if(gROOT->FindObject("h2VH_Tg")) delete gROOT->FindObject("h2VH_Tg");
	//TH2F *h2VH_Tg=new TH2F("h2VH_Tg","-X0_tr:-Y0_tr ; -Y_{0}^{tr} (mm); -X_{0}^{tr} (mm)",
	//	  120,-12.0,12.0,120,-10.0,12.0);
	track0->Draw("-X0_tr:-Y0_tr >> h2VH_Tg",theCut,"contz"); 
	TH2* h2VH_Tg=(TH2*)gROOT->FindObject("h2VH_Tg");
	//h2VH_Tg->SetTitle("-X0_tr:-Y0_tr (cut)");

	pCan->cd(4);
	gPad->SetRightMargin(0.15);
	if(gROOT->FindObject("h2TP_Tg")) delete gROOT->FindObject("h2TP_Tg");
	//TH2F *h2TP_Tg=new TH2F("h2TP_Tg","Theta0_tr:Phi0_tr ; #theta_{0}^{tr} (rad); #phi_{0}^{tr} (rad)",
	//	  30,-0.03,0.03,700,-0.25,0.1);
	track0->Draw("Theta0_tr:Phi0_tr >> h2TP_Tg",theCut,"contz"); 
	TH2* h2TP_Tg=(TH2*)gROOT->FindObject("h2TP_Tg");
	//h2TP_Tg->SetTitle("Theta0_tr:Phi0_tr (cut)");

	pCan->Modified();
	pCan->Update();
}


/////////////////////////////////////////////////////////////////
void DoEvent(Int_t event, Int_t x, Int_t y, TObject *selected)
{
	//do not response to mouse movement
	if(event>20 || event<=10) return;

	double px = gPad->AbsPixeltoX(x);
	gX0 = gPad->PadtoX(px);
	double py = gPad->AbsPixeltoY(y);
	gY0 = gPad->PadtoY(py);

	//printf("event=%d, x=%d, y=%d, selected_type=%s, selected_name=%s, p2x=%f, p2y=%f\n", event, x, y, selected->IsA()->GetName(),selected->GetName(),px,py);

	TVirtualPad *padsav = 0;

	int holeindex=-1;
	if(strcmp("h2VH_VB",selected->GetName())==0) 
	{
		padsav = gPad;
		if(event==11) DrawElipse(x,y);
		else if(event==12) DrawElipse(x,y,3.0,6.0);
		holeindex=GetHoleIndex(px,py,gX0,gY0);
		
		ReDrawMainCan(padsav->GetCanvas()); 
		if(holeindex>0) VarInHole(holeindex);
	}
	else if(strcmp("CUTGVBXY",selected->GetName())==0) 
	{
		padsav = gPad;
		holeindex=GetHoleIndex(px,py,gX0,gY0);
		
		ReDrawMainCan(padsav->GetCanvas()); 
		if(holeindex>0) VarInHole(holeindex);
	}

	if(padsav) 
	{
		//cout<<"cd "<<padsav->GetName()<<endl; 
		padsav->cd();
	}

}


/////////////////////////////////////////////////////////////////

void DoFieldEffect(TCut &theCut,int holeindex)
{
	TTree *track0=(TTree*)gROOT->FindObject("track0");

	TCanvas *pCanF=new TCanvas("pCanF","Field effect in Target Plane",650,30,1350,860);
	pCanF->Clear();
	pCanF->Divide(4,4,0.001,0.001);

	TH1 *h1, *h12; h1=0; h12=0;
	TH2 *h2, *h22; h2=0; h22=0;
	char hdest[255],hname[255],target[255],htitle[255];
	double pMean=0,pRMS=0;
	char pHistoRange[255];

	const char *strKey[5]={"X","Y","Theta","Phi","Delta"};
	const char *strKeyTitle[5]={"X","Y","#theta","#phi","#delta"};
	const char *strUnit[5]={"(mm)","(mm)","(mrad)","(mrad)","(x 10^{-4})"};

	bool bPlotTargetPlane=true;
	////////////////////////////////////

	int firstxbin, lastxbin;
	double xxmin, xxmax;

	pCanF->cd(1);
	gPad->SetRightMargin(0.15); 

	if(holeindex>0 && holeindex<78)
	{
		sprintf(pHistoRange,"(12,%.1f,%.1f,18,%.1f,%.1f)",gX0-4.5,gX0+4.5,gY0-13.3,gY0+13.3 + 13.3);
		sprintf(hname,"h2VH_vb_%02d",holeindex);
		sprintf(target,"-Xvb_tr:-Yvb_tr");
		sprintf(htitle,"Sieve Hole %02d: Hits at Sieve; -Y_{vb}^{tr} (mm) ; -X_{vb}^{tr} (mm)",holeindex);
	}
	else
	{
		sprintf(hname,"h1VH_vb_%02d",holeindex);
		sprintf(target,"-Xvb_tr:-Yvb_tr");
		sprintf(htitle,"%s; -Y_{vb}^{tr} (mm) ; -X_{vb}^{tr} (mm) ",gKey);
	}
		
	h2=DrawH2_RMS(hname,htitle,target,theCut,pHistoRange,"contz",1,4);
	

	TText *text=0;
	double xx=1.0-gPad->GetRightMargin()-0.01; 
	double yy=1.0-gPad->GetTopMargin()-0.01; 
	TPaveText *pt = new TPaveText(gPad->GetLeftMargin()+0.01,yy-0.3,xx,yy,"brNDC");
	pt->SetBorderSize(0);
	pt->SetFillColor(0);
	pt->SetTextColor(2);
 
	//if(holeindex>0 && holeindex<78) text=pt->AddText("Sieve Hole %02d",holeindex);
	if(UseHelmField)
	{
		text=pt->AddText(Form("TargetFieldRatio=%0.1f",HelmCurrentRatio));
	}
	else
	{
		text=pt->AddText("NO Target Field");
	}
	text=pt->AddText(Form("Beam = %.3f",Beam));
	if(!bIsCombinedTree)
	{  
		text=pt->AddText(Form("Target_A = %.0f",TargetAtomicNumber));
	} 
	else
	{  
		text=pt->AddText("Multi-Targets");	
	}
	pt->Draw("same");


	////////////////////////////////////
	pCanF->cd(2);

	sprintf(hname,"h2DTheta_tgplaneVsPvb_%02d",holeindex);
	sprintf(target,"(Thetavb_tr-Theta0_tr)*1000:Pvb");
	sprintf(htitle,"#theta_{proj2tg}^{tr} - #theta_{tg}^{tr} Vs P_{vb}; P_{vb} (GeV) ; #theta_{proj2tg}^{tr} - #theta_{tg}^{tr} (mrad)");

	h2=DrawH2_RMS(hname,htitle,target,theCut,pHistoRange,"",0,0,3.0);

	GetValidBinIndexX(h2,0.15,firstxbin,lastxbin);
	xxmin=h2->GetXaxis()->GetBinCenter(firstxbin);
	xxmax=h2->GetXaxis()->GetBinCenter(lastxbin);

	h2=DrawH2(hname,htitle,target,theCut,pHistoRange,"prof",4,20);
	TF1 f1_3("f1_3","[0]/x",xxmin,xxmax);
	f1_3.SetLineColor(6); 
	h2->Fit(&f1_3,"R","",xxmin,xxmax);
	//h2->Fit("pol3","R","",h2->GetXaxis()->GetBinCenter(10),h2->GetXaxis()->GetBinCenter(h2->GetNbinsX()-10));

	
	if(gStyle->GetOptFit()!=0)
		pt = new TPaveText(0.5,0.42,0.9,0.55,"brNDC");
	else 
		pt = new TPaveText(0.5,0.34,0.9,0.6,"brNDC");
	
	pt->SetBorderSize(0);
	pt->SetFillColor(0);
	pt->SetTextColor(2);
	pt->AddText("f(x)=p_{0}/x");
	if(gStyle->GetOptFit()==0)
	{
		pt->AddText(Form("P0 = %g",f1_3.GetParameter(0)));
	}
	pt->Draw("same");

	////////////////////////////////////
	pCanF->cd(3);

	sprintf(hname,"h2DX_tgplaneVsPvb_%02d",holeindex);
	sprintf(target,"X_proj2tg_tr-Xtg_tr:Pvb");
	sprintf(htitle,"X_{proj2tg}^{tr} - X_{tg}^{tr} Vs P_{vb}; P_{vb} (GeV) ; X_{proj2tg}^{tr} - X_{tg}^{tr} (mm)");

	h2=DrawH2_RMS(hname,htitle,target,theCut,pHistoRange,"",0,0,3.0);
	
	//using the previous result in pad 2
	//GetValidBinIndex(h2,0.15,firstxbin,lastxbin);
	//xxmin=h2->GetXaxis()->GetBinCenter(firstxbin);
	//xxmax=h2->GetXaxis()->GetBinCenter(lastxbin);

	TF1 f1_2("f1_2","[0]/x+[1]",xxmin,xxmax);	
	f1_2.SetLineColor(6);
	
	h2=DrawH2(hname,htitle,target,theCut,pHistoRange,"prof",4,20);
	h2->Fit(&f1_2,"R","",xxmin,xxmax);

	if(gStyle->GetOptFit()!=0)
		pt = new TPaveText(0.5,0.42,0.9,0.55,"brNDC");
	else 
		pt = new TPaveText(0.5,0.30,0.9,0.6,"brNDC");
	
	pt->SetBorderSize(0);
	pt->SetFillColor(0);
	pt->SetTextColor(2);
	pt->AddText("f(x)=p_{0}/x+p_{1}");
	if(gStyle->GetOptFit()==0)
	{
		pt->AddText(Form("P0 = %g",f1_2.GetParameter(0)));
		pt->AddText(Form("P1 = %g",f1_2.GetParameter(1))); 
	}
	pt->Draw("same");

	////////////////////////////////////
	pCanF->cd(4);

	sprintf(hname,"h2DPVsP0_%02d",holeindex);
	sprintf(target,"(P0-P_rec)*1000:P0");
	sprintf(htitle,"P_{0}-P{rec} Vs P_{0}; P_{0} (GeV) ; P_{0}-P{rec} (MeV)");

	h2=DrawH2_RMS(hname,htitle,target,theCut,pHistoRange,"prof",4,20,3.0);

	xxmin=h2->GetXaxis()->GetBinCenter(10);
	xxmax=h2->GetXaxis()->GetBinCenter(h2->GetNbinsX()-10);
	

	////////////////////////////////////
	//pCanF->cd(5);

	for(int i=0; i<4;i++)
	{
		pCanF->cd(i+5);

		sprintf(htitle,"%s_{tg}^{tr}(black), %s_{proj2tg}^{tr}(blue) %s",
			strKeyTitle[i],strKeyTitle[i],strUnit[i]);

		//this part to get the range
		sprintf(hname,"h1%s_proj2tg_tr",strKey[i]);
		if(i<2) sprintf(target,"%s_proj2tg_tr",strKey[i]);
		else sprintf(target,"%svb_tr*1000",strKey[i]);
		if(i==2) h1=DrawH1_RMS(hname,htitle,target,theCut,pHistoRange,"",4,7.0);
		else h1=DrawH1_RMS(hname,htitle,target,theCut,pHistoRange,"",4,3.0);

		sprintf(hname,"h1%stg_tr",strKey[i]);
		if(i<2) sprintf(target,"%stg_tr",strKey[i]);
		else sprintf(target,"%s0_tr*1000",strKey[i]);
		h12=DrawH1(hname,htitle,target,theCut,pHistoRange,"same",1);
		if(h12->GetMaximum() > 1.2*h1->GetMaximum())
		{
			h12->Draw();
			h1->Draw("same");
		}
	} 

	////////////////////////////////////
if(true)
{
	pCanF->cd(9);
	
	for(int i=0; i<4;i++)
	{
		pCanF->cd(i+9);

		sprintf(htitle,"%s_{proj2tg}^{tr} - %s_{tg}^{tr} %s",
			strKeyTitle[i],strKeyTitle[i],strUnit[i]);

		//this part to get the range
		sprintf(hname,"h1D%s_tg_tr_%02d",strKey[i],holeindex);
		if(i<2)
		{
			sprintf(target,"%s_proj2tg_tr-%stg_tr",strKey[i],strKey[i]);
		}
		else
		{
			sprintf(target,"(%svb_tr-%s0_tr)*1000",strKey[i],strKey[i]);
		}
		h1=DrawH1_RMS(hname,htitle,target,theCut,pHistoRange,"",1,3.0);
		pRMS=h1->GetRMS();
		if(!(i<2 && pRMS<0.01)) FitGaus(h1);
	} 
	

	////////////////////////////////////

	////////////////////////////////////
	//Just plot the target plate
	////////////////////////////////////
	
	//pCanF->cd(13);

	for(int i=0; i<4;i++)
	{
		pCanF->cd(i+13);

		sprintf(htitle,"%s_{proj2tg}^{tr} - %s_{drift2tg}^{tr} %s",
			strKeyTitle[i],strKeyTitle[i],strUnit[i]);

		//this part to get the range
		sprintf(hname,"h1D%s_tgplane_tr_%02d",strKey[i],holeindex);
		if(i<2)
		{
			sprintf(target,"%s_proj2tg_tr-%stg_rec_tr",strKey[i],strKey[i]);
		}
		else
		{
			sprintf(target,"(%svb_tr-%stg_rec_tr)*1000",strKey[i],strKey[i]);
		}
		h1=DrawH1_RMS(hname,htitle,target,theCut,pHistoRange,"",4,3.0);
		pRMS=h1->GetRMS();
		if(!(i<2 && pRMS<0.01)) FitGaus(h1);
	} 
}
	////////////////////////////////////

	pCanF->cd();
	system("mkdir -p graph");
	//cout<<"gKey="<<gKey<<endl;
	if(!bIsCombinedTree)
	{  
		pCanF->SaveAs(Form("graph/FieldEffect_SieveHole%02d_E%.3f_Helm%.0fdeg_R%.1f_A%.0f_%s.png",holeindex,
			Beam,HelmRotAngle1*180/3.141593,HelmCurrentRatio,TargetAtomicNumber,gKey));	
	} 
	else
	{ 
		pCanF->SaveAs(Form("graph/FieldEffect_SieveHole%02d_E%.3f_Helm%.0fdeg_R%.1f_combined_%s.png",holeindex,
			Beam,HelmRotAngle1*180/3.141593,HelmCurrentRatio,gKey));
	}
}


void CheckRec(char *filename="",const char *key="")
{
	// using signal/slot in TCanvas/TPad to get feedback about processed events. 

	SetThisStyle();
	if(strlen(filename)>5)
	{
		cout<<"trying to open root file "<<filename<<endl;
		TFile::Open(filename);
	}

	if(gROOT->GetListOfFiles()->GetEntries()<1) {
		cout<<"no root file opened yet, quit ... \n"; 
		return;
	}

	if(!pCan) pCan=new TCanvas("pCan","Sieve Slit",630,800);
	pCan->cd();
	ReadConfig();

	TTree* track0=(TTree*)gROOT->FindObject("track0");

	TString file=gROOT->GetFile()->GetName();
	if(strlen(key)<2)
	{
		if(file.Last('/')>=0) file.Remove(0,file.Last('/')+1);
		file.Remove(file.Length()-5,5);
		key=file.Data();
	}
	sprintf(gKey,"%s",key);

	pCan->Clear();
	pCan->Divide(2,2,0.001,0.001);

	TH2F *h2VH_VB=0,*h2TP_VB=0;

	pCan->cd(1);
	gPad->SetRightMargin(0.15);
	if(gROOT->FindObject("h2VH_VB")) delete gROOT->FindObject("h2VH_VB");
	h2VH_VB=new TH2F("h2VH_VB", Form("%s; -Yvb_tr (mm); -Xvb_tr (mm)",key),
		50,-24.0,13.5,42,-43.6,49.5);
	track0->Draw("-Xvb_tr:-Yvb_tr >> h2VH_VB","TrackClass>5","contz"); 

	pCan->cd(2);
	gPad->SetRightMargin(0.15);
	if(gROOT->FindObject("h2TP_VB")) delete gROOT->FindObject("h2TP_VB");
	h2TP_VB=new TH2F("h2TP_VB","Thetavb_tr:Phivb_tr; #phi_{vb}^{tr} (rad); #theta_{vb}^{tr} (rad)",
		30,-0.03,0.03,60,-0.06,0.06);
	track0->Draw("Thetavb_tr:Phivb_tr >> h2TP_VB","TrackClass>5","contz"); 

	pCan->cd(3);
	gPad->SetRightMargin(0.15);
	//TH2F *h2VH_Tg=new TH2F("h2VH_Tg","-X0_tr:-Y0_tr ; -Y_{0}^{tr} (mm); -X_{0}^{tr} (mm)",
	//	  120,-12.0,12.0,120,-10.0,12.0);
	track0->Draw("-X0_tr:-Y0_tr >> h2VH_Tg","TrackClass>5","contz"); 

	pCan->cd(4);
	gPad->SetRightMargin(0.15);
	//TH2F *h2TP_Tg=new TH2F("h2TP_Tg","Theta0_tr:Phi0_tr ; #theta_{0}^{tr} (rad); #phi_{0}^{tr} (rad)",
	//	  30,-0.03,0.03,700,-0.25,0.1);
	track0->Draw("Theta0_tr:Phi0_tr >> h2TP_Tg","TrackClass>5","contz"); 


	pCan->ToggleToolBar();
	pCan->cd(1);
	pCan->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,
		"DoEvent(Int_t,Int_t,Int_t,TObject*)");
}


void ViewRecAll(char *filename="",const char *key="")
{
	// using signal/slot in TCanvas/TPad to get feedback about processed events. 

	SetThisStyle();
	if(strlen(filename)>5)
	{
		cout<<"trying to open root file "<<filename<<endl;
		TFile::Open(filename);
	}

	if(gROOT->GetListOfFiles()->GetEntries()<1) {
		cout<<"no root file opened yet, quit ... \n"; 
		return;
	}

	ReadConfig();

	TString file=gROOT->GetFile()->GetName();
	if(strlen(key)<2)
	{
		if(file.Last('/')>=0) file.Remove(0,file.Last('/')+1);
		file.Remove(file.Length()-5,5);
		key=file.Data();
	}
	sprintf(gKey,"%s",key);

	
	PlotRec();
}


void FieldEffect(char* extraCut="1>0")
{
	TCut all="TrackClass>5 && abs(Xfp_tr)<800 && abs(Yvb)<100";
	TCut theCut=all;

	if(gROOT->FindObject("CUTGVBXY")) theCut+="CUTGVBXY";	
	if(gROOT->FindObject("CUTG")) theCut+="CUTG";
	theCut+=extraCut;
		
	int holeindex=0;
	DoFieldEffect(theCut,holeindex);
}
