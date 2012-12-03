#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TF1.h"
#include "TPad.h"
#include "TStyle.h"

#include "TText.h"
#include "TPaveText.h"

#include "HRSRecUseDB.hh"

#include "HRSTransport.hh"
#include "HRSTransform_TCSNHCS.hh"

using namespace std;

//flat random number generator between [0,1)
double fRand();
double fRandGaus(double m=0.0, double s=1.0);
//return random number in [low,High) following a*x+c prob density
double fLinearRand(double a=1.0,double c=0.0,double low=0.0,double high=1.0);


typedef double (*func)(double);

double fPol1(double x)
{
	double a=1.0,b=0.0;
	return a*x+b;
}

double fPol2(double x)
{
	double a=1.0,b=0.0,c=0.0;
	return a*x*x+b*x+c;
}

//generate a non-uniform x using the constrain(distribution) of f(x)
double RandomOfFunction(func f, double low,double high)
{
	//input:
	//	low,high	  lower,higher limit of x
	double x,y;
	double ylow=f(low),yhigh=f(high);
	if(ylow<yhigh) 
	{
		double tmp=ylow;
		ylow=yhigh;
		yhigh=tmp;
	}

	do{
		x=low+(high-low)*((double)rand())/(double(RAND_MAX));
		y=ylow+(yhigh-ylow)*((double)rand())/(double(RAND_MAX));

		if(f(x)>y) break;
	}while(true);
	return x;
}


bool SNAKEThruHRS(int pIsLeftArm, double pEndPlaneAngle, double pXtg_BPM_tr, 
				  double pHRSMomentum, int iFieldRotation, int iExperiment,
				  double* pV5tg_tr, double* pV5fp_tr);

void SetStyle()
{
	/*
	gStyle->SetPalette(1);
	gStyle->SetFillColor(0);

	double pPadRightMargin=0.05;
	gStyle->SetPadRightMargin(pPadRightMargin); 
	gStyle->SetStatX(1.0-gStyle->GetPadRightMargin());
	gStyle->SetStatY(1.0-gStyle->GetPadTopMargin());

	gStyle->SetOptStat(0);	
	gStyle->SetOptFit(011);
	gStyle->SetStatH(0.07);
	gStyle->SetStatW(0.25);
	gStyle->SetStatStyle(4000);

	gStyle->SetTitleH(0.085);
	gStyle->SetTitleBorderSize(0);
	gStyle->SetTitleX(0.15);
	*/


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


TF1* DoFitGaus(TH1* h1, double range_in_sigma=1.0)
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


TF1* DoFitGaus(char* hName, double range_in_sigma=1.0)
{
	TH1* h1=(TH1 *)gROOT->FindObject(hName);
	return DoFitGaus(h1,range_in_sigma);
}


TH1* DrawTH1_RMS(TTree *tree, char *hname, char *title, char *target, TCut &cut, 
				 char* hrange, const char *option="", int color=1, double RMS_factor=3.0)
{
	char hdest[512];
	TH1 *h1=0;
	//if the option contains 'same' or 'SAME', do not search the range
	if(!(strstr(option,"same") || strstr(option,"SAME")))  
	{
		double pMean,pRMS;
		TH1 *h1 = (TH1*) (gROOT->FindObject(hname));
		if(h1)  {delete h1;}
		sprintf(hdest,"%s >> %s",target,hname);
		tree->Draw(hdest,cut,"");
		h1 = (TH1*) (gROOT->FindObject(hname));
		pMean=h1->GetMean();
		pRMS=h1->GetRMS();
		sprintf(hrange,"(60,%.3f,%.3f)",pMean-RMS_factor*pRMS,pMean+RMS_factor*pRMS);
	}
	//cout<<"DrawTH1(): "<<tree->GetName()<<"->Draw(\""<<target<<" >> "<<hname<<hrange<<"\");"
	//	<<" title=\""<<title<<"\"\n";

	sprintf(hdest,"%s >> %s %s",target,hname,hrange);
	h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	tree->Draw(hdest,cut,option);
	h1 = (TH1*) (gROOT->FindObject(hname));
	h1->SetTitle(title);
	h1->SetLineColor(color);
	if(!(strstr(option,"same") || strstr(option,"SAME")))
		h1->GetYaxis()->SetRangeUser(h1->GetMinimum(),1.2*h1->GetMaximum());  
	return h1;
}


TH1* DrawTH1(TTree *tree, char *hname, char *title, char *target, TCut &cut, 
			 char* hrange,const char *option="", int color=4)
{
	char hdest[512];
	sprintf(hdest,"%s >> %s %s",target,hname,hrange);
	TH1* h1 = (TH1*) (gROOT->FindObject(hname));
	if(h1)  {delete h1;}
	tree->Draw(hdest,cut,option);
	h1 = (TH1*) (gROOT->FindObject(hname));
	h1->SetTitle(title);
	h1->SetLineColor(color);
	if(!(strstr(option,"same") || strstr(option,"SAME")))   
		h1->GetYaxis()->SetRangeUser(h1->GetMinimum(),1.2*h1->GetMaximum());  
	return h1;
}


void PlotDelta(TTree *s,char *grname="Delta.png")
{
	TCanvas *pCan=new TCanvas("pCan","",1350,850);
	pCan->Divide(5,5,0.001,0.001);
	TCut pCut="";
	TPad *pPad=0;
	TH1 *h1,*h12; h1=0;h12=0;
	TH2 *h2,*h22; h2=0;h22=0;

	char hName[255],strRange[255],strTitle[255],strTg[255];
	const char *strKey[5]={"X","Theta","Y","Phi","Delta"};
	const char *strKeyTitle[5]={"X","#theta","Y","#phi","#delta"};
	const char *strUnit[5]={"(mm)","(mrad)","(mm)","(mrad)","(x 10^{-4})"};

	int i=0;
	for(i=0; i<5;i++)
	{
		pPad=(TPad*)(pCan->cd(i+1));
		if(i==4)
		{
			sprintf(hName,"h1d%s",strKey[i]);
			sprintf(strTg,"(%srec-%s)*10000",strKey[i],strKey[i]);
			sprintf(strTitle,"%s_{rec}-%s_{tg} %s",strKeyTitle[i],strKeyTitle[i],strUnit[i]);
		}
		else
		{
			sprintf(hName,"h1d%s_tr",strKey[i]);
			sprintf(strTg,"(%srec_tr-%stg_tr)*1000",strKey[i],strKey[i]);
			sprintf(strTitle,"%s_{rec}^{tr}-%s_{tg}^{tr} %s",strKeyTitle[i],strKeyTitle[i],strUnit[i]);
		}

		h1=DrawTH1_RMS(s,hName,strTitle,strTg,pCut,strRange,"",1,3.0);
		DoFitGaus(h1,1.0);
	} 

	for(i=0; i<4;i++)
	{
		pPad=(TPad*)(pCan->cd(i+6));
		sprintf(hName,"h1d%s_lab",strKey[i]);
		sprintf(strTg,"(%srec_lab-%stg_lab)*1000",strKey[i],strKey[i]);
		sprintf(strTitle,"%s_{rec}^{lab}-%s_{tg}^{lab} %s",
			strKeyTitle[i],strKeyTitle[i],strUnit[i]);

		h1=DrawTH1_RMS(s,hName,strTitle,strTg,pCut,strRange,"",1,3.0);
		DoFitGaus(h1,1.0);
	} 
	i=4;
	pCan->cd(i+6);
	sprintf(hName,"h1dX0");
	sprintf(strTg,"(Xtg_BPM_tr-Xtg_tr)*1000");
	sprintf(strTitle,"X_{tg}^{tr}_BPM-X_{tg}^{tr} (mm)");
	h1=DrawTH1_RMS(s,hName,strTitle,strTg,pCut,strRange,"",1,3.0);
	DoFitGaus(h1,1.0);


	for(i=0; i<4;i++)
	{
		pPad=(TPad*)(pCan->cd(i+11));

		sprintf(strTitle,"%s_{tg}^{tr}(black), %s_{rec}^{tr}(blue) %s",
			strKeyTitle[i],strKeyTitle[i],strUnit[i]);

		//this part to get the range
		sprintf(hName,"h1%srec_tr",strKey[i]);
		sprintf(strTg,"%srec_tr*1000",strKey[i]);
		h1=DrawTH1_RMS(s,hName,strTitle,strTg,pCut,strRange,"",4,3.0);

		sprintf(hName,"h1%stg_tr",strKey[i]);
		sprintf(strTg,"%stg_tr*1000",strKey[i]);
		h12=DrawTH1(s,hName,strTitle,strTg,pCut,strRange,"same",1);
		if(h12->GetMaximum() > 1.2*h1->GetMaximum())
		{
			h12->Draw();
			h1->Draw("same");
		}
	} 

	i=4;
	pPad=(TPad*)(pCan->cd(i+11));
	sprintf(strTitle,"%s(black), %s_{rec}(blue) %s",
		strKeyTitle[i],strKeyTitle[i],strUnit[i]);
	//this part will get the range
	sprintf(hName,"h1%srec",strKey[i]);
	sprintf(strTg,"%srec*10000",strKey[i]);
	h1=DrawTH1_RMS(s,hName,strTitle,strTg,pCut,strRange,"",4,3.0);

	sprintf(hName,"h1%s",strKey[i]);
	sprintf(strTg,"%s*10000",strKey[i]);
	h12=DrawTH1(s,hName,strTitle,strTg,pCut,strRange,"same",1);	
	if(h12->GetMaximum() > 1.2*h1->GetMaximum())
	{
		h12->Draw();
		h1->Draw("same");
	}

	for(i=0; i<4;i++)
	{
		pPad=(TPad*)(pCan->cd(i+16));

		sprintf(strTitle,"%s_{tg}^{lab}(black), %s_{rec}^{lab}(blue) %s",
			strKeyTitle[i],strKeyTitle[i],strUnit[i]);

		//this part will get the range
		sprintf(hName,"h1%srec_lab",strKey[i]);
		sprintf(strTg,"%srec_lab*1000",strKey[i]);
		h1=DrawTH1_RMS(s,hName,strTitle,strTg,pCut,strRange,"",4,3.0);

		sprintf(hName,"h1%s_lab",strKey[i]);
		sprintf(strTg,"%stg_lab*1000",strKey[i]);
		h12=DrawTH1(s,hName,strTitle,strTg,pCut,strRange,"same",1);
		if(h12->GetMaximum() > 1.2*h1->GetMaximum())
		{
			h12->Draw();
			h1->Draw("same");
		}
	} 

	
	for(i=0; i<4;i++)
	{
		pPad=(TPad*)(pCan->cd(i+21));

		sprintf(strTitle,"%s_{tg}^{tr}(black), %s_{rec_db}^{tr}(blue) %s",
			strKeyTitle[i],strKeyTitle[i],strUnit[i]);

		//this part to get the range
		sprintf(hName,"h1%srec_db_tr",strKey[i]);
		sprintf(strTg,"%srec_db_tr*1000",strKey[i]);
		h1=DrawTH1_RMS(s,hName,strTitle,strTg,pCut,strRange,"",4,3.0);

		sprintf(hName,"h1%stg_tr_",strKey[i]);
		sprintf(strTg,"%stg_tr*1000",strKey[i]);
		h12=DrawTH1(s,hName,strTitle,strTg,pCut,strRange,"same",1);
		if(h12->GetMaximum() > 1.2*h1->GetMaximum())
		{
			h12->Draw();
			h1->Draw("same");
		}
	} 

	cout<<"grname="<<grname<<endl;
	pCan->SaveAs(grname);

}

//dependence on BPM resolution for x0_tr,y0_tr,theta0_tr,phi0_tr,delta
void _old_DrawBPMDependence(char *treename="s",char *grname="BPMDependence.png")
{
	bool DoFitPol1=true;
	int OldOptFitStyle=gStyle->GetOptFit();
	double OldStatY=gStyle->GetStatY();
	if(DoFitPol1)
	{
		gStyle->SetOptFit(0010);
		if(gStyle->GetOptFit()==8 || gStyle->GetOptFit()==9) gStyle->SetStatH(0.08);
		else gStyle->SetStatH(0.16);
		gStyle->SetStatY(0.4);
	}

	TTree *s=(TTree*)gROOT->FindObject(treename); 
	TCut theCut="";
	system("mkdir -p graph");

	const int nitemx=6;
	const int nitemy=5;
	char **ylist,**ynamelist;
	char *ylist_lab[nitemy]={"(Xrec_lab-Xtg_lab)*1000","(Yrec_lab-Ytg_lab)*1000",
		"(Thetarec_lab-Thetatg_lab)*1000","(Phirec_lab-Phitg_lab)*1000",
		"(Deltarec-Delta)*10000"};
	char *ynamelist_lab[nitemy]={"dX_lab","dY_lab","dTheta_lab","dPhi_lab","dDelta"};

	char *ylist_tr[nitemy]={"(Xrec_tr-Xtg_tr)*1000","(Yrec_tr-Ytg_tr)*1000",
		"(Thetarec_tr-Thetatg_tr)*1000","(Phirec_tr-Phitg_tr)*1000",
		"(Deltarec-Delta)*10000"};
	char *ynamelist_tr[nitemy]={"dX_tr","dY_tr","dTheta_tr","dPhi_tr","dDelta"};

	char *yunit[nitemy]={"(mm)","(mm)","(mrad)","(mrad)","(x 10^{-4})"};

	char *xlist[nitemx]={"BPMRes*1000","Xtg_tr*1000","Ytg_tr*1000","Thetatg_tr*1000",
		"Phitg_tr*1000","Delta*10000"};
	char *xnamelist[nitemx]={"BPMRes","Xtg_tr","Ytg_tr","Thetatg_tr","Phitg_tr","Delta"};
	char *xunit[nitemx]={"(mm)","(mm)","(mm)","(mrad)","(mrad)","(x 10^{-4})"};

	int iDrawContent=3;  //can be 1(sigma only), 2(mean+sigma),3(2-D histo+mean+sigma)
	
	int ncol=iDrawContent;
	int nrow=nitemy;
	if(iDrawContent<3 || nrow>5)
	{
		if((nrow%3)==0) 
		{
			ncol*=3;
			nrow=nrow/3;
		}
		else if((nrow%2)==0) 
		{
			ncol*=2;
			nrow=nrow/2;
		}
	}
	TCanvas *pCanR=new TCanvas("pCanR","Resolution",320*ncol,240*nrow);
	pCanR->Divide(ncol,nrow,0.001,0.001);

	char hname[255],hdest[500],htitle[255];

	TH1 *h1mean,*h1sigma;
	TH2 *h2=0;

	for(int k=0;k<2;k++)
	{
		if(k==0)
		{
			ylist=ylist_lab;
			ynamelist=ynamelist_lab;
		}
		if(k==1)
		{
			ylist=ylist_tr;
			ynamelist=ynamelist_tr;
		}

		for (int i=0;i<nitemx;i++)
		{
			for (int j=0;j<nitemy;j++)
			{
				pCanR->cd(j*iDrawContent+1);
				sprintf(hname,"h2%sVs%s",ynamelist[j],xnamelist[i]);
				sprintf(htitle,"%s Vs %s; %s %s; %s %s",ynamelist[j],xnamelist[i],
					xnamelist[i],xunit[i],ynamelist[j],yunit[j]);
				sprintf(hdest,"%s:%s >> %s",ylist[j],xlist[i],hname);
				h2=(TH2*) gROOT->FindObject(hname);
				if(h2) delete h2;

				//cout<<"Draw 2-D histo "<<hname<<"(\""<<ylist[i]<<":"<<xlist[i]<<"\")"
				//	<<" and title=\""<<htitle<<"\"\n"; 

				s->Draw(hdest,theCut,"");	
				h2=(TH2*) gROOT->FindObject(hname);
				h2->SetTitle(htitle);

				//search for start bin index and end bin index
				TH1* h1projx=(TH1*) (h2->ProjectionX(Form("%s_projx",hname)));
				int firstxbin=0;
				int lastxbin=h2->GetNbinsX(); 
				double maxentries=h1projx->GetMaximum(); 
				for(int ib=0;ib<lastxbin;ib++)
				{
					if(h1projx->GetBinContent(ib)>=0.1*maxentries) 
					{
						firstxbin=ib;
						break;
					}
				}	
				for(int ib=lastxbin;ib>firstxbin;ib--)
				{
					if(h1projx->GetBinContent(ib)>=0.1*maxentries) 
					{
						lastxbin=ib;
						break;
					}
				}
				delete h1projx;

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
					pCanR->cd(j*1+1); h1sigma->Draw();	
				}
				else if(iDrawContent==2)
				{
					pCanR->cd(j*2+1); h1mean->Draw();
					pCanR->cd(j*2+2); h1sigma->Draw();	
				}
				else if(iDrawContent==3) 
				{
					pCanR->cd(j*3+1); gPad->SetRightMargin(0.15);h2->Draw("contz");
					pCanR->cd(j*3+2); h1mean->Draw();
					pCanR->cd(j*3+3); h1sigma->Draw();	
				}
			}
			pCanR->cd();

			cout<<"grname="<<Form("%s_%s_%s",xnamelist[i],((k==0)?"lab":"tr"),grname)<<endl;
			pCanR->SaveAs(Form("%s_%s_%s",xnamelist[i],((k==0)?"lab":"tr"),grname));	
		}
	}

	if(DoFitPol1)
	{
		gStyle->SetOptFit(OldOptFitStyle);
		gStyle->SetStatY(OldStatY);
	}
}

//This routine will polt YList's dependence on each Xlist item
//
//plot dependence on XX for x0_tr,y0_tr,theta0_tr,phi0_tr,delta (plot1) and 
//x0_lab,y0_lab,theta0_lab,phi0_lab,delta (plot2)
//int iDrawContent=3;  //can be 1(sigma only), 2(mean+sigma),3(2-D histo+mean+sigma)
//One can provide XX as char** arrays (No more than 6 XX, please)
//if _nitemx==0, will use all 6 XX listed in default char** arrays, which are
//"BPMRes","Xtg_tr","Ytg_tr","Thetatg_tr","Phitg_tr","Delta"
//if _nitemx==[-1,-6], will plot for just one item spesified by abs(_nitemx)
//if _nitemx>0, one should provide these arrays:
//char **_xlist, char **_xnamelist, char **_xunit,
//if _nitemy==0, will use all 5 YY listed in default arrays, 
//if _nitemy>0, will only plot with provided arrays, no more than 10 YY items
void DrawDependenceCore(const char *treename="s",const char *cut="",
						const char *_grname="BPMDependence.png", int iDrawContent=3,
						int _nitemx=0, char **_xlist=0, char **_xnamelist=0,char **_xunit=0,
						int _nitemy=0, char **_ylist=0, char **_ynamelist=0,char **_yunit=0)
{	
	bool DoFitPol1=true;
	int OldOptFitStyle=gStyle->GetOptFit();
	double OldStatY=gStyle->GetStatY();
	if(DoFitPol1)
	{
		gStyle->SetOptFit(0010);
		if(gStyle->GetOptFit()==8 || gStyle->GetOptFit()==9) gStyle->SetStatH(0.08);
		else gStyle->SetStatH(0.16);
		gStyle->SetStatY(0.4);
	}

	//int iDrawContent=3;  //can be 1(sigma only), 2(mean+sigma),3(2-D histo+mean+sigma)
	TTree *s=(TTree*)gROOT->FindObject(treename); 
	TCut theCut=cut;
	char grname[512];

	const int nitemy=10;
	char *ylist[nitemy]={"(Xrec_lab-Xtg_lab)*1000","(Yrec_lab-Ytg_lab)*1000",
		"(Thetarec_lab-Thetatg_lab)*1000","(Phirec_lab-Phitg_lab)*1000",
		"(Deltarec-Delta)*10000",
		"(Xrec_tr-Xtg_tr)*1000","(Yrec_tr-Ytg_tr)*1000",
		"(Thetarec_tr-Thetatg_tr)*1000","(Phirec_tr-Phitg_tr)*1000",
		"(Deltarec-Delta)*10000"};
	char *ynamelist[nitemy]={"dX_lab","dY_lab","dTheta_lab","dPhi_lab","dDelta",
		"dX_tr","dY_tr","dTheta_tr","dPhi_tr","dDelta"};

	char *yunit[nitemy]={"(mm)","(mm)","(mrad)","(mrad)","(x 10^{-4})",
		"(mm)","(mm)","(mrad)","(mrad)","(x 10^{-4})"};

	const int nitemx=6;
	char *xlist[nitemx]={"BPMRes*1000","Xtg_tr*1000","Ytg_tr*1000","Thetatg_tr*1000",
		"Phitg_tr*1000","Delta*10000"};
	char *xnamelist[nitemx]={"BPMRes","Xtg_tr","Ytg_tr","Thetatg_tr","Phitg_tr","Delta"};
	char *xunit[nitemx]={"(mm)","(mm)","(mm)","(mrad)","(mrad)","(x 10^{-4})"};

	int nitemy_start=0, nitemy_end=8;
	int nitemx_start=0, nitemx_end=nitemx-1;
	if(_nitemx>0) 
	{
		if(_nitemx>6)
		{
			cout<<"Warning: DrawDependenceCore() can plot no more tham 6 plots at a time."<<endl;
			_nitemx=6;
		}
		nitemx_end=_nitemx-1;
		for(int i=0;i<_nitemx;i++)
		{
			xlist[i]=_xlist[i];
			xnamelist[i]=_xnamelist[i];
			xunit[i]=_xunit[i];
		}
	}
	else if(_nitemx<=-1 && _nitemx>=-6)
	{
		nitemx_start=nitemx_end=-_nitemx;
	}

	if(_nitemy>0) 
	{
		if(_nitemy>10)
		{
			cout<<"Warning: DrawDependenceCore() can plot no more tham 10 plots at a time."<<endl;
			_nitemx=10;
		}
		nitemy_end=_nitemy-1;
		for(int i=0;i<_nitemy;i++)
		{
			ylist[i]=_ylist[i];
			ynamelist[i]=_ynamelist[i];
			yunit[i]=_yunit[i];
		}
	}

	int ncol=iDrawContent;
	int nrow=nitemy_end-nitemy_start+1;
	if(iDrawContent<3 || nrow>5)
	{
		if((nrow%3)==0) 
		{
			ncol*=3;
			nrow=nrow/3;
		}
		else if((nrow%2)==0) 
		{
			ncol*=2;
			nrow=nrow/2;
		}
	}
	TCanvas *pCanR=new TCanvas("pCanR","Resolution",320*ncol,240*nrow);
	pCanR->Divide(ncol,nrow,0.001,0.001);

	char hname[255],hdest[500],htitle[255];

	TH1 *h1mean,*h1sigma;
	TH2 *h2=0;
	for (int i=nitemx_start;i<=nitemx_end;i++)
	{
		for (int j=nitemy_start;j<=nitemy_end;j++)
		{
			pCanR->cd(j*iDrawContent+1);
			sprintf(hname,"h2%sVs%s",ynamelist[j],xnamelist[i]);
			sprintf(htitle,"%s Vs %s; %s %s; %s %s",ynamelist[j],xnamelist[i],
				xnamelist[i],xunit[i],ynamelist[j],yunit[j]);
			sprintf(hdest,"%s:%s >> %s",ylist[j],xlist[i],hname);
			h2=(TH2*) gROOT->FindObject(hname);
			if(h2) delete h2;

			//cout<<"Draw 2-D histo "<<hname<<"(\""<<ylist[i]<<":"<<xlist[i]<<"\")"
			//	<<" and title=\""<<htitle<<"\"\n"; 

			s->Draw(hdest,theCut,"");	
			h2=(TH2*) gROOT->FindObject(hname);
			h2->SetTitle(htitle);

			//search for start bin index and end bin index
			TH1* h1projx=(TH1*) (h2->ProjectionX(Form("%s_projx",hname)));
			int firstxbin=0;
			int lastxbin=h2->GetNbinsX(); 
			double maxentries=h1projx->GetMaximum(); 
			for(int ib=0;ib<lastxbin;ib++)
			{
				if(h1projx->GetBinContent(ib)>=0.15*maxentries) 
				{
					firstxbin=ib;
					break;
				}
			}	
			for(int ib=lastxbin;ib>firstxbin;ib--)
			{
				if(h1projx->GetBinContent(ib)>=0.15*maxentries) 
				{
					lastxbin=ib;
					break;
				}
			}
			delete h1projx;

			h2->FitSlicesY(0,firstxbin,lastxbin);
			h1mean  = (TH1*)gDirectory->Get(Form("%s_1",hname));
			h1mean->SetTitle(Form("%s Mean Vs %s; %s %s; %s Mean %s",
				ynamelist[j],xnamelist[i],xnamelist[i],xunit[i],ynamelist[j],yunit[j]));
			h1mean->SetMarkerStyle(20); 

			h1sigma = (TH1*)gDirectory->Get(Form("%s_2",hname));
			h1sigma->SetMarkerStyle(20); 
			//h1sigma->SetTitle(Form("%s Resolution Vs %s; %s %s; %s Resolution %s",
			h1sigma->SetTitle(Form("%s Res. Vs %s; %s %s; %s Res. %s",
				ynamelist[j],xnamelist[i],xnamelist[i],xunit[i],ynamelist[j],yunit[j]));

			if(DoFitPol1) h1sigma->Fit("pol1","Q");

			if(iDrawContent==1) 
			{
				pCanR->cd(j*1+1); h1sigma->Draw();	
			}
			else if(iDrawContent==2)
			{
				pCanR->cd(j*2+1); h1mean->Draw();
				pCanR->cd(j*2+2); h1sigma->Draw();	
			}
			else if(iDrawContent==3) 
			{
				pCanR->cd(j*3+1); gPad->SetRightMargin(0.15);h2->Draw("contz");
				pCanR->cd(j*3+2); h1mean->Draw();
				pCanR->cd(j*3+3); h1sigma->Draw();	
			}
		}
		pCanR->cd();

		sprintf(grname,"%s_%s",xnamelist[i],_grname);
		cout<<"grname="<<grname<<endl;
		pCanR->SaveAs(grname);	
	}

	if(DoFitPol1)
	{
		gStyle->SetOptFit(OldOptFitStyle);
		gStyle->SetStatY(OldStatY);
	}
}

//dependence on BPM resolution for x0_tr,y0_tr,theta0_tr,phi0_tr,delta
void DrawDependence(TTree *s,char *_grname="BPMDependence.png",int plotBPM=0)
{

	const int nitemy=5;

	char *ylist_lab[nitemy]={"(Xrec_lab-Xtg_lab)*1000","(Yrec_lab-Ytg_lab)*1000",
		"(Thetarec_lab-Thetatg_lab)*1000","(Phirec_lab-Phitg_lab)*1000",
		"(Deltarec-Delta)*10000"};
	char *ynamelist_lab[nitemy]={"dX_lab","dY_lab","dTheta_lab","dPhi_lab","dDelta"};

	char *ylist_tr[nitemy]={"(Xrec_tr-Xtg_tr)*1000","(Yrec_tr-Ytg_tr)*1000",
		"(Thetarec_tr-Thetatg_tr)*1000","(Phirec_tr-Phitg_tr)*1000",
		"(Deltarec-Delta)*10000"};
	char *ynamelist_tr[nitemy]={"dX_tr","dY_tr","dTheta_tr","dPhi_tr","dDelta"};

	char *yunit[nitemy]={"(mm)","(mm)","(mrad)","(mrad)","(x 10^{-4})"};

	const int nitemx=5;
	char *xlist[nitemx]={"Xtg_tr*1000","Ytg_tr*1000","Thetatg_tr*1000",
		"Phitg_tr*1000","Delta*10000"};
	char *xnamelist[nitemx]={"Xtg_tr","Ytg_tr","Thetatg_tr","Phitg_tr","Delta"};
	char *xunit[nitemx]={"(mm)","(mm)","(mrad)","(mrad)","(x 10^{-4})"};

	char *xlist_BPM[]={"BPMRes*1000"};
	char *xnamelist_BPM[]={"BPMRes"};
	char *xunit_BPM[]={"(mm)"};

	char grname[512];

	int iDrawConternt=3;
	//only plot for BPMRes
	if(plotBPM)
	{
		sprintf(grname,"lab_%s",_grname);
		//DrawDependenceCore(s->GetName(),"",grname,iDrawConternt,
		//	-1,0,0,0,nitemy,ylist_lab,ynamelist_lab,yunit); 
		DrawDependenceCore(s->GetName(),"",grname,iDrawConternt,
			1,xlist_BPM,xnamelist_BPM,xunit_BPM,nitemy,ylist_lab,ynamelist_lab,yunit); 

		sprintf(grname,"tr_%s",_grname);
		//DrawDependenceCore(s->GetName(),"",grname,iDrawConternt,
		//	-1,0,0,0,nitemy,ylist_tr,ynamelist_tr,yunit);
		DrawDependenceCore(s->GetName(),"",grname,iDrawConternt,
			1,xlist_BPM,xnamelist_BPM,xunit_BPM,nitemy,ylist_tr,ynamelist_tr,yunit);
	}
	else
	{
		sprintf(grname,"lab_%s",_grname);
		//DrawDependenceCore(s->GetName(),"",grname,iDrawConternt,
		//	0,0,0,0,nitemy,ylist_lab,ynamelist_lab,yunit); 
		DrawDependenceCore(s->GetName(),"",grname,iDrawConternt,
			5,xlist,xnamelist,xunit,nitemy,ylist_lab,ynamelist_lab,yunit); 

		sprintf(grname,"tr_%s",_grname);
		//DrawDependenceCore(s->GetName(),"",grname,iDrawConternt,
		//	0,0,0,0,nitemy,ylist_tr,ynamelist_tr,yunit);
		DrawDependenceCore(s->GetName(),"",grname,iDrawConternt,
			5,xlist,xnamelist,xunit,nitemy,ylist_tr,ynamelist_tr,yunit);
	}

	//default core: plot 9_YY vs XX, there are 5 XX in total
	DrawDependenceCore(s->GetName(),"",grname,1);
}


//input:
//iExperiment=10: g2p normal, 11: g2p 484816_shim; 12 g2p 403216_shim; 13 g2p 400016_shim 
//iSmearX0=0,1,2,3 means Xtg_BPM_tr=0, Xtg_BPM_tr, 1mm_gaussian_smeared_XBMP_tr, 0-2.5mm_gaussian_smeared_XBMP_tr
//iSourceDistr=0,1,2 means pV5tg_tr[5] is in delta, flat, gaussian distribution
//iArm = 0,1,2 means random, left, right
void TestSNAKE(int iNEvent, int iExperiment, int iSmearX0, int iSourceDistr, int iArm)
{	
	srand(time(0));
	SetStyle();
	

	const double deg=acos(0.0)/90;

	int    Index=0;
	int    iLeftArm=1;
	double pHRSAngle=5.69*deg;  //6deg
	double pXtg_BPM_tr=0;
	double pHRSMomentum=2.251;
	int    iFieldRotation=90;  //90 deg	
	double pV5fp_tr[5]={0,0,0,0,0};
	double pV5tg_tr[5]={0,0,0,0,0};
	double pV5rec_tr[5]={0,0,0,0,0};

	double pV5tg_lab[5]={0,0,0,0,0};
	double pV5rec_lab[5]={0,0,0,0,0};
	double pBPMRes=0.001;  //BPM resolution, will try 0 - 2.5 if iSmaeX0>=3 

	TFile* pFile= new TFile(Form("snake_Exp%02d_X%d_Distr%d.root",
		iExperiment,iSmearX0,iSourceDistr),"recreate");
	TTree* pSnake=new TTree("s","snake reconstruction");
	pSnake->Branch("Index",&Index,"Index/I");

	pSnake->Branch("BPMRes",&pBPMRes,"BPMRes/D");
	pSnake->Branch("Left",&iLeftArm,"Left/I");
	pSnake->Branch("FieldRotation",&iFieldRotation,"FieldRotation/I");
	pSnake->Branch("HRSAngle",&pHRSAngle,"HRSAngle/D");
	pSnake->Branch("P0",&pHRSMomentum,"P0/D");
	pSnake->Branch("Xtg_BPM_tr",&pXtg_BPM_tr,"Xtg_BPM_tr/D");

	pSnake->Branch("Xfp_tr",&pV5fp_tr[0],"Xfp_tr/D");
	pSnake->Branch("Thetafp_tr",&pV5fp_tr[1],"Thetafp_tr/D");
	pSnake->Branch("Yfp_tr",&pV5fp_tr[2],"Yfp_tr/D");
	pSnake->Branch("Phifp_tr",&pV5fp_tr[3],"Phifp_tr/D");

	pSnake->Branch("Xtg_tr",&pV5tg_tr[0],"Xtg_tr/D");
	pSnake->Branch("Thetatg_tr",&pV5tg_tr[1],"Thetatg_tr/D");
	pSnake->Branch("Ytg_tr",&pV5tg_tr[2],"Ytg_tr/D");
	pSnake->Branch("Phitg_tr",&pV5tg_tr[3],"Phitg_tr/D");
	pSnake->Branch("Delta",&pV5tg_tr[4],"Delta/D");

	pSnake->Branch("Xrec_tr",&pV5rec_tr[0],"Xrec_tr/D");
	pSnake->Branch("Thetarec_tr",&pV5rec_tr[1],"Thetarec_tr/D");
	pSnake->Branch("Yrec_tr",&pV5rec_tr[2],"Yrec_tr/D");
	pSnake->Branch("Phirec_tr",&pV5rec_tr[3],"Phirec_tr/D");
	pSnake->Branch("Deltarec",&pV5rec_tr[4],"Deltarec/D");

	pSnake->Branch("Xtg_lab",&pV5tg_lab[0],"Xtg_lab/D");
	pSnake->Branch("Thetatg_lab",&pV5tg_lab[1],"Thetatg_lab/D");
	pSnake->Branch("Ytg_lab",&pV5tg_lab[2],"Ytg_lab/D");
	pSnake->Branch("Phitg_lab",&pV5tg_lab[3],"Phitg_lab/D");
	pSnake->Branch("Ztg_lab",&pV5tg_lab[4],"Ztg_lab/D");

	pSnake->Branch("Xrec_lab",&pV5rec_lab[0],"Xrec_lab/D");
	pSnake->Branch("Thetarec_lab",&pV5rec_lab[1],"Thetarec_lab/D");
	pSnake->Branch("Yrec_lab",&pV5rec_lab[2],"Yrec_lab/D");
	pSnake->Branch("Phirec_lab",&pV5rec_lab[3],"Phirec_lab/D");
	pSnake->Branch("Zrec_lab",&pV5tg_lab[4],"Zrec_lab/D");


	//rec by DB
	double pV5rec_db_tr[5];
	pSnake->Branch("Xrec_db_tr",&pV5rec_db_tr[0],"Xrec_db_tr/D");
	pSnake->Branch("Thetarec_db_tr",&pV5rec_db_tr[1],"Thetarec_db_tr/D");
	pSnake->Branch("Yrec_db_tr",&pV5rec_db_tr[2],"Yrec_db_tr/D");
	pSnake->Branch("Phirec_db_tr",&pV5rec_db_tr[3],"Phirec_db_tr/D");
	pSnake->Branch("Deltarec_db",&pV5rec_db_tr[4],"Deltarec_db/D");

	HRSRecUseDB *pRecDBL=new HRSRecUseDB("L","db_L.vdc.dat");
	HRSRecUseDB *pRecDBR=new HRSRecUseDB("R","db_R.vdc.dat");
	HRSRecUseDB *pRecDB=0;

	//do smear: 0: no smear; 1: flat; 2 gaus 
	int NThrown=0;
	Index=0;
	while(Index<iNEvent)
	{
		if(iSourceDistr==1)
		{
			pV5tg_tr[0] = 0.010 * 2 * (fRand()-0.5);
			pV5tg_tr[1] = 0.070 * 2 * (fRand()-0.5);
			pV5tg_tr[2] = 0.010 * 2 * (fRand()-0.5);
			pV5tg_tr[3] = 0.035 * 2 * (fRand()-0.5);
			pV5tg_tr[4] = 0.050 * 2 * (fRand()-0.5);
		}
		else if(iSourceDistr>=2)
		{
			pV5tg_tr[0] = fRandGaus(0,0.00010/2);
			pV5tg_tr[1] = fRandGaus(0,0.010/2);
			pV5tg_tr[2] = fRandGaus(0,0.00010/2);
			pV5tg_tr[3] = fRandGaus(0,0.010/2);
			pV5tg_tr[4] = fRandGaus(0,0.030/2);
		}

		if(iSmearX0<=0) pXtg_BPM_tr=0;
		else if(iSmearX0==1) pXtg_BPM_tr=pV5tg_tr[0];
		if(iSmearX0>=2)
		{
			if(iSmearX0>=3) 
			{
				if(Index<iNEvent/4) pBPMRes=(int(11.0*fRand()))*0.00025;
				else pBPMRes=(int(fLinearRand(1.0,0.0,0.0,11.0)))*0.00025;
			}
			else pBPMRes=0.001;
			//BPM reconstruction uncertianty 0.001 m
			pXtg_BPM_tr=pV5tg_tr[0]+fRandGaus(0,pBPMRes);
		}

		//iArm = 0,1,2 means random, left, right
		if(iArm==1) iLeftArm=1;
		else if(iArm==2) iLeftArm=0;
		else iLeftArm = (rand()>0.5*RAND_MAX)?0:1;
		//specify the reconstruction  package
		pRecDB=(iLeftArm)?pRecDBL:pRecDBR;

		for(int j=0;j<5;j++) pV5rec_tr[j]=pV5tg_tr[j]; 
		bool bGoodParticle=SNAKEThruHRS(iLeftArm, pHRSAngle, pXtg_BPM_tr, 
			pHRSMomentum,iFieldRotation, iExperiment, pV5rec_tr, pV5fp_tr);

		//throw away this event if delta_rec>=1.0
		if(bGoodParticle && pV5rec_tr[4]<1.0) 
		{
			Transform::X_TCS2HCS(pV5tg_tr[0],pV5tg_tr[2],0.0,pHRSAngle,pV5tg_lab[0],pV5tg_lab[2],pV5tg_lab[4]);
			Transform::P_TCS2HCS(pV5tg_tr[1],pV5tg_tr[3],pHRSAngle,pV5tg_lab[1],pV5tg_lab[3]);

			Transform::X_TCS2HCS(pV5rec_tr[0],pV5rec_tr[2],0.0,pHRSAngle,pV5rec_lab[0],pV5rec_lab[2],pV5rec_lab[4]);
			Transform::P_TCS2HCS(pV5rec_tr[1],pV5rec_tr[3],pHRSAngle,pV5rec_lab[1],pV5rec_lab[3]);

   
			//Reconstruct use optics database
			pRecDB->CalcTargetCoords(pV5fp_tr,pV5rec_db_tr);

			pSnake->Fill();
			Index++;
		}
		NThrown++;
	}
	pFile->Write();
	cout<<"\n"<<Index<<"/"<<NThrown<<" events saved into ntuple"<<endl;

	char grname[255];
	sprintf(grname,"Delta_Exp%02d_X%d_Distr%d.png",iExperiment,iSmearX0,iSourceDistr);
	PlotDelta(pSnake,grname);
	if(iSourceDistr>0)
	{
		if(iSmearX0==2)
		{
			sprintf(grname,"Dependence_Exp%02d_Distr%d.png",iExperiment,iSourceDistr);
			DrawDependence(pSnake,grname,0);
		}
		else if(iSmearX0>=3)
		{
			//with BPMRes
			sprintf(grname,"Dependence_Exp%02d_Distr%d.png",iExperiment,iSourceDistr);
			DrawDependence(pSnake,grname,1);
		}
	}
	pFile->Delete();
}


void PlotSNAKE( int iExperiment, int iSmearX0, int iSourceDistr)
{	
	SetStyle();

	TFile* pFile= new TFile(Form("snake_Exp%02d_X%d_Distr%d.root",
		iExperiment,iSmearX0,iSourceDistr));
	TTree* pSnake=(TTree*)gROOT->FindObject("s");
	
	char grname[255];
	sprintf(grname,"Delta_Exp%02d_X%d_Distr%d.png",iExperiment,iSmearX0,iSourceDistr);
	PlotDelta(pSnake,grname);
	if(iSourceDistr>0)
	{
		if(iSmearX0==2)
		{
			sprintf(grname,"Dependence_Exp%02d_Distr%d.png",iExperiment,iSourceDistr);
			DrawDependence(pSnake,grname,0);
		}
		else if(iSmearX0>=3)
		{
			//with BPMRes
			sprintf(grname,"Dependence_Exp%02d_Distr%d.png",iExperiment,iSourceDistr);
			DrawDependence(pSnake,grname,1);
		}
	}
	pFile->Delete();
}

