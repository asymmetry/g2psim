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

void SetStyle4000()
{
	gStyle->SetPalette(1);
	gStyle->SetFillColor(0);

	double pPadRightMargin=0.05;
	gStyle->SetPadRightMargin(pPadRightMargin); 
	gStyle->SetStatX(1.0-gStyle->GetPadRightMargin());
	gStyle->SetStatY(1.0-gStyle->GetPadTopMargin());

	gStyle->SetOptStat(0);	
	gStyle->SetOptFit(011);
	gStyle->SetStatH(0.08);
	gStyle->SetStatW(0.25);
	gStyle->SetStatStyle(4000);

	gStyle->SetTitleH(0.085);
	gStyle->SetTitleBorderSize(0);
	gStyle->SetTitleX(0.15);
	
}

void SetStyle()
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
		TPaveText *pt = new TPaveText(xx,0.25,xx+0.45,0.50,"brNDC");
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

TH1* DrawTH1_RMS(TTree *tree, char *hname, char *title, char *target, TCut &cut, char* hrange, 
				 char *option="", int color=1, double RMS_factor=3.0)
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

TH1* DrawTH1(TTree *tree, char *hname, char *title, char *target, TCut &cut, char* hrange, 
			 char *option="", int color=4)
{
	char hdest[512];
	TTree *track0=(TTree*)gROOT->FindObject("track0");

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


void PlotDelta(char *treename="s",char *grname="Delta.png")
{
	SetStyle();
	TTree *s=(TTree*)gROOT->FindObject(treename);

	TCanvas *pCan=new TCanvas("pCan","",1350,850);
	pCan->Divide(5,4,0.001,0.001);
	TCut pCut="";
	TPad *pPad=0;
	TH1 *h1; h1=0;
	TH2 *h2; h2=0;

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
		h1=DrawTH1_RMS(s,hName,strTitle,strTg,pCut,strRange,"",1,3.0);
		delete h1;

		sprintf(hName,"h1%stg_tr",strKey[i]);
		sprintf(strTg,"%stg_tr*1000",strKey[i]);
		h1=DrawTH1(s,hName,strTitle,strTg,pCut,strRange,"",1);
		sprintf(hName,"h1%srec_tr",strKey[i]);
		sprintf(strTg,"%srec_tr*1000",strKey[i]);
		h1=DrawTH1(s,hName,strTitle,strTg,pCut,strRange,"same",4);
	} 

	i=4;
	pPad=(TPad*)(pCan->cd(i+11));
	sprintf(strTitle,"%s(black), %s_{rec}(blue) %s",
		strKeyTitle[i],strKeyTitle[i],strUnit[i]);
	//this part to get the range
	sprintf(hName,"h1%srec",strKey[i]);
	sprintf(strTg,"%srec*10000",strKey[i]);
	h1=DrawTH1_RMS(s,hName,strTitle,strTg,pCut,strRange,"",1,3.0);
	delete h1;

	sprintf(hName,"h1%s",strKey[i]);
	sprintf(strTg,"%s*10000",strKey[i]);
	h1=DrawTH1(s,hName,strTitle,strTg,pCut,strRange,"",1);
	sprintf(hName,"h1%srec",strKey[i]);
	sprintf(strTg,"%srec*10000",strKey[i]);
	h1=DrawTH1(s,hName,strTitle,strTg,pCut,strRange,"same",4);

	for(i=0; i<4;i++)
	{
		pPad=(TPad*)(pCan->cd(i+16));

		sprintf(strTitle,"%s_{tg}^{lab}(black), %s_{rec}^{lab}(blue) %s",
			strKeyTitle[i],strKeyTitle[i],strUnit[i]);

		//this part to get the range
		sprintf(hName,"h1%srec_lab",strKey[i]);
		sprintf(strTg,"%srec_lab*1000",strKey[i]);
		h1=DrawTH1_RMS(s,hName,strTitle,strTg,pCut,strRange,"",1,3.0);
		delete h1;

		sprintf(hName,"h1%s_lab",strKey[i]);
		sprintf(strTg,"%stg_lab*1000",strKey[i]);
		h1=DrawTH1(s,hName,strTitle,strTg,pCut,strRange,"",1);
		sprintf(hName,"h1%srec_lab",strKey[i]);
		sprintf(strTg,"%srec_lab*1000",strKey[i]);
		h1=DrawTH1(s,hName,strTitle,strTg,pCut,strRange,"same",4);
	} 

	cout<<"grname="<<grname<<endl;
	pCan->SaveAs(grname);

}


//dependence on BPM resolution for x0_tr,y0_tr,theta0_tr,phi0_tr,delta
void DrawBPMDependence(char *treename="s",int iDrawContent=3, 
					   char *grname="BPMDependence.png")
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
	TCut theCut="";
	system("mkdir -p graph");

	const int nitemx=1;
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

	char *xlist[nitemx]={"BPMRes*1000"};
	char *xnamelist[nitemx]={"BPMRes"};
	char *xunit[nitemx]={"(mm)"};

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

