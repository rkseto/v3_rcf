#include "TROOT.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TRandom.h"
#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TBox.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TPolyLine.h"
#include "TLegend.h"
#include "TPaveText.h"

void mpcstatwpp_nice_newplt (bool iscentral = false, int nstp=0, double frac_typeA = 0.10,double frac_typeApp = 0.10, double frac_sysC = 0.03, double frac_sysB1 = 0.072, double frac_sysB2 = 0.007, double frac_sysB4 = 0.01, double frac_sysB7 = 0.01, float ffrag=0.6, float rg_s = 1.34,float rgpp_s=1.34,
		      char* Nucleus="Au",
		      int pA_lumi=450,
		      int pp_lumi=50
) {

  //************* FOR DEBUGGING ONLY****************
  //  cout<<"****************** setting nstp=5 for DEBUGGING ***************"<<endl;
  //  nstp=5;
  // END FOR DEBUGGING


  //  the last 3 in the argument are just for labeling
  // THIS VERSION ALSO MAKES NICE PLOTS ALA BUP 2013
  // THIS VERSION ADDS ERROR IN PP
  cout<<" iscentral="<< iscentral<<endl; 
  cout<<" nstp="<<  nstp<<endl; 
  cout<<" frac_typeA="<<  frac_typeA<<endl; 
  cout<<" frac_typeApp="<<  frac_typeApp<<endl; 
  cout<<" frac_sysC="<<  frac_sysC<<endl; 
  cout<<" frac_sysB1 (1/r corr)="<<  frac_sysB1<<endl; 
  cout<<" frac_sysB2 (inc) ="<<  frac_sysB2<<endl; 
  cout<<" frac_sysB4 (1/r_dau) ="<<  frac_sysB4<<endl; 
  cout<<" frac_sysB7 (1/r pp)="<<  frac_sysB7<<endl; 
  cout<<" ffrag="<<  ffrag<<endl; 
  cout<<" rg_s="<<   rg_s <<endl; 
  cout<<" rgpp_s="<<  rgpp_s<<endl; 

  
  bool override_typeA = true; 
  bool rand_points = false; 
  double frac_sysB = 0.0;
  int iseed = 7;
  
  int nstepc=10;   // these should end up being ~ 60 or so
  int nstepb1=10;
  int nstepb4=10;
  int nstepb7=10;
  int nstepb2=10;  

  if(nstp>0){
    nstepc=nstp;
    nstepb1=nstp;
    nstepb4=nstp;
    nstepb7=nstp;
    nstepb2=nstp;  
  }

  /*
    B1 - sys error on 1/r correlated dau and pp  0.06
    B4 - sys error on dau 1/r   0.02
    B7 - sys error on pp 1/r    0.02
    B2 -  sys error on i'/i     0.02
    C - sys error on ncoll      0.03
  */  


  // the following numbers come from sarah's full simulation
  // we dont have them for all rapidity bins so just assume all the same

  // first for dAu
  float inclusive_s=2000.;  // this will be reset using rg
  float pi0_s=1500.;
  float gamma_pi0_calc_s=1.03;
  
  // then for pp
  float inclusivepp_s=2000.; // this will be reset using rgpp
  float pi0pp_s=1556.;
  float gamma_pi0_calcpp_s=1.03;

  // set inclusive_s and inclusivepp_s using r_gamma 
  inclusive_s=pi0_s*gamma_pi0_calc_s*rg_s;
  inclusivepp_s=pi0pp_s*gamma_pi0_calcpp_s*rgpp_s;

  
  bool printout = true;

  cout << " " << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "J.L. Nagle -- compiled macro, modified from code from R. Seto 05-17-2010 modified again by R. Seto 02-06-2012,  modified chisq based on Annote 563" << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << " " << endl;
  cout << "Override Staterr = " << override_typeA << " frac stat = " << frac_typeA << " sysB = " << frac_sysB 
       << " sysC = " << frac_sysC << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << " " << endl;
  
  //Reset ROOT
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //gStyle->SetOptDate(1);
  //  gROOT -> ProcessLine (".X init.C");

  //  gStyle -> SetTitleXSize (0.05);
  //  gStyle -> SetTitleYSize (0.05);
  gStyle -> SetTitleYSize (0.045);
  //  gROOT -> ForceStyle();


  
  float r_gamma_sav=0;
  float r_gammapp_sav=0;
  float RdAu_sav=0;
  
  // --- INPUT FILES ---
  char dat1[50]="anaphpythia_ana10.root";
  cout << "== Input file with R_dAu versus eta (name = " << dat1 << ")" << endl;
  TFile *fc = new TFile(dat1,"read");
  
  //  char dat2[50]="focal_5.root";
  char dat2[50]="epsgeom.root"; 
  cout << "== Input file with actual R_G values from EPS09 (name = " << dat2 << ")" << endl;
  TFile *fb = new TFile(dat2,"read");
  
  // note that the index numbers for the calculations are the same (i.e. perfect
  // correspondence) between the two files
  
  //-----------------------------------------------------------
  // --- GET HISTOGRAMS OF INTEREST ---
  TH1D *etafG[40];
  TH1D *etafG2[34];
  for (int i=0;i<32;i++) {
    char tempname[256];
    char tempname2[256];
    sprintf(tempname,"etafG_p4mbdi_w%1d",i);
    if(iscentral)sprintf(tempname,"etafG_p4cedi_w%1d",i);
    etafG[i] = (TH1D*) fc->Get(tempname);
    if (etafG[i] == NULL) cout << "Problem with etafG get " << i <<" "<<tempname<< endl;

    sprintf(tempname2,"etafG_p4mbfr_w%1d",i);    
    if(iscentral)sprintf(tempname2,"etafG_p4cefr_w%1d",i);

    etafG2[i] = (TH1D*) fc->Get(tempname2);
    if (etafG2[i] == NULL) cout << "Problem with etafG2 get " << i <<" "<<tempname2<< endl;
    int nbinsx=etafG[i]->GetNbinsX();
    for(int j=1; j<=nbinsx; j++) {

      int j1=int(etafG[i]->GetBinContent(j));
      int j2=int(etafG2[i]->GetBinContent(j));
      
      //      cout<<" dir frag "<<j1<<" "<<j2<<endl;
      float xx=j1+j2*ffrag;
      
      //      float bb=etafG[i]->Integral()/etafG2[i]->Integral()*(ffrag/(1-ffrag));
      //      cout<<" bb="<<bb<<endl;
      //     float xx=j1+j2*bb;
      
      int j3= int(xx);
      etafG[i]->SetBinContent(j,j3);
    }
    
  }
  
  
  TH1D *G_x2[34];
  for (int i=1;i<32;i++) {
    char tempname[256];
    //    sprintf(tempname,"Gmbq02_%1d",i);
    sprintf(tempname,"Gmbq09_%1d",i);
    if(iscentral)sprintf(tempname,"Gcentq02_%1d",i);
    G_x2[i] = (TH1D*) fb->Get(tempname);
    if (G_x2[i] == NULL) cout << "Problem with G2_x get " << i << endl;
  }
  

  int jjj=1; // regular

  //  jjj=16; // upper
  //  jjj=17; // lower


  TH1D *save = new TH1D("save","save", 1, 0., 1.0);
  etafG[1]->Copy(*save);
  etafG[jjj]->Copy(*etafG[1]);
  save->Copy(*etafG[jjj]);
  etafG2[1]->Copy(*save);
  etafG2[jjj]->Copy(*etafG2[1]);
  save->Copy(*etafG2[jjj]);
  G_x2[1]->Copy(*save);
  G_x2[jjj]->Copy(*G_x2[1]);
  save->Copy(*G_x2[jjj]);


  // rebinning
  int rebin_factor =1;
  for (int i=0;i<32;i++) {
    etafG[i]->Rebin(rebin_factor);
    //    if (i>0) etafG[i]->Divide(etafG[0]);
    if (i>1) etafG[i]->Divide(etafG[0]); // don't divide out the pp for the central value. Do it during the calculation
  }
  
  // clone 34-36 to have a place for inclusive, pi0, gammatopi0calc==========
  etafG[34]= (TH1D*) etafG[1]->Clone("etafGinclusive");
  etafG[35]= (TH1D*) etafG[1]->Clone("etafGpi0");
  etafG[36]= (TH1D*) etafG[1]->Clone("etafGgampi0calc");
  etafG[37]= (TH1D*) etafG[1]->Clone("etafGinclusivepp");
  etafG[38]= (TH1D*) etafG[1]->Clone("etafGpi0pp");
  etafG[39]= (TH1D*) etafG[1]->Clone("etafGgampi0calcpp");

  for(int j=1; j<=etafG[1]->GetNbinsX(); j++) {
    // the following numbers come from sarah's full simulation
    // we dont have them for all rapidity bins so just assume all the same
    // first for dAu
    float inclusive=inclusive_s;
    float pi0=pi0_s;
    float gamma_pi0_calc=gamma_pi0_calc_s;
    
    // then for pp
    float inclusivepp=inclusivepp_s;
    float pi0pp=pi0pp_s;
    float gamma_pi0_calcpp=gamma_pi0_calcpp_s;

    float gamma_pi0_meas=inclusive/pi0;
    float r_gamma=gamma_pi0_meas/gamma_pi0_calc;

    r_gamma_sav=r_gamma;
    float ndirgamma_meas=inclusive*(1.-1./r_gamma);
    float alpha=etafG[1]->GetBinContent(j)/ndirgamma_meas; // this is a scale factor between etafG and the number of events given by sarah. The normalization is taken to be arbitrary. Put the statistical errors in by hand as frag_typeA 

    inclusive*=alpha;
    pi0*=alpha;

    etafG[34]->SetBinContent(j,inclusive);
    etafG[34]->SetBinError(j,0.0);
    etafG[35]->SetBinContent(j,pi0);
    etafG[35]->SetBinError(j,0.0);
    etafG[36]->SetBinContent(j,gamma_pi0_calc);
    etafG[36]->SetBinError(j,0.0);

    float gamma_pi0_measpp=inclusivepp/pi0pp;
    float r_gammapp=gamma_pi0_measpp/gamma_pi0_calcpp;
    r_gammapp_sav=r_gammapp;
    float ndirgamma_measpp=inclusivepp*(1.-1./r_gammapp);
    float alphapp=etafG[0]->GetBinContent(j)/ndirgamma_measpp;
    inclusivepp*=alphapp;
    pi0pp*=alphapp;
    etafG[37]->SetBinContent(j,inclusivepp);
    etafG[37]->SetBinError(j,0.0);
    etafG[38]->SetBinContent(j,pi0pp);
    etafG[38]->SetBinError(j,0.0);
    etafG[39]->SetBinContent(j,gamma_pi0_calcpp);
    etafG[39]->SetBinError(j,0.0);

    RdAu_sav=ndirgamma_meas/ndirgamma_measpp;
  } // j
  
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // to do - figure out etafG[1,0] using r_gamma and make sure you get the original values back. it has to be dont 1-j i.e. for each bin
  // now get back etafG[1] and etafG[0]
  for(int j=1; j<=etafG[1]->GetNbinsX(); j++) {
    
    // first for dAu
    float inclusive=etafG[34]->GetBinContent(j);
    float pi0=etafG[35]->GetBinContent(j);
    float gamma_pi0_calc=etafG[36]->GetBinContent(j);
    float gamma_pi0_meas=inclusive/pi0;
    float r_gamma=gamma_pi0_meas/gamma_pi0_calc;
    float ndirgamma_meas=inclusive*(1.-1./r_gamma);

    // then for pp
    float inclusivepp=etafG[37]->GetBinContent(j);
    float pi0pp=etafG[38]->GetBinContent(j);
    float gamma_pi0_calcpp=    etafG[39]->GetBinContent(j);
    float gamma_pi0_measpp=inclusivepp/pi0pp;
    float r_gammapp=gamma_pi0_measpp/gamma_pi0_calcpp;
    float ndirgamma_measpp=inclusivepp*(1.-1./r_gammapp);
    // now do the comparison - they should match 
    //    cout<<" >> j dau etafg pp etafg "<<j<<" "<<ndirgamma_meas<<" "<< etafG[1]->GetBinContent(j)<<" "<< ndirgamma_measpp<<" "<<etafG[0]->GetBinContent(j)<<endl;
    
    // OK WE USE THE R_gamma method to find rdau whichc is etafg[1]
    float RdAu=ndirgamma_meas/ndirgamma_measpp;
    etafG[1]->SetBinContent(j,RdAu);
    etafG[1]->SetBinError(j,0.0);
    //    cout<<" ************************Rdau ="<<RdAu<<endl;
  } // j
 
  //===========================================================
  
  // zero out bins which are out of new configuration of MPCEX
  // the real "fake" data is stored in this one etafG[1]
  double etamin=3.1;
  double etamax=3.8;
  cout << "New MPCEX configuration with eta range:  " << etamin << " to " << etamax << endl;
  int nbinsx=etafG[1]->GetNbinsX();
  for(int j=1; j<=nbinsx; j++) {
    if ((etafG[1]->GetBinLowEdge(j) < etamin) || (etafG[1]->GetBinLowEdge(j) > etamax) ) {
      etafG[1]->SetBinContent(j,0.0);
      etafG[1]->SetBinError(j,0.0);
    } else {
      // if selected -> randomly throw points based on statistical uncertainty
      TF1 *fgaus = new TF1("fgaus","TMath::Exp(-x*x/(2.*1.0*1.0))",-4.0,4.0);
      fgaus->SetNpx(10000);
      volatile double foobar;
      for (int ii=0; ii<iseed;ii++) foobar = fgaus->GetRandom(); // how to plant a new seed number here ?
      if (override_typeA) {
	etafG[1]->SetBinError(j, etafG[1]->GetBinContent(j) * 
			      sqrt(frac_typeA*frac_typeA+frac_typeApp*frac_typeApp)
);
	if (rand_points) etafG[1]->SetBinContent(j, etafG[1]->GetBinContent(j)*(1.0+sqrt(frac_typeA*frac_typeA+frac_typeApp*frac_typeApp)*fgaus->GetRandom()));
      }
    }
  } //nbinsx
  
  // calculate modified chi^2 between data (stored in etafG[1]) and various theories
  // and store in an array
  double mod_chi2_total[32];
  double mod_best_typeCoffset[32];
  double mod_best_typeB1offset[32];
  double mod_best_itiltB1[32];
  double mod_best_typeB2offset[32];
  double mod_best_itiltB2[32];
  double mod_best_typeB4offset[32];
  double mod_best_itiltB4[32];
  double mod_best_typeB7offset[32];
  double mod_best_itiltB7[32];


  for (int itest=0;itest<32;itest++) mod_chi2_total[itest] = 9999999.0;
  for (int itest=0;itest<32;itest++) mod_best_typeCoffset[itest] = 9999999.0;

  for (int itest=0;itest<32;itest++) mod_best_typeB1offset[itest] = 9999999.0;
  for (int itest=0;itest<32;itest++) mod_best_itiltB1[itest] = 9999999.0;
  for (int itest=0;itest<32;itest++) mod_best_typeB2offset[itest] = 9999999.0;
  for (int itest=0;itest<32;itest++) mod_best_itiltB2[itest] = 9999999.0;
  for (int itest=0;itest<32;itest++) mod_best_typeB4offset[itest] = 9999999.0;
  for (int itest=0;itest<32;itest++) mod_best_itiltB4[itest] = 9999999.0;
  for (int itest=0;itest<32;itest++) mod_best_typeB7offset[itest] = 9999999.0;
  for (int itest=0;itest<32;itest++) mod_best_itiltB7[itest] = 9999999.0;


  // loop over different theories
  for (int itest=1;itest<32;itest++) { 

    cout<<" theory step "<<itest<<" of 31"<<endl;
    
    // need to consider type B and C separated, and allow B to tilt or not tilt (exclusive to one case for now)
    
    
    for (int isysC=0; isysC<nstepc; isysC++) {
      
      for (int isysB2=0; isysB2<nstepb1; isysB2++) {
	for (int itiltB2=0; itiltB2<2; itiltB2++) {
	  
	  for (int isysB1=0; isysB1<nstepb1; isysB1++) {
	    for (int itiltB1=0; itiltB1<2; itiltB1++) {
	      
	      for (int isysB4=0; isysB4<nstepb4; isysB4++) {
		for (int itiltB4=0; itiltB4<2; itiltB4++) {
		  
		  for (int isysB7=0; isysB7<nstepb7; isysB7++) {
		    for (int itiltB7=0; itiltB7<2; itiltB7++) {
		      
		      double mod_chi2_total_temp = 0.0;
		      double nsigma_sysC = (-3.0 + (((double) isysC)*(6.0/float(nstepc))));
		      double nsigma_sysB2 = (-3.0 + (((double) isysB2)*(6.0/float(nstepb2))));
		      double nsigma_sysB1 = (-3.0 + (((double) isysB1)*(6.0/float(nstepb1))));
		      double nsigma_sysB4 = (-3.0 + (((double) isysB4)*(6.0/float(nstepb4))));
		      double nsigma_sysB7 = (-3.0 + (((double) isysB7)*(6.0/float(nstepb7))));
		      
		      double mult_sysC = 1.0 + (nsigma_sysC * frac_sysC); // multiplicative systematic factor
		      double mult_sysB2 = 1.0 + (nsigma_sysB2 * frac_sysB2); // multiplicative systematic factor
		      double mult_sysB1 = 1.0 + (nsigma_sysB1 * frac_sysB1); // multiplicative systematic factor
		      double mult_sysB4 = 1.0 + (nsigma_sysB4 * frac_sysB4); // multiplicative systematic factor
		      double mult_sysB7 = 1.0 + (nsigma_sysB7 * frac_sysB7); // multiplicative systematic factor
		      
		      // loop over all points
		      for (int ipt=1; ipt<=etafG[1]->GetNbinsX(); ipt++) {
			if (etafG[1]->GetBinContent(ipt)==0.0) continue; // no data point here
			
			int j=ipt;
			float inclusive=etafG[34]->GetBinContent(j);
			float pi0=etafG[35]->GetBinContent(j);		      
			float gamma_pi0_calc=etafG[36]->GetBinContent(j);
			float inclusivepp=etafG[37]->GetBinContent(j);
			float pi0pp=etafG[38]->GetBinContent(j);
			float gamma_pi0_calcpp=    etafG[39]->GetBinContent(j);
			
			
			double tiltminrange = 3.1;
			double tiltmaxrange = 3.8;
			double tiltmidpoint = (tiltmaxrange + tiltminrange)/2.0;
			// this should range from +1 to -1 across the range....
			double tiltfactor = 1.0 + nsigma_sysB1 * frac_sysB1 * 
			  ((etafG[1]->GetBinCenter(ipt) - tiltmidpoint) / ((tiltmaxrange - tiltminrange)/2.0));
			double tiltfactor2 = 1.0 + nsigma_sysB2 * frac_sysB2 * 
			  ((etafG[1]->GetBinCenter(ipt) - tiltmidpoint) / ((tiltmaxrange - tiltminrange)/2.0));
			double tiltfactor4 = 1.0 + nsigma_sysB4 * frac_sysB4 * 
			  ((etafG[1]->GetBinCenter(ipt) - tiltmidpoint) / ((tiltmaxrange - tiltminrange)/2.0));
			double tiltfactor7 = 1.0 + nsigma_sysB7 * frac_sysB7 * 
			  ((etafG[1]->GetBinCenter(ipt) - tiltmidpoint) / ((tiltmaxrange - tiltminrange)/2.0));
			
			float ffB4=mult_sysB4;
			float ffB1=mult_sysB1;
			float ffB7=mult_sysB7;
			float ffB4pp=mult_sysB4;
			float ffB1pp=mult_sysB1;
			float ffB7pp=mult_sysB7;

			float ffB2=mult_sysB2;			// this is error on all RdAu
			float ffB2pp=mult_sysB2;			

			if (itiltB4 > 0)ffB4=tiltfactor4;
			if (itiltB1 > 0)ffB1=tiltfactor;
			if (itiltB7 > 0)ffB7=tiltfactor7;
			if (itiltB4 > 0)ffB4pp=tiltfactor4;
			if (itiltB1 > 0)ffB1pp=tiltfactor;
			if (itiltB7 > 0)ffB7pp=tiltfactor7;

			if (itiltB2 > 0)ffB1pp=tiltfactor2;
			if (itiltB2 > 0)ffB1pp=tiltfactor2;

			float gamma_pi0_meas=inclusive/pi0;
			float r_gamma=gamma_pi0_meas/gamma_pi0_calc;
			r_gamma*=ffB1*ffB4;
			float ndirgamma_meas=inclusive*(1.-1./r_gamma);
			// then for pp
			float gamma_pi0_measpp=inclusivepp/pi0pp;
			float r_gammapp=gamma_pi0_measpp/gamma_pi0_calcpp;
			r_gammapp*=ffB1*ffB7;
			float ndirgamma_measpp=inclusivepp*(1.-1./r_gammapp);
			// OK WE USE THE R_gamma method to fund rdau whichc is etafg[1]
			float RdAu=ndirgamma_meas/ndirgamma_measpp;
			RdAu*=ffB2;
			float etafG1_errors=RdAu;  // THIS IS THE BIG STATEMENT

			// this mod chisq assume that B1,2,4,7 are the same for dAu and pp
			mod_chi2_total_temp += pow( (( etafG1_errors - etafG[itest]->GetBinContent(ipt)) /
						     etafG[1]->GetBinError(ipt)*mult_sysC*ffB1*ffB2*ffB4*ffB7) , 2.0);
		      } // end loop over all points
		      
		      // then add systematic once
		    mod_chi2_total_temp += pow(nsigma_sysC,2.0) + pow(nsigma_sysB1,2.0)+ pow(nsigma_sysB2,2.0)+ pow(nsigma_sysB4,2.0)+ pow(nsigma_sysB7,2.0);
		    
		    //	    if (itest==17) printf("nsigB = %2.1f nsigC %2.1f tilt %1d modchi = %3.2f\n",nsigma_sysB,nsigma_sysC,itiltB,mod_chi2_total_temp);
		    
		    if (mod_chi2_total_temp < mod_chi2_total[itest]) {
		      mod_chi2_total[itest] = mod_chi2_total_temp;
		      mod_best_typeB1offset[itest] = nsigma_sysB1;
		      mod_best_itiltB1[itest] = itiltB1;
		      mod_best_typeB4offset[itest] = nsigma_sysB4;
		      mod_best_itiltB4[itest] = itiltB4;
		      mod_best_typeB7offset[itest] = nsigma_sysB7;
		      mod_best_itiltB7[itest] = itiltB7;
		      mod_best_typeCoffset[itest] = nsigma_sysC;
		    }
		    
		    } //sys B7 tilt	    
		  } //sys B7
		} //sys B4 tilt	    
	      } //sys B4
	    } //sys B1 tilt	    
	  } //sys B1
	} //sys B2 tilt	    
      } //sys B2
    } // sysC end loop over possible systematic uncertainties
    
  } // end loop over all itest (all theories)
  
  // now determine the minimum, all theories within minimum + 1, all theories within minimum + 1.6^2 (? 90% C.L.)
  bool within_onesigma[32];
  bool within_16sigma[32];
  
  double chi2_best = 9999999;
  // do not include itest=1 since that was thrown from and now modified !!
  for (int itest=2; itest<32; itest++) if (mod_chi2_total[itest] < chi2_best) chi2_best = mod_chi2_total[itest];
  
  cout << "Best modified chi^2 = " << chi2_best << endl;
  
  for (int itest=1; itest<32; itest++) {
    
    printf("ITEST=%2d (RdAu(eta=2.8)=%3.2f)::: mchi2=%3.2f  bestBshift=%2.1f (tilt=%1d) bestCshift=%2.1f\n",
	   itest,etafG[itest]->GetBinContent(etafG[itest]->FindBin(2.7)),
	   mod_chi2_total[itest], mod_best_typeB1offset[itest], (int) mod_best_itiltB1[itest], mod_best_typeCoffset[itest]);
    
    if (mod_chi2_total[itest] < (chi2_best +1.0)) {
      within_onesigma[itest] = true;
    } else {
      within_onesigma[itest] = false;
    }
    if (mod_chi2_total[itest] < (chi2_best + pow(1.6,2.))) {
      within_16sigma[itest] = true;
    } else {
      within_16sigma[itest] = false;
    }
  }
      
  //=====================================================================  
  // clone to [1] to [32] and [33] to set up histograms
  etafG[32]= (TH1D*) etafG[1]->Clone("etafGmin");
  etafG[33]= (TH1D*) etafG[1]->Clone("etafGmax");
  G_x2[32]= (TH1D*) G_x2[1]->Clone("G_x2min");
  G_x2[33]= (TH1D*) G_x2[1]->Clone("G_x2max");
  
  for (int j=0; j<etafG[1]->GetNbinsX(); j++){
    float deltamin=0;
    float deltamax=0;
      for (int i=1;i<32;i++) {
	if(i==jjj)continue;
	float delt=etafG[i]->GetBinContent(j)-etafG[jjj]->GetBinContent(j);
	if(delt<0)deltamin+=delt*delt;
	if(delt>0)deltamax+=delt*delt;
      }
      deltamin=sqrt(deltamin);
      deltamax=sqrt(deltamax);
      float valmin=etafG[jjj]->GetBinContent(j)-deltamin;
      if(valmin<0)valmin=0;
      float valmax=etafG[jjj]->GetBinContent(j)+deltamax;
      etafG[32]->SetBinContent(j,valmin);
      etafG[33]->SetBinContent(j,valmax);

  }

  cout<<" NUMBER OF BINS "<<G_x2[1]->GetNbinsX()<<endl;

  for (int j=0; j<G_x2[1]->GetNbinsX(); j++){

    float deltamin=0;
    float deltamax=0;
    float dmin=9999.;
    float dmax=-9999.;
      for (int i=1;i<32;i++) {
	if(i==jjj)continue;
	float delt=G_x2[i]->GetBinContent(j)-G_x2[jjj]->GetBinContent(j);
	if(G_x2[jjj]->GetBinContent(j)<dmin)dmin=G_x2[jjj]->GetBinContent(j);
	if(G_x2[jjj]->GetBinContent(j)>dmax)dmax=G_x2[jjj]->GetBinContent(j);
	if(delt<0)deltamin+=delt*delt;
	if(delt>0)deltamax+=delt*delt;
      }
      deltamin=sqrt(deltamin);
      deltamax=sqrt(deltamax);
      float valmin=G_x2[jjj]->GetBinContent(j)-deltamin;
      if(valmin<0)valmin=0;
      float valmax=G_x2[jjj]->GetBinContent(j)+deltamax;
      G_x2[32]->SetBinContent(j,valmin);
      G_x2[33]->SetBinContent(j,valmax);

 }

//=====================================================================  

  //--------------------------------------------
  // --- DRAW STUFF ---

  TCanvas *c1 = new TCanvas("c1","c1",0,0,1400,700);
  c1->Divide(2,1);
  c1->cd(1);  

  double x1=+3.0;
  double x2=+4.0;
  double y1=+0.0;
  double y2=+1.0;
  TH1F *hr1 = new TH1F("hr1","hr1",100,x1,x2);
  hr1->SetMinimum(y1);
  hr1->SetMaximum(y2);
  //  TH1F *hr1 = c1->DrawFrame(x1,y1,x2,y2);
  hr1->SetXTitle("Photon Pseudorapidity");
  hr1->SetYTitle("R_{dAu}");
  hr1->DrawCopy();

  // TLegend box with labels
  TLegend *tl = new TLegend(0.11,0.10,0.7,0.25,"","brNDC");
  tl->SetTextSize(0.02);
  char tempname2[256];
  sprintf(tempname2,"Assume TypeA=%3.2f TypeApp=%3.2f TypeB2=%3.2f TypeC=%3.2f",frac_typeA,frac_typeApp,frac_sysB2,frac_sysC);
  tl->AddEntry("",tempname2,"");
  sprintf(tempname2,"B1(1/r corr)=%3.2f B4(1/r dAu)=%3.2f B7(1/r pp)=%3.2f",frac_sysB1,frac_sysB4,frac_sysB7);
  tl->AddEntry("",tempname2,"");
  sprintf(tempname2,"r_gamma=%3.2f r_gammapp=%3.2f R_dAu=%3.2f",r_gamma_sav,r_gammapp_sav,RdAu_sav);
  tl->AddEntry("",tempname2,"");
  if (iscentral) {
  sprintf(tempname2,"Central    ffrag=%3.2f",ffrag);
  tl->AddEntry("",tempname2,"");
  } else {
  sprintf(tempname2,"MinBias    ffrag=%3.2f",ffrag);
  tl->AddEntry("",tempname2,"");
  }
  tl->AddEntry("","Blue=1 RMS Band, Magenta=90% C.L. Band","");
  tl->Draw("same");

  etafG[1]->SetLineColor(2);
  etafG[1]->SetLineWidth(5);
  etafG[1]->SetMarkerStyle(20);
  etafG[1]->SetMarkerColor(2);
  etafG[1]->DrawCopy("p,e,SAME");   

  for (int i=2;i<32;i++) {
    etafG[i]->SetLineColor(1);
    etafG[i]->SetLineWidth(3);
    if (within_16sigma[i]) etafG[i]->SetLineColor(6);
    if (within_16sigma[i]) etafG[i]->SetLineWidth(3);
    if (within_onesigma[i]) etafG[i]->SetLineColor(4);
    etafG[i]->DrawCopy("HIST SAME");
  }

  etafG[1]->DrawCopy("p,e,SAME");

  etafG[32]->SetLineColor(7);
  etafG[33]->SetLineColor(7);
  etafG[32]->SetLineWidth(3);
  etafG[33]->SetLineWidth(3);
  etafG[32]->DrawCopy("HIST SAME");
  etafG[33]->DrawCopy("HIST SAME");


  c1->cd(2);

  x1=-4.0;
  x2=+0.0;
  y1=+0.0;
  y2=+1.5;
  TH1F *hr3 = new TH1F("hr3","hr3",100,x1,x2);
  hr3->SetMinimum(y1);
  hr3->SetMaximum(y2);
  hr3->SetXTitle("log10(x)");
  //  hr3->SetYTitle("EPS09 Nuclear Mod. R_{G}   [Q^{2}=1.69 GeV^{2}]");
  hr3->SetYTitle("EPS09 Nuclear Mod. R_{G}   [Q^{2}=9 GeV^{2}]");
  hr3->DrawCopy();
  
  // box corresponds to the x-range of sensitivity of the measurement (approx).
  TBox *b = new TBox(-2.6,0.05,-1.9,1.45);
  b->SetFillColor(5);
  //  b->Draw("same");

  for (int i=1;i<32;i++) {
    G_x2[i]->SetLineWidth(3);
    G_x2[i]->SetLineColor(1);
    if (within_16sigma[i]) G_x2[i]->SetLineColor(6);
    if (within_16sigma[i]) G_x2[i]->SetLineWidth(3);
    if (within_onesigma[i]) G_x2[i]->SetLineColor(4);
    G_x2[i]->DrawCopy("SAME");
  }

  G_x2[32]->SetLineColor(7);
  G_x2[33]->SetLineColor(7);
  G_x2[32]->SetLineWidth(5);
  G_x2[33]->SetLineWidth(5);
  G_x2[32]->DrawCopy("HIST SAME");
  G_x2[33]->DrawCopy("HIST SAME");


  if(nstp<10)sprintf(tempname2,"minbias_nstp%1.0f_A%3.2f_App%3.2f_C%3.2f_1B%3.2f_2B%3.2f_4B%3.2f_7B%3.2f_ffrag%3.2f_rg%3.2f_rgpp%3.2f.gif",
		     float(nstp), frac_typeA, frac_typeApp, frac_sysC, frac_sysB1, frac_sysB2, frac_sysB4, frac_sysB7, ffrag, rg_s,rgpp_s);
  if(nstp>9)sprintf(tempname2,"minbias_nstp%2.0f_A%3.2f_App%3.2f_C%3.2f_1B%3.2f_2B%3.2f_4B%3.2f_7B%3.2f_ffrag%3.2f_rg%3.2f_rgpp%3.2f.gif",
     float(nstp), frac_typeA, frac_typeApp,frac_sysC, frac_sysB1, frac_sysB2, frac_sysB4, frac_sysB7, ffrag, rg_s,rgpp_s);
  if(nstp<10 && iscentral)sprintf(tempname2,"central_nstp%1.0f_A%3.2f_App%3.2f_C%3.2f_1B%3.2f_2B%3.2f_4B%3.2f_7B%3.2f_ffrag%3.2f_rg%3.2f_rgpp%3.2f.gif",
     float(nstp), frac_typeA, frac_typeApp, frac_sysC, frac_sysB1, frac_sysB2, frac_sysB4, frac_sysB7, ffrag, rg_s,rgpp_s);
  if(nstp>9 && iscentral)sprintf(tempname2,"central_nstp%2.0f_A%3.2f_App%3.2f_C%3.2f_1B%3.2f_2B%3.2f_4B%3.2f_7B%3.2f_ffrag%3.2f_rg%3.2f_rgpp%3.2f.gif",
     float(nstp), frac_typeA, frac_typeApp, frac_sysC, frac_sysB1, frac_sysB2, frac_sysB4, frac_sysB7, ffrag, rg_s,rgpp_s);

    cout<<tempname2<<endl;

  c1->Update();
  //  if (printout) c1->Print("plot.gif");
  //  if (printout) c1->Print(tempname2);




  //--------------------------------------------
  // --- DRAW STUFF2 ---

  double x[7];
  double xu[7];
  double xl[7];
  double y[7];

  double y1l[7];
  double y1u[7];
  double y2l[7];
  double y2u[7];
  double y3l[7];
  double y3u[7];
  double y4l[7];
  double y4u[7];
  double y5l[7];
  double y5u[7];

  int itestmax=0;
  int itestmin=0;
  int itest1sigmamax=0;
  int itest1sigmamin=0;
  int itest16sigmamax=0;
  int itest16sigmamin=0;

  float testmax=-1000;
  float testmin=1000;
  float test1sigmamax=-1000;
  float test1sigmamin=1000;
  float test16sigmamax=-1000;
  float test16sigmamin=1000;

  for (int itest=1; itest<32; itest++) {

    if(etafG[itest]->GetBinContent(4)>testmax){
      testmax=etafG[itest]->GetBinContent(4);
      itestmax=itest;
    }

    if(etafG[itest]->GetBinContent(4)<testmin){
      testmin=etafG[itest]->GetBinContent(4);
      itestmin=itest;
    }

    if(within_onesigma[itest] && etafG[itest]->GetBinContent(4)>test1sigmamax){
      test1sigmamax=etafG[itest]->GetBinContent(4);
      itest1sigmamax=itest;
    }

    if(within_onesigma[itest] && etafG[itest]->GetBinContent(4)<test1sigmamin){
      test1sigmamin=etafG[itest]->GetBinContent(4);
      itest1sigmamin=itest;
    }

    if(within_16sigma[itest] && etafG[itest]->GetBinContent(4)>test16sigmamax){
      test16sigmamax=etafG[itest]->GetBinContent(4);
      itest16sigmamax=itest;
    }

    if(within_16sigma[itest] && etafG[itest]->GetBinContent(4)<test16sigmamin){
      test16sigmamin=etafG[itest]->GetBinContent(4);
      itest16sigmamin=itest;
    }

  }

  cout<<" test min max "<<itestmin<<" "<<itestmax<<endl;

    for(int j=2; j<=9; j++) {
      x[j-2]=etafG[1]->GetBinCenter(j)+.05;
      y[j-2]=etafG[1]->GetBinContent(j);

      xu[j-1]=0.01;
      xl[j-1]=0.01;

      // here the width of the line
      y1u[j-2]=0.001;
      y1l[j-2]=0.001;

      // here 1 sigma
      y2u[j-2]=etafG[itest1sigmamax]->GetBinContent(j)-etafG[1]->GetBinContent(j);
      y2l[j-2]=-(etafG[itest1sigmamin]->GetBinContent(j)-etafG[1]->GetBinContent(j));

      // here 90% CL
      y3u[j-2]=etafG[itest16sigmamax]->GetBinContent(j)-etafG[1]->GetBinContent(j);
      y3l[j-2]=-(etafG[itest16sigmamin]->GetBinContent(j)-etafG[1]->GetBinContent(j));

      // here is the rest
      y4u[j-2]=etafG[itestmax]->GetBinContent(j)-etafG[1]->GetBinContent(j);
      y4l[j-2]=-(etafG[itestmin]->GetBinContent(j)-etafG[1]->GetBinContent(j));

      // here the outer extent
      y5u[j-2]=etafG[33]->GetBinContent(j)-etafG[1]->GetBinContent(j);
      y5l[j-2]=-(etafG[32]->GetBinContent(j)-etafG[1]->GetBinContent(j));

      cout<<y4l[j-2]<<" "<<y4u[j-2]<<endl;
    }


  TCanvas *c2 = new TCanvas("c2","c2",0,0,700,700);

  TH1F *hFrame2=  c2->DrawFrame(3.1,0.,3.8,1.0);
  hFrame2->GetXaxis()->SetTitle("#eta");
  hFrame2->GetYaxis()->SetTitle("R_{dAu}");

  TGraphAsymmErrors *gr1 = new TGraphAsymmErrors(7,x,y,xl,xu,y1l,y1u);
  gr1->SetFillColor(kRed);

  TGraphAsymmErrors *gr2 = new TGraphAsymmErrors(7,x,y,xl,xu,y2l,y2u);
  gr2->SetFillColor(kBlue);
  //  gr2->SetFillColor(kBlue+2);
  //  gr2->SetFillColor(kPink);
  //  gr2->SetFillColor(kRed);
  //  gr2->SetFillStyle(3023);

  TGraphAsymmErrors *gr3 = new TGraphAsymmErrors(7,x,y,xl,xu,y3l,y3u);
  gr3->SetFillColor(kCyan);
  //  gr3->SetFillColor(kBlue);
  //  gr3->SetFillColor(kPink+1);
  //  gr3->SetFillColor(kMagenta);
  //  gr3->SetFillColor(kRed);
  //  gr3->SetFillStyle(3002);

  TGraphAsymmErrors *gr4 = new TGraphAsymmErrors(7,x,y,xl,xu,y4l,y4u);
  //  gr4->SetFillColor(kBlack);
  //  gr4->SetFillColor(kYellow);
  gr4->SetFillColor(kBlack);
  gr4->SetFillStyle(3005);

  TGraphAsymmErrors *gr5 = new TGraphAsymmErrors(7,x,y,xl,xu,y5l,y5u);
  //  gr5->SetFillColor(kCyan);
  //  gr5->SetFillColor(kYellow);
  gr5->SetFillColor(kBlack);
  gr5->SetFillStyle(3004);

  //  gr5->Draw("a4");
  gr5->Draw("3");
  gr4->Draw("3");
  gr3->Draw("3");
  gr2->Draw("3");
  gr1->Draw("3");

  c2->Update();
  if (printout) c2->Print("plot2wppnice_newplt.png");

  //====================================================================================
  //====================================================================================

  double gx[1000];
  double gxu[1000];
  double gxl[1000];
  double gy[1000];

  double gy1l[1000];
  double gy1u[1000];
  double gy2l[1000];
  double gy2u[1000];
  double gy3l[1000];
  double gy3u[1000];
  double gy4l[1000];
  double gy4u[1000];
  double gy5l[1000];
  double gy5u[1000];

  for(int j=0; j<1000; j++) {
    
    gx[j]=G_x2[1]->GetBinCenter(j);
    gy[j]=G_x2[1]->GetBinContent(j);
    
    gxu[j]=0.01;
    gxl[j]=0.01;
    
    // here the width of the line
    gy1u[j]=0.001;
    gy1l[j]=0.001;
    
    // here 1 sigma
    gy2u[j]=G_x2[itest1sigmamax]->GetBinContent(j)-G_x2[1]->GetBinContent(j);
    gy2l[j]=-(G_x2[itest1sigmamin]->GetBinContent(j)-G_x2[1]->GetBinContent(j));
    
    // here 90% CL
    gy3u[j]=G_x2[itest16sigmamax]->GetBinContent(j)-G_x2[1]->GetBinContent(j);
    gy3l[j]=-(G_x2[itest16sigmamin]->GetBinContent(j)-G_x2[1]->GetBinContent(j));
    
    // here is the rest
    gy4u[j]=G_x2[itestmax]->GetBinContent(j)-G_x2[1]->GetBinContent(j);
    gy4l[j]=-(G_x2[itestmin]->GetBinContent(j)-G_x2[1]->GetBinContent(j));
    
    // here the outer extent
    gy5u[j]=G_x2[33]->GetBinContent(j)-G_x2[1]->GetBinContent(j);
    gy5l[j]=-(G_x2[32]->GetBinContent(j)-G_x2[1]->GetBinContent(j));
    
  }
  
  
  TCanvas *c3 = new TCanvas("c3","c3",0,0,700,700);

  TH1F *hFrame3=c3->DrawFrame(-3.,0.,0.,1.5);
  hFrame3->GetXaxis()->SetTitle("Log_{10}(x_{gluon})");
  //  hFrame3->GetYaxis()->SetTitle("EPS09 Nuclear Mod. R_{G} [Q^{2}=1.69GeV^{2}]");
  hFrame3->GetYaxis()->SetTitle("EPS09 Nuclear Mod. R_{G} [Q^{2}=9GeV^{2}]");


  // draw box between log(x2) of -3 and -2
   Double_t xpoly[5] = {-3.,-3.,-2.,-2.0,-3.};
   Double_t ypoly[5] = {0.,1.5,1.5,0.,0.};
   TPolyLine *pline = new TPolyLine(5,xpoly,ypoly);
   pline->SetFillColor(18);
   pline->SetLineColor(1);
   pline->SetLineWidth(1);
   pline->Draw("f");
   pline->Draw();


   // THE 568 here keep the limits from being plotted for log(x2)>-1.7
   //  TGraphAsymmErrors *ggr1 = new TGraphAsymmErrors(1000,gx,gy,gxl,gxu,gy1l,gy1u);
  TGraphAsymmErrors *ggr1 = new TGraphAsymmErrors(568,gx,gy,gxl,gxu,gy1l,gy1u);
  ggr1->SetFillColor(kRed);

  //  TGraphAsymmErrors *ggr2 = new TGraphAsymmErrors(1000,gx,gy,gxl,gxu,gy2l,gy2u);
  TGraphAsymmErrors *ggr2 = new TGraphAsymmErrors(568,gx,gy,gxl,gxu,gy2l,gy2u);
  ggr2->SetFillColor(kBlue);
  //  ggr2->SetFillColor(kBlue+2);
  //   ggr2->SetFillColor(kPink);

  //  TGraphAsymmErrors *ggr3 = new TGraphAsymmErrors(1000,gx,gy,gxl,gxu,gy3l,gy3u);
  TGraphAsymmErrors *ggr3 = new TGraphAsymmErrors(568,gx,gy,gxl,gxu,gy3l,gy3u);
 ggr3->SetFillColor(kCyan);
 // ggr3->SetFillColor(kBlue);
  //  ggr3->SetFillColor(kPink+1);

 //  TGraphAsymmErrors *ggr4 = new TGraphAsymmErrors(1000,gx,gy,gxl,gxu,gy4l,gy4u);
 TGraphAsymmErrors *ggr4 = new TGraphAsymmErrors(568,gx,gy,gxl,gxu,gy4l,gy4u);
  ggr4->SetFillColor(kBlack);
  ggr4->SetFillStyle(3005);

  //  TGraphAsymmErrors *ggr5 = new TGraphAsymmErrors(1000,gx,gy,gxl,gxu,gy5l,gy5u);
  TGraphAsymmErrors *ggr5 = new TGraphAsymmErrors(568,gx,gy,gxl,gxu,gy5l,gy5u);
  ggr5->SetFillColor(kBlack);
  ggr5->SetFillStyle(3004);

  ggr5->Draw("3");
  ggr4->Draw("3");
  ggr3->Draw("3");
  ggr2->Draw("3");
  ggr1->Draw("3");



  for (int i=1;i<32;i++) {
    G_x2[i]->SetLineWidth(0);
    G_x2[i]->SetLineColor(1);
    G_x2[i]->DrawCopy("SAME");
  }

  G_x2[32]->SetLineWidth(0);
  G_x2[33]->SetLineWidth(0);
  G_x2[32]->SetLineColor(1);
  G_x2[33]->SetLineColor(1);
  G_x2[32]->DrawCopy("HIST SAME");
  G_x2[33]->DrawCopy("HIST SAME");

  char buffer[1000];
  //  char Nucleus[4]="Au";

  TPaveText* Tpl3 = new TPaveText (0.4, 0.12, 0.8, 0.28, "ndc");
  sprintf(buffer,"#sqrt{s} = 200 GeV");
  Tpl3 -> AddText (buffer);
  Tpl3 -> SetTextFont (62);
  Tpl3 -> SetTextSize (0.03);
  Tpl3 -> SetBorderSize (0);
  Tpl3 -> SetFillStyle (0);
  Tpl3 -> Draw();

  TPaveText* Tpl6 = new TPaveText (0.4, 0.10, 0.8, 0.20, "ndc");
  sprintf(buffer," L_{d%s}dt=%d nb^{-1}",Nucleus,pA_lumi);
  Tpl6 -> AddText (buffer);
  Tpl6 -> SetTextFont (62);
  Tpl6 -> SetTextSize (0.03);
  Tpl6 -> SetBorderSize (0);
  Tpl6 -> SetFillStyle (0);
  Tpl6 -> Draw();


  TPaveText* Tpl4 = new TPaveText (0.15, 0.78, 0.55, 0.78, "ndc");
  //  sprintf(buffer,"d+%s",Nucleus," R_{CP}=0-20/60-88"); // we will want this one
  sprintf(buffer,"d+%s",Nucleus);
  //  Tpl4 -> AddText ("p+Au");
  Tpl4 -> AddText (buffer);
  Tpl4 -> SetTextFont (62);
  Tpl4 -> SetTextSize (0.05);
  Tpl4 -> SetBorderSize (0);
  Tpl4 -> SetFillStyle (0);
  Tpl4 -> Draw();


  TPaveText* Tpl5 = new TPaveText (0.15, 0.7, 0.55, 0.7, "ndc");
  //  sprintf(buffer,"d+%s",Nucleus," R_{CP}=0-20/60-88"); // we will want this one
  sprintf(buffer,"R_{MB-P} (0-88%%/60-88%%)");
  if(iscentral)sprintf(buffer,"R_{CP} (0-20%%/60-88%%)");
  //  Tpl4 -> AddText ("p+Au");
  Tpl5 -> AddText (buffer);
  Tpl5 -> SetTextFont (62);
  Tpl5 -> SetTextSize (0.04);
  Tpl5 -> SetBorderSize (0);
  Tpl5 -> SetFillStyle (0);
  Tpl5 -> Draw();

  c3->Update();
  if(iscentral){
  if (printout) c3->Print("plot3wppnice_central_newplt.png");
  if (printout) c3->Print("plot3wppnice_central_newplt.pdf");
  }else{
  if (printout) c3->Print("plot3wppnice_newplt.png");
  if (printout) c3->Print("plot3wppnice_newplt.pdf");
  }
  //====================================================================================
  //====================================================================================


  //-----------------------------------------------------------

  if (printout) {
    cout << "Save to output file" << endl;
    TDirectory *dir = gDirectory;
    TFile*  hfile = new TFile("mpcstat.root","RECREATE","ROOT file with histograms");
    dir->GetList()->Write();
    hfile->Write();
  }

}  

		 
