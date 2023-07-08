#include<stdio.h>
#include <iostream>
#include <fstream>
#include <new>
#include <stdlib.h>
#include <vector>

// load ROOT headers
#include "TChain.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TVector3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TLeaf.h"
#include "TComplex.h"
#include "TRandom3.h"
#include <iostream>

//star library
#include <StEpdGeom.h>


using namespace std;

const int norder=6;
const int n =1;
const int Nbins=1;

const double Mass_Proton = 0.938272;
const double Mass_Pion = 0.139568;
const double Mass_Kaon = 0.493646;
const double  midrapidity = -1.0450222;
const double rapidityLow = 0.0;
const double rapidityHigh = 3.0;
const double ptLow = 0.0;
const double ptHigh = 6.0;

vector<TH1*> hist_list;

TH1D* heta=NULL;
TH1D* hv2_pi_p = NULL;
TH1D* hM=NULL;
TH1D* hm=NULL;
TH2D* hpartId=NULL;
TH2D* hcumA1_1_M = NULL;

TChain* chain = 0;


int cent16=-1;
int nfiles=0;

bool set_chain(char* f_input, int nfiles); 

void init_hist();
void write_out_hist();

void cal();

TComplex get_CorA1_1(TComplex* Qn);
TComplex get_CorB1_1(TComplex* Qn,TComplex* qn);

int main(int argc, char* argv[]){
 if(argc !=4){
    cout<<"Need the input file, number of files to run, the centrality bin"<<endl;
    return 0;
  }

 // int start = atoi(argv[1]);
 // int end = atoi(argv[2]);
  //cent should from 0 to 15
  nfiles = atoi(argv[2]);
  cent16 = atoi(argv[3]);

  set_chain(argv[1],nfiles);
  init_hist();
  cal();
  write_out_hist();
  return 0;
}

TComplex get_CorA1_1(TComplex* Qn){
  if(!Qn){
    cout<<"invalid Qn!"<<endl;
    return TComplex(0,0);
  }
  
  TComplex M = Qn[0];
  TComplex res = (Qn[1]*TComplex::Conjugate(Qn[1])-M)/M/(M-1.);
  //cout<<" Qn[0] "<<Qn[0]<<" Qn[1] "<<Qn[1]<<" res "<<res<<endl;
  return res;
}

TComplex get_CorB1_1(TComplex* Qn,TComplex* qn){
  if(!qn){
    cout<<"invalid qn !"<<endl;
  }

  TComplex M = Qn[0];
  TComplex m = qn[0];

  TComplex res = qn[1]*TComplex::Conjugate(Qn[1])-m;
  res = res/m/(M-1.);

  return res;
}

bool set_chain(char* f_input, int nfiles){
  char infile [ 200 ]; 
      
  ifstream in_txt(f_input);
  string f;
  vector<string> files_list;
  while(getline(in_txt,f)){
//    cout<<"check in file "<<f<<endl;
    files_list.push_back(f);
  }
  cout<<"There are "<<files_list.size()<<" files "<<endl;
  chain = new TChain ( "Autree" );
   
  for(int s=0;s<files_list.size();s++){
    if(s>nfiles) break;
    chain->Add(files_list[s].c_str());
  }
  Long64_t nentries = chain->GetEntries();
  cout<<"Number total entries: "<<nentries<<endl;
}



void init_hist(){
  heta = new TH1D("heta","Pseudorapidity",Nbins,rapidityLow,rapidityHigh);
  heta->GetXaxis()->SetTitle("Eta");
  hist_list.push_back(heta);

  hv2_pi_p = new TH1D("hv2_pi_p","Pi Plus",Nbins,rapidityLow,rapidityHigh);
  hv2_pi_p->GetXaxis()->SetTitle("Eta");
  hist_list.push_back(hv2_pi_p);
  
  hM = new TH1D("hM","total number of tracks",200,-0.5,199.5);
  hM->GetXaxis()->SetTitle("Number of Tracks");
  hist_list.push_back(hM);

  hm = new TH1D("hm","total number of tracks interested",200,-0.5,199.5);
  hm->GetXaxis()->SetTitle("Number of Tracks");
  hist_list.push_back(hm);
  
  hpartId = new TH2D("hpartId","Particle Id mRapidity>0",6,-0.5,5.5,200,-0.5,199.5); 
  hpartId->GetXaxis()->SetTitle("Particle Id");
  hpartId->GetYaxis()->SetTitle("Number of Tracks");
  hist_list.push_back(hpartId);
 
  hcumA1_1_M = new TH2D("hcumA1_1_M","cumA1_1 M",200,-0.5,199.5,200,-1,1);
  hcumA1_1_M->GetXaxis()->SetTitle("M");
  hcumA1_1_M->GetYaxis()->SetTitle("cumA1_1"); 
  hist_list.push_back(hcumA1_1_M);
}

void write_out_hist(){
  cout<<"Write out hist !"<<endl;
  TFile* ofile = new TFile(Form("output_cent%d_nfiles%d.root",cent16,nfiles),"RECREATE");
  for(int i=0;i<hist_list.size();i++) hist_list[i]->Write();
  //ofile->Close();
}

void cal(){
  
  TComplex sumCumA1_1 = TComplex(0,0);
  TComplex sumCumB1_1[Nbins];
  double sumWA = 0;
  double sumWB[Nbins];

  TComplex sumCorA1_1 = 0;
  TComplex sumCorA1 = 0;
  double sumWA1_1 = 0;
  double sumWA1 = 0;

  TComplex sumCorB1_1[Nbins];
  TComplex sumCorB1[Nbins];
  double sumWB1_1[Nbins];
  double sumWB1[Nbins];
  
  for(int ibin=0;ibin<Nbins;ibin++){
    sumWB[ibin]=0;
    sumCumB1_1[ibin]=TComplex(0,0);
    sumCorB1_1[ibin]=TComplex(0,0);
    sumCorB1[ibin]=TComplex(0,0);
    sumWB1_1[ibin]=0;
    sumWB1[ibin]=0;
  }
  
  Long64_t nentries = chain->GetEntries(); 
  
  for ( Long64_t ievent = 0 ; ievent < nentries ; ievent ++ ){
    if ( ievent + 1 > nentries ) break; 
		
    if ( ( ievent + 1 ) % 1000 == 0 ) std::cout << "Analyzed event " << ievent + 1 << std::endl << std::endl; 
    chain -> GetEntry ( ievent );

    //info access
    TLeaf* leaf_runId = (TLeaf*)chain->GetLeaf("runId");
    Int_t runId = (Int_t)leaf_runId->GetValue(0);

    TLeaf* leaf_eventId = (TLeaf*)chain->GetLeaf("eventId");
    Int_t eventId = (Int_t)leaf_eventId->GetValue(0);

    TLeaf* leaf_bField = (TLeaf*)chain->GetLeaf("bField");
    Float_t bField = (Float_t)leaf_bField->GetValue(0);

    TLeaf* leaf_Vx = (TLeaf*)chain->GetLeaf("Vx");
    TLeaf* leaf_Vy = (TLeaf*)chain->GetLeaf("Vy");
    TLeaf* leaf_Vz = (TLeaf*)chain->GetLeaf("Vz");
    TVector3 V((Float_t)leaf_Vx->GetValue(0),(Float_t)leaf_Vy->GetValue(0),(Float_t)leaf_Vz->GetValue(0));

    TLeaf* leaf_centrality = (TLeaf*)chain->GetLeaf("centrality");
    Int_t centrality = (Int_t)leaf_centrality->GetValue(0);

    TLeaf* leaf_tracknumber = (TLeaf*)chain->GetLeaf("tracknumber");
    Int_t numberOfInputTracks = (Int_t)leaf_tracknumber->GetValue(0);

    TLeaf* leaf_PID = (TLeaf*)chain->GetLeaf("PID");
    TLeaf* leaf_Charge = (TLeaf*)chain->GetLeaf("Charge");
    TLeaf* leaf_Px = (TLeaf*)chain->GetLeaf("Px");
    TLeaf* leaf_Py = (TLeaf*)chain->GetLeaf("Py");
    TLeaf* leaf_Pz = (TLeaf*)chain->GetLeaf("Pz");


    TLeaf* leaf_nEPDhits = (TLeaf*)chain->GetLeaf("nEPDhits");
    Int_t nEPDhits = (UShort_t)leaf_nEPDhits->GetValue(0);
    TLeaf* leaf_EPDid = (TLeaf*)chain->GetLeaf("EPDid");
    TLeaf* leaf_EPDnMip = (TLeaf*)chain->GetLeaf("EPDnMip");

    Double_t reaction_plane = 0.0;
    
    if(numberOfInputTracks <= 4 || numberOfInputTracks > 195) continue;    // event cut on track numbers to prevent segmentation violation
    if((centrality != cent16) && (cent16>=0) ) continue;    // select one numberOfInputTracks
    
    //define Q in one event
    TComplex t_Qn[norder];
    TComplex t_qn[Nbins][norder];
    
    for(int s=0;s<norder;s++){
      t_Qn[s]=TComplex(0.,0.);
      for(int t=0;t<Nbins;t++){
        t_qn[t][s] = TComplex(0.,0.);
      } 
    }

    double ct_pid[6]={0,0,0,0,0,0}; 
    
    for(int iter=0;iter<200;iter++){ 
    for ( Int_t itrack = 0 ; itrack < numberOfInputTracks ; itrack ++ ){
      Short_t pid = (Short_t)leaf_PID->GetValue(itrack);
      Short_t charge = (Short_t)leaf_Charge->GetValue(itrack);
      Float_t px = (Float_t)leaf_Px->GetValue(itrack);
      Float_t py = (Float_t)leaf_Py->GetValue(itrack);
      Float_t pz = (Float_t)leaf_Pz->GetValue(itrack);
      Double_t pt = TMath::Sqrt(px*px + py*py);
      Double_t p  = TMath::Sqrt(px*px + py*py + pz*pz);
      TVector3 pMom(px,py,pz);
      Double_t eta = pMom.Eta();
      Double_t phi = pMom.Phi();
      if ( phi < 0.0                ) phi += 2.0 * TMath::Pi ();
      if ( phi > 2.0 * TMath::Pi () ) phi -= 2.0 * TMath::Pi ();
      Double_t energy = 0.0;
      if(pid == 0) energy = TMath::Sqrt(p*p+Mass_Proton*Mass_Proton);
      if(pid == 1) energy = TMath::Sqrt(p*p+Mass_Pion*Mass_Pion);
      if(pid == 2) energy = TMath::Sqrt(p*p+Mass_Kaon*Mass_Kaon);
      Double_t y = 0.5*TMath::Log((energy+pz)/(energy-pz));
      Double_t mRapidity = y - midrapidity;
      
      if(mRapidity>0.0 && charge>0) ct_pid[pid]++;
      if(mRapidity>0.0 && charge<0) ct_pid[pid+3]++;


      if(mRapidity>0){
        //Q calculation for all rapidity
        for(int s=0;s<norder;s++){
	  t_Qn[s]+=TComplex(TMath::Cos(n*s*phi),TMath::Sin(n*s*phi));
	}
      }

      if(pid == 1 && charge > 0 && mRapidity > 0.0){
        Int_t bin = -1;
        bin = heta->FindBin(mRapidity) - 1;
        if(bin >= 0 && bin < Nbins){
	  for(int s=0;s<norder;s++){
	    t_qn[bin][s]+=TComplex(TMath::Cos(n*s*phi),TMath::Sin(n*s*phi));
	  }
	}
      }
    } // loop tracks
    }
    for(int s=0;s<6;++s){
      hpartId->Fill(s*1.0,ct_pid[s]);
    } 

//    cout<<"M: "<<t_Qn[0]<<endl; 
    //one particle correlation
    TComplex corA1;
    TComplex corB1[Nbins];
    
    hM->Fill(t_Qn[0].Re());
    //two particle correlation: phi1-phii2
    TComplex corA1_1;
    TComplex corB1_1[Nbins];
    
    
    if(t_Qn[0].Re()>1){
      corA1 = t_Qn[1]/t_Qn[0]; //one particle
      corA1_1 = get_CorA1_1(t_Qn); //two particle
     
      sumCorA1+=corA1*t_Qn[0].Re();
      sumWA1+=t_Qn[0].Re();
      sumCorA1_1+=corA1_1*t_Qn[0].Re()*(t_Qn[0].Re()-1);
      sumWA1_1+=t_Qn[0].Re()*(t_Qn[0].Re()-1);
    }
    else{
      corA1 = TComplex(0,0);
      corA1_1 = TComplex(0,0);
    }

    for(int ibin=0;ibin<Nbins;ibin++){
      hm->Fill(t_qn[ibin][0].Re());
      if(t_qn[ibin][0].Re()>1){
        corB1[ibin] = t_qn[ibin][1]/t_qn[ibin][0];
        corB1_1[ibin] = get_CorB1_1(t_Qn,t_qn[ibin]);

        sumCorB1[ibin]+=corB1[ibin]*t_qn[ibin][0].Re();
        sumWB1[ibin]+=t_qn[ibin][0].Re();
        sumCorB1_1[ibin]+=corB1_1[ibin]*t_qn[ibin][0].Re();
        sumWB1_1[ibin]+=t_qn[ibin][0].Re()*(t_qn[ibin][0].Re()-1);
      }
      else{
        corB1[ibin] = 0;
        corB1_1[ibin] = 0;
      }
    }
    
    

    //2 particle cumulent
    //should based on large number of particles
    TComplex cumA1_1=TComplex(0,0);
    TComplex cumB1_1[Nbins];
    for(int ibin=0;ibin<Nbins;ibin++){
      cumB1_1[ibin] = TComplex(0,0);
    }
    
    if(t_Qn[0].Re()>1){
    //number of tracks > 1 
      cumA1_1 = corA1_1 - corA1*TComplex::Conjugate(corA1);
      double wA = t_Qn[0].Re()*(t_Qn[0].Re()-1); 
//      cout<<" corA1_1: "<<corA1_1<<" corA1 "<<corA1<<" "
//          <<"cumA1_1 "<<cumA1_1<<" wA "<<wA<<endl;
      sumCumA1_1+=cumA1_1*wA;
      sumWA+=wA;
      hcumA1_1_M->Fill(t_Qn[0].Re(),cumA1_1);
      for(int ibin=0;ibin<Nbins;ibin++){
        if(t_qn[ibin][0].Re()>1){
	//number of tracks > 1
          cumB1_1[ibin] = corB1_1[ibin] -corB1[ibin]*TComplex::Conjugate(corA1); 
          double wB =  t_qn[ibin][0].Re()*(t_Qn[0].Re()-1);
          sumCumB1_1[ibin]+=cumB1_1[ibin]*wB;
          sumWB[ibin]+=wB;
//          cout<<"ibin "<<ibin<<" "
//              <<" corB1_1 "<<corB1_1[ibin]<<" corB1 "<<corB1[ibin]<<" "
//              <<" cumB1_1 "<< cumB1_1[ibin]<<" wB "<<wB<<endl;
	}
      }
    }
  }//loop events
  
//  chain->Reset(); 
 
  cout<<"flow calculation !"<<endl; 
  TComplex avgCorA1_1 = sumCorA1_1/sumWA1_1;
  TComplex avgCorA1 = sumCorA1/sumWA1;
  TComplex fCumA1_1 = avgCorA1_1-avgCorA1*TComplex::Conjugate(avgCorA1);
  cout<<"final cal versoin: "<<fCumA1_1<<endl;
    
  //flow calculation
  TComplex avgCumA1_1 = 0;
  if(sumWA>0) avgCumA1_1 = sumCumA1_1/sumWA;

  for(int ibin=0;ibin<Nbins;ibin++){
    cout<<"ibin "<<ibin<<" "<<sumWB[ibin]<<endl;
    if(sumWB[ibin]>0 && sumWA>0){
      cout<<"calculate bin "<<ibin<<endl;
      TComplex avgCumB1_1=sumCumB1_1[ibin]/sumWB[ibin];
      double v2 = avgCumB1_1.Re()/TMath::Sqrt(avgCumA1_1.Re());
      cout<<"CumA1_1: "<<avgCumA1_1.Re()<<" "
          <<"CumB1_1: "<<avgCumB1_1.Re()<<" "
          <<"v2: "<<v2<<endl;
      hv2_pi_p->SetBinContent(ibin+1,v2);
    }  
  }
  
  //sum then cal
  for(int ibin=0;ibin<Nbins;ibin++){
    cout<<"ibin "<<ibin<<endl;
    if(sumWB1_1[ibin]>0 && sumWB1[ibin]>0){
      TComplex avgCorB1_1 = sumCorB1_1[ibin]/sumWB1_1[ibin];
      TComplex avgCorB1 = sumCorB1[ibin]/sumWB1[ibin];
      
      TComplex fCumB1_1 = avgCorB1_1 - avgCorB1*TComplex::Conjugate(avgCorA1);
      double f_v2 = fCumB1_1.Re()/TMath::Sqrt(fCumA1_1.Re());
      cout<<"fCumA1_1: "<<fCumA1_1.Re()<<" "
          <<"fCumB1_1: "<<fCumB1_1.Re()<<" "
	  <<"f_v2: "<<f_v2<<endl;
    }
  }
}
