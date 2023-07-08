// Preprocessor to avoid loading hearder file multiple times
#ifndef ywTreeMaker_def
#define ywTreeMaker_def

// Load C/C++ header files
#include <stdio.h>
#include <iostream>
#include <vector>

// Load ROOT header files
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"

// Load STARLibrary header files
#include "StMaker.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoHelix.h"
#include "StRoot/StPicoEvent/StPicoBbcHit.h"
#include "StRoot/StPicoEvent/StPicoEpdHit.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoTrackCovMatrix.h"

class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;

const Double_t PI = TMath::Pi();
const Int_t Ncentralities =  16;
const Int_t N_track       = 195;

class ywTreeMaker : public StMaker {
  private:
    StPicoDstMaker*       mPicoDstMaker;
    TString                   JobIdName;
    Double_t                    cutTest;
    TFile*                   outputFile;
    Double_t                  Mass_Pion;
    Double_t                  Mass_Kaon;
    Double_t                Mass_Proton;
    Double_t                midrapidity;
    
    TTree*                      FxtTree;
    Int_t               tree_runId;
    Int_t               tree_eventId;
    Float_t             tree_bField;
    Float_t             tree_Vx;
    Float_t             tree_Vy;
    Float_t             tree_Vz;
    UShort_t            tree_centrality;
    UShort_t            tree_tracknumber;
    UShort_t            tree_PID[N_track];
    Short_t             tree_Charge[N_track];
    Float_t             tree_Px[N_track];
    Float_t             tree_Py[N_track];
    Float_t             tree_Pz[N_track];
    UShort_t            tree_BBCadc[32];
    UShort_t            tree_nEPDhits;
    Short_t             tree_EPDid[744];
    Float_t             tree_EPDnMip[744];
    Int_t               tree_EPDadc[744];
    Float_t             tree_ZdcSmdEastHorizontal[8];
    Float_t             tree_ZdcSmdEastVertical[8];
    
    TH1D*                hist_eventcuts;
    TH1D*                hist_trackcuts;
    TH1D*                hist_centralities;
    TH1D*                hist_Vz;
    TH1D*                hist_Vr;
    TH2D*                hist_VxVy;
    TH1D*                hist_trackmult;
    TH1D*                hist_refmult;
    TH1D*                hist_tofmult;
    TH2D*                hist_trackmult_refmult;
    TH2D*                hist_trackmult_tofmult;
    TH2D*                hist_refmult_tofmult;
    
    TH1D*                hist_ndedx;
    TH1D*                hist_nhits;
    TH1D*                hist_ratio;
    TH1D*                hist_DCA;
    TH1D*                hist_pt;
    TH1D*                hist_eta;
    TH1D*                hist_phi;
    TH2D*                hist_dedx;
    TH2D*                hist_beta;
    TH2D*                hist_mass;
    
    TH1D*                hist_mult_proton;
    TH1D*                hist_pt_proton;
    TH1D*                hist_y_proton;
    TH1D*                hist_eta_proton;
    TH1D*                hist_phi_proton;
    TH2D*                hist_pt_y_proton;
    TH2D*                hist_pt_eta_proton;
    TH2D*                hist_dedx_proton;
    TH2D*                hist_beta_proton;
    TH2D*                hist_mass_proton;
    
    TH1D*                hist_mult_piPlus;
    TH1D*                hist_pt_piPlus;
    TH1D*                hist_y_piPlus;
    TH1D*                hist_eta_piPlus;
    TH1D*                hist_phi_piPlus;
    TH2D*                hist_pt_y_piPlus;
    TH2D*                hist_pt_eta_piPlus;
    TH2D*                hist_dedx_piPlus;
    TH2D*                hist_beta_piPlus;
    TH2D*                hist_mass_piPlus;
    
    TH1D*                hist_mult_piMinus;
    TH1D*                hist_pt_piMinus;
    TH1D*                hist_y_piMinus;
    TH1D*                hist_eta_piMinus;
    TH1D*                hist_phi_piMinus;
    TH2D*                hist_pt_y_piMinus;
    TH2D*                hist_pt_eta_piMinus;
    TH2D*                hist_dedx_piMinus;
    TH2D*                hist_beta_piMinus;
    TH2D*                hist_mass_piMinus;
    
    TH1D*                hist_mult_Kplus;
    TH1D*                hist_pt_Kplus;
    TH1D*                hist_y_Kplus;
    TH1D*                hist_eta_Kplus;
    TH1D*                hist_phi_Kplus;
    TH2D*                hist_pt_y_Kplus;
    TH2D*                hist_pt_eta_Kplus;
    TH2D*                hist_dedx_Kplus;
    TH2D*                hist_beta_Kplus;
    TH2D*                hist_mass_Kplus;
    
    TH1D*                hist_mult_Kminus;
    TH1D*                hist_pt_Kminus;
    TH1D*                hist_y_Kminus;
    TH1D*                hist_eta_Kminus;
    TH1D*                hist_phi_Kminus;
    TH2D*                hist_pt_y_Kminus;
    TH2D*                hist_pt_eta_Kminus;
    TH2D*                hist_dedx_Kminus;
    TH2D*                hist_beta_Kminus;
    TH2D*                hist_mass_Kminus;
  public:
    ywTreeMaker(StPicoDstMaker* Maker, TString JobId, Int_t EventsNumber, Double_t inputParameter);
    virtual              ~ywTreeMaker();
    Int_t                        Init();
    Int_t                        Make();
    Int_t                      Finish();
    Bool_t   IsGoodRun(Int_t runnumber);
    ClassDef(ywTreeMaker,1) // Class title
};

// End of preprocessor
#endif
