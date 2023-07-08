// Load header files
#include "ywTreeMaker.h"

// Class implementation in CINT
ClassImp(ywTreeMaker)

// Constructor
ywTreeMaker::ywTreeMaker(StPicoDstMaker* Maker, TString JobId, Int_t EventsNumber, Double_t inputParameter) : StMaker() {
    // Initialize parameters
    Mass_Proton = 0.938272; Mass_Pion = 0.139568; Mass_Kaon = 0.493646; midrapidity = -1.0450222;
    mPicoDstMaker = Maker; JobIdName = JobId; cutTest = inputParameter;
    JobIdName.Append(".root"); // Name output file by assigned Job ID
}

// Destructor
ywTreeMaker::~ywTreeMaker() {}

Int_t ywTreeMaker::Init() {
    outputFile = new TFile(JobIdName,"recreate");
    // ywTree
    FxtTree = new TTree("Autree","TTree to hold FXT events and tracks");
    FxtTree->Branch("runId",&tree_runId,"runId/I");
    FxtTree->Branch("eventId",&tree_eventId,"eventId/I");
    FxtTree->Branch("bField",&tree_bField,"bField/F");
    FxtTree->Branch("Vx",&tree_Vx,"Vx/F");
    FxtTree->Branch("Vy",&tree_Vy,"Vy/F");
    FxtTree->Branch("Vz",&tree_Vz,"Vz/F");
    FxtTree->Branch("centrality",&tree_centrality,"centrality/s");
    FxtTree->Branch("tracknumber",&tree_tracknumber,"tracknumber/s");
    FxtTree->Branch("PID",tree_PID,"PID[tracknumber]/s");
    FxtTree->Branch("Charge",tree_Charge,"Charge[tracknumber]/S");
    FxtTree->Branch("Px",tree_Px,"Px[tracknumber]/F");
    FxtTree->Branch("Py",tree_Py,"Py[tracknumber]/F");
    FxtTree->Branch("Pz",tree_Pz,"Pz[tracknumber]/F");
    FxtTree->Branch("BBCadc",tree_BBCadc,"BBCadc[32]/s");
    FxtTree->Branch("nEPDhits",&tree_nEPDhits,"nEPDhits/s");
    FxtTree->Branch("EPDid",tree_EPDid,"EPDid[nEPDhits]/S");
    FxtTree->Branch("EPDnMip",tree_EPDnMip,"EPDnMip[nEPDhits]/F");
    FxtTree->Branch("EPDadc",tree_EPDadc,"EPDadc[nEPDhits]/F");
    FxtTree->Branch("ZdcSmdEastHorizontal",tree_ZdcSmdEastHorizontal,"ZdcSmdEastHorizontal[8]/F");
    FxtTree->Branch("ZdcSmdEastVertical",tree_ZdcSmdEastVertical,"ZdcSmdEastVertical[8]/F");
    // QA plots
    hist_eventcuts = new TH1D("eventCuts","# of Selected Events after Each Cut",7,0.5,7.5);
    hist_eventcuts->GetXaxis()->SetBinLabel(1,"Before Cuts");
    hist_eventcuts->GetXaxis()->SetBinLabel(2,"After Run Number Cut");
    hist_eventcuts->GetXaxis()->SetBinLabel(3,"After Trigger ID Cut");
    hist_eventcuts->GetXaxis()->SetBinLabel(4,"After V_{Z} Cut");
    hist_eventcuts->GetXaxis()->SetBinLabel(5,"After V_{R} Cut");
    hist_eventcuts->GetXaxis()->SetBinLabel(6,"After Min-bias Cut");
    hist_eventcuts->GetYaxis()->SetTitle("# of Events");
    
    hist_trackcuts = new TH1D("trackCuts","# of Primary Tracks after Each Cut",3,0.5,3.5);
    hist_trackcuts->GetXaxis()->SetBinLabel(1,"Before Cuts");
    hist_trackcuts->GetXaxis()->SetBinLabel(2,"After Cuts");
    hist_eventcuts->GetYaxis()->SetTitle("# of Tracks");
    
    hist_centralities = new TH1D("hist_centralities","Centrality",Ncentralities,-0.5,Ncentralities+1+0.5);
    hist_centralities->GetXaxis()->SetTitle("Centrality bin");
    hist_centralities->GetYaxis()->SetTitle("# of events");
    
    hist_Vz = new TH1D("hist_Vz","V_{Z} [cm]",500,180.0,220.0);
    hist_Vz->GetXaxis()->SetTitle("V_{Z} [cm]");
    hist_Vz->GetYaxis()->SetTitle("# of events");
    
    hist_Vr = new TH1D("hist_Vr","V_{R} [cm]",500,0.0,20.0);
    hist_Vr->GetXaxis()->SetTitle("V_{R} [cm]");
    hist_Vr->GetYaxis()->SetTitle("# of events");
    
    hist_VxVy = new TH2D("hist_VxVy","V_{Y} [cm] vs. V_{X} [cm]",500,-5.0,5.0,500,-5.0,5.0);
    hist_VxVy->GetXaxis()->SetTitle("V_{X} [cm]");
    hist_VxVy->GetYaxis()->SetTitle("V_{Y} [cm]");
    
    hist_refmult = new TH1D("hist_refmult","Reference multiplicity",1001,-0.5,1000.5);
    hist_refmult->GetXaxis()->SetTitle("RefMult");
    hist_refmult->GetXaxis()->SetTitle("# of events");
    
    hist_tofmult = new TH1D("hist_tofmult","TOF multiplicity",1001,-0.5,1000.5);
    hist_tofmult->GetXaxis()->SetTitle("TofMult");
    hist_tofmult->GetXaxis()->SetTitle("# of events");
    
    hist_trackmult = new TH1D("hist_trackmult","Actual track multiplicity",1501,-0.5,1500.5);
    hist_trackmult->GetXaxis()->SetTitle("TrackMult");
    hist_trackmult->GetXaxis()->SetTitle("# of events");
    
    hist_trackmult_refmult = new TH2D("hist_trackmult_refmult","Actual track multiplicity vs. RefMult",1501,-0.5,1500.5,1001,-0.5,1000.5);
    hist_trackmult_refmult->GetXaxis()->SetTitle("TrackMult");
    hist_trackmult_refmult->GetXaxis()->SetTitle("RefMult");
    
    hist_trackmult_tofmult = new TH2D("hist_trackmult_tofmult","Actual track multiplicity vs. TofMult",1501,-0.5,1500.5,1001,-0.5,1000.5);
    hist_trackmult_tofmult->GetXaxis()->SetTitle("TrackMult");
    hist_trackmult_tofmult->GetXaxis()->SetTitle("TofMult");
    
    hist_refmult_tofmult = new TH2D("hist_refmult_tofmult","RefMult vs. TofMult",1001,-0.5,1000.5,1001,-0.5,1000.5);
    hist_refmult_tofmult->GetXaxis()->SetTitle("TofMult");
    hist_refmult_tofmult->GetYaxis()->SetTitle("RefMult");
    
    hist_ndedx = new TH1D("hist_ndedx","ndedx",101,-0.5,100.5);
    hist_ndedx->GetXaxis()->SetTitle("ndedx");
    hist_ndedx->GetYaxis()->SetTitle("# of tracks");
    
    hist_nhits = new TH1D("hist_nhits","nhits",101,-0.5,100.5);
    hist_nhits->GetXaxis()->SetTitle("nhits");
    hist_nhits->GetYaxis()->SetTitle("# of tracks");
    
    hist_ratio = new TH1D("hist_ratio","nhitsFit/nhitsPoss",200,0.0,2.0);
    hist_ratio->GetXaxis()->SetTitle("nhitsFit/nhitsPoss");
    hist_ratio->GetYaxis()->SetTitle("# of tracks");
    
    hist_DCA = new TH1D("hist_DCA","DCA [cm]",100,0.0,10.0);
    hist_DCA->GetXaxis()->SetTitle("DCA [cm]");
    hist_DCA->GetYaxis()->SetTitle("# of tracks");
    
    hist_pt = new TH1D("hist_pt","p_{T} [GeV/c]",1000,0.0,10.0);
    hist_pt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hist_pt->GetYaxis()->SetTitle("# of tracks");
    
    hist_eta = new TH1D("hist_eta","#eta",500,-5.0,5.0);
    hist_eta->GetXaxis()->SetTitle("#eta");
    hist_eta->GetYaxis()->SetTitle("# of tracks");
    
    hist_phi = new TH1D("hist_phi","#phi [Radian]",1000,-1.5*PI,1.5*PI);
    hist_phi->GetXaxis()->SetTitle("#phi [Radian]");
    hist_phi->GetYaxis()->SetTitle("# of tracks");
    
    hist_dedx = new TH2D("hist_dedx","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
    hist_dedx->GetXaxis()->SetTitle("q*|p| (GeV/c)");
    hist_dedx->GetYaxis()->SetTitle("dE/dx (keV/cm)");
    
    hist_beta = new TH2D("hist_beta","1/#beta vs q*|p|",1000,-5.0,5.0,500,0.0,5.0);
    hist_beta->GetXaxis()->SetTitle("q*|p| (GeV/c)");
    hist_beta->GetYaxis()->SetTitle("1/#beta");
    
    hist_mass = new TH2D("hist_mass","m^{2} vs q*|p|",1000,-5.0,5.0,1000,-0.6,4.0);
    hist_mass->GetXaxis()->SetTitle("q*|p| (GeV/c)");
    hist_mass->GetYaxis()->SetTitle("m^{2} (GeV/c^{2})^{2}");
    
    hist_pt_proton = new TH1D("hist_pt_proton","p_{T} [GeV/c]",1000,0.0,5.0);
    hist_pt_proton->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hist_pt_proton->GetYaxis()->SetTitle("# of tracks");
    
    hist_eta_proton = new TH1D("hist_eta_proton","#eta",500,-5.0,5.0);
    hist_eta_proton->GetXaxis()->SetTitle("#eta");
    hist_eta_proton->GetYaxis()->SetTitle("# of tracks");
    
    hist_y_proton = new TH1D("hist_y_proton","y",500,-5.0,5.0);
    hist_y_proton->GetXaxis()->SetTitle("Rapidity y");
    hist_y_proton->GetYaxis()->SetTitle("# of tracks");
    
    hist_pt_y_proton = new TH2D("hist_pt_y_proton","p_{T} [GeV/c] vs. y",700,-3.0,3.0,500,0.0,3.5);
    hist_pt_y_proton->GetXaxis()->SetTitle("y");
    hist_pt_y_proton->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    
    hist_pt_eta_proton = new TH2D("hist_pt_eta_proton","p_{T} [GeV/c] vs. #eta",700,-3.0,3.0,500,0.0,3.5);
    hist_pt_eta_proton->GetXaxis()->SetTitle("#eta");
    hist_pt_eta_proton->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    
    hist_phi_proton = new TH1D("hist_phi_proton","#phi [Radian]",1000,-1.5*PI,1.5*PI);
    hist_phi_proton->GetXaxis()->SetTitle("#phi [Radian]");
    hist_phi_proton->GetYaxis()->SetTitle("# of tracks");
    
    hist_mult_proton = new TH1D("hist_mult_proton","Proton track multiplicity",1001,-0.5,1000.5);
    hist_mult_proton->GetXaxis()->SetTitle("K^{#plus} #");
    hist_mult_proton->GetXaxis()->SetTitle("# of events");
    
    hist_dedx_proton = new TH2D("hist_dedx_proton","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
    hist_dedx_proton->GetXaxis()->SetTitle("q*|p| (GeV/c)");
    hist_dedx_proton->GetYaxis()->SetTitle("dE/dx (keV/cm)");
    
    hist_beta_proton = new TH2D("hist_beta_proton","1/#beta vs q*|p|",1000,-5.0,5.0,500,0.0,5.0);
    hist_beta_proton->GetXaxis()->SetTitle("q*|p| (GeV/c)");
    hist_beta_proton->GetYaxis()->SetTitle("1/#beta");
    
    hist_mass_proton = new TH2D("hist_mass_proton","m^{2} vs q*|p|",1000,-5.0,5.0,1000,-0.6,4.0);
    hist_mass_proton->GetXaxis()->SetTitle("q*|p| (GeV/c)");
    hist_mass_proton->GetYaxis()->SetTitle("m^{2} (GeV/c^{2})^{2}");
    
    hist_pt_piPlus = new TH1D("hist_pt_piPlus","p_{T} [GeV/c]",1000,0.0,5.0);
    hist_pt_piPlus->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hist_pt_piPlus->GetYaxis()->SetTitle("# of tracks");
    
    hist_eta_piPlus = new TH1D("hist_eta_piPlus","#eta",500,-5.0,5.0);
    hist_eta_piPlus->GetXaxis()->SetTitle("#eta");
    hist_eta_piPlus->GetYaxis()->SetTitle("# of tracks");
    
    hist_y_piPlus = new TH1D("hist_y_piPlus","y",500,-5.0,5.0);
    hist_y_piPlus->GetXaxis()->SetTitle("Rapidity y");
    hist_y_piPlus->GetYaxis()->SetTitle("# of tracks");
    
    hist_pt_y_piPlus = new TH2D("hist_pt_y_piPlus","p_{T} [GeV/c] vs. y",700,-3.0,3.0,500,0.0,3.5);
    hist_pt_y_piPlus->GetXaxis()->SetTitle("y");
    hist_pt_y_piPlus->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    
    hist_pt_eta_piPlus = new TH2D("hist_pt_eta_piPlus","p_{T} [GeV/c] vs. #eta",700,-3.0,3.0,500,0.0,3.5);
    hist_pt_eta_piPlus->GetXaxis()->SetTitle("#eta");
    hist_pt_eta_piPlus->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    
    hist_phi_piPlus = new TH1D("hist_phi_piPlus","#phi [Radian]",1000,-1.5*PI,1.5*PI);
    hist_phi_piPlus->GetXaxis()->SetTitle("#phi [Radian]");
    hist_phi_piPlus->GetYaxis()->SetTitle("# of tracks");
    
    hist_mult_piPlus = new TH1D("hist_mult_piPlus","#pi^{#plus} track multiplicity",1001,-0.5,1000.5);
    hist_mult_piPlus->GetXaxis()->SetTitle("K^{#plus} #");
    hist_mult_piPlus->GetXaxis()->SetTitle("# of events");
    
    hist_dedx_piPlus = new TH2D("hist_dedx_piPlus","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
    hist_dedx_piPlus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
    hist_dedx_piPlus->GetYaxis()->SetTitle("dE/dx (keV/cm)");
    
    hist_beta_piPlus = new TH2D("hist_beta_piPlus","1/#beta vs q*|p|",1000,-5.0,5.0,500,0.0,5.0);
    hist_beta_piPlus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
    hist_beta_piPlus->GetYaxis()->SetTitle("1/#beta");
    
    hist_mass_piPlus = new TH2D("hist_mass_piPlus","m^{2} vs q*|p|",1000,-5.0,5.0,1000,-0.6,4.0);
    hist_mass_piPlus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
    hist_mass_piPlus->GetYaxis()->SetTitle("m^{2} (GeV/c^{2})^{2}");
    
    hist_pt_piMinus = new TH1D("hist_pt_piMinus","p_{T} [GeV/c]",1000,0.0,5.0);
    hist_pt_piMinus->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hist_pt_piMinus->GetYaxis()->SetTitle("# of tracks");
    
    hist_eta_piMinus = new TH1D("hist_eta_piMinus","#eta",500,-5.0,5.0);
    hist_eta_piMinus->GetXaxis()->SetTitle("#eta");
    hist_eta_piMinus->GetYaxis()->SetTitle("# of tracks");
    
    hist_y_piMinus = new TH1D("hist_y_piMinus","y",500,-5.0,5.0);
    hist_y_piMinus->GetXaxis()->SetTitle("Rapidity y");
    hist_y_piMinus->GetYaxis()->SetTitle("# of tracks");
    
    hist_pt_y_piMinus = new TH2D("hist_pt_y_piMinus","p_{T} [GeV/c] vs. y",700,-3.0,3.0,500,0.0,3.5);
    hist_pt_y_piMinus->GetXaxis()->SetTitle("y");
    hist_pt_y_piMinus->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    
    hist_pt_eta_piMinus = new TH2D("hist_pt_eta_piMinus","p_{T} [GeV/c] vs. #eta",700,-3.0,3.0,500,0.0,3.5);
    hist_pt_eta_piMinus->GetXaxis()->SetTitle("#eta");
    hist_pt_eta_piMinus->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    
    hist_phi_piMinus = new TH1D("hist_phi_piMinus","#phi [Radian]",1000,-1.5*PI,1.5*PI);
    hist_phi_piMinus->GetXaxis()->SetTitle("#phi [Radian]");
    hist_phi_piMinus->GetYaxis()->SetTitle("# of tracks");
    
    hist_mult_piMinus = new TH1D("hist_mult_piMinus","#pi^{#minus} track multiplicity",1001,-0.5,1000.5);
    hist_mult_piMinus->GetXaxis()->SetTitle("K^{#minus} #");
    hist_mult_piMinus->GetXaxis()->SetTitle("# of events");
    
    hist_dedx_piMinus = new TH2D("hist_dedx_piMinus","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
    hist_dedx_piMinus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
    hist_dedx_piMinus->GetYaxis()->SetTitle("dE/dx (keV/cm)");
    
    hist_beta_piMinus = new TH2D("hist_beta_piMinus","1/#beta vs q*|p|",1000,-5.0,5.0,500,0.0,5.0);
    hist_beta_piMinus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
    hist_beta_piMinus->GetYaxis()->SetTitle("1/#beta");
    
    hist_mass_piMinus = new TH2D("hist_mass_piMinus","m^{2} vs q*|p|",1000,-5.0,5.0,1000,-0.6,4.0);
    hist_mass_piMinus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
    hist_mass_piMinus->GetYaxis()->SetTitle("m^{2} (GeV/c^{2})^{2}");
    
    hist_pt_Kplus = new TH1D("hist_pt_Kplus","p_{T} [GeV/c]",1000,0.0,5.0);
    hist_pt_Kplus->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hist_pt_Kplus->GetYaxis()->SetTitle("# of tracks");
    
    hist_eta_Kplus = new TH1D("hist_eta_Kplus","#eta",500,-5.0,5.0);
    hist_eta_Kplus->GetXaxis()->SetTitle("#eta");
    hist_eta_Kplus->GetYaxis()->SetTitle("# of tracks");
    
    hist_y_Kplus = new TH1D("hist_y_Kplus","y",500,-5.0,5.0);
    hist_y_Kplus->GetXaxis()->SetTitle("Rapidity y");
    hist_y_Kplus->GetYaxis()->SetTitle("# of tracks");
    
    hist_pt_y_Kplus = new TH2D("hist_pt_y_Kplus","p_{T} [GeV/c] vs. y",700,-3.0,3.0,500,0.0,3.5);
    hist_pt_y_Kplus->GetXaxis()->SetTitle("y");
    hist_pt_y_Kplus->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    
    hist_pt_eta_Kplus = new TH2D("hist_pt_eta_Kplus","p_{T} [GeV/c] vs. #eta",700,-3.0,3.0,500,0.0,3.5);
    hist_pt_eta_Kplus->GetXaxis()->SetTitle("#eta");
    hist_pt_eta_Kplus->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    
    hist_phi_Kplus = new TH1D("hist_phi_Kplus","#phi [Radian]",1000,-1.5*PI,1.5*PI);
    hist_phi_Kplus->GetXaxis()->SetTitle("#phi [Radian]");
    hist_phi_Kplus->GetYaxis()->SetTitle("# of tracks");
    
    hist_mult_Kplus = new TH1D("hist_mult_Kplus","K^{#plus} track multiplicity",1001,-0.5,1000.5);
    hist_mult_Kplus->GetXaxis()->SetTitle("K^{#plus} #");
    hist_mult_Kplus->GetXaxis()->SetTitle("# of events");
    
    hist_dedx_Kplus = new TH2D("hist_dedx_Kplus","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
    hist_dedx_Kplus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
    hist_dedx_Kplus->GetYaxis()->SetTitle("dE/dx (keV/cm)");
    
    hist_beta_Kplus = new TH2D("hist_beta_Kplus","1/#beta vs q*|p|",1000,-5.0,5.0,500,0.0,5.0);
    hist_beta_Kplus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
    hist_beta_Kplus->GetYaxis()->SetTitle("1/#beta");
    
    hist_mass_Kplus = new TH2D("hist_mass_Kplus","m^{2} vs q*|p|",1000,-5.0,5.0,1000,-0.6,4.0);
    hist_mass_Kplus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
    hist_mass_Kplus->GetYaxis()->SetTitle("m^{2} (GeV/c^{2})^{2}");
    
    hist_pt_Kminus = new TH1D("hist_pt_Kminus","p_{T} [GeV/c]",1000,0.0,5.0);
    hist_pt_Kminus->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hist_pt_Kminus->GetYaxis()->SetTitle("# of tracks");
    
    hist_eta_Kminus = new TH1D("hist_eta_Kminus","#eta",500,-5.0,5.0);
    hist_eta_Kminus->GetXaxis()->SetTitle("#eta");
    hist_eta_Kminus->GetYaxis()->SetTitle("# of tracks");
    
    hist_y_Kminus = new TH1D("hist_y_Kminus","y",500,-5.0,5.0);
    hist_y_Kminus->GetXaxis()->SetTitle("Rapidity y");
    hist_y_Kminus->GetYaxis()->SetTitle("# of tracks");
    
    hist_pt_y_Kminus = new TH2D("hist_pt_y_Kminus","p_{T} [GeV/c] vs. y",700,-3.0,3.0,500,0.0,3.5);
    hist_pt_y_Kminus->GetXaxis()->SetTitle("y");
    hist_pt_y_Kminus->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    
    hist_pt_eta_Kminus = new TH2D("hist_pt_eta_Kminus","p_{T} [GeV/c] vs. #eta",700,-3.0,3.0,500,0.0,3.5);
    hist_pt_eta_Kminus->GetXaxis()->SetTitle("#eta");
    hist_pt_eta_Kminus->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    
    hist_phi_Kminus = new TH1D("hist_phi_Kminus","#phi [Radian]",1000,-1.5*PI,1.5*PI);
    hist_phi_Kminus->GetXaxis()->SetTitle("#phi [Radian]");
    hist_phi_Kminus->GetYaxis()->SetTitle("# of tracks");
    
    hist_mult_Kminus = new TH1D("hist_mult_Kminus","K^{#minus} track multiplicity",1001,-0.5,1000.5);
    hist_mult_Kminus->GetXaxis()->SetTitle("K^{#minus} #");
    hist_mult_Kminus->GetXaxis()->SetTitle("# of events");
    
    hist_dedx_Kminus = new TH2D("hist_dedx_Kminus","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
    hist_dedx_Kminus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
    hist_dedx_Kminus->GetYaxis()->SetTitle("dE/dx (keV/cm)");
    
    hist_beta_Kminus = new TH2D("hist_beta_Kminus","1/#beta vs q*|p|",1000,-5.0,5.0,500,0.0,5.0);
    hist_beta_Kminus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
    hist_beta_Kminus->GetYaxis()->SetTitle("1/#beta");
    
    hist_mass_Kminus = new TH2D("hist_mass_Kminus","m^{2} vs q*|p|",1000,-5.0,5.0,1000,-0.6,4.0);
    hist_mass_Kminus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
    hist_mass_Kminus->GetYaxis()->SetTitle("m^{2} (GeV/c^{2})^{2}");
    return kStOK;
}

Int_t ywTreeMaker::Finish() {
    outputFile->Write();
    return kStOK;
}

Int_t ywTreeMaker::Make() {
    hist_eventcuts->Fill(1); // Count # of events before any cuts
    StPicoEvent* picoEvent = mPicoDstMaker->picoDst()->event(); // Get Event pointer
    if( picoEvent ) { // Ensure event pointer is not NULL
        if( IsGoodRun(picoEvent->runId()) ) { // Select only good runs
            hist_eventcuts->Fill(2); // Count # of good run events
            Int_t triggerFlag = 0;
            std::vector<unsigned int> triggersCollection = picoEvent->triggerIds();
            for(unsigned int i=0; i<triggersCollection.size(); i++) {
                if(triggersCollection.at(i) == 620052) triggerFlag++;
            }
            if( triggerFlag > 0 ) { // Trigger id cut
                hist_eventcuts->Fill(3); // Count # of events after trigger id cut
                TVector3 V = picoEvent->primaryVertex();
                hist_Vz->Fill(V.Z()); hist_VxVy->Fill(V.X(),V.Y());
                if( V.Z() > 198.0 && V.Z() < 202.0 ) { // Vz cut
                    hist_eventcuts->Fill(4);    // Count # of events after Vz cut
                    hist_Vr->Fill(TMath::Sqrt(TMath::Power(V.X(),2.0) + TMath::Power(V.Y()+2.0,2.0)));
                    if( TMath::Power(V.X(),2.0) + TMath::Power(V.Y()+2.0,2.0) < TMath::Power(2.0,2.0) ) { // Vr cut
                        hist_eventcuts->Fill(5); // Count # of events after Vr cut
                        Int_t tracksNumber = mPicoDstMaker->picoDst()->numberOfTracks(); // Get the number of Primary Tracks
                        Int_t Ntracks = 0, N_PiPlus = 0, N_PiMinus = 0, N_Proton = 0, N_Kplus = 0, N_Kminus = 0;
                        // 1st Primary Tracks loop to get centrality starts
                        for( Int_t itr = 0; itr < tracksNumber; itr++ ) {
                            StPicoTrack* track = mPicoDstMaker->picoDst()->track(itr); // Get Track pointer
                            if(!track->isPrimary()) continue; // Get Primary Tracks
                            Ntracks++;
                            hist_trackcuts->Fill(1);
                            hist_nhits->Fill(track->nHits());
                            hist_ndedx->Fill(track->nHitsDedx());
                            hist_ratio->Fill((Double_t)track->nHitsFit()/track->nHitsPoss());
                            hist_DCA->Fill(track->gDCA(V).Mag());
                        } // 1st Primary tracks loop ends
                        hist_refmult->Fill(picoEvent->refMultHalfEast());
                        hist_trackmult->Fill(Ntracks);
                        hist_tofmult->Fill(picoEvent->nBTOFMatch());
                        hist_trackmult_refmult->Fill(Ntracks,picoEvent->refMultHalfEast());
                        hist_trackmult_tofmult->Fill(Ntracks,picoEvent->nBTOFMatch());
                        hist_refmult_tofmult->Fill(picoEvent->refMultHalfEast(),picoEvent->nBTOFMatch());
                        // Get centrality
                        Int_t centrality = Ncentralities;
                        if( Ntracks >=   3 && Ntracks <=   4 ) centrality =  0; // 75-80
                        if( Ntracks >=   5 && Ntracks <=   6 ) centrality =  1; // 70-75
                        if( Ntracks >=   7 && Ntracks <=   9 ) centrality =  2; // 65-70
                        if( Ntracks >=  10 && Ntracks <=  13 ) centrality =  3; // 60-65
                        if( Ntracks >=  14 && Ntracks <=  17 ) centrality =  4; // 55-60
                        if( Ntracks >=  18 && Ntracks <=  23 ) centrality =  5; // 50-55
                        if( Ntracks >=  24 && Ntracks <=  29 ) centrality =  6; // 45-50
                        if( Ntracks >=  30 && Ntracks <=  37 ) centrality =  7; // 40-45
                        if( Ntracks >=  38 && Ntracks <=  46 ) centrality =  8; // 35-40
                        if( Ntracks >=  47 && Ntracks <=  56 ) centrality =  9; // 30-35
                        if( Ntracks >=  57 && Ntracks <=  68 ) centrality = 10; // 25-30
                        if( Ntracks >=  69 && Ntracks <=  82 ) centrality = 11; // 20-25
                        if( Ntracks >=  83 && Ntracks <=  98 ) centrality = 12; // 15-20
                        if( Ntracks >=  99 && Ntracks <= 117 ) centrality = 13; // 10-15
                        if( Ntracks >= 118 && Ntracks <= 140 ) centrality = 14; // 5-10
                        if( Ntracks >= 141 && Ntracks <= 195 ) centrality = 15; // 0-5
                        // Select only min-bias events
                        if( centrality != Ncentralities ) {
                            hist_eventcuts->Fill(6); // Count # of min-bias events after all cuts
                            hist_centralities->Fill(centrality+1);
                            // Prepare ywTree parameters
                            tree_Vx = V.X(); tree_Vy = V.Y(); tree_Vz = V.Z();
                            tree_runId = picoEvent->runId(); tree_eventId = picoEvent->eventId(); tree_bField = picoEvent->bField();
                            tree_centrality = centrality; /*tree_tracknumber = Ntracks;*/ Int_t realTrackIndex = 0;
                            for(Int_t i = 0;i < 16;i++) {
                                tree_BBCadc[i] = picoEvent->bbcAdcEast(i); // BBC east: 0-15
                                tree_BBCadc[i+16] = picoEvent->bbcAdcWest(i); // BBC west: 16-31
                            }
                            for(Int_t i = 0;i < 8;i++) {
                                tree_ZdcSmdEastHorizontal[i] = picoEvent->ZdcSmdEastHorizontal(i);
                                tree_ZdcSmdEastVertical[i] = picoEvent->ZdcSmdEastVertical(i);
                            }
                            tree_nEPDhits = mPicoDstMaker->picoDst()->numberOfEpdHits();
                            if(tree_nEPDhits) {
                                for(Int_t hit=0;hit<tree_nEPDhits;hit++) {
                                    StPicoEpdHit *pEpdHit = mPicoDstMaker->picoDst()->epdHit(hit);
                                    if(!pEpdHit) continue;
                                    if(!pEpdHit->isGood()) continue;
                                    tree_EPDid[hit]   = pEpdHit->id();
                                    tree_EPDnMip[hit] = pEpdHit->nMIP();
                                    tree_EPDadc[hit]  = pEpdHit->adc();
                                }
                            }
                            // 2nd Loop through tracks to fill histograms
                            for( Int_t itr = 0; itr < tracksNumber; itr++ ) {
                                StPicoTrack* track = mPicoDstMaker->picoDst()->track(itr); // Get Track pointer
                                if(!track->isPrimary()) continue;
                                if( track->nHitsDedx() >= 10
                                    && (Double_t)track->nHitsFit()/track->nHitsPoss() >= 0.51
                                    && track->pPt() > 0.2 //&& track->pt() < 2.0
                                    && track->gDCA(V).Mag() <= 3.0
                                    && track->nHits() >= 15
                                   ) { // Track cuts
                                    hist_trackcuts->Fill(2);
                                    Double_t pt = track->pPt(), pz = track->pMom().Z(), eta = track->pMom().Eta(), phi = track->pMom().Phi(),
                                    trackP = track->pMom().Mag();
                                    Short_t charge = track->charge();
                                    hist_phi->Fill(phi); hist_eta->Fill(eta); hist_pt->Fill(pt);
                                    hist_dedx->Fill(charge*trackP,track->dEdx());
                                    Double_t Beta = -999.0, mass2 = 0.0;
                                    Int_t trackTofIndex = track->bTofPidTraitsIndex();
                                    if(trackTofIndex >= 0) Beta = mPicoDstMaker->picoDst()->btofPidTraits(trackTofIndex)->btofBeta();
                                    if( Beta != -999.0 /*&& TMath::Abs(Beta) < 1.0*/ ) {
                                        mass2 = trackP*trackP * ( ( 1.0 / ( Beta*Beta ) ) - 1.0 );
                                        hist_beta->Fill(charge*trackP,1.0/Beta);
                                        hist_mass->Fill(charge*trackP,mass2);
                                    }
                                    Double_t energy_proton = TMath::Sqrt(trackP*trackP+Mass_Proton*Mass_Proton);
                                    Double_t energy_pion   = TMath::Sqrt(trackP*trackP+Mass_Pion  *Mass_Pion);
                                    Double_t energy_kaon   = TMath::Sqrt(trackP*trackP+Mass_Kaon  *Mass_Kaon);
                                    Double_t rap_Proton = 0.5*TMath::Log((energy_proton+pz)/(energy_proton-pz));
                                    Double_t rap_Pion   = 0.5*TMath::Log((energy_pion  +pz)/(energy_pion  -pz));
                                    Double_t rap_Kaon   = 0.5*TMath::Log((energy_kaon  +pz)/(energy_kaon  -pz));
                                    Double_t mRapidity = 999.0;
                                    Bool_t isProton = kFALSE, isPion = kFALSE, isKaon = kFALSE;
                                    if( TMath::Abs( track->nSigmaProton() ) < 2.0 && ( Beta != -999.0 && mass2 > 0.6 && mass2 < 1.1 )
                                        && track->charge() > 0
                                       ) { // PID Proton
                                        isProton = kTRUE; mRapidity = rap_Proton; N_Proton++;
                                        hist_eta_proton->Fill(eta);
                                        hist_phi_proton->Fill(phi);
                                        hist_pt_proton->Fill(pt);
                                        hist_y_proton->Fill(mRapidity);
                                        hist_pt_y_proton->Fill(mRapidity,pt);
                                        hist_pt_eta_proton->Fill(eta,pt);
                                        hist_dedx_proton->Fill(charge*trackP,track->dEdx());
                                        hist_beta_proton->Fill(charge*trackP,1.0/Beta);
                                        hist_mass_proton->Fill(charge*trackP,mass2);
                                    } // PID Proton ends
                                    if( TMath::Abs( track->nSigmaPion() ) < 2.0
                                       && (( Beta != -999.0 && mass2 > -0.1 && mass2 < 0.12
                                             //&& ( 0.4 <= trackP || ( ( 0.0 < trackP && trackP < 0.4 ) && ( fabs(mass2) > 0.002 ) ) ) // for piM
                                             && ( 0.25 <= trackP || ( ( 0.0 < trackP && trackP < 0.25 ) && ( fabs(mass2) > 0.00476 ) ) ) // for piP
                                            )
                                           )
                                       && (!isProton)
                                       ) { // PID Pions
                                        isPion = true; mRapidity = rap_Pion;
                                        if(charge > 0) {
                                            N_PiPlus++;
                                            hist_eta_piPlus->Fill(eta);
                                            hist_phi_piPlus->Fill(phi);
                                            hist_pt_piPlus->Fill(pt);
                                            hist_y_piPlus->Fill(mRapidity);
                                            hist_pt_y_piPlus->Fill(mRapidity,pt);
                                            hist_pt_eta_piPlus->Fill(eta,pt);
                                            hist_dedx_piPlus->Fill(charge*trackP,track->dEdx());
                                            hist_beta_piPlus->Fill(charge*trackP,1.0/Beta);
                                            hist_mass_piPlus->Fill(charge*trackP,mass2);
                                        }
                                        if(charge < 0) {
                                            N_PiMinus++;
                                            hist_eta_piMinus->Fill(eta);
                                            hist_phi_piMinus->Fill(phi);
                                            hist_pt_piMinus->Fill(pt);
                                            hist_y_piMinus->Fill(mRapidity);
                                            hist_pt_y_piMinus->Fill(mRapidity,pt);
                                            hist_pt_eta_piMinus->Fill(eta,pt);
                                            hist_dedx_piMinus->Fill(charge*trackP,track->dEdx());
                                            hist_beta_piMinus->Fill(charge*trackP,1.0/Beta);
                                            hist_mass_piMinus->Fill(charge*trackP,mass2);
                                        }
                                    } // PID Pion ends
                                    if( TMath::Abs( track->nSigmaKaon() ) < 2.0 && ( Beta != -999.0 && mass2 > 0.21 && mass2 < 0.28 )
                                       && (!( isPion || isProton ))
                                       ) { // PID Kaons
                                        isKaon = true;
                                        mRapidity = rap_Kaon;
                                        if(charge > 0) {
                                            N_Kplus++;
                                            hist_eta_Kplus->Fill(eta);
                                            hist_phi_Kplus->Fill(phi);
                                            hist_pt_Kplus->Fill(pt);
                                            hist_y_Kplus->Fill(mRapidity);
                                            hist_pt_y_Kplus->Fill(mRapidity,pt);
                                            hist_pt_eta_Kplus->Fill(eta,pt);
                                            hist_dedx_Kplus->Fill(charge*trackP,track->dEdx());
                                            hist_beta_Kplus->Fill(charge*trackP,1.0/Beta);
                                            hist_mass_Kplus->Fill(charge*trackP,mass2);
                                        }
                                        if(charge < 0) {
                                            N_Kminus++;
                                            hist_eta_Kminus->Fill(eta);
                                            hist_phi_Kminus->Fill(phi);
                                            hist_pt_Kminus->Fill(pt);
                                            hist_y_Kminus->Fill(mRapidity);
                                            hist_pt_y_Kminus->Fill(mRapidity,pt);
                                            hist_pt_eta_Kminus->Fill(eta,pt);
                                            hist_dedx_Kminus->Fill(charge*trackP,track->dEdx());
                                            hist_beta_Kminus->Fill(charge*trackP,1.0/Beta);
                                            hist_mass_Kminus->Fill(charge*trackP,mass2);
                                        }
                                    } // PID Kaon ends
                                    if(isProton || isPion || isKaon) {
                                        tree_Charge[realTrackIndex] = charge;
                                        if(isProton) tree_PID[realTrackIndex] = 0;
                                        if(isPion) tree_PID[realTrackIndex] = 1;
                                        if(isKaon) tree_PID[realTrackIndex] = 2;
                                        tree_Px[realTrackIndex] = track->pMom().X();
                                        tree_Py[realTrackIndex] = track->pMom().Y();
                                        tree_Pz[realTrackIndex] = pz;
                                        realTrackIndex += 1;
                                    } // Tree loop ends
                                } // track cuts end
                            } // 2nd track loop ends
                            hist_mult_proton->Fill(N_Proton);
                            hist_mult_piPlus->Fill(N_PiPlus);
                            hist_mult_piMinus->Fill(N_PiMinus);
                            hist_mult_Kplus->Fill(N_Kplus);
                            hist_mult_Kminus->Fill(N_Kminus);
                            tree_tracknumber = N_Proton + N_PiPlus + N_PiMinus + N_Kplus + N_Kminus;
                            FxtTree->Fill();
                        } // Min-bias event selection ends
                    } // Vr cut ends
                } // Vz cut ends
            } // Trigger id cut ends
        } // Good run cut ends
    } // Event pointer not NULL cut ends
    return kStOK;
}

Bool_t ywTreeMaker::IsGoodRun(Int_t runnumber){
    Bool_t goodRunFlag = kFALSE;
    Int_t goodRunList[170] = {
        19151031, 19151034, 19151036, 19151039, 19151041, 19151043, 19151044, 19151045,
        19151046, 19151047, 19151048, 19151049, 19151050, 19151052, 19151053, 19151054,
        19151055, 19151056, 19151066, 19151067, 19151068, 19151069, 19151070, 19151071,
        19151072, 19151082, 19151083, 19151084, 19152002, 19152003, 19152008, 19152009,
        19152010, 19152014, 19152016, 19152021, 19152023, 19152024, 19152025, 19152027,
        19152028, 19152029, 19152030, 19152031, 19152032, 19152033, 19152034, 19152035,
        19152036, 19152037, 19152038, 19152039, 19152040, 19152041, 19152042, 19152043,
        19152044, 19152045, 19152046, 19152048, 19152051, 19152052, 19152053, 19152054,
        19152055, 19152071, 19152073, 19152074, 19152075, 19152076, 19152081, 19153001,
        19153002, 19153003, 19153004, 19153007, 19153009, 19153010, 19153011, 19153012,
        19153013, 19153014, 19153015, 19153016, 19153017, 19153018, 19153019, 19153020,
        19153021, 19153022, 19153024, 19153025, 19153027, 19153028, 19153029, 19153031,
        19153033, 19153034, 19153035, 19153036, 19153037, 19153042, 19153043, 19153044,
        19153050, 19153051, 19153052, 19153053, 19153054, 19153055, 19153056, 19153057,
        19153058, 19153059, 19153061, 19153062, 19153063, 19153064, 19153066, 19154001,
        19154002, 19154005, 19154007, 19154027, 19154028, 19154029, 19154030, 19154031,
        19154032, 19154036, 19154037, 19154038, 19154039, 19154040, 19154041, 19154044,
        19154045, 19154046, 19154047, 19154048, 19154049, 19154052, 19154053, 19154054,
        19154055, 19154056, 19154057, 19154058, 19154061, 19154063, 19154064, 19154065,
        19154066, 19154067, 19155001, 19155003, 19155004, 19155005, 19155006, 19155008,
        19155009, 19155010, 19155011, 19155016, 19155017, 19155018, 19155019, 19155020,
        19155021, 19155022
    };
    for(Int_t runnumberIndex = 0;runnumberIndex < 170;runnumberIndex++) {
        if(goodRunList[runnumberIndex] == runnumber) {
            goodRunFlag = kTRUE;
            break;
        }
    }
    return goodRunFlag;
}
