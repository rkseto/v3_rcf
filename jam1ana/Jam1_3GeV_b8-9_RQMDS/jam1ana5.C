// for looking at v3
// load c++ and c headers
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include "TSystem.h"
#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TFrame.h"
#include "TROOT.h"
#include "TLeaf.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TFormula.h"
#include "TPaveLabel.h"
#include "TFile.h"
#include "TMath.h"
#include "TLeaf.h"
#include "TComplex.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TLorentzVector.h"

#define PI 3.14159
#define TWOPI 6.28318

std::ifstream file;
TChain *chain;
Int_t nentries=0;

double ybeam=0.;
const double m_proton = 0.938;
// charged
const double m_pion = 0.140;
const double m_kaon = 0.494;

int eventnum=0;
int testdummy=0;
int number_particles=0;
int refmultfxt=0;
double b_impact=0;
int npart_proj=0;
int npart_targ=0;
int nelas_proj=0;
int ninelas_proj=0;
int nelas_targ=0;
int ninelas_targ=0;
float dummy2=0;
double rot=0;

int part_id_in=0;
double px_in=0;
double py_in=0;
double pz_in=0;
double mass_in=0;
double x_in=0;
double y_in=0;
double z_in=0;
double t_in=0;
double rapidity_in=0;

// -- new stuff from smash --
double p0_in;
int so_id, charge_in;

double form_time_in,xsecfac_in,t_last_coll_in;
int ncoll_in,proc_id_origin_in,proc_type_origin_in,pdg_mother1_in,pdg_mother2_in;
//-----


int nev=0;

  FILE* fInputFile;
//  char read[200];

// TLeaves that store track parameters
/*
TLeaf * leaf_PID; 
TLeaf * leaf_Charge;
TLeaf * leaf_Px;
TLeaf * leaf_Py;
TLeaf * leaf_Pz;
TLeaf * leaf_Pt;
TLeaf * leaf_Pmag;
TLeaf * leaf_Eta;
TLeaf * leaf_Phi;
TLeaf * leaf_E;
TLeaf * leaf_Rapidity;      
*/


// for my thrown event
int nparticles=0;
double phitrack[1200000];
double pttrack[1200000];
double ytrack[1200000];
int nparticles_previous=0;
double phitrack_previous[1200000];
double pttrack_previous[1200000];
double ytrack_previous[1200000];
const int nparticlesmin=6;


bool readit(int itype, int ido, int itrack);
//define a function with parameters
Double_t fitf(Double_t *xx,Double_t *par);
int throwevent(int ievent);

void jam1ana5(int itype=5,int ntodo=-1, int irot=0, int i_inject=0, double Ecm=7.7, double bmin=0., double bmax=14.) {
    
  if(itype<=0){
    std::cout<<" aruments are: itype, 1=AMPT(make link to data) 2=RQMD 3=my generator(don't double inject flow) 4=smash 5=jam1; ntodo=num events (-1 runs all); irot = 0 no rotation of reaction plane 1=rotate; i_inject=1  to inject flow; Ecm for calculations of eta regions and cuts "<<std::endl;
    exit(1);
  }
  std::cout<<" Starting"<<std::endl;
  std::cout<<"itype="<<itype<<" ntodo="<<ntodo<<" irot="<<irot<<" i_inject="<<i_inject<<" Ecm="<<Ecm<<" bmin="<<bmin<<" bmax="<<bmax<<std::endl;

  double Ebeam=Ecm/2.;
  double pzbeam=sqrt(Ebeam*Ebeam-m_proton*m_proton);
  double ybeam=0.5*log((Ebeam+pzbeam)/(Ebeam-pzbeam));
  std::cout<<" ybeam="<<ybeam<<std::endl;
  
  double v_injected[6]={0.1,0.2,0.15,0.05,0.,0.};  // injected flow

  // initialize random numbers note: you don't seem to have to do this for RQMD
  TRandom3 * ran0 = new TRandom3 ( 0 );
 
  std::vector<int> eventnum_v;
  std::vector<int> number_particles_v;
  std::vector<int> refmultfxt_v;
  std::vector<double> b_impact_v;
  std::vector<int> npart_proj_v;
  std::vector<int> npart_targ_v;
  std::vector<int> nelas_proj_v;
  std::vector<int> ninelas_proj_v;
  std::vector<int> nelas_targ_v;
  std::vector<int> ninelas_targ_v;
  std::vector <double> rot_v;
  //  std::vector<std::vector <int>> part_id_v;
  std::vector<std::vector <int> > part_id_v;
  std::vector<std::vector <double> > px_v;
  std::vector<std::vector <double> > py_v;
  std::vector<std::vector <double> > pz_v;
  std::vector<std::vector <double> > mass_v;
  std::vector<std::vector <double> > x_v;
  std::vector<std::vector <double> > y_v;
  std::vector<std::vector <double> > z_v;
  std::vector<std::vector <double> > t_v;
  std::vector<std::vector <double> > energy_v;
  std::vector<std::vector <double> > pt_v;
  std::vector<std::vector <double> > rapidity_v;
  std::vector<std::vector <double> > phi_v;
  std::vector<std::vector <double> > eta_v;

  std::vector<std::vector <double> > form_time_v;
  std::vector<std::vector <double> > t_last_coll_v;
  std::vector<std::vector <double> > ncollpart_v;
    
  // we will define our histograms here
  // event level
  TH1D *number_particlesh = new TH1D("number_particlesh","number_particlesh",1000,0.,10000);
  TH1D *b_impactallh = new TH1D("impactallh","impactallh before b cut",150,0.,15);
  TH1D *b_impacth = new TH1D("impacth","impacth",150,0.,15);
  TH1D *bsq_impacth = new TH1D("bsqimpacth","bsqimpacth",200,0.,200);
  TH1D *n_participantsh = new TH1D("n_participantsh","n_participantsh",400,0.,400.);

  TH1D *fxtmulth = new TH1D("fxtmulth","fxtmulth",500,0.,500);
  TH2D *number_particles_bsqh = new TH2D("number_particles_bsqh","number_particles_bsqh",500,0.,1000., 200,0.,200);
  TH2D *fxtmult_bsqh = new TH2D("fxtmult_bsqh","fxtmult_bsqh",500,0.,500., 200,0.,200);
  TH2D *nparticipants_bsqh = new TH2D("nparticipants_bsqh","nparticipants_bsqh",500,0.,500., 200,0.,200);
  
  // particle level
  TH1D *infoh = new TH1D("infoh","1-itype ",10,0,10);
  TH1D *part_idh = new TH1D("part_idh","part_id",2000,-10000.,10000);
  TH1D *part_id2h = new TH1D("part_id2h","part_id for clusters",1000,0.,2000000000);

  TH1D *pxh = new TH1D("pxh","pxh",1000,0.,5);
  TH1D *pyh = new TH1D("pyh","pyh",1000,0.,5);
  TH1D *pzh = new TH1D("pzh","pzh",100,0.,Ecm);

  TH1D *xh = new TH1D("xh","log xh",1000,0.,6);
  TH1D *yh = new TH1D("yh","log yh",1000,0.,6);
  TH1D *zh = new TH1D("zh","log zh",1000,0.,6);
  TH1D *th = new TH1D("th","log th",1000,0.,6);

  TH1D *ncollparth = new TH1D("ncollparth","ncollparth",50,0.,50.);
  TH1D *formtimeh = new TH1D("formtimeh","formtimeh",800,-10.,190.);
  TH1D *lastcolltimeh = new TH1D("lastcolltimeh","lastcolltimeh",800,0.,200.);
  
  TH2D *rapidityptrawh = new TH2D("rapidityptrawh","rapidityptrawh",480,-6.,6.,500,0.0,5.0);
  
  TH2D *rapiditypth = new TH2D("rapiditypth","rapiditypth",120,-6.,6.,100,0.0,5.0);

  TH2D *rapiditypt_piph = new TH2D("rapiditypt_piph","rapiditypt_piph",120,-6.,6.,100,0.0,5.0);
  TH2D *rapiditypt_pimh = new TH2D("rapiditypt_pimh","rapiditypt_pimh",120,-6.,6.,100,0.0,5.0);
  TH2D *rapiditypt_kph = new TH2D("rapiditypt_kph","rapiditypt_kph",120,-6.,6.,100,0.0,5.0);
  TH2D *rapiditypt_kmh = new TH2D("rapiditypt_kmh","rapiditypt_kmh",120,-6.,6.,100,0.0,5.0);
  TH2D *rapiditypt_prh = new TH2D("rapiditypt_prh","rapiditypt_prh",120,-6.,6.,100,0.0,5.0);
  
  TH2D *ptpzh = new TH2D("ptpzh","ptpzh all participant nucleons",100,0.0,5.0,200.0,-Ecm,Ecm);
  TH2D *etarapidityh = new TH2D("etarapidityh","etarapidityh all participant nucleons",100, -8.0, 8.0,100,-8.0,8.0);
  TH2D *xvsz_allnh = new TH2D("xvsz_allnh","xvszh all participant nucleons",200,-200.,200.,200,-200.,200.);
  TH2D *xvsy_allnh = new TH2D("xvsy_allnh","xvsyh all participant nucleons",200,-200.,200.,200,-200.,200.);
  TH2D *xvst_allnh = new TH2D("xvst_allnh","xvsth all participant nucleons",200,-200.,200.,245,1.,50.);
  TH2D *ptvst_allnh = new TH2D("ptvst_allnh","pt vs t all participant nucleons",100,0.,5.,245,1.,50.);
  TH2D *etavst_allnh = new TH2D("etavst_allnh","eta vs t participant all nucleons",100,-10.,10.,245,1.,50.);
  TH2D *rapidityvst_allnh = new TH2D("rapidityvst_allnh","rapidity vs t all participant nucleons",100,-10.,10.,245,1.,50.);
  TH2D *xvsrapidity_allnucleonsh = new TH2D("xvsrapidity_allnucleonsh","xvsetah all nucleons",800,-4.,4.,200,-50.,50.);

  TH2D *pxvst_allnh = new TH2D("pxvst_allnh","pxvsth all participant  nucleons",200,-1.,1.,245,1.,50.);
  TH2D *pyvst_allnh = new TH2D("pyvst_allnh","pyvsth all participant nucleons",200,-1.,1.,245,1.,50.);
  TH2D *pzvst_allnh = new TH2D("pzvst_allnh","pzvsth all participant nucleons",200,-Ecm,Ecm,245,1.,50.);  

  TH2D *pxvstposrap_allnh = new TH2D("pxvstposrap_allnh","pxvstposraph all participant  nucleons",200,-1.,1.,245,1.,50.);
  TH2D *pyvstposrap_allnh = new TH2D("pyvstposrap_allnh","pyvstposraph all participant nucleons",200,-1.,1.,245,1.,50.);
  TH2D *pzvstposrap_allnh = new TH2D("pzvstposrap_allnh","pzvstposraph all participant nucleons",200,-Ecm,Ecm,245,1.,50.);  

  TH2D *pxvstybeam_allnh = new TH2D("pxvstybeam_allnh","pxvstybeamh all participant  nucleons",200,-1.,1.,245,1.,50.);
  TH2D *pyvstybeam_allnh = new TH2D("pyvstybeam_allnh","pyvstybeamh all participant nucleons",200,-1.,1.,245,1.,50.);
  TH2D *pzvstybeam_allnh = new TH2D("pzvstybeam_allnh","pzvstybeamh all participant nucleons",200,-Ecm,Ecm,245,1.,50.);  

  TProfile *cth1vst_allnh = new TProfile("cth1vst_allnh","cth1vsth all participant nucleons",245,1.,50.,-1.,1.);  
  TProfile *cth1vstposrap_allnh = new TProfile("cth1vstposrap_allnh","cth1vstposraph all participant nucleons",245,1.,50.,-1.,1.);  
  TProfile *cth1vstybeam_allnh = new TProfile("cth1vstybeam_allnh","cth1vstybeamh all participant nucleons",245,1.,50.,-1.,1.);  

  TProfile *cth2vst_allnh = new TProfile("cth2vst_allnh","cth2vsth all participant nucleons",245,1.,50.,-1.,1.);  
  TProfile *cth2vstposrap_allnh = new TProfile("cth2vstposrap_allnh","cth2vstposraph all participant nucleons",245,1.,50.,-1.,1.);  
  TProfile *cth2vstybeam_allnh = new TProfile("cth2vstybeam_allnh","cth2vstybeamh all participant nucleons",245,1.,50.,-1.,1.);  

  TProfile *cth3vst_allnh = new TProfile("cth3vst_allnh","cth3vsth all participant nucleons",245,1.,50.,-1.,1.);  
  TProfile *cth3vstposrap_allnh = new TProfile("cth3vstposrap_allnh","cth3vstposraph all participant nucleons",245,1.,50.,-1.,1.);  
  TProfile *cth3vstybeam_allnh = new TProfile("cth3vstybeam_allnh","cth3vstybeamh all participant nucleons",245,1.,50.,-1.,1.);  
  
  TH1D *pth = new TH1D("pth","pt",100,0.,5);
  TH2D *xvsyh = new TH2D("xvsyh","xvsyh",1000,-200.,200.,1000,-200.,200.);
  TH2D *zvsth = new TH2D("zvsth","zvsth",200,-200.,200.,245,1.,50.);

  TH2D *xvsth = new TH2D("xvsth","xvsth",200,-200.,200.,245,1.,50.);
  TH2D *yvsth = new TH2D("yvsth","yvsth",200,-200.,200.,245,1.,50.);
  
  TH2D *rvsth = new TH2D("rvsth","rvsth",300,0.,300.,245,1.,50.);
  TH1D *rapidityh = new TH1D("rapidityh","rapidityh",120,-6.,6.);
  TH1D *phih = new TH1D("phih","phi",100,-TWOPI,TWOPI);
  TH1D *etah = new TH1D("etah","etah",100,-5.,5.);  

  TH1D *phi_p_cuth = new TH1D("phi_p_cuth","phi_p_cuth",100,-TWOPI,TWOPI);
  TH1D *phi_p_cut2h = new TH1D("phi_p_cut2h","phi_p_cut2h",100,-TWOPI,TWOPI);
  TH2D *xvsy_p_cuth = new TH2D("xvsy_p_cuth","xvsy_p_cuth",1000,-200.,200.,1000,-200.,200.);
  TH2D *xvsy_p_cut2h = new TH2D("xvsy_p_cut2h","xvsy_p_cut2h",1000,-200.,200.,1000,-200.,200.);
  TH2D *xvszh = new TH2D("xvszh","xvszh all part nucleons eta>etamin",150,-150.,150.,200,-200.,200.);

  TH2D *xvsy_ptsmall_ybeamh = new TH2D("xvsy_ptsmall_ybeamh","xvsy_ptsmall_ybeamh",1000,-200.,200.,1000,-200.,200.);
  TH2D *xvsy_not_ptsmall_ybeamh = new TH2D("xvsy_not_ptsmall_ybeamh","xvsy_not_ptsmall_ybeamh",1000,-200.,200.,1000,-200.,200.);  

  TH2D *xvsyybeamh = new TH2D("xvsyybeamh","xvsyybeamh",1000,-200.,200.,1000,-200.,200.);

  TH2D *xvsy_form_time_h = new TH2D("xvsy_form_time_h","xvsy_form_time_h",1000,-200.,200.,1000,-200.,200.);
  TH2D *xvsy_t_last_coll_h = new TH2D("xvsy_t_last_coll_h","xvsy_t_last_coll_h",1000,-200.,200.,1000,-200.,200.);
  TH2D *xvsy_ncollpart_h = new TH2D("xvsy_ncollpart_h","xvsy_ncollpart_h",1000,-200.,200.,1000,-200.,200.);
  
  TH2D *xvszybeamh = new TH2D("xvszybeamh","xvszybeamh all part nucleons eta>ybeam",150,-150.,150.,200,-200.,200.);
  TH1D *part_idybeamh = new TH1D("part_idybeamh","part_id eta>ybeam",2500,0.,2500);

  TH1D *phi_p_posraph = new TH1D("phi_p_posraph","phi_p_posraph",100,-TWOPI,TWOPI);
  TH1D *phi_p_negraph = new TH1D("phi_p_negraph","phi_p_negraph",100,-TWOPI,TWOPI);
  TH2D *xvsy_p_posraph = new TH2D("xvsy_p_posraph","xvsy_p_posraph",1000,-200.,200.,1000,-200.,200.);
  TH2D *xvsy_p_negraph = new TH2D("xvsy_p_negraph","xvsy_p_negraph",1000,-200.,200.,1000,-200.,200.);

  TH2D *xvsy_p_binm4 = new TH2D("xvsy_p_binm4","xvsy_p_binm4",1000,-200.,200.,1000,-200.,200.);
  TH2D *xvsy_p_binm3 = new TH2D("xvsy_p_binm3","xvsy_p_binm3",1000,-200.,200.,1000,-200.,200.);
  TH2D *xvsy_p_binm2 = new TH2D("xvsy_p_binm2","xvsy_p_binm2",1000,-200.,200.,1000,-200.,200.);
  TH2D *xvsy_p_binm1 = new TH2D("xvsy_p_binm1","xvsy_p_binm1",1000,-200.,200.,1000,-200.,200.);
  TH2D *xvsy_p_bin0 = new TH2D("xvsy_p_bin0","xvsy_p_bin0",1000,-200.,200.,1000,-200.,200.);
  TH2D *xvsy_p_binp1 = new TH2D("xvsy_p_binp1","xvsy_p_binp1",1000,-200.,200.,1000,-200.,200.);
  TH2D *xvsy_p_binp2 = new TH2D("xvsy_p_binp2","xvsy_p_binp2",1000,-200.,200.,1000,-200.,200.);
  TH2D *xvsy_p_binp3 = new TH2D("xvsy_p_binp3","xvsy_p_binp3",1000,-200.,200.,1000,-200.,200.);
  TH2D *xvsy_p_binp4 = new TH2D("xvsy_p_binp4","xvsy_p_binp4",1000,-200.,200.,1000,-200.,200.);

  TH2D *xvsy_p_pt0_binm4 = new TH2D("xvsy_p_pt0_binm4","xvsy_p_pt0_binm4",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt0_binm3 = new TH2D("xvsy_p_pt0_binm3","xvsy_p_pt0_binm3",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt0_binm2 = new TH2D("xvsy_p_pt0_binm2","xvsy_p_pt0_binm2",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt0_binm1 = new TH2D("xvsy_p_pt0_binm1","xvsy_p_pt0_binm1",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt0_bin0 = new TH2D("xvsy_p_pt0_bin0","xvsy_p_pt0_bin0",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt0_binp1 = new TH2D("xvsy_p_pt0_binp1","xvsy_p_pt0_binp1",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt0_binp2 = new TH2D("xvsy_p_pt0_binp2","xvsy_p_pt0_binp2",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt0_binp3 = new TH2D("xvsy_p_pt0_binp3","xvsy_p_pt0_binp3",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt0_binp4 = new TH2D("xvsy_p_pt0_binp4","xvsy_p_pt0_binp4",500,-100.,100.,500,-100.,100.);      
  
  TH2D *xvsy_p_pt1_binm4 = new TH2D("xvsy_p_pt1_binm4","xvsy_p_pt1_binm4",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt1_binm3 = new TH2D("xvsy_p_pt1_binm3","xvsy_p_pt1_binm3",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt1_binm2 = new TH2D("xvsy_p_pt1_binm2","xvsy_p_pt1_binm2",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt1_binm1 = new TH2D("xvsy_p_pt1_binm1","xvsy_p_pt1_binm1",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt1_bin0 = new TH2D("xvsy_p_pt1_bin0","xvsy_p_pt1_bin0",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt1_binp1 = new TH2D("xvsy_p_pt1_binp1","xvsy_p_pt1_binp1",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt1_binp2 = new TH2D("xvsy_p_pt1_binp2","xvsy_p_pt1_binp2",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt1_binp3 = new TH2D("xvsy_p_pt1_binp3","xvsy_p_pt1_binp3",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt1_binp4 = new TH2D("xvsy_p_pt1_binp4","xvsy_p_pt1_binp4",500,-100.,100.,500,-100.,100.);      
  
  TH2D *xvsy_p_pt2_binm4 = new TH2D("xvsy_p_pt2_binm4","xvsy_p_pt2_binm4",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt2_binm3 = new TH2D("xvsy_p_pt2_binm3","xvsy_p_pt2_binm3",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt2_binm2 = new TH2D("xvsy_p_pt2_binm2","xvsy_p_pt2_binm2",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt2_binm1 = new TH2D("xvsy_p_pt2_binm1","xvsy_p_pt2_binm1",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt2_bin0 = new TH2D("xvsy_p_pt2_bin0","xvsy_p_pt2_bin0",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt2_binp1 = new TH2D("xvsy_p_pt2_binp1","xvsy_p_pt2_binp1",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt2_binp2 = new TH2D("xvsy_p_pt2_binp2","xvsy_p_pt2_binp2",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt2_binp3 = new TH2D("xvsy_p_pt2_binp3","xvsy_p_pt2_binp3",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt2_binp4 = new TH2D("xvsy_p_pt2_binp4","xvsy_p_pt2_binp4",500,-100.,100.,500,-100.,100.);      

  TH2D *xvsy_p_pt3_binm4 = new TH2D("xvsy_p_pt3_binm4","xvsy_p_pt3_binm4",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt3_binm3 = new TH2D("xvsy_p_pt3_binm3","xvsy_p_pt3_binm3",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt3_binm2 = new TH2D("xvsy_p_pt3_binm2","xvsy_p_pt3_binm2",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt3_binm1 = new TH2D("xvsy_p_pt3_binm1","xvsy_p_pt3_binm1",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt3_bin0 = new TH2D("xvsy_p_pt3_bin0","xvsy_p_pt3_bin0",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt3_binp1 = new TH2D("xvsy_p_pt3_binp1","xvsy_p_pt3_binp1",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt3_binp2 = new TH2D("xvsy_p_pt3_binp2","xvsy_p_pt3_binp2",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt3_binp3 = new TH2D("xvsy_p_pt3_binp3","xvsy_p_pt3_binp3",500,-100.,100.,500,-100.,100.);
  TH2D *xvsy_p_pt3_binp4 = new TH2D("xvsy_p_pt3_binp4","xvsy_p_pt3_binp4",500,-100.,100.,500,-100.,100.);      
  

  
  TH2D *xvsyallh = new TH2D("xvsyallh","xvsyallh",1000,-200.,200.,1000,-200.,200.);
  TH2D *zvstallh = new TH2D("zvstallh","zvstallh",200,-200.,200.,245,1.,50.);
  TH2D *xvsyspech = new TH2D("xvsyspech","xvsyspech",1000,-200.,200.,1000,-200.,200.);
  TH2D *xvszspech = new TH2D("xvszspech","xvszspech",150,-150.,150.,200,-200.,200.);
  TH2D *zvstspech = new TH2D("zvstspech","zvstspech",200,-200.,200.,245,1.,50.);

  TH1D *ncollpartspech = new TH1D("ncollpartspech","ncollpartspech",50,0.,50.);
  
  TH1D *reactionplanerotationh = new TH1D("reactionplanerotationh","reactionplanerotation",180,-PI,PI);
  TH1D *PSI1h = new TH1D("PSI1h","PSI1h",180,-PI,PI);
  TH1D *PSI2h = new TH1D("PSI2h","PSI2h",180,-PI,PI);
  TH1D *PSI3h = new TH1D("PSI3h","PSI3h",180,-PI,PI);
  TH1D *PSI4h = new TH1D("PSI4h","PSI4h",180,-PI,PI);

  TH2D *x0y0h = new TH2D("x0y0h","x0y0h",1000,-200.,200.,1000,-200.,200.);
  TH2D *x0y0posraph = new TH2D("x0y0posraph","x0y0posraph",1000,-200.,200.,1000,-200.,200.);
  TH2D *x0y0negraph = new TH2D("x0y0negraph","x0y0negraph",1000,-200.,200.,1000,-200.,200.);

  TH1D *phixyh = new TH1D("phixyh","phixyh",180,-TWOPI,TWOPI);
  TH1D *phixyposraph = new TH1D("phixyposraph","phixyposraph",180,-TWOPI,TWOPI);
  TH1D *phixynegraph = new TH1D("phixynegraph","phixynegraph",180,-TWOPI,TWOPI);

  TH1D *rsqh = new TH1D("rsqh","rsqh",100,0.,200.);
  TH1D *rsqposraph = new TH1D("rsqposraph","rsqposraph",100,0.,200.);
  TH1D *rsqnegraph = new TH1D("rsqnegraph","rsqnegraph",100,0.,200.);

  TH1D *rsqavh = new TH1D("rsqavh","rsqavh",100,0.,200.);
  TH1D *rsqavposraph = new TH1D("rsqavposraph","rsqavposraph",100,0.,200.);
  TH1D *rsqavnegraph = new TH1D("rsqavnegraph","rsqavnegraph",100,0.,200.);

  TH1D *eps1h = new TH1D("eps1h","eps1h p/n particpants",100,0.,1.);
  TH1D *eps2h = new TH1D("eps2h","eps2h p/n particpants",100,0.,1.);
  TH1D *eps3h = new TH1D("eps3h","eps3h p/n particpants",100,0.,1.);
  TH1D *eps4h = new TH1D("eps4h","eps4h p/n particpants",100,0.,1.);
  TH1D *eps5h = new TH1D("eps5h","eps5h p/n particpants",100,0.,1.);
  TH1D *eps6h = new TH1D("eps6h","eps6h p/n particpants",100,0.,1.);

  TH1D *eps1posraph = new TH1D("eps1posraph","eps1posraph eta>ybeam/2 p/n particpants",100,0.,1.);
  TH1D *eps2posraph = new TH1D("eps2posraph","eps2posraph eta>ybeam/2p/n particpants",100,0.,1.);
  TH1D *eps3posraph = new TH1D("eps3posraph","eps3posraph eta>ybeam/2 p/n particpants not recentered",100,0.,1.);

  TH1D *eps4posraph = new TH1D("eps4posraph","eps4posraph eta>ybeam/2 p/n particpants",100,0.,1.);
  TH1D *eps5posraph = new TH1D("eps5posraph","eps5posraph eta>ybeam/2 p/n particpants",100,0.,1.);
  TH1D *eps6posraph = new TH1D("eps6posraph","eps6posraph eta>ybeam/2 p/n particpants",100,0.,1.);

  TH1D *eps1negraph = new TH1D("eps1negraph","eps1negraph eta<-ybeam/2 p/n particpants",100,0.,1.);
  TH1D *eps2negraph = new TH1D("eps2negraph","eps2negraph eta<-ybeam/2 p/n particpants",100,0.,1.);
  TH1D *eps3negraph = new TH1D("eps3negraph","eps3negraph eta<-ybeam/2 p/n particpants not recentered",100,0.,1.);

  TH1D *eps4negraph = new TH1D("eps4negraph","eps4negraph eta<-ybeam/2 p/n particpants",100,0.,1.);
  TH1D *eps5negraph = new TH1D("eps5negraph","eps5negraph eta<-ybeam/2 p/n particpants",100,0.,1.);
  TH1D *eps6negraph = new TH1D("eps6negraph","eps6negraph eta<-ybeam/2 p/n particpants",100,0.,1.);

  TH1D *PSI1posraph = new TH1D("PSI1posraph","PSI1posraph",180,-PI,PI);
  TH1D *PSI1negraph = new TH1D("PSI1negraph","PSI1negraph",180,-PI,PI);
  TH1D *PSI2posraph = new TH1D("PSI2posraph","PSI2posraph",180,-PI,PI);
  TH1D *PSI2negraph = new TH1D("PSI2negraph","PSI2negraph",180,-PI,PI);  
  TH1D *PSI3posraph = new TH1D("PSI3posraph","PSI3posraph",180,-PI,PI);
  TH1D *PSI3negraph = new TH1D("PSI3negraph","PSI3negraph",180,-PI,PI);
  
  TH2D *PSI1_roth = new TH2D("PSI1_roth","PSI1_roth",180,-PI,PI,180,-PI,PI);
  TH2D *PSI2_roth = new TH2D("PSI2_roth","PSI2_roth",180,-PI,PI,180,-PI,PI);
  TH2D *PSI3_roth = new TH2D("PSI3_roth","PSI3_roth",180,-PI,PI,180,-PI,PI);
  TH2D *PSI4_roth = new TH2D("PSI4_roth","PSI4_roth",180,-PI,PI,180,-PI,PI);
  
  TH1D *v1_RP1h = new TH1D("v1_RP1h","v1_RP1",100,-0.5,0.5);
  TH1D *v2_RP1h = new TH1D("v2_RP1h","v2_RP1",100,-0.5,0.5);
  TH1D *v3_RP1h = new TH1D("v3_RP1h","v3_RP1",100,-0.5,0.5);
  TH1D *v4_RP1h = new TH1D("v4_RP1h","v4_RP1",100,-0.5,0.5);
  TH1D *v5_RP1h = new TH1D("v5_RP1h","v5_RP1",100,-0.5,0.5);
  TH1D *v6_RP1h = new TH1D("v6_RP1h","v6_RP1",100,-0.5,0.5);  
  TH1D *v2_RP2h = new TH1D("v2_RP2h","v2_RP2",100,-0.5,0.5);
  TH1D *v4_RP2h = new TH1D("v4_RP2h","v4_RP2",100,-0.5,0.5);
  TH1D *v6_RP2h = new TH1D("v6_RP2h","v6_RP2",100,-0.5,0.5);  
  TH1D *v3_RP3h = new TH1D("v3_RP3h","v3_RP3",100,-0.5,0.5);
  TH1D *v6_RP3h = new TH1D("v6_RP3h","v6_RP3",100,-0.5,0.5);


  TProfile *cth_1_rapidity_all_h = new TProfile("cth_1_rapidity_all_h","cth_1_rapidity_all_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_1_pt_all_h = new TProfile("cth_1_pt_all_h","cth_1_pt_all_h",20,0.,2.,-1.,1.);
  TProfile *cth_1_b_all_h = new TProfile("cth_1_b_all_h","cth_1_b_all_h",80,0.,200.,-1.,1.);  
  TProfile *cth_1_rapidity_pip_h = new TProfile("cth_1_rapidity_pip_h","cth_1_rapidity_pip_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_1_pt_pip_h = new TProfile("cth_1_pt_pip_h","cth_1_pt_pip_h",20,0.,2.,-1.,1.);
  TProfile *cth_1_b_pip_h = new TProfile("cth_1_b_pip_h","cth_1_b_pip_h",80,0.,200.,-1.,1.);  
  TProfile *cth_1_rapidity_pim_h = new TProfile("cth_1_rapidity_pim_h","cth_1_rapidity_pim_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_1_pt_pim_h = new TProfile("cth_1_pt_pim_h","cth_1_pt_pim_h",20,0.,2.,-1.,1.);
  TProfile *cth_1_b_pim_h = new TProfile("cth_1_b_pim_h","cth_1_b_pim_h",80,0.,200.,-1.,1.);  
  TProfile *cth_1_rapidity_kp_h = new TProfile("cth_1_rapidity_kp_h","cth_1_rapidity_kp_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_1_pt_kp_h = new TProfile("cth_1_pt_kp_h","cth_1_pt_kp_h",20,0.,2.,-1.,1.);
  TProfile *cth_1_b_kp_h = new TProfile("cth_1_b_kp_h","cth_1_b_kp_h",80,0.,200.,-1.,1.);  
  TProfile *cth_1_rapidity_km_h = new TProfile("cth_1_rapidity_km_h","cth_1_rapidity_km_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_1_pt_km_h = new TProfile("cth_1_pt_km_h","cth_1_pt_km_h",20,0.,2.,-1.,1.);
  TProfile *cth_1_b_km_h = new TProfile("cth_1_b_km_h","cth_1_b_km_h",80,0.,200.,-1.,1.);  
  TProfile *cth_1_rapidity_pr_h = new TProfile("cth_1_rapidity_pr_h","cth_1_rapidity_pr_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_1_pt_pr_h = new TProfile("cth_1_pt_pr_h","cth_1_pt_pr_h",20,0.,2.,-1.,1.);
  TProfile *cth_1_b_pr_h = new TProfile("cth_1_b_pr_h","cth_1_b_pr_h",80,0.,200.,-1.,1.);  
  TProfile *cth_1_rapidity_deut_h = new TProfile("cth_1_rapidity_deut_h","cth_1_rapidity_deut_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_1_pt_deut_h = new TProfile("cth_1_pt_deut_h","cth_1_pt_deut_h",20,0.,2.,-1.,1.);
  TProfile *cth_1_b_deut_h = new TProfile("cth_1_b_deut_h","cth_1_b_deut_h",80,0.,200.,-1.,1.);  
  
  TProfile *cth_2_rapidity_all_h = new TProfile("cth_2_rapidity_all_h","cth_2_rapidity_all_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_2_pt_all_h = new TProfile("cth_2_pt_all_h","cth_2_pt_all_h",20,0.,2.,-1.,1.);
  TProfile *cth_2_b_all_h = new TProfile("cth_2_b_all_h","cth_2_b_all_h",80,0.,200.,-1.,1.);  
  TProfile *cth_2_rapidity_pip_h = new TProfile("cth_2_rapidity_pip_h","cth_2_rapidity_pip_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_2_pt_pip_h = new TProfile("cth_2_pt_pip_h","cth_2_pt_pip_h",20,0.,2.,-1.,1.);
  TProfile *cth_2_b_pip_h = new TProfile("cth_2_b_pip_h","cth_2_b_pip_h",80,0.,200.,-1.,1.);  
  TProfile *cth_2_rapidity_pim_h = new TProfile("cth_2_rapidity_pim_h","cth_2_rapidity_pim_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_2_pt_pim_h = new TProfile("cth_2_pt_pim_h","cth_2_pt_pim_h",20,0.,2.,-1.,1.);
  TProfile *cth_2_b_pim_h = new TProfile("cth_2_b_pim_h","cth_2_b_pim_h",80,0.,200.,-1.,1.);  
  TProfile *cth_2_rapidity_kp_h = new TProfile("cth_2_rapidity_kp_h","cth_2_rapidity_kp_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_2_pt_kp_h = new TProfile("cth_2_pt_kp_h","cth_2_pt_kp_h",20,0.,2.,-1.,1.);
  TProfile *cth_2_b_kp_h = new TProfile("cth_2_b_kp_h","cth_2_b_kp_h",80,0.,200.,-1.,1.);  
  TProfile *cth_2_rapidity_km_h = new TProfile("cth_2_rapidity_km_h","cth_2_rapidity_km_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_2_pt_km_h = new TProfile("cth_2_pt_km_h","cth_2_pt_km_h",20,0.,2.,-1.,1.);
  TProfile *cth_2_b_km_h = new TProfile("cth_2_b_km_h","cth_2_b_km_h",80,0.,200.,-1.,1.);  
  TProfile *cth_2_rapidity_pr_h = new TProfile("cth_2_rapidity_pr_h","cth_2_rapidity_pr_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_2_pt_pr_h = new TProfile("cth_2_pt_pr_h","cth_2_pt_pr_h",20,0.,2.,-1.,1.);
  TProfile *cth_2_b_pr_h = new TProfile("cth_2_b_pr_h","cth_2_b_pr_h",80,0.,200.,-1.,1.);  
  TProfile *cth_2_rapidity_deut_h = new TProfile("cth_2_rapidity_deut_h","cth_2_rapidity_deut_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_2_pt_deut_h = new TProfile("cth_2_pt_deut_h","cth_2_pt_deut_h",20,0.,2.,-1.,1.);
  TProfile *cth_2_b_deut_h = new TProfile("cth_2_b_deut_h","cth_2_b_deut_h",80,0.,200.,-1.,1.);  
  
  TProfile *cth_3_rapidity_all_h = new TProfile("cth_3_rapidity_all_h","cth_3_rapidity_all_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_pt_all_h = new TProfile("cth_3_pt_all_h","cth_3_pt_all_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_b_all_h = new TProfile("cth_3_b_all_h","cth_3_b_all_h",80,0.,200.,-1.,1.);  
  TProfile *cth_3_rapidity_pip_h = new TProfile("cth_3_rapidity_pip_h","cth_3_rapidity_pip_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_pt_pip_h = new TProfile("cth_3_pt_pip_h","cth_3_pt_pip_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_b_pip_h = new TProfile("cth_3_b_pip_h","cth_3_b_pip_h",80,0.,200.,-1.,1.);  
  TProfile *cth_3_rapidity_pim_h = new TProfile("cth_3_rapidity_pim_h","cth_3_rapidity_pim_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_pt_pim_h = new TProfile("cth_3_pt_pim_h","cth_3_pt_pim_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_b_pim_h = new TProfile("cth_3_b_pim_h","cth_3_b_pim_h",80,0.,200.,-1.,1.);  
  TProfile *cth_3_rapidity_kp_h = new TProfile("cth_3_rapidity_kp_h","cth_3_rapidity_kp_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_pt_kp_h = new TProfile("cth_3_pt_kp_h","cth_3_pt_kp_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_b_kp_h = new TProfile("cth_3_b_kp_h","cth_3_b_kp_h",80,0.,200.,-1.,1.);  
  TProfile *cth_3_rapidity_km_h = new TProfile("cth_3_rapidity_km_h","cth_3_rapidity_km_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_pt_km_h = new TProfile("cth_3_pt_km_h","cth_3_pt_km_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_b_km_h = new TProfile("cth_3_b_km_h","cth_3_b_km_h",80,0.,200.,-1.,1.);  
  TProfile *cth_3_rapidity_pr_h = new TProfile("cth_3_rapidity_pr_h","cth_3_rapidity_pr_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_pt_pr_h = new TProfile("cth_3_pt_pr_h","cth_3_pt_pr_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_b_pr_h = new TProfile("cth_3_b_pr_h","cth_3_b_pr_h",80,0.,200.,-1.,1.);  
  TProfile *cth_3_rapidity_deut_h = new TProfile("cth_3_rapidity_deut_h","cth_3_rapidity_deut_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_pt_deut_h = new TProfile("cth_3_pt_deut_h","cth_3_pt_deut_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_b_deut_h = new TProfile("cth_3_b_deut_h","cth_3_b_deut_h",80,0.,200.,-1.,1.);  

  //-----------------------------------------
  TProfile *cth_3_rapidity_ptcut1_cent_all_h = new TProfile("cth_3_rapidity_ptcut1_cent_all_h","cth_3_rapidity_ptcut1_cent_all_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_rapidity_ptcut1_cent_pip_h = new TProfile("cth_3_rapidity_ptcut1_cent_pip_h","cth_3_rapidity_ptcut1_cent_pip_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_rapidity_ptcut1_cent_pim_h = new TProfile("cth_3_rapidity_ptcut1_cent_pim_h","cth_3_rapidity_ptcut1_cent_pim_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_rapidity_ptcut1_cent_kp_h = new TProfile("cth_3_rapidity_ptcut1_cent_kp_h","cth_3_rapidity_ptcut1_cent_kp_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_rapidity_ptcut1_cent_km_h = new TProfile("cth_3_rapidity_ptcut1_cent_km_h","cth_3_rapidity_ptcut1_cent_km_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_rapidity_ptcut1_cent_pr_h = new TProfile("cth_3_rapidity_ptcut1_cent_pr_h","cth_3_rapidity_ptcut1_cent_pr_h",120,-6.,6.,-1.,1.);
  TProfile *cth_3_rapidity_ptcut2_cent_pr_h = new TProfile("cth_3_rapidity_ptcut2_cent_pr_h","cth_3_rapidity_ptcut2_cent_pr_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_rapidity_ptcut1_cent_deut_h = new TProfile("cth_3_rapidity_ptcut1_cent_deut_h","cth_3_rapidity_ptcut1_cent_deut_h",120,-6.,6.,-1.,1.);
  TProfile *cth_3_rapidity_ptcut2_cent_deut_h = new TProfile("cth_3_rapidity_ptcut2_cent_deut_h","cth_3_rapidity_ptcut2_cent_deut_h",120,-6.,6.,-1.,1.);  

  TProfile *cth_3_rapidity_ptcut1_mid_all_h = new TProfile("cth_3_rapidity_ptcut1_mid_all_h","cth_3_rapidity_ptcut1_mid_all_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_rapidity_ptcut1_mid_pip_h = new TProfile("cth_3_rapidity_ptcut1_mid_pip_h","cth_3_rapidity_ptcut1_mid_pip_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_rapidity_ptcut1_mid_pim_h = new TProfile("cth_3_rapidity_ptcut1_mid_pim_h","cth_3_rapidity_ptcut1_mid_pim_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_rapidity_ptcut1_mid_kp_h = new TProfile("cth_3_rapidity_ptcut1_mid_kp_h","cth_3_rapidity_ptcut1_mid_kp_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_rapidity_ptcut1_mid_km_h = new TProfile("cth_3_rapidity_ptcut1_mid_km_h","cth_3_rapidity_ptcut1_mid_km_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_rapidity_ptcut1_mid_pr_h = new TProfile("cth_3_rapidity_ptcut1_mid_pr_h","cth_3_rapidity_ptcut1_mid_pr_h",120,-6.,6.,-1.,1.);
  TProfile *cth_3_rapidity_ptcut2_mid_pr_h = new TProfile("cth_3_rapidity_ptcut2_mid_pr_h","cth_3_rapidity_ptcut2_mid_pr_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_rapidity_ptcut1_mid_deut_h = new TProfile("cth_3_rapidity_ptcut1_mid_deut_h","cth_3_rapidity_ptcut1_mid_deut_h",120,-6.,6.,-1.,1.);
  TProfile *cth_3_rapidity_ptcut2_mid_deut_h = new TProfile("cth_3_rapidity_ptcut2_mid_deut_h","cth_3_rapidity_ptcut2_mid_deut_h",120,-6.,6.,-1.,1.);  

  TProfile *cth_3_rapidity_ptcut1_periph_all_h = new TProfile("cth_3_rapidity_ptcut1_periph_all_h","cth_3_rapidity_ptcut1_periph_all_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_rapidity_ptcut1_periph_pip_h = new TProfile("cth_3_rapidity_ptcut1_periph_pip_h","cth_3_rapidity_ptcut1_periph_pip_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_rapidity_ptcut1_periph_pim_h = new TProfile("cth_3_rapidity_ptcut1_periph_pim_h","cth_3_rapidity_ptcut1_periph_pim_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_rapidity_ptcut1_periph_kp_h = new TProfile("cth_3_rapidity_ptcut1_periph_kp_h","cth_3_rapidity_ptcut1_periph_kp_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_rapidity_ptcut1_periph_km_h = new TProfile("cth_3_rapidity_ptcut1_periph_km_h","cth_3_rapidity_ptcut1_periph_km_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_rapidity_ptcut1_periph_pr_h = new TProfile("cth_3_rapidity_ptcut1_periph_pr_h","cth_3_rapidity_ptcut1_periph_pr_h",120,-6.,6.,-1.,1.);
  TProfile *cth_3_rapidity_ptcut2_periph_pr_h = new TProfile("cth_3_rapidity_ptcut2_periph_pr_h","cth_3_rapidity_ptcut2_periph_pr_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_rapidity_ptcut1_periph_deut_h = new TProfile("cth_3_rapidity_ptcut1_periph_deut_h","cth_3_rapidity_ptcut1_periph_deut_h",120,-6.,6.,-1.,1.);
  TProfile *cth_3_rapidity_ptcut2_periph_deut_h = new TProfile("cth_3_rapidity_ptcut2_periph_deut_h","cth_3_rapidity_ptcut2_periph_deut_h",120,-6.,6.,-1.,1.);  

  TProfile *cth_3_pt_ycut1_cent_all_h = new TProfile("cth_3_pt_ycut1_cent_all_h","cth_3_pt_ycut1_cent_all_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_pt_ycut1_cent_pip_h = new TProfile("cth_3_pt_ycut1_cent_pip_h","cth_3_pt_ycut1_cent_pip_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_pt_ycut1_cent_pim_h = new TProfile("cth_3_pt_ycut1_cent_pim_h","cth_3_pt_ycut1_cent_pim_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_pt_ycut1_cent_kp_h = new TProfile("cth_3_pt_ycut1_cent_kp_h","cth_3_pt_ycut1_cent_kp_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_pt_ycut1_cent_km_h = new TProfile("cth_3_pt_ycut1_cent_km_h","cth_3_pt_ycut1_cent_km_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_pt_ycut1_cent_pr_h = new TProfile("cth_3_pt_ycut1_cent_pr_h","cth_3_pt_ycut1_cent_pr_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_pt_ycut1_cent_deut_h = new TProfile("cth_3_pt_ycut1_cent_deut_h","cth_3_pt_ycut1_cent_deut_h",20,0.,2.,-1.,1.);

  TProfile *cth_3_pt_ycut1_mid_all_h = new TProfile("cth_3_pt_ycut1_mid_all_h","cth_3_pt_ycut1_mid_all_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_pt_ycut1_mid_pip_h = new TProfile("cth_3_pt_ycut1_mid_pip_h","cth_3_pt_ycut1_mid_pip_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_pt_ycut1_mid_pim_h = new TProfile("cth_3_pt_ycut1_mid_pim_h","cth_3_pt_ycut1_mid_pim_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_pt_ycut1_mid_kp_h = new TProfile("cth_3_pt_ycut1_mid_kp_h","cth_3_pt_ycut1_mid_kp_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_pt_ycut1_mid_km_h = new TProfile("cth_3_pt_ycut1_mid_km_h","cth_3_pt_ycut1_mid_km_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_pt_ycut1_mid_pr_h = new TProfile("cth_3_pt_ycut1_mid_pr_h","cth_3_pt_ycut1_mid_pr_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_pt_ycut1_mid_deut_h = new TProfile("cth_3_pt_ycut1_mid_deut_h","cth_3_pt_ycut1_mid_deut_h",20,0.,2.,-1.,1.);

  
  TProfile *cth_3_pt_ycut1_periph_all_h = new TProfile("cth_3_pt_ycut1_periph_all_h","cth_3_pt_ycut1_periph_all_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_pt_ycut1_periph_pip_h = new TProfile("cth_3_pt_ycut1_periph_pip_h","cth_3_pt_ycut1_periph_pip_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_pt_ycut1_periph_pim_h = new TProfile("cth_3_pt_ycut1_periph_pim_h","cth_3_pt_ycut1_periph_pim_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_pt_ycut1_periph_kp_h = new TProfile("cth_3_pt_ycut1_periph_kp_h","cth_3_pt_ycut1_periph_kp_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_pt_ycut1_periph_km_h = new TProfile("cth_3_pt_ycut1_periph_km_h","cth_3_pt_ycut1_periph_km_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_pt_ycut1_periph_pr_h = new TProfile("cth_3_pt_ycut1_periph_pr_h","cth_3_pt_ycut1_periph_pr_h",20,0.,2.,-1.,1.);
  TProfile *cth_3_pt_ycut1_periph_deut_h = new TProfile("cth_3_pt_ycut1_periph_deut_h","cth_3_pt_ycut1_periph_deut_h",20,0.,2.,-1.,1.);

  TProfile *cth_3_b_ycut1_ptcut1_all_h = new TProfile("cth_3_b_ycut1_ptcut1_all_h","cth_3_b_ycut1_ptcut1_all_h",80,0.,200.,-1.,1.);  
  TProfile *cth_3_b_ycut1_ptcut1_pip_h = new TProfile("cth_3_b_ycut1_ptcut1_pip_h","cth_3_b_ycut1_ptcut1_pip_h",80,0.,200.,-1.,1.);  
  TProfile *cth_3_b_ycut1_ptcut1_pim_h = new TProfile("cth_3_b_ycut1_ptcut1_pim_h","cth_3_b_ycut1_ptcut1_pim_h",80,0.,200.,-1.,1.);  
  TProfile *cth_3_b_ycut1_ptcut1_kp_h = new TProfile("cth_3_b_ycut1_ptcut1_kp_h","cth_3_b_ycut1_ptcut1_kp_h",80,0.,200.,-1.,1.);  
  TProfile *cth_3_b_ycut1_ptcut1_km_h = new TProfile("cth_3_b_ycut1_ptcut1_km_h","cth_3_b_ycut1_ptcut1_km_h",80,0.,200.,-1.,1.);  
  TProfile *cth_3_b_ycut1_ptcut1_pr_h = new TProfile("cth_3_b_ycut1_ptcut1_pr_h","cth_3_b_ycut1_ptcut1_pr_h",80,0.,200.,-1.,1.);
  TProfile *cth_3_b_ycut1_ptcut2_pr_h = new TProfile("cth_3_b_ycut1_ptcut2_pr_h","cth_3_b_ycut1_ptcut2_pr_h",80,0.,200.,-1.,1.);  
  TProfile *cth_3_b_ycut1_ptcut1_deut_h = new TProfile("cth_3_b_ycut1_ptcut1_deut_h","cth_3_b_ycut1_ptcut1_deut_h",80,0.,200.,-1.,1.);
  TProfile *cth_3_b_ycut1_ptcut2_deut_h = new TProfile("cth_3_b_ycut1_ptcut2_deut_h","cth_3_b_ycut1_ptcut2_deut_h",80,0.,200.,-1.,1.);  

// eps_x
double epsxymin=-100000.;
double epsxymax=100000.;

// note eps_3=sqrt(eps_x^2+eps_y^2)/rsq1 - vs rapidity, pt, b
 
  TProfile *eps_x_rapidity_ptcut1_cent_pr_h = new TProfile("eps_x_rapidity_ptcut1_cent_pr_h","eps_x_rapidity_ptcut1_cent_pr_h",120,-6.,6.,epsxymin,epsxymin);
  TProfile *eps_x_rapidity_ptcut2_cent_pr_h = new TProfile("eps_x_rapidity_ptcut2_cent_pr_h","eps_x_rapidity_ptcut2_cent_pr_h",120,-6.,6.,epsxymin,epsxymin);  
  TProfile *eps_x_rapidity_ptcut1_mid_pr_h = new TProfile("eps_x_rapidity_ptcut1_mid_pr_h","eps_x_rapidity_ptcut1_mid_pr_h",120,-6.,6.,epsxymin,epsxymin);
  TProfile *eps_x_rapidity_ptcut2_mid_pr_h = new TProfile("eps_x_rapidity_ptcut2_mid_pr_h","eps_x_rapidity_ptcut2_mid_pr_h",120,-6.,6.,epsxymin,epsxymin);  
  TProfile *eps_x_rapidity_ptcut1_periph_pr_h = new TProfile("eps_x_rapidity_ptcut1_periph_pr_h","eps_x_rapidity_ptcut1_periph_pr_h",120,-6.,6.,epsxymin,epsxymin);
  TProfile *eps_x_rapidity_ptcut2_periph_pr_h = new TProfile("eps_x_rapidity_ptcut2_periph_pr_h","eps_x_rapidity_ptcut2_periph_pr_h",120,-6.,6.,epsxymin,epsxymin);  

  TProfile *eps_x_pt_ycut1_cent_pr_h = new TProfile("eps_x_pt_ycut1_cent_pr_h","eps_x_pt_ycut1_cent_pr_h",20,0.,2.,epsxymin,epsxymin);
  TProfile *eps_x_pt_ycut1_mid_pr_h = new TProfile("eps_x_pt_ycut1_mid_pr_h","eps_x_pt_ycut1_mid_pr_h",20,0.,2.,epsxymin,epsxymin);
  TProfile *eps_x_pt_ycut1_periph_pr_h = new TProfile("eps_x_pt_ycut1_periph_pr_h","eps_x_pt_ycut1_periph_pr_h",20,0.,2.,epsxymin,epsxymin);

  TProfile *eps_x_b_ycut1_ptcut1_pr_h = new TProfile("eps_x_b_ycut1_ptcut1_pr_h","eps_x_b_ycut1_ptcut1_pr_h",80,0.,200.,epsxymin,epsxymin);
  TProfile *eps_x_b_ycut1_ptcut2_pr_h = new TProfile("eps_x_b_ycut1_ptcut2_pr_h","eps_x_b_ycut1_ptcut2_pr_h",80,0.,200.,epsxymin,epsxymin);  

//eps_y
  TProfile *eps_y_rapidity_ptcut1_cent_pr_h = new TProfile("eps_y_rapidity_ptcut1_cent_pr_h","eps_y_rapidity_ptcut1_cent_pr_h",120,-6.,6.,epsxymin,epsxymin);
  TProfile *eps_y_rapidity_ptcut2_cent_pr_h = new TProfile("eps_y_rapidity_ptcut2_cent_pr_h","eps_y_rapidity_ptcut2_cent_pr_h",120,-6.,6.,epsxymin,epsxymin);  
  TProfile *eps_y_rapidity_ptcut1_mid_pr_h = new TProfile("eps_y_rapidity_ptcut1_mid_pr_h","eps_y_rapidity_ptcut1_mid_pr_h",120,-6.,6.,epsxymin,epsxymin);
  TProfile *eps_y_rapidity_ptcut2_mid_pr_h = new TProfile("eps_y_rapidity_ptcut2_mid_pr_h","eps_y_rapidity_ptcut2_mid_pr_h",120,-6.,6.,epsxymin,epsxymin);  
  TProfile *eps_y_rapidity_ptcut1_periph_pr_h = new TProfile("eps_y_rapidity_ptcut1_periph_pr_h","eps_y_rapidity_ptcut1_periph_pr_h",120,-6.,6.,epsxymin,epsxymin);
  TProfile *eps_y_rapidity_ptcut2_periph_pr_h = new TProfile("eps_y_rapidity_ptcut2_periph_pr_h","eps_y_rapidity_ptcut2_periph_pr_h",120,-6.,6.,epsxymin,epsxymin);  

  TProfile *eps_y_pt_ycut1_cent_pr_h = new TProfile("eps_y_pt_ycut1_cent_pr_h","eps_y_pt_ycut1_cent_pr_h",20,0.,2.,epsxymin,epsxymin);
  TProfile *eps_y_pt_ycut1_mid_pr_h = new TProfile("eps_y_pt_ycut1_mid_pr_h","eps_y_pt_ycut1_mid_pr_h",20,0.,2.,epsxymin,epsxymin);
  TProfile *eps_y_pt_ycut1_periph_pr_h = new TProfile("eps_y_pt_ycut1_periph_pr_h","eps_y_pt_ycut1_periph_pr_h",20,0.,2.,epsxymin,epsxymin);

  TProfile *eps_y_b_ycut1_ptcut1_pr_h = new TProfile("eps_y_b_ycut1_ptcut1_pr_h","eps_y_b_ycut1_ptcut1_pr_h",80,0.,200.,epsxymin,epsxymin);
  TProfile *eps_y_b_ycut1_ptcut2_pr_h = new TProfile("eps_y_b_ycut1_ptcut2_pr_h","eps_y_b_ycut1_ptcut2_pr_h",80,0.,200.,epsxymin,epsxymin);  

//rsq
double rsqmin=0.;
double rsqmax=100000.;
  TProfile *rsq_rapidity_ptcut1_cent_pr_h = new TProfile("rsq_rapidity_ptcut1_cent_pr_h","rsq_rapidity_ptcut1_cent_pr_h",120,-6.,6.,rsqmin,rsqmax);
  TProfile *rsq_rapidity_ptcut2_cent_pr_h = new TProfile("rsq_rapidity_ptcut2_cent_pr_h","rsq_rapidity_ptcut2_cent_pr_h",120,-6.,6.,rsqmin,rsqmax);  
  TProfile *rsq_rapidity_ptcut1_mid_pr_h = new TProfile("rsq_rapidity_ptcut1_mid_pr_h","rsq_rapidity_ptcut1_mid_pr_h",120,-6.,6.,rsqmin,rsqmax);
  TProfile *rsq_rapidity_ptcut2_mid_pr_h = new TProfile("rsq_rapidity_ptcut2_mid_pr_h","rsq_rapidity_ptcut2_mid_pr_h",120,-6.,6.,rsqmin,rsqmax);  
  TProfile *rsq_rapidity_ptcut1_periph_pr_h = new TProfile("rsq_rapidity_ptcut1_periph_pr_h","rsq_rapidity_ptcut1_periph_pr_h",120,-6.,6.,rsqmin,rsqmax);
  TProfile *rsq_rapidity_ptcut2_periph_pr_h = new TProfile("rsq_rapidity_ptcut2_periph_pr_h","rsq_rapidity_ptcut2_periph_pr_h",120,-6.,6.,rsqmin,rsqmax);  

  TProfile *rsq_pt_ycut1_cent_pr_h = new TProfile("rsq_pt_ycut1_cent_pr_h","rsq_pt_ycut1_cent_pr_h",20,0.,2.,rsqmin,rsqmax);
  TProfile *rsq_pt_ycut1_mid_pr_h = new TProfile("rsq_pt_ycut1_mid_pr_h","rsq_pt_ycut1_mid_pr_h",20,0.,2.,rsqmin,rsqmax);
  TProfile *rsq_pt_ycut1_periph_pr_h = new TProfile("rsq_pt_ycut1_periph_pr_h","rsq_pt_ycut1_periph_pr_h",20,0.,2.,rsqmin,rsqmax);

  TProfile *rsq_b_ycut1_ptcut1_pr_h = new TProfile("rsq_b_ycut1_ptcut1_pr_h","rsq_b_ycut1_ptcut1_pr_h",80,0.,200.,rsqmin,rsqmax);
  TProfile *rsq_b_ycut1_ptcut2_pr_h = new TProfile("rsq_b_ycut1_ptcut2_pr_h","rsq_b_ycut1_ptcut2_pr_h",80,0.,200.,rsqmin,rsqmax);  


  //-------------------------------------------
  
  TProfile *cth_1_eta_all_h = new TProfile("cth_1_eta_all_h","cth_1_eta_all_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_1_eta_pip_h = new TProfile("cth_1_eta_pip_h","cth_1_eta_pip_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_1_eta_pim_h = new TProfile("cth_1_eta_pim_h","cth_1_eta_pim_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_1_eta_kp_h = new TProfile("cth_1_eta_kp_h","cth_1_eta_kp_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_1_eta_km_h = new TProfile("cth_1_eta_km_h","cth_1_eta_km_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_1_eta_pr_h = new TProfile("cth_1_eta_pr_h","cth_1_eta_pr_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_1_eta_deut_h = new TProfile("cth_1_eta_deut_h","cth_1_eta_deut_h",120,-6.,6.,-1.,1.);  
  
  TProfile *cth_2_eta_all_h = new TProfile("cth_2_eta_all_h","cth_2_eta_all_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_2_eta_pip_h = new TProfile("cth_2_eta_pip_h","cth_2_eta_pip_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_2_eta_pim_h = new TProfile("cth_2_eta_pim_h","cth_2_eta_pim_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_2_eta_kp_h = new TProfile("cth_2_eta_kp_h","cth_2_eta_kp_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_2_eta_km_h = new TProfile("cth_2_eta_km_h","cth_2_eta_km_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_2_eta_pr_h = new TProfile("cth_2_eta_pr_h","cth_2_eta_pr_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_2_eta_deut_h = new TProfile("cth_2_eta_deut_h","cth_2_eta_deut_h",120,-6.,6.,-1.,1.);  

  TProfile *cth_3_eta_all_h = new TProfile("cth_3_eta_all_h","cth_3_eta_all_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_eta_pip_h = new TProfile("cth_3_eta_pip_h","cth_3_eta_pip_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_eta_pim_h = new TProfile("cth_3_eta_pim_h","cth_3_eta_pim_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_eta_kp_h = new TProfile("cth_3_eta_kp_h","cth_3_eta_kp_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_eta_km_h = new TProfile("cth_3_eta_km_h","cth_3_eta_km_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_eta_pr_h = new TProfile("cth_3_eta_pr_h","cth_3_eta_pr_h",120,-6.,6.,-1.,1.);  
  TProfile *cth_3_eta_deut_h = new TProfile("cth_3_eta_deut_h","cth_3_eta_deut_h",120,-6.,6.,-1.,1.);
  
  TH1D *costermh = new TH1D("costermh","costerm",100,-2,2);

  double ylimitextended= ybeam*2.0;
  if(ybeam>5.)ylimitextended= ybeam*1.5;
  TH1D *rapidityext_all_h = new TH1D("rapidityext_allh","rapidityext_allh",120,-ylimitextended,ylimitextended);
  TH1D *etaext_all_h = new TH1D("etaext_allh","etaext_allh",100,-ylimitextended,ylimitextended);
  TProfile *cth_1_rapidityext_all_h = new TProfile("cth_1_rapidityext_all_h","cth_1_rapidityext_all_h",120,-ylimitextended,ylimitextended,-1.,1.);
  TProfile *cth_1_etaext_all_h = new TProfile("cth_1_etaext_all_h","cth_1_etaext_all_h",120,-ylimitextended,ylimitextended,-1.,1.);
  TH1D *rapidityext_pn_h = new TH1D("rapidityext_pnh","rapidityext_pnh",120,-ylimitextended,ylimitextended);
  TH1D *etaext_pn_h = new TH1D("etaext_pnh","etaext_pnh",100,-ylimitextended,ylimitextended);
  TProfile *cth_1_rapidityext_pn_h = new TProfile("cth_1_rapidityext_pn_h","cth_1_rapidityext_pn_h",120,-ylimitextended,ylimitextended,-1.,1.);
  TProfile *cth_1_etaext_pn_h = new TProfile("cth_1_etaext_pn_h","cth_1_etaext_pn_h",120,-ylimitextended,ylimitextended,-1.,1.);
  TH1D *rapidityext_pnbar_h = new TH1D("rapidityext_pnbarh","rapidityext_pnbarh",120,-ylimitextended,ylimitextended);
  TH1D *etaext_pnbar_h = new TH1D("etaext_pnbarh","etaext_pnbarh",100,-ylimitextended,ylimitextended);
  TProfile *cth_1_rapidityext_pnbar_h = new TProfile("cth_1_rapidityext_pnbar_h","cth_1_rapidityext_pnbar_h",120,-ylimitextended,ylimitextended,-1.,1.);
  TProfile *cth_1_etaext_pnbar_h = new TProfile("cth_1_etaext_pnbar_h","cth_1_etaext_pnbar_h",120,-ylimitextended,ylimitextended,-1.,1.);
  
  infoh->Fill(1,itype);
  
  int ncount=1000;
  //    int ncount=100;
  //  int ncount=1;

  //  int itype=1; // ampt
  //int itype=2; // rqmd
    //int itype=4; // smash
      //int itype=5; // jam1

  //=================================1=======================================
  // read in event loop

  // load events

  // for JAM1

  TDatabasePDG *db = new TDatabasePDG();
  
  char infile [ 200 ]; 
  sprintf ( infile , "jam1data.root" );
  Int_t nfile = 0; TChain * chain = new TChain ( "jam" ); nfile += chain -> Add( infile );

  Int_t nentries = 0;
  if(itype==5){
  nentries = chain -> GetEntries ();
  std::cout << std::endl << "Added " << nfile << " files, " << "# of events is " << nentries << std::endl << std::endl;
  }
  
  if(!readit(itype,0,0)){
    std::cout<<" file data not found"<<std::endl;
    exit(1);
  }

  if(ntodo>1200000){
    ntodo=1200000;
    std::cout<<" because of array size max number of events is 1.2M"<<std::endl;
  }
  Int_t number_particles = 0;
  Int_t centrality = 0; 
  Int_t npart = 0; 
	//  std::ifstream file("amptdata.dat");
  int imax=9999999;
  if(itype==5)imax=nentries;
  //  std::cout<<" imax="<<imax<<std::endl;
  for(int i=0; i<imax; i++){   //begin event loop read in
    if(itype==5){
      //      std::cout<<" i="<<i<<std::endl;
      chain -> GetEntry ( i );
      // read track multiplicity
      TLeaf * leaf_refmult = chain -> GetLeaf ( "mul" );
      number_particles = leaf_refmult -> GetValue ( 0 );
      //      std::cout<<" number of tracks = "<< number_particles<<std::endl;
      // get event centrality bin
      TLeaf * leaf_centrality = chain -> GetLeaf ( "b" );
      b_impact = leaf_centrality -> GetValue ( 0 );
      //      std::cout<<" b_impact= "<<b_impact<<std::endl;
      // get event npart
      TLeaf * leaf_npart = chain -> GetLeaf ( "Npart" );
      npart_proj = leaf_npart -> GetValue ( 0 );
      //      std::cout<<" nparticipants (for jam1 labeled npart_prog, npart_targ=0)="<<npart_proj<<std::endl;

      /*
        // read TLeaves that store track parameters
        TLeaf * leaf_PID      = chain -> GetLeaf ( "pid" );
        TLeaf * leaf_Px       = chain -> GetLeaf ( "px" );
        TLeaf * leaf_Py       = chain -> GetLeaf ( "py" );
        TLeaf * leaf_Pz       = chain -> GetLeaf ( "pz" );
	TLeaf * leaf_ks   = chain -> GetLeaf ( "ks" );
	TLeaf * leaf_x    = chain -> GetLeaf ( "x" );
	TLeaf * leaf_y    = chain -> GetLeaf ( "y" );
	TLeaf * leaf_z    = chain -> GetLeaf ( "z" );
	TLeaf * leaf_t     = chain -> GetLeaf ( "t" );
	TLeaf * leaf_frx    = chain -> GetLeaf ( "frx" );
	TLeaf * leaf_fry    = chain -> GetLeaf ( "fry" );
	TLeaf * leaf_frz    = chain -> GetLeaf ( "frz" );
	TLeaf * leaf_frt     = chain -> GetLeaf ( "frt" );
      */
      
	
      //      chain->Show();
    }

        // read TLeaves that store track parameters
        TLeaf * leaf_PID      = chain -> GetLeaf ( "pid" );
        TLeaf * leaf_Px       = chain -> GetLeaf ( "px" );
        TLeaf * leaf_Py       = chain -> GetLeaf ( "py" );
        TLeaf * leaf_Pz       = chain -> GetLeaf ( "pz" );
	TLeaf * leaf_ks   = chain -> GetLeaf ( "ks" );
	TLeaf * leaf_x    = chain -> GetLeaf ( "x" );
	TLeaf * leaf_y    = chain -> GetLeaf ( "y" );
	TLeaf * leaf_z    = chain -> GetLeaf ( "z" );
	TLeaf * leaf_t     = chain -> GetLeaf ( "t" );
	TLeaf * leaf_frx    = chain -> GetLeaf ( "frx" );
	TLeaf * leaf_fry    = chain -> GetLeaf ( "fry" );
	TLeaf * leaf_frz    = chain -> GetLeaf ( "frz" );
	TLeaf * leaf_frt     = chain -> GetLeaf ( "frt" );

	//    std::cout<<" 1 "<<std::endl;


    
    if(!readit(itype,1,0))break; // read next event
    //    std::cout<<"  test b_impact bmin bmax "<<b_impact<<" "<<bmin<<" "<<bmax<<" i="<<i<<std::endl;
    //    if(b_impact<bmin || b_impact>bmax)continue; // this does funny things
    
    nev++;
    //    std::cout<<" 3 nentries="<<nentries<<std::endl;
    //    if ( (itype==2) && (nev + 1 > nentries) ) break;     
    //    std::cout<<" 3.1 "<<std::endl;
    if(nev%ncount==0)std::cout<<" events read="<<nev<<std::endl;
    if(ntodo>0 && nev > ntodo)break;

    //    std::cout<<" num particles="<<number_particles<<std::endl;

    //    std::cout<<" 2 "<<std::endl;
    
    std::vector <int> part_id;
    std::vector <double> px;
    std::vector <double> py;
    std::vector <double> pz;
    std::vector <double> mass;
    std::vector <double> x;
    std::vector <double> y;
    std::vector <double> z;
    std::vector <double> t;
    std::vector <double> energy;
    std::vector <double> pt;
    std::vector <double> rapidity;
    std::vector <double> phi;
    std::vector <double> eta;

    std::vector <double> form_time;    
    std::vector <double> t_last_coll;
    std::vector <double> ncollpart;

    // generate uniformly randomly distributed "reaction plane"
    if(irot>0)rot = ran0->Rndm()*2.0 * PI-PI;
    rot_v.push_back(rot);

    //    std::cout<<" number_particles="<<number_particles<<std::endl;    

    //    std::cout<<" 3 "<<std::endl;

    refmultfxt=0;
    for (int i=0; i< number_particles; i++){
      Int_t itrack=i;
      Int_t       part_id_in = -1000.; 
	Int_t charge_in=0;
	Double_t pmag_in =  -1000.; 
	Double_t energy_in=-10000;
	Double_t theta_in=-1000.;
	Double_t eta_in = -1000.; 
	Double_t phi_in = -1000.;
	Double_t pt_in = -1000.; 
	Double_t rapidity_in = -1000.; 

	//	std::cout<<" read particle number = "<<itrack<<std::endl;
      if(itype==5){

	part_id_in = leaf_PID -> GetValue ( itrack );

	if(part_id_in>10000 && part_id_in !=1001001000 )cout<<"   debug 1b  - 2 part_id="<<part_id_in<<endl; // debug 1b
	//	cout<<" pid="<<part_id_in<<endl;
	if(part_id_in==1001001000){      //jam1 pid for deuteron
	  part_id_in = 6201;  // redo id of deuteron to 6201
	  mass_in=1.87705;    //deuteron mass GeV
	}else if(part_id_in==30333) {
	  // phi(1680)
	  part_id_in = 6202;  // redo id of phi 1680  to 6202
	  mass_in=1.680;    // mass GeV
	}else if(part_id_in==50223) {
	  // f1(1510)
	  part_id_in = 6203;  // redo id of phi 1680  to 6202
	  mass_in=1.518;    // mass GeV
	}else{
	  TParticlePDG *part = db->GetParticle(part_id_in);
	  mass_in = part->Mass();
	}

	//	std::cout<<" part_id_in="<<part_id_in<<" pdgmass="<<mass_in<<std::endl;
	Int_t charge_in=0;
	if(part_id_in>0)charge_in=1;
	if(part_id_in<0)charge_in=-1;
	
	px_in = leaf_Px -> GetValue ( itrack );
	py_in = leaf_Py -> GetValue ( itrack );
	pz_in = leaf_Pz -> GetValue ( itrack );
	pt_in = sqrt(px_in*px_in+py_in*py_in);
	pmag_in =  sqrt(px_in*px_in+py_in*py_in+pz_in*pz_in);
	energy_in=-10000;
	if(mass_in>-1000.) energy_in = sqrt(mass_in*mass_in+pmag_in*pmag_in);
	TLorentzVector particle(px_in,py_in,pz_in,energy_in);
	theta_in=particle.Theta();
	eta_in = particle.PseudoRapidity();
	phi_in = particle.Phi();
	rapidity_in = particle.Rapidity();
	
	Double_t     ax = leaf_x -> GetValue ( itrack );
	Double_t     ay = leaf_y -> GetValue ( itrack );
	Double_t     az = leaf_z -> GetValue ( itrack );
	Double_t     at = leaf_t -> GetValue ( itrack );
	
	Double_t     frx = leaf_frx -> GetValue ( itrack );
	Double_t     fry = leaf_fry -> GetValue ( itrack );
	Double_t     frz = leaf_frz -> GetValue ( itrack );
	Double_t     frt = leaf_frt -> GetValue ( itrack );
	
	x_in=ax;
	y_in=ay;
	z_in=az;
	t_in=at;


	//count refmultfxt
	double rapidity_lab=rapidity_in+ybeam/2.;
	double eta_lab=rapidity_lab;
	double etaminfxt=-2.;
	double etamaxfxt=0.;
	double ptminfxt=0.200;
	if( pt_in > ptminfxt && eta_lab>etaminfxt && eta_lab<etamaxfxt)refmultfxt++; 
	//	if(pt_in<0.01 && energy_in>0.)cout<<" x_in="<<x_in<<" pt_in="<<pt_in<<" emergy_in="<<energy_in<<" b_impact="<<b_impact<<endl;
      } // end jam1 track read
      //     cout<<" refmultfxt="<<refmultfxt<<"nbaryons="<<nbaryons<<" b_impact="<<b_impact<<" number_particles="<<number_particles<<endl;
      //      file>>part_id_in>>px_in>>py_in>>pz_in>>mass_in>>x_in>>y_in>>z_in>>t_in;
      //      std::cout<<" i th particle="<<i<<std::endl;
      if(!readit(itype,2,i)){     //read one particle
	std::cout<<" something went wrong. reached end of file in middle of event, exiting"<<std::endl;
	exit(1);
      }
      //      std::cout<<" event="<<nev<<" i="<<i<<" part_id_in="<<part_id_in<<" number_partices="<<number_particles<<std::endl;

      // rotate the event around z axis if irot = 1

      if(irot>0){
	double px_inp = px_in*cos(rot)-py_in*sin(rot);
	double py_inp = px_in*sin(rot)+py_in*cos(rot);
	px_in=px_inp;
	py_in=py_inp;
	double x_inp = x_in*cos(rot)-y_in*sin(rot);
	double y_inp = x_in*sin(rot)+y_in*cos(rot);
	x_in=x_inp;
	y_in=y_inp;
	reactionplanerotationh->Fill(rot);  
      }

      // injecting flow into system     
      if(i_inject>0){
	// inject flow
	double PSI_injected_direction=0.;
	double rot_flow_injected=0;
	double phi_orig=atan2(py_in,px_in);
	for(int nn_injected=0; nn_injected<6; nn_injected++){
	  int n_injected=nn_injected+1;
	  rot_flow_injected+= -2/n_injected*v_injected[nn_injected]*sin(n_injected*(phi_orig-PSI_injected_direction));
	}
	double px_inp = px_in*cos(rot_flow_injected)-py_in*sin(rot_flow_injected);
	double py_inp = px_in*sin(rot_flow_injected)+py_in*cos(rot_flow_injected);
	px_in=px_inp;
	py_in=py_inp;
	double x_inp = x_in*cos(rot_flow_injected)-y_in*sin(rot_flow_injected);
	double y_inp = x_in*sin(rot_flow_injected)+y_in*cos(rot_flow_injected);
	x_in=x_inp;
	y_in=y_inp;
      }
      // end - injecting flow

	      
      part_id.push_back(part_id_in);
      px.push_back(px_in);
      py.push_back(py_in);
      pz.push_back(pz_in);
      mass.push_back(mass_in);
      x.push_back(x_in);
      y.push_back(y_in);
      z.push_back(z_in);
      t.push_back(t_in);
      
      energy.push_back(sqrt(px_in*px_in+py_in*py_in+pz_in*pz_in+mass_in*mass_in));
      pt.push_back(sqrt(px_in*px_in+py_in*py_in));
      if(itype==5){
	rapidity.push_back(rapidity_in);
      }else{
	rapidity.push_back(0.5*log((energy[i]+pz_in)/(energy[i]-pz_in)));
      }
      phi.push_back(atan2(py_in,px_in));
      double theta=acos(pz_in/sqrt(px_in*px_in+py_in*py_in+pz_in*pz_in));
      eta.push_back(-log(tan(theta/2.)));
      //debug      std::cout<<" going to next particles i="<<i<<std::endl;

      // when you inject flow we ignore t_last_coll, ncoll, formation time !!!!

      form_time.push_back(form_time_in);
      t_last_coll.push_back(t_last_coll_in);
      ncollpart.push_back(ncoll_in);


      //      std::cout<<" event="<<nev<<" tracknum="<<i<<" p="<<px_in<<" "<<py_in<<" "<<pz_in<<" xyz="<<x_in<<" "<<y_in<<" "<<z_in<<" t_in="<<t_in<<" b="<<b_impact<<" mass="<<mass_in<<std::endl;
      
    } // end particle read in loop
        
    eventnum_v.push_back(eventnum);
    number_particles_v.push_back(number_particles);
    refmultfxt_v.push_back(refmultfxt);
    b_impact_v.push_back(b_impact);
    npart_proj_v.push_back(npart_proj);
    npart_targ_v.push_back(npart_targ);
    nelas_proj_v.push_back(nelas_proj);
    ninelas_proj_v.push_back(ninelas_proj);
    nelas_targ_v.push_back(nelas_targ);
    ninelas_targ_v.push_back(ninelas_targ);
    part_id_v.push_back(part_id);
    px_v.push_back(px);
    py_v.push_back(py);
    pz_v.push_back(pz);
    mass_v.push_back(mass);
    x_v.push_back(x);
    y_v.push_back(y);
    z_v.push_back(z);
    t_v.push_back(t);
    energy_v.push_back(energy);
    pt_v.push_back(pt);
    rapidity_v.push_back(rapidity);
    phi_v.push_back(phi);
    eta_v.push_back(eta);        

    form_time_v.push_back(form_time);
    t_last_coll_v.push_back(t_last_coll);
    ncollpart_v.push_back(ncollpart);

      if(!readit(itype,3,i)){std::cout<<" something went wrong reading b_impact"<<std::endl; break;}     //at end of event for SMASH read one line to gete b_impact

  } // event loop read in
  readit(itype,-1,0);
  std::cout<<" FINISHED READING IN FILE, number of events="<<nev<<std::endl;
  // nev might be 1 more than it should be!!

  nev--;  
  
  //===============================2=================================== 
  // define a bunch of stuff outside the loops

  double PSI[6]={0}; // reaction plane
  
  std::vector <int> part_id;
  std::vector <double> px;
  std::vector <double> py;
  std::vector <double> pz;
  std::vector <double> mass;
  std::vector <double> x;
  std::vector <double> y;
  std::vector <double> z;
  std::vector <double> t;
  std::vector <double> energy;
  std::vector <double> pt;
  std::vector <double> rapidity;
  std::vector <double> phi;
  std::vector <double> eta;

  std::vector <double> form_time;
  std::vector <double> t_last_coll;
  std::vector <double> ncollpart;

  int nevbmin=0;
  int nevbmax=0;
  
  // event loop 2   plotting stuff and finding reaction plane

  //  std::cout<<" 2 nev="<<nev<<" ncount="<<ncount<<std::endl;

  for(int iev=0; iev<nev; iev++){
    if(iev%ncount==0)std::cout<<" 2 eventnum="<<iev<<std::endl;
    if(ntodo>0 && iev > ntodo)break;

    eventnum=eventnum_v[iev];
    number_particles=number_particles_v[iev];
    refmultfxt=refmultfxt_v[iev];
    b_impact=b_impact_v[iev];
    npart_proj=npart_proj_v[iev];
    npart_targ=npart_targ_v[iev];
    nelas_proj=nelas_proj_v[iev];
    ninelas_proj=ninelas_proj_v[iev];
    nelas_targ=nelas_targ_v[iev];
    ninelas_targ=ninelas_targ_v[iev];
    rot = rot_v[iev];
    part_id=part_id_v[iev];
    px=px_v[iev];
    py=py_v[iev];
    pz=pz_v[iev];
    mass=mass_v[iev];
    x=x_v[iev];
    y=y_v[iev];
    z=z_v[iev];
    t=t_v[iev];
    energy=energy_v[iev];
    pt=pt_v[iev];
    rapidity=rapidity_v[iev];
    phi=phi_v[iev];
    eta=eta_v[iev];

    form_time=form_time_v[iev];
    t_last_coll=t_last_coll_v[iev];
    ncollpart=ncollpart_v[iev];
    
    if(b_impact<bmin)nevbmin++;
    if(b_impact<bmax)nevbmax++;

    //    std::cout<<" filling b_impact="<<b_impact<<std::endl;

    b_impactallh->Fill(b_impact);
    
    if(b_impact<bmin || b_impact>bmax)continue; 
    
    number_particlesh->Fill(number_particles);
    b_impacth->Fill(b_impact);
    bsq_impacth->Fill(b_impact*b_impact);
    n_participantsh->Fill(npart_proj+npart_targ);

    double ptmin= 0.001;// no cut // 0.0001; // a cut of p=0 is the same as setting tmin=0.5
    double tmin=0.5; // to cut spectators set tmin=0.5  - don't make cuts on t fro RQMD
    double tmax=9999999; // rqmd seems to set this to 200 which is set in the inputfile 
   //    double bmin=0;
    //    double bmax=14.;

    //    bmin=3;
    //    bmax=5;

    //    double etamin=0.5;
    double etamin=ybeam/2.;
    //    double etamin=ybeam;

    int useonlynucleons = 1e9;  // to use only nucleons set to 1e9. to use any particles set to -1e9   
    int nbaryons=0;
  
    // plot stuff   
    for (int i=0; i< number_particles; i++){

      double cth1=cos(phi[i]);
      double cth2=cos(2*phi[i]);
      double cth3=cos(3*phi[i]);
      
      rapidityptrawh->Fill(rapidity[i],pt[i]);
  
      xvsyallh->Fill(x[i],y[i]);
      zvstallh->Fill(z[i],t_last_coll[i]);

      if(pt[i]<0.315 && abs(rapidity[i])>0.95*ybeam && abs(rapidity[i])<1.05*ybeam){      
	xvsy_ptsmall_ybeamh->Fill(x[i],y[i]);
      }else{
	xvsy_not_ptsmall_ybeamh->Fill(x[i],y[i]); 
      }
      
      if(pt[i]<0.3 && abs(rapidity[i])>0.8*ybeam){
	xvsyspech->Fill(x[i],y[i]);
	xvszspech->Fill(x[i],z[i]);
	zvstspech->Fill(z[i],t_last_coll[i]);
	ncollpartspech->Fill(ncollpart[i]);
      }
      if(abs(eta[i])>ybeam)part_idybeamh->Fill(part_id[i]);

      double ddy=2.*ybeam/9.;
      int iybin=rapidity[i]/ddy;
      int iptbin=-1;
      if(pt[i]>0.)iptbin=0;
      if(pt[i]>0.100)iptbin=1;
      if(pt[i]>0.250)iptbin=2;
      if(pt[i]>1.000)iptbin=3;

      //      cout<<" ddy= "<<ddy<<" y= "<<rapidity[i]<<" iybin= "<<iybin<<endl;

      //      double ptmin2=ptmin;
      double ptmin2=0.2; //to cut spectators

      if(abs(part_id[i])==2212 || abs(part_id[i])==2112  ){

      if(pt[i]>ptmin2){                   // cut to plot participants - do I want a y> 0.5*ybeam or a pid=proton/neutron cut?
	if(iybin==-4)xvsy_p_binm4->Fill(x[i],y[i]);
	if(iybin==-3)xvsy_p_binm3->Fill(x[i],y[i]);
	if(iybin==-4)xvsy_p_binm2->Fill(x[i],y[i]);
	if(iybin==-1)xvsy_p_binm1->Fill(x[i],y[i]);
	if(iybin==0)xvsy_p_bin0->Fill(x[i],y[i]);
	if(iybin==1)xvsy_p_binp1->Fill(x[i],y[i]);
	if(iybin==2)xvsy_p_binp2->Fill(x[i],y[i]);
	if(iybin==3)xvsy_p_binp3->Fill(x[i],y[i]);
	if(iybin==4)xvsy_p_binp4->Fill(x[i],y[i]);
      }

      xvsrapidity_allnucleonsh->Fill(rapidity[i],x[i]);
      if(iybin==-4 && iptbin==0)xvsy_p_pt0_binm4->Fill(x[i],y[i]);
      if(iybin==-3 && iptbin==0)xvsy_p_pt0_binm3->Fill(x[i],y[i]);
      if(iybin==-2 && iptbin==0)xvsy_p_pt0_binm2->Fill(x[i],y[i]);
      if(iybin==-1 && iptbin==0)xvsy_p_pt0_binm1->Fill(x[i],y[i]);
      if(iybin==0 && iptbin==0)xvsy_p_pt0_bin0->Fill(x[i],y[i]);
      if(iybin==1 && iptbin==0)xvsy_p_pt0_binp1->Fill(x[i],y[i]);
      if(iybin==2 && iptbin==0)xvsy_p_pt0_binp2->Fill(x[i],y[i]);
      if(iybin==3 && iptbin==0)xvsy_p_pt0_binp3->Fill(x[i],y[i]);
      if(iybin==4 && iptbin==0)xvsy_p_pt0_binp4->Fill(x[i],y[i]);
      
      if(iybin==-4 && iptbin==1)xvsy_p_pt1_binm4->Fill(x[i],y[i]);
      if(iybin==-3 && iptbin==1)xvsy_p_pt1_binm3->Fill(x[i],y[i]);
      if(iybin==-2 && iptbin==1)xvsy_p_pt1_binm2->Fill(x[i],y[i]);
      if(iybin==-1 && iptbin==1)xvsy_p_pt1_binm1->Fill(x[i],y[i]);
      if(iybin==0 && iptbin==1)xvsy_p_pt1_bin0->Fill(x[i],y[i]);
      if(iybin==1 && iptbin==1)xvsy_p_pt1_binp1->Fill(x[i],y[i]);
      if(iybin==2 && iptbin==1)xvsy_p_pt1_binp2->Fill(x[i],y[i]);
      if(iybin==3 && iptbin==1)xvsy_p_pt1_binp3->Fill(x[i],y[i]);
      if(iybin==4 && iptbin==1)xvsy_p_pt1_binp4->Fill(x[i],y[i]);
      
      if(iybin==-4 && iptbin==2)xvsy_p_pt2_binm4->Fill(x[i],y[i]);
      if(iybin==-3 && iptbin==2)xvsy_p_pt2_binm3->Fill(x[i],y[i]);
      if(iybin==-2 && iptbin==2)xvsy_p_pt2_binm2->Fill(x[i],y[i]);
      if(iybin==-1 && iptbin==2)xvsy_p_pt2_binm1->Fill(x[i],y[i]);
      if(iybin==0 && iptbin==2)xvsy_p_pt2_bin0->Fill(x[i],y[i]);
      if(iybin==1 && iptbin==2)xvsy_p_pt2_binp1->Fill(x[i],y[i]);
      if(iybin==2 && iptbin==2)xvsy_p_pt2_binp2->Fill(x[i],y[i]);
      if(iybin==3 && iptbin==2)xvsy_p_pt2_binp3->Fill(x[i],y[i]);
      if(iybin==4 && iptbin==2)xvsy_p_pt2_binp4->Fill(x[i],y[i]);
      
      if(iybin==-4 && iptbin==3)xvsy_p_pt3_binm4->Fill(x[i],y[i]);
      if(iybin==-3 && iptbin==3)xvsy_p_pt3_binm3->Fill(x[i],y[i]);
      if(iybin==-2 && iptbin==3)xvsy_p_pt3_binm2->Fill(x[i],y[i]);
      if(iybin==-1 && iptbin==3)xvsy_p_pt3_binm1->Fill(x[i],y[i]);
      if(iybin==0 && iptbin==3)xvsy_p_pt3_bin0->Fill(x[i],y[i]);
      if(iybin==1 && iptbin==3)xvsy_p_pt3_binp1->Fill(x[i],y[i]);
      if(iybin==2 && iptbin==3)xvsy_p_pt3_binp2->Fill(x[i],y[i]);
      if(iybin==3 && iptbin==3)xvsy_p_pt3_binp3->Fill(x[i],y[i]);
      if(iybin==4 && iptbin==3)xvsy_p_pt3_binp4->Fill(x[i],y[i]);

      }

      if(abs(part_id[i])==2212 ||abs(part_id[i])==2112 ||abs(part_id[i])==3122 ||abs(part_id[i])==2214 ||abs(part_id[i])==2114 ||abs(part_id[i])==1114 ||abs(part_id[i])==2224 )nbaryons++;
      //*******************************************
      if(pt[i]<ptmin)continue; // kill pt=0 spectators
      if(pt[i]<0.3 && abs(rapidity[i])>0.8*ybeam)continue; // kill spectators
      if(abs(rapidity[i])>0.5*ybeam)continue; // kill particles not in central rapidity

      
      rapiditypth->Fill(rapidity[i],pt[i]);      
      if(part_id[i]==211)rapiditypt_piph->Fill(rapidity[i],pt[i]);      
      if(part_id[i]==-211)rapiditypt_pimh->Fill(rapidity[i],pt[i]);      
      if(part_id[i]==321)rapiditypt_kph->Fill(rapidity[i],pt[i]);      
      if(part_id[i]==-321)rapiditypt_kmh->Fill(rapidity[i],pt[i]);      
      if(part_id[i]==2212)rapiditypt_prh->Fill(rapidity[i],pt[i]);      

      rapidityext_all_h->Fill(rapidity[i]);
      etaext_all_h->Fill(eta[i]);
      if(part_id[i]==2212 ||part_id[i]==2112)rapidityext_pn_h->Fill(rapidity[i]);
      if(part_id[i]==2212 ||part_id[i]==2112)etaext_pn_h->Fill(eta[i]);
      if(part_id[i]==-2212 ||part_id[i]==-2112)rapidityext_pnbar_h->Fill(rapidity[i]);
      if(part_id[i]==-2212 ||part_id[i]==-2112)etaext_pnbar_h->Fill(eta[i]);
      
      if(part_id[i]==2112||part_id[i]==2212 || part_id[i]>useonlynucleons) {

	// study all particles incl spectators      
	// pt=0 then it gets rid of anything with t=0.2
	//	if(pt[i]<0.0001)std::cout<<" pt="<<pt[i]<<" t="<<t[i]<<" iev="<<iev<<" i="<<i<<std::endl;
	//	if(t[i]==0.2)std::cout<<" pt="<<pt[i]<<" t="<<t[i]<<" iev="<<iev<<" i="<<i<<std::endl;
	if(pt[i]<ptmin)continue;
	if(t[i]>tmax)continue; // AMPT only follows stuff 400*0.2 fm =80 fm
	ptpzh->Fill(pt[i],pz[i] );    
	etarapidityh->Fill(eta[i],rapidity[i]);
	xvsz_allnh->Fill(x[i],z[i]);
	xvsy_allnh->Fill(x[i],y[i]);     
	xvst_allnh->Fill(x[i],t_last_coll[i]);     
	ptvst_allnh->Fill(pt[i],t_last_coll[i]); 
	etavst_allnh->Fill(eta[i],t_last_coll[i]); 
	rapidityvst_allnh->Fill(rapidity[i],t_last_coll[i]); 
	pxvst_allnh->Fill(px[i],t_last_coll[i]);
	pyvst_allnh->Fill(py[i],t_last_coll[i]);
	pzvst_allnh->Fill(pz[i],t_last_coll[i]);
	cth1vst_allnh->Fill(t[i],cth1);
	cth2vst_allnh->Fill(t[i],cth2);
	cth3vst_allnh->Fill(t[i],cth3);

	if(eta[i]>etamin){
	  pxvstposrap_allnh->Fill(px[i],t_last_coll[i]);
	  pyvstposrap_allnh->Fill(py[i],t_last_coll[i]);
	  pzvstposrap_allnh->Fill(pz[i],t_last_coll[i]);
	  cth1vstposrap_allnh->Fill(t[i],cth1);
	  cth2vstposrap_allnh->Fill(t[i],cth2);
	  cth3vstposrap_allnh->Fill(t[i],cth3);
	  
	  if(eta[i]>ybeam){
	    pxvstybeam_allnh->Fill(px[i],t_last_coll[i]);
	    pyvstybeam_allnh->Fill(py[i],t_last_coll[i]);
	    pzvstybeam_allnh->Fill(pz[i],t_last_coll[i]);
	    cth1vstybeam_allnh->Fill(t[i],cth1);
	    cth2vstybeam_allnh->Fill(t[i],cth2);
	    cth3vstybeam_allnh->Fill(t[i],cth3);
	  }
	}
      }

      if(form_time[i]<0.)xvsy_form_time_h->Fill(x[i],y[i]);
      if(t_last_coll[i]<0.4)xvsy_t_last_coll_h->Fill(x[i],y[i]);
      if(ncollpart[i]<0.5)xvsy_ncollpart_h->Fill(x[i],y[i]);

      if(pt[i]<ptmin)continue; // only plot participants
      if(t[i]>tmax)continue; // AMPT only follows stuff 400*0.2 fm =80 fm
      if(abs(rapidity[i])>0.5*ybeam)continue; // URQMD - really kill spectators
      
      //      if(abs(eta[i])>5.)continue;
      //      if(phi[i]==0.)std::cout<<" phi="<<i<<" "<<phi[i]<<" "<<px[i]<<" "<<py[i]<<" "<<pz[i]<<" eta="<<eta[i]<<std::endl;
      xvsyh->Fill(x[i],y[i]);
      
      if(eta[i]>etamin && (part_id[i]==2112||part_id[i]==2212|| part_id[i]>useonlynucleons) && b_impact>bmin && b_impact<bmax) xvsy_p_cuth->Fill(x[i],y[i]);
      if(eta[i]<-etamin && (part_id[i]==2112||part_id[i]==2212|| part_id[i]>useonlynucleons) && b_impact>bmin && b_impact<bmax) xvsy_p_cut2h->Fill(x[i],y[i]);  
      if(abs(eta[i])>etamin && (part_id[i]==2112||part_id[i]==2212|| part_id[i]>useonlynucleons))xvszh->Fill(x[i],z[i]);
      if(abs(eta[i])>ybeam && (part_id[i]==2112||part_id[i]==2212|| part_id[i]>useonlynucleons))xvszybeamh->Fill(x[i],z[i]);
      if(abs(eta[i])>ybeam && (part_id[i]==2112||part_id[i]==2212|| part_id[i]>useonlynucleons))xvsyybeamh->Fill(x[i],y[i]);
      if(eta[i]>etamin && (part_id[i]==2112||part_id[i]==2212|| part_id[i]>useonlynucleons) && b_impact>bmin && b_impact<bmax)phi_p_cuth->Fill(phi[i]);
      if(eta[i]<-etamin && (part_id[i]==2112||part_id[i]==2212|| part_id[i]>useonlynucleons) && b_impact>bmin && b_impact<bmax)phi_p_cut2h->Fill(phi[i]);

      if(eta[i]>0.)phi_p_posraph->Fill(phi[i]);  // AMPT uses etamin
      if(eta[i]<-0.)phi_p_negraph->Fill(phi[i]); // AMPT uses etamin
      if(eta[i]>0.)xvsy_p_posraph->Fill(x[i],y[i]); // AMPT uses etamin
      if(eta[i]<-0.)xvsy_p_negraph->Fill(x[i],y[i]); // AMPT uses etamin
      
      zvsth->Fill(z[i],t_last_coll[i]);
      xvsth->Fill(x[i],t_last_coll[i]);
      yvsth->Fill(y[i],t_last_coll[i]);
      rvsth->Fill(sqrt(x[i]*x[i]+y[i]*y[i]),t_last_coll[i]);
      part_idh->Fill(part_id[i]);
      part_id2h->Fill(part_id[i]);
      pth->Fill(pt[i]);
      rapidityh->Fill(rapidity[i]);
      phih->Fill(phi[i]);
      pxh->Fill(px[i]);
      pyh->Fill(py[i]);
      pzh->Fill(pz[i]);
      xh->Fill(log10(x[i]));
      yh->Fill(log10(y[i]));
      zh->Fill(log10(z[i]));
      th->Fill(log10(t[i]));
      etah->Fill(eta[i]);	

      ncollparth->Fill(ncollpart[i]);
      formtimeh->Fill(form_time[i]);
      lastcolltimeh->Fill(t_last_coll[i]);

    } // end particle loop

    fxtmulth->Fill(refmultfxt);        
    number_particles_bsqh->Fill(number_particles,b_impact*b_impact);
    fxtmult_bsqh->Fill(refmultfxt,b_impact*b_impact);    
    nparticipants_bsqh->Fill(npart_proj+npart_targ,b_impact*b_impact); 
    
    //      cout<<" refmultfxt="<<refmultfxt<<" nbaryons="<<nbaryons<<" b_impact="<<b_impact<<" number_particles="<<number_particles<<endl;
    //-------- now calculate <cth> -------------------------------------

    for (int i=0; i< number_particles; i++){
      
      double cth1=cos(phi[i]);
      double cth2=cos(2*phi[i]);
      double cth3=cos(3*phi[i]);
      double sth3=sin(3*phi[i]);

      double ddy=2.*ybeam/9.;
      int iybin=rapidity[i]/ddy;
      int iptbin=-1;
      if(pt[i]>0.)iptbin=0;
      if(pt[i]>0.100)iptbin=1;
      if(pt[i]>0.250)iptbin=2;
      if(pt[i]>1.000)iptbin=3;

      // recenter (from xy plots)
      double xx0=0;
      double yy0=0;
      if(iybin==-4 && iptbin==0){ xx0=-6.355;  yy0=0.;}
      if(iybin==-3 && iptbin==0){ xx0=-6.241;  yy0=0.;}
      if(iybin==-2 && iptbin==0){ xx0=-3.480;  yy0=0.;}
      if(iybin==-1 && iptbin==0){ xx0=-1.422;  yy0=0.;}
      if(iybin==0 && iptbin==0){ xx0=0.;  yy0=0.;}
      if(iybin==1 && iptbin==0){ xx0=1.422;  yy0=0.;}
      if(iybin==2 && iptbin==0){ xx0=3.480;  yy0=0.;}
      if(iybin==3 && iptbin==0){ xx0=6.241;  yy0=0.;}
      if(iybin==4 && iptbin==0){ xx0=6.355;  yy0=0.;}				      

      if(iybin==-4 && iptbin==1){ xx0=-6.521;  yy0=0.;}
      if(iybin==-3 && iptbin==1){ xx0=-6.361;  yy0=0.;}
      if(iybin==-2 && iptbin==1){ xx0=-3.675;  yy0=0.;}
      if(iybin==-1 && iptbin==1){ xx0=-1.754;  yy0=0.;}
      if(iybin==0 && iptbin==1){ xx0=0.;  yy0=0.;}
      if(iybin==1 && iptbin==1){ xx0=1.754;  yy0=0.;}
      if(iybin==2 && iptbin==1){ xx0=3.675;  yy0=0.;}
      if(iybin==3 && iptbin==1){ xx0=6.361;  yy0=0.;}
      if(iybin==4 && iptbin==1){ xx0=6.521;  yy0=0.;}				      

      if(iybin==-4 && iptbin==2){ xx0=-7.991;  yy0=0.;}
      if(iybin==-3 && iptbin==2){ xx0=-7.814;  yy0=0.;}
      if(iybin==-2 && iptbin==2){ xx0=-5.786;  yy0=0.;}
      if(iybin==-1 && iptbin==2){ xx0=-3.332;  yy0=0.;}
      if(iybin==0 && iptbin==2){ xx0=0.;  yy0=0.;}
      if(iybin==1 && iptbin==2){ xx0=3.332;  yy0=0.;}
      if(iybin==2 && iptbin==2){ xx0=5.786;  yy0=0.;}
      if(iybin==3 && iptbin==2){ xx0=7.814;  yy0=0.;}
      if(iybin==4 && iptbin==2){ xx0=7.991;  yy0=0.;}				      

      if(iybin==-4 && iptbin==3){ xx0=-13.617;  yy0=0.;}
      if(iybin==-3 && iptbin==3){ xx0=-11.830;  yy0=0.;}
      if(iybin==-2 && iptbin==3){ xx0=-8.495;  yy0=0.;}
      if(iybin==-1 && iptbin==3){ xx0=-4.700;  yy0=0.;}
      if(iybin==0 && iptbin==3){ xx0=0.;  yy0=0.;}
      if(iybin==1 && iptbin==3){ xx0=4.700;  yy0=0.;}
      if(iybin==2 && iptbin==3){ xx0=8.495;  yy0=0.;}
      if(iybin==3 && iptbin==3){ xx0=11.830;  yy0=0.;}
      if(iybin==4 && iptbin==3){ xx0=13.617;  yy0=0.;}				      
      double rsqmax=0;
      if(iybin==-4 && iptbin==0){ rsqmax=6.075 ;} 
      if(iybin==-3 && iptbin==0){ rsqmax=6.178 ;}   
      if(iybin==-2 && iptbin==0){ rsqmax=7.638 ;}   
      if(iybin==-1 && iptbin==0){ rsqmax=7.475 ;}   
      if(iybin==0 && iptbin==0){ rsqmax=7.556 ;}    
      if(iybin==1 && iptbin==0){ rsqmax=7.475 ;}    
      if(iybin==2 && iptbin==0){ rsqmax=7.638 ;}    
      if(iybin==3 && iptbin==0){ rsqmax=6.178 ;}    
      if(iybin==4 && iptbin==0){ rsqmax=6.075 ;}    
      
      if(iybin==-4 && iptbin==1){ rsqmax=6.394 ;}   
      if(iybin==-3 && iptbin==1){ rsqmax=7.241 ;}   
      if(iybin==-2 && iptbin==1){ rsqmax=11.782 ;}   
      if(iybin==-1 && iptbin==1){ rsqmax=13.167 ;}   
      if(iybin==0 && iptbin==1){ rsqmax=13.822 ;}    
      if(iybin==1 && iptbin==1){ rsqmax=13.167 ;}    
      if(iybin==2 && iptbin==1){ rsqmax=11.782 ;}    
      if(iybin==3 && iptbin==1){ rsqmax=7.241 ;}    
      if(iybin==4 && iptbin==1){ rsqmax=6.394 ;}    

      if(iybin==-4 && iptbin==2){ rsqmax=11.842 ;}   
      if(iybin==-3 && iptbin==2){ rsqmax=18.870 ;}   
      if(iybin==-2 && iptbin==2){ rsqmax=26.792 ;}   
      if(iybin==-1 && iptbin==2){ rsqmax=31.493 ;}   
      if(iybin==0 && iptbin==2){ rsqmax=34.015 ;}    
      if(iybin==1 && iptbin==2){ rsqmax=31.493 ;}    
      if(iybin==2 && iptbin==2){ rsqmax=26.792 ;}    
      if(iybin==3 && iptbin==2){ rsqmax=18.870 ;}    
      if(iybin==4 && iptbin==2){ rsqmax=11.842 ;}    

      if(iybin==-4 && iptbin==3){ rsqmax=25.471 ;}   
      if(iybin==-3 && iptbin==3){ rsqmax=33.599 ;}   
      if(iybin==-2 && iptbin==3){ rsqmax=41.697 ;}   
      if(iybin==-1 && iptbin==3){ rsqmax=47.842 ;}   
      if(iybin==0 && iptbin==3){ rsqmax=51.307 ;}    
      if(iybin==1 && iptbin==3){ rsqmax=47.482 ;}    
      if(iybin==2 && iptbin==3){ rsqmax=41.697 ;}    
      if(iybin==3 && iptbin==3){ rsqmax=33.599 ;}    
      if(iybin==4 && iptbin==3){ rsqmax=25.471 ;}    

      //      rsqmax*=1.5;

      rsqmax=rsqmax*rsqmax;		 
      
      
      double rsq1=(x[i]-xx0)*(x[i]-xx0)+(y[i]-yy0)*(y[i]-yy0);
      
      if(pt[i]<ptmin)continue; // kill pt=0 participants
      
      double etamaxtpc=1000; // 1.5
      double etamintpc=-1000.; //-1.5
      double rapiditymaxcut=6.0;
      double rapiditymincut=-6.0;
      double ptmincut=0.;//  0.050;    // to simulate what particles we use
	double ptmaxcut=999999.;//    2.0;      // to simulate what particles we use

      double rapfrac=0.75;
      double rapiditymaxcut_2=ybeam*rapfrac;    // for pt distributions
      double rapiditymincut_2=-ybeam*rapfrac;   // for pt distributions
      
      if(pt[i]>ptmincut && pt[i]<ptmaxcut ){
	cth_1_rapidityext_all_h->Fill(rapidity[i],cth1);
	cth_1_etaext_all_h->Fill(eta[i],cth1);
	if(part_id[i]==2212 ||part_id[i]==2112)cth_1_rapidityext_pn_h->Fill(rapidity[i],cth1);
	if(part_id[i]==2212 ||part_id[i]==2112)cth_1_etaext_pn_h->Fill(eta[i],cth1);
	if(part_id[i]==-2212 ||part_id[i]==-2112)cth_1_rapidityext_pnbar_h->Fill(rapidity[i],cth1);
	if(part_id[i]==-2212 ||part_id[i]==-2112)cth_1_etaext_pnbar_h->Fill(eta[i],cth1);
      }
     
      if(eta[i]> etamintpc && eta[i]<etamaxtpc && rapidity[i] > rapiditymincut && rapidity[i] < rapiditymaxcut && pt[i]>ptmincut && pt[i]<ptmaxcut  ){
	
	cth_1_rapidity_all_h->Fill(rapidity[i],cth1);
	if(rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2)cth_1_pt_all_h->Fill(pt[i],cth1*eta[i]/abs(eta[i]));
	if(rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2)cth_1_b_all_h->Fill(b_impact*b_impact,cth1*eta[i]/abs(eta[i]));
	if(part_id[i]==211) cth_1_rapidity_pip_h->Fill(rapidity[i],cth1);
	if(part_id[i]==211 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_1_pt_pip_h->Fill(pt[i],cth1*eta[i]/abs(eta[i]));
	if(part_id[i]==211 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_1_b_pip_h->Fill(b_impact*b_impact,cth1*eta[i]/abs(eta[i]));
	if(part_id[i]==-211) cth_1_rapidity_pim_h->Fill(rapidity[i],cth1);
	if(part_id[i]==-211 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_1_pt_pim_h->Fill(pt[i],cth1*eta[i]/abs(eta[i]));
	if(part_id[i]==-211 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_1_b_pim_h->Fill(b_impact*b_impact,cth1*eta[i]/abs(eta[i]));
	if(part_id[i]==321) cth_1_rapidity_kp_h->Fill(rapidity[i],cth1);
	if(part_id[i]==321 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_1_pt_kp_h->Fill(pt[i],cth1*eta[i]/abs(eta[i]));
	if(part_id[i]==321 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_1_b_kp_h->Fill(b_impact*b_impact,cth1*eta[i]/abs(eta[i]));
	if(part_id[i]==-321) cth_1_rapidity_km_h->Fill(rapidity[i],cth1);
	if(part_id[i]==-321 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_1_pt_km_h->Fill(pt[i],cth1*eta[i]/abs(eta[i]));
	if(part_id[i]==-321 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_1_b_km_h->Fill(b_impact*b_impact,cth1*eta[i]/abs(eta[i]));
	if(part_id[i]==2212) cth_1_rapidity_pr_h->Fill(rapidity[i],cth1);
	if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_1_pt_pr_h->Fill(pt[i],cth1*eta[i]/abs(eta[i]));
	if(part_id[i]==2212) cth_1_b_pr_h->Fill(b_impact*b_impact,cth1*eta[i]/abs(eta[i]));
	if(part_id[i]==6201) cth_1_rapidity_deut_h->Fill(rapidity[i],cth1);
	if(part_id[i]==6201 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_1_pt_deut_h->Fill(pt[i],cth1*eta[i]/abs(eta[i]));
	if(part_id[i]==6201) cth_1_b_deut_h->Fill(b_impact*b_impact,cth1*eta[i]/abs(eta[i]));
	
	cth_2_rapidity_all_h->Fill(rapidity[i],cth2);
	if(rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2)cth_2_pt_all_h->Fill(pt[i],cth2);
	if(rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2)cth_2_b_all_h->Fill(b_impact*b_impact,cth2);
	if(part_id[i]==211) cth_2_rapidity_pip_h->Fill(rapidity[i],cth2);
	if(part_id[i]==211 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_2_pt_pip_h->Fill(pt[i],cth2);
	if(part_id[i]==211 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_2_b_pip_h->Fill(b_impact*b_impact,cth2);
	if(part_id[i]==-211) cth_2_rapidity_pim_h->Fill(rapidity[i],cth2);
	if(part_id[i]==-211 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_2_pt_pim_h->Fill(pt[i],cth2);
	if(part_id[i]==-211 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_2_b_pim_h->Fill(b_impact*b_impact,cth2);
	if(part_id[i]==321) cth_2_rapidity_kp_h->Fill(rapidity[i],cth2);
	if(part_id[i]==321 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_2_pt_kp_h->Fill(pt[i],cth2);
	if(part_id[i]==321 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_2_b_kp_h->Fill(b_impact*b_impact,cth2);
	if(part_id[i]==-321) cth_2_rapidity_km_h->Fill(rapidity[i],cth2);
	if(part_id[i]==-321 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_2_pt_km_h->Fill(pt[i],cth2);
	if(part_id[i]==-321 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_2_b_km_h->Fill(b_impact*b_impact,cth2);
	if(part_id[i]==2212) cth_2_rapidity_pr_h->Fill(rapidity[i],cth2);
	if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_2_pt_pr_h->Fill(pt[i],cth2);
	if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_2_b_pr_h->Fill(b_impact*b_impact,cth2);
	if(part_id[i]==6201) cth_2_rapidity_deut_h->Fill(rapidity[i],cth2);
	if(part_id[i]==6201 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_2_pt_deut_h->Fill(pt[i],cth2);
	if(part_id[i]==6201 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_2_b_deut_h->Fill(b_impact*b_impact,cth2);
	
	cth_3_rapidity_all_h->Fill(rapidity[i],cth3);
	if(rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2)cth_3_pt_all_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	if(rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2)cth_3_b_all_h->Fill(b_impact*b_impact,cth3*eta[i]/abs(eta[i]));
	if(part_id[i]==211) cth_3_rapidity_pip_h->Fill(rapidity[i],cth3);
	if(part_id[i]==211 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_pip_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	if(part_id[i]==211 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_b_pip_h->Fill(b_impact*b_impact,cth3*eta[i]/abs(eta[i]));
	if(part_id[i]==-211) cth_3_rapidity_pim_h->Fill(rapidity[i],cth3);
	if(part_id[i]==-211 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_pim_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	if(part_id[i]==-211 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_b_pim_h->Fill(b_impact*b_impact,cth3*eta[i]/abs(eta[i]));
	if(part_id[i]==321) cth_3_rapidity_kp_h->Fill(rapidity[i],cth3);
	if(part_id[i]==321 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_kp_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	if(part_id[i]==321 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_b_kp_h->Fill(b_impact*b_impact,cth3*eta[i]/abs(eta[i]));
	if(part_id[i]==-321) cth_3_rapidity_km_h->Fill(rapidity[i],cth3);
	if(part_id[i]==-321 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_km_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	if(part_id[i]==-321 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_b_km_h->Fill(b_impact*b_impact,cth3*eta[i]/abs(eta[i]));
	if(part_id[i]==2212) cth_3_rapidity_pr_h->Fill(rapidity[i],cth3);
	if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_pr_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_b_pr_h->Fill(b_impact*b_impact,cth3*eta[i]/abs(eta[i]));
	if(part_id[i]==6201) cth_3_rapidity_deut_h->Fill(rapidity[i],cth3);
	if(part_id[i]==6201 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_deut_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	if(part_id[i]==6201 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_b_deut_h->Fill(b_impact*b_impact,cth3*eta[i]/abs(eta[i]));

//-----------------------------------------------------------
	double centralitybin=0.;
	bool cent=false;
	bool mid=false;
	bool periph=false;
	//the max b^2 is 250 or b=15.8
	// update centrality bins for JAM July 12 2022
	// max b is 13.4  b^2 max=180
	//	if(b_impact>=0. && b_impact<5.){  //0-10%
	if(b_impact>=0. && b_impact<4.24){  //0-10%
	  cent=true;
	  centralitybin=1;
	}
	//	if(b_impact>=5. && b_impact<9.){ //10-40%
	if(b_impact>=4.24 && b_impact<8.5){ //10-40%
	  mid=true;
	  centralitybin=2;
	}
	//	if(b_impact>=9. && b_impact<12.25){ //40-60%
	if(b_impact>=8.5 && b_impact<10.4){ //40-60%
	  periph=true;
	  centralitybin=3;
	}
	
	double ycut1 = 0.5;
	double ptcut1min_pi=0.18;
	double ptcut1max_pi=1.6;
	double ptcut1min_k=0.18;
	double ptcut1max_k=1.6;	
	double ptcut1min_pr=0.4;
	double ptcut1max_pr=2.0;
	double ptcut2min_pr=1.0;
	double ptcut2max_pr=2.5;
	double ptcut1min_deut=0.4*2;
	double ptcut1max_deut=2.0*2;
	double ptcut2min_deut=1.0*2;
	double ptcut2max_deut=2.5*2;

	if(cent){
	  if(pt[i]>ptcut1min_pi && pt[i]<ptcut1max_pi ){
	    cth_3_rapidity_ptcut1_cent_all_h->Fill(rapidity[i],cth3);
	    if(part_id[i]==211) cth_3_rapidity_ptcut1_cent_pip_h->Fill(rapidity[i],cth3);
	    if(part_id[i]==-211) cth_3_rapidity_ptcut1_cent_pim_h->Fill(rapidity[i],cth3);

	  }
	  if(pt[i]>ptcut1min_k && pt[i]<ptcut1max_k ){
	    if(part_id[i]==321) cth_3_rapidity_ptcut1_cent_kp_h->Fill(rapidity[i],cth3);
	    if(part_id[i]==-321) cth_3_rapidity_ptcut1_cent_km_h->Fill(rapidity[i],cth3);
	  }
	  if(pt[i]>ptcut1min_pr && pt[i]<ptcut1max_pr ){
	    if(part_id[i]==2212) cth_3_rapidity_ptcut1_cent_pr_h->Fill(rapidity[i],cth3);
	    if(part_id[i]==2212 && rsq1<rsqmax)eps_x_rapidity_ptcut1_cent_pr_h->Fill(rapidity[i],rsq1*cth3);
	    if(part_id[i]==2212&& rsq1<rsqmax)eps_y_rapidity_ptcut1_cent_pr_h->Fill(rapidity[i],rsq1*sth3);
	    if(part_id[i]==2212&& rsq1<rsqmax)rsq_rapidity_ptcut1_cent_pr_h->Fill(rapidity[i],rsq1);
	  }
	  if(pt[i]>ptcut2min_pr && pt[i]<ptcut2max_pr ){
	    if(part_id[i]==2212) cth_3_rapidity_ptcut2_cent_pr_h->Fill(rapidity[i],cth3);
	    if(part_id[i]==2212&& rsq1<rsqmax)eps_x_rapidity_ptcut2_cent_pr_h->Fill(rapidity[i],rsq1*cth3);
	    if(part_id[i]==2212&& rsq1<rsqmax)eps_y_rapidity_ptcut2_cent_pr_h->Fill(rapidity[i],rsq1*sth3);
	    if(part_id[i]==2212&& rsq1<rsqmax)rsq_rapidity_ptcut2_cent_pr_h->Fill(rapidity[i],rsq1);
	  }

	  if(pt[i]>ptcut1min_deut && pt[i]<ptcut1max_deut ){
	    if(part_id[i]==6201) cth_3_rapidity_ptcut1_cent_deut_h->Fill(rapidity[i],cth3);
	  }
	  if(pt[i]>ptcut2min_deut && pt[i]<ptcut2max_deut ){
	    if(part_id[i]==6201) cth_3_rapidity_ptcut2_cent_deut_h->Fill(rapidity[i],cth3);
	  }
	}
	if(mid){
	  if(pt[i]>ptcut1min_pi && pt[i]<ptcut1max_pi ){
	    cth_3_rapidity_ptcut1_mid_all_h->Fill(rapidity[i],cth3);
	    if(part_id[i]==211) cth_3_rapidity_ptcut1_mid_pip_h->Fill(rapidity[i],cth3);
	    if(part_id[i]==-211) cth_3_rapidity_ptcut1_mid_pim_h->Fill(rapidity[i],cth3);
	  }
	  if(pt[i]>ptcut1min_k && pt[i]<ptcut1max_k ){
	    if(part_id[i]==321) cth_3_rapidity_ptcut1_mid_kp_h->Fill(rapidity[i],cth3);
	    if(part_id[i]==-321) cth_3_rapidity_ptcut1_mid_km_h->Fill(rapidity[i],cth3);
	  }
	  if(pt[i]>ptcut1min_pr && pt[i]<ptcut1max_pr ){
	    if(part_id[i]==2212) cth_3_rapidity_ptcut1_mid_pr_h->Fill(rapidity[i],cth3);
	    if(part_id[i]==2212&& rsq1<rsqmax)eps_x_rapidity_ptcut1_mid_pr_h->Fill(rapidity[i],rsq1*cth3);
	    if(part_id[i]==2212&& rsq1<rsqmax)eps_y_rapidity_ptcut1_mid_pr_h->Fill(rapidity[i],rsq1*sth3);
	    if(part_id[i]==2212&& rsq1<rsqmax)rsq_rapidity_ptcut1_mid_pr_h->Fill(rapidity[i],rsq1); 
	  }
	  if(pt[i]>ptcut2min_pr && pt[i]<ptcut2max_pr ){
	    if(part_id[i]==2212) cth_3_rapidity_ptcut2_mid_pr_h->Fill(rapidity[i],cth3);
	    if(part_id[i]==2212&& rsq1<rsqmax)eps_x_rapidity_ptcut2_mid_pr_h->Fill(rapidity[i],rsq1*cth3);
	    if(part_id[i]==2212&& rsq1<rsqmax)eps_y_rapidity_ptcut2_mid_pr_h->Fill(rapidity[i],rsq1*sth3);
	    if(part_id[i]==2212&& rsq1<rsqmax)rsq_rapidity_ptcut2_mid_pr_h->Fill(rapidity[i],rsq1); 
	  }
	  if(pt[i]>ptcut1min_deut && pt[i]<ptcut1max_deut ){
	    if(part_id[i]==6201) cth_3_rapidity_ptcut1_mid_deut_h->Fill(rapidity[i],cth3);
	  }
	  if(pt[i]>ptcut2min_deut && pt[i]<ptcut2max_deut ){
	    if(part_id[i]==6201) cth_3_rapidity_ptcut2_mid_deut_h->Fill(rapidity[i],cth3);
	  }
	}
	if(periph){
	  if(pt[i]>ptcut1min_pi && pt[i]<ptcut1max_pi ){
	    cth_3_rapidity_ptcut1_periph_all_h->Fill(rapidity[i],cth3);
	    if(part_id[i]==211) cth_3_rapidity_ptcut1_periph_pip_h->Fill(rapidity[i],cth3);
	    if(part_id[i]==-211) cth_3_rapidity_ptcut1_periph_pim_h->Fill(rapidity[i],cth3);
	  }
	  if(pt[i]>ptcut1min_k && pt[i]<ptcut1max_k ){
	    if(part_id[i]==321) cth_3_rapidity_ptcut1_periph_kp_h->Fill(rapidity[i],cth3);
	    if(part_id[i]==-321) cth_3_rapidity_ptcut1_periph_km_h->Fill(rapidity[i],cth3);
	  }
	  if(pt[i]>ptcut1min_pr && pt[i]<ptcut1max_pr ){
	    if(part_id[i]==2212) cth_3_rapidity_ptcut1_periph_pr_h->Fill(rapidity[i],cth3);
	    if(part_id[i]==2212&& rsq1<rsqmax)eps_x_rapidity_ptcut1_periph_pr_h->Fill(rapidity[i],rsq1*cth3);
	    if(part_id[i]==2212&& rsq1<rsqmax)eps_y_rapidity_ptcut1_periph_pr_h->Fill(rapidity[i],rsq1*sth3);
	    if(part_id[i]==2212&& rsq1<rsqmax)rsq_rapidity_ptcut1_periph_pr_h->Fill(rapidity[i],rsq1);
	  }
	  if(pt[i]>ptcut2min_pr && pt[i]<ptcut2max_pr ){
	    if(part_id[i]==2212) cth_3_rapidity_ptcut2_periph_pr_h->Fill(rapidity[i],cth3);
	    if(part_id[i]==2212&& rsq1<rsqmax)eps_x_rapidity_ptcut2_periph_pr_h->Fill(rapidity[i],rsq1*cth3);
	    if(part_id[i]==2212&& rsq1<rsqmax)eps_y_rapidity_ptcut2_periph_pr_h->Fill(rapidity[i],rsq1*sth3);
	    if(part_id[i]==2212&& rsq1<rsqmax)rsq_rapidity_ptcut2_periph_pr_h->Fill(rapidity[i],rsq1);
	  }
	  if(pt[i]>ptcut1min_deut && pt[i]<ptcut1max_deut ){
	    if(part_id[i]==6201) cth_3_rapidity_ptcut1_periph_deut_h->Fill(rapidity[i],cth3);
	  }
	  if(pt[i]>ptcut2min_deut && pt[i]<ptcut2max_deut ){
	    if(part_id[i]==6201) cth_3_rapidity_ptcut2_periph_deut_h->Fill(rapidity[i],cth3);
	  }
	}

	if(abs(rapidity[i])<ycut1){
	  if(cent){
	    if(rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2)cth_3_pt_ycut1_cent_all_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==211 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_ycut1_cent_pip_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==-211 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_ycut1_cent_pim_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==321 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_ycut1_cent_kp_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==-321 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_ycut1_cent_km_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_ycut1_cent_pr_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==6201 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_ycut1_cent_deut_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2&& rsq1<rsqmax)eps_x_pt_ycut1_cent_pr_h->Fill(pt[i],rsq1*cth3*eta[i]/abs(eta[i])); 
	    if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2&& rsq1<rsqmax)  eps_y_pt_ycut1_cent_pr_h->Fill(pt[i],rsq1*sth3);
	    if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2&& rsq1<rsqmax)  rsq_pt_ycut1_cent_pr_h->Fill(pt[i],rsq1); 

	  }	  
	  if(mid){
	    if(rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2)cth_3_pt_ycut1_mid_all_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==211 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_ycut1_mid_pip_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==-211 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_ycut1_mid_pim_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==321 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_ycut1_mid_kp_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==-321 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_ycut1_mid_km_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_ycut1_mid_pr_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==6201 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_ycut1_mid_deut_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2&& rsq1<rsqmax)eps_x_pt_ycut1_mid_pr_h->Fill(pt[i],rsq1*cth3*eta[i]/abs(eta[i]));  
	    if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2&& rsq1<rsqmax)  eps_y_pt_ycut1_mid_pr_h->Fill(pt[i],rsq1*sth3);
	    if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2&& rsq1<rsqmax)  rsq_pt_ycut1_mid_pr_h->Fill(pt[i],rsq1);   

	  }	  
	  if(periph){
	    if(rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2)cth_3_pt_ycut1_periph_all_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==211 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_ycut1_periph_pip_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==-211 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_ycut1_periph_pim_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==321 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_ycut1_periph_kp_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==-321 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_ycut1_periph_km_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_ycut1_periph_pr_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==6201 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_pt_ycut1_periph_deut_h->Fill(pt[i],cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2&& rsq1<rsqmax)eps_x_pt_ycut1_periph_pr_h->Fill(pt[i],rsq1*cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2&& rsq1<rsqmax)  eps_y_pt_ycut1_periph_pr_h->Fill(pt[i],rsq1*sth3);
	    if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2&& rsq1<rsqmax)  rsq_pt_ycut1_periph_pr_h->Fill(pt[i],rsq1);

	  }
	}
	
	if(abs(rapidity[i])<ycut1){
	  if(pt[i]>ptcut1min_pi && pt[i]<ptcut1max_pi ){
	    if(rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2)cth_3_b_ycut1_ptcut1_all_h->Fill(b_impact*b_impact,cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==211 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_b_ycut1_ptcut1_pip_h->Fill(b_impact*b_impact,cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==-211 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_b_ycut1_ptcut1_pim_h->Fill(b_impact*b_impact,cth3*eta[i]/abs(eta[i]));
	  }
	  if(pt[i]>ptcut1min_k && pt[i]<ptcut1max_k ){
	    if(part_id[i]==321 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_b_ycut1_ptcut1_kp_h->Fill(b_impact*b_impact,cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==-321 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_b_ycut1_ptcut1_km_h->Fill(b_impact*b_impact,cth3*eta[i]/abs(eta[i]));
	  }
	  if(pt[i]>ptcut1min_pr && pt[i]<ptcut1max_pr ){
	    if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_b_ycut1_ptcut1_pr_h->Fill(b_impact*b_impact,cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==6201 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_b_ycut1_ptcut1_deut_h->Fill(b_impact*b_impact,cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2&& rsq1<rsqmax)eps_x_b_ycut1_ptcut1_pr_h->Fill(b_impact*b_impact,rsq1*cth3*eta[i]/abs(eta[i])); 
	    if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2&& rsq1<rsqmax)  eps_y_b_ycut1_ptcut1_pr_h->Fill(b_impact*b_impact,rsq1*sth3);
	    if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2&& rsq1<rsqmax)  rsq_b_ycut1_ptcut1_pr_h->Fill(b_impact*b_impact,rsq1);
	  }
	  if(pt[i]>ptcut2min_pr && pt[i]<ptcut2max_pr ){
	    if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_b_ycut1_ptcut2_pr_h->Fill(b_impact*b_impact,cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==6201 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2) cth_3_b_ycut1_ptcut2_deut_h->Fill(b_impact*b_impact,cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2&& rsq1<rsqmax)eps_x_b_ycut1_ptcut2_pr_h->Fill(b_impact*b_impact,rsq1*cth3*eta[i]/abs(eta[i]));
	    if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2&& rsq1<rsqmax)  eps_y_b_ycut1_ptcut2_pr_h->Fill(b_impact*b_impact,rsq1*sth3); 
	    if(part_id[i]==2212 && rapidity[i] > rapiditymincut_2 && rapidity[i] < rapiditymaxcut_2&& rsq1<rsqmax)  rsq_b_ycut1_ptcut2_pr_h->Fill(b_impact*b_impact,rsq1); 
	  }
	}
	

 //-----------------------------------------------------------
	
	cth_1_eta_all_h->Fill(eta[i],cth1);
	if(part_id[i]==211) cth_1_eta_pip_h->Fill(eta[i],cth1);
	if(part_id[i]==-211) cth_1_eta_pim_h->Fill(eta[i],cth1);
	if(part_id[i]==321) cth_1_eta_kp_h->Fill(eta[i],cth1);
	if(part_id[i]==-321) cth_1_eta_km_h->Fill(eta[i],cth1);
	if(part_id[i]==2212) cth_1_eta_pr_h->Fill(eta[i],cth1);
	if(part_id[i]==6201) cth_1_eta_deut_h->Fill(eta[i],cth1);
	
	cth_2_eta_all_h->Fill(eta[i],cth2);
	if(part_id[i]==211) cth_2_eta_pip_h->Fill(eta[i],cth2);
	if(part_id[i]==-211) cth_2_eta_pim_h->Fill(eta[i],cth2);
	if(part_id[i]==321) cth_2_eta_kp_h->Fill(eta[i],cth2);
	if(part_id[i]==-321) cth_2_eta_km_h->Fill(eta[i],cth2);
	if(part_id[i]==2212) cth_2_eta_pr_h->Fill(eta[i],cth2);
	if(part_id[i]==6201) cth_2_eta_deut_h->Fill(eta[i],cth2);

	cth_3_eta_all_h->Fill(eta[i],cth3);
	if(part_id[i]==211) cth_3_eta_pip_h->Fill(eta[i],cth3);
	if(part_id[i]==-211) cth_3_eta_pim_h->Fill(eta[i],cth3);
	if(part_id[i]==321) cth_3_eta_kp_h->Fill(eta[i],cth3);
	if(part_id[i]==-321) cth_3_eta_km_h->Fill(eta[i],cth3);
	if(part_id[i]==2212) cth_3_eta_pr_h->Fill(eta[i],cth3);
	if(part_id[i]==6201) cth_3_eta_deut_h->Fill(eta[i],cth3);
	
      }
      
    } // end particle loop
    
    
//----------------------------------------- finding eccentricities ---------------------

// right now we are NOT moving the x and y coord system to the center of grav
    double eps[6]={0};      // eccentricity
    double eps_x[6]={0};      
    double eps_y[6]={0};
    double epsposrap[6]={0};      // eccentricity
    double eps_xposrap[6]={0};      
    double eps_yposrap[6]={0};
    double epsnegrap[6]={0};      // eccentricity
    double eps_xnegrap[6]={0};      
    double eps_ynegrap[6]={0};
    double  PSI_eps[6]={0};   // reaction plane from x and y coordidates


    double x0=0;
    double y0=0;
    double x0posrap=0;
    double y0posrap=0;
    double x0negrap=0;
    double y0negrap=0;
    int numusedparticles=0;	
    int numusedparticlesposrap=0;	
    int numusedparticlesnegrap=0;	
    for (int i=0; i< number_particles; i++){
      // These cuts must depend on energy?
      if(abs(eta[i])>5.)continue;
     if(pt[i]<ptmin)continue; // only plot participants
      if(t[i]>tmax)continue; // AMPT only follows stuff 400*0.2 fm =80 fm
      if(!(part_id[i]==2112||part_id[i]==2212|| part_id[i]>useonlynucleons))continue; // only look at p and n
      //      if(abs(eta[i])<etamin)continue;  // use eta=[1.5,5] for RP (RPD)
      numusedparticles++;
      x0+=x[i];
      y0+=y[i];
      if(eta[i]>0.){    // used etamin for AMPT
	numusedparticlesposrap++;
	x0posrap+=x[i];
	y0posrap+=y[i];
      }
      if(eta[i]<-0.){   // used etamin for AMPT
	numusedparticlesnegrap++;
	x0negrap+=x[i];
	y0negrap+=y[i];
      }
    } // end particle loop  i 
    x0/=double(numusedparticles);
    y0/=double(numusedparticles);
    x0posrap/=double(numusedparticlesposrap);    // check for zero
    y0posrap/=double(numusedparticlesposrap);    // check for zero
    x0negrap/=double(numusedparticlesnegrap);   // check for zero
    y0negrap/=double(numusedparticlesnegrap);   // check for zero
    //    std::cout<<" x y="<<x0<<" "<<y0<<" "<<x0posrap<<" "<<y0posrap<<" "<<x0negrap<<" "<<y0negrap<<std::endl;
  x0y0h->Fill(x0,y0);  
  x0y0posraph->Fill(x0posrap,y0posrap);
  x0y0negraph->Fill(x0negrap,y0negrap);

  // set x0 and y0 - I don't think you recenter RKS Jan 26, 2022
  x0=0.;
  y0=0.;
  //  x0posrap=4.;
  x0posrap=0.;
  y0posrap=0.;
  //  x0negrap=-4.;
  x0negrap=0.;
  y0negrap=0.;


    double phixy=0;
    double rsq=0;	
    double rsqav=0;	

    double phixyposrap=0;
    double rsqposrap=0;	
    double rsqavposrap=0;	

    double phixynegrap=0;
    double rsqnegrap=0;	
    double rsqavnegrap=0;	

    //    int numusedparticles=0;	
    //    int numusedparticlesposrap=0;	
    //    int numusedparticlesnegrap=0;	
    numusedparticles=0;	
    numusedparticlesposrap=0;	
    numusedparticlesnegrap=0;	
    for (int i=0; i< number_particles; i++){
      // These cuts must depend on energy?
      if(abs(eta[i])>5.)continue;
      if(pt[i]<ptmin)continue; // only plot participants
      if(t[i]>tmax)continue; // AMPT only follows stuff 400*0.2 fm =80 fm
      if(!(part_id[i]==2112||part_id[i]==2212|| part_id[i]>useonlynucleons))continue; // only look at p and n
      //      if(abs(eta[i])<etamin)continue;  // use eta=[1.5,5] for RP (RPD)
      phixy=atan2(y[i],x[i]);
      rsq=x[i]*x[i]+y[i]*y[i];
      rsqav+=rsq;
      numusedparticles++;
      if(eta[i]>etamin){
	phixyposrap=atan2(y[i]-y0posrap,x[i]-x0posrap);
	rsqposrap=(x[i]-x0posrap)*(x[i]-x0posrap)+(y[i]-y0posrap)*(y[i]-y0posrap);
	rsqavposrap+=rsqposrap;
	numusedparticlesposrap++;
      }
      if(eta[i]<-etamin){
	phixynegrap=atan2(y[i]-y0negrap,x[i]-x0negrap);
	rsqnegrap=(x[i]-x0negrap)*(x[i]-x0negrap)+(y[i]-y0negrap)*(y[i]-y0negrap);
	rsqavnegrap+=rsqnegrap;
	numusedparticlesnegrap++;
      }
      
      phixyh->Fill(phixy); 
      phixyposraph->Fill(phixyposrap);
      phixynegraph->Fill(phixynegrap);

      rsqh->Fill(rsq); 
      rsqposraph->Fill(rsqposrap); 
      rsqnegraph->Fill(rsqnegrap); 

      for(int nn=0; nn<6 ; nn++){
	int nharmonic=nn+1;
	eps_x[nn]+=rsq*cos(nharmonic*phixy);
	eps_y[nn]+=rsq*sin(nharmonic*phixy);    
	if(eta[i]>etamin){
	  eps_xposrap[nn]+=rsqposrap*cos(nharmonic*phixyposrap);
	  eps_yposrap[nn]+=rsqposrap*sin(nharmonic*phixyposrap);
	}
	if(eta[i]<-etamin){
	  eps_xnegrap[nn]+=rsqnegrap*cos(nharmonic*phixynegrap);
	  eps_ynegrap[nn]+=rsqnegrap*sin(nharmonic*phixynegrap);	
	}
      } //nn
    } // end particle loop  i 

    rsqav/=double(numusedparticles);
    rsqavposrap/=double(numusedparticlesposrap);    // check for zero
    rsqavnegrap/=double(numusedparticlesnegrap);   // check for zero

    rsqavh->Fill(rsqav);
    rsqavposraph->Fill(rsqavposrap); 
    rsqavnegraph->Fill(rsqavnegrap); 

    for(int nn=0; nn<6 ; nn++){
      int nharmonic=nn+1;
      PSI_eps[nn]=atan2(eps_y[nn],eps_x[nn])/nharmonic;
      
      eps_x[nn]/=double(numusedparticles); // check for zero
      eps_y[nn]/=double(numusedparticles); // check for zero
      eps[nn]=sqrt(eps_x[nn]*eps_x[nn]+eps_y[nn]*eps_y[nn])/rsqav;

      eps_xposrap[nn]/=double(numusedparticlesposrap); // check for zero
      eps_yposrap[nn]/=double(numusedparticlesposrap); // check for zero
      epsposrap[nn]=sqrt(eps_xposrap[nn]*eps_xposrap[nn]+eps_yposrap[nn]*eps_yposrap[nn])/rsqavposrap;
     eps_xnegrap[nn]/=double(numusedparticlesnegrap); // check for zero
      eps_ynegrap[nn]/=double(numusedparticlesnegrap); // check for zero
      epsnegrap[nn]=sqrt(eps_xnegrap[nn]*eps_xnegrap[nn]+eps_ynegrap[nn]*eps_ynegrap[nn])/rsqavnegrap;

    }
    //    std::cout<<" eps 1 2 3 ="<< eps[0]<<" "<<eps[1]<<" "<<eps[2]<<" "<<std::endl;
    eps1h->Fill(eps[1-1]);
    eps2h->Fill(eps[2-1]);
    eps3h->Fill(eps[3-1]);
    eps4h->Fill(eps[4-1]);
    eps5h->Fill(eps[5-1]);
    eps6h->Fill(eps[6-1]);

  eps1posraph->Fill(epsposrap[1-1]); 
  eps2posraph->Fill(epsposrap[2-1]); 
  eps3posraph->Fill(epsposrap[3-1]); 
  eps4posraph->Fill(epsposrap[4-1]); 
  eps5posraph->Fill(epsposrap[5-1]); 
  eps6posraph->Fill(epsposrap[6-1]); 

  eps1negraph->Fill(epsnegrap[1-1]); 
  eps2negraph->Fill(epsnegrap[2-1]); 
  eps3negraph->Fill(epsnegrap[3-1]); 
  eps4negraph->Fill(epsnegrap[4-1]); 
  eps5negraph->Fill(epsnegrap[5-1]); 
  eps6negraph->Fill(epsnegrap[6-1]); 
 
    /*
    double PSI_eps1posrap=atan2(eps_yposrap[0],eps_xposrap[0])/1;
    double PSI_eps1negrap=atan2(eps_ynegrap[0],eps_xnegrap[0])/1;
    double PSI_eps2posrap=atan2(eps_yposrap[1],eps_xposrap[1])/2;
    double PSI_eps2negrap=atan2(eps_ynegrap[1],eps_xnegrap[1])/2;
    double PSI_eps3posrap=atan2(eps_yposrap[2],eps_xposrap[2])/3;
    double PSI_eps3negrap=atan2(eps_ynegrap[2],eps_xnegrap[2])/3;
    */    
/*
    PSI_eps1h->Fill(PSI_eps[1-1]);
    PSI_eps2h->Fill(PSI_eps[2-1]);
    PSI_eps3h->Fill(PSI_eps[3-1]);
    PSI_eps4h->Fill(PSI_eps[4-1]);

    PSI_eps1_roth->Fill(PSI_eps[1-1],rot);
    PSI_eps2_roth->Fill(PSI_eps[2-1],rot);
    PSI_eps3_roth->Fill(PSI_eps[3-1],rot);
    PSI_eps4_roth->Fill(PSI_eps[4-1],rot);

    PSI_eps1posraph->Fill(PSI_eps1posrap);
    PSI_eps1negraph->Fill(PSI_eps1negrap);
    PSI_eps2posraph->Fill(PSI_eps2posrap);
    PSI_eps2negraph->Fill(PSI_eps2negrap);
    PSI_eps3posraph->Fill(PSI_eps3posrap);
    PSI_eps3negraph->Fill(PSI_eps3negrap);
*/

//----------------------------------------- end finding eccentricities ---------------------

    //Finding reaction plane
    double Qx[6]={0};      // used to find the reaction plane
    double Qy[6]={0};
    double Qxposrap[6]={0};      // used to find the reaction plane
    double Qyposrap[6]={0};
    double Qxnegrap[6]={0};      // used to find the reaction plane
    double Qynegrap[6]={0};
  
    double Weight=1;
    for (int i=0; i< number_particles; i++){
      //      if(t[i]>tmax)continue; // AMPT only follows stuff 400*0.2 fm =80 fm - leave out cut - allow decays
      // These cuts must depend on energy?
      if(abs(eta[i])>5)continue;
      //      if(abs(eta[i])<1.5)continue;  // use eta=[1.5,5] for RP (RPD)
      for(int nn=0; nn<6 ; nn++){
	int nharmonic=nn+1;
	Weight=1;
	if(nharmonic==1 || nharmonic==3){
	  Weight=eta[i];
	  //	  Weight=1;
	}
	Qx[nn]+=Weight*cos(nharmonic*phi[i]);
	Qy[nn]+=Weight*sin(nharmonic*phi[i]);
	Weight=1;
	if(eta[i]>0.)Qxposrap[nn]+=Weight*cos(nharmonic*phi[i]);
	if(eta[i]>0.)Qyposrap[nn]+=Weight*sin(nharmonic*phi[i]);
	if(eta[i]<0.)Qxnegrap[nn]+=Weight*cos(nharmonic*phi[i]);
	if(eta[i]<0.)Qynegrap[nn]+=Weight*sin(nharmonic*phi[i]);

      }
    } // end particle loop
    for(int nn=0; nn<6 ; nn++){
      int nharmonic=nn+1;
      PSI[nn]=atan2(Qy[nn],Qx[nn])/nharmonic;
    }
    double PSI1posrap=atan2(Qyposrap[0],Qxposrap[0])/1;
    double PSI1negrap=atan2(Qynegrap[0],Qxnegrap[0])/1;
    double PSI2posrap=atan2(Qyposrap[1],Qxposrap[1])/2;
    double PSI2negrap=atan2(Qynegrap[1],Qxnegrap[1])/2;
    double PSI3posrap=atan2(Qyposrap[2],Qxposrap[2])/3;
    double PSI3negrap=atan2(Qynegrap[2],Qxnegrap[2])/3;
    
    PSI1h->Fill(PSI[1-1]);
    PSI2h->Fill(PSI[2-1]);
    PSI3h->Fill(PSI[3-1]);
    PSI4h->Fill(PSI[4-1]);

    PSI1_roth->Fill(PSI[1-1],rot);
    PSI2_roth->Fill(PSI[2-1],rot);
    PSI3_roth->Fill(PSI[3-1],rot);
    PSI4_roth->Fill(PSI[4-1],rot);

    PSI1posraph->Fill(PSI1posrap);
    PSI1negraph->Fill(PSI1negrap);
    PSI2posraph->Fill(PSI2posrap);
    PSI2negraph->Fill(PSI2negrap);
    PSI3posraph->Fill(PSI3posrap);
    PSI3negraph->Fill(PSI3negrap);

  } // end event loop 2

  std::cout<<" bmin frac="<<float(nevbmin)/float(nev)<<" bmax frac="<<float(nevbmax)/float(nev)<<std::endl;
    
  //===============================3==============================
  // event loop 3 Find vn's using reaction plane
  double vsumtot_a[6][6]={0};   // sum over all events, all particles [harmonic,RP]
  int nvsumtot_a=0; // count over all events, all partices
  for(int iev=0; iev<nev; iev++){
    if(iev%ncount==0)std::cout<<" 3 eventnum="<<iev<<std::endl;
    if(ntodo>0 && iev > ntodo)break;

    eventnum=eventnum_v[iev];
    number_particles=number_particles_v[iev];
    refmultfxt=refmultfxt_v[iev];
    b_impact=b_impact_v[iev];
    npart_proj=npart_proj_v[iev];
    npart_targ=npart_targ_v[iev];
    nelas_proj=nelas_proj_v[iev];
    ninelas_proj=ninelas_proj_v[iev];
    nelas_targ=nelas_targ_v[iev];
    ninelas_targ=ninelas_targ_v[iev];
    part_id=part_id_v[iev];
    px=px_v[iev];
    py=py_v[iev];
    pz=pz_v[iev];
    mass=mass_v[iev];
    x=x_v[iev];
    y=y_v[iev];
    z=z_v[iev];
    t=t_v[iev];
    energy=energy_v[iev];
    pt=pt_v[iev];
    rapidity=rapidity_v[iev];
    phi=phi_v[iev];
    eta=eta_v[iev];

    if(b_impact<bmin || b_impact>bmax)continue;
    
    double vsum_a[6][6]={0};  // sum over all particles [harmonic,RP]
    int nvsum_a=0; // count over all partices  

    for (int i=0; i< number_particles; i++){
      if(abs(eta[i])>5)continue;
      if(abs(eta[i])>1.)continue; // only calculate v's for eta=[-1,1]
     
      for(int nn=0; nn<6 ; nn++){
	int nharmonic=nn+1;	
	for(int mm=0; mm<6 ; mm++){
	  int mRP=mm+1;
	  if(nharmonic%mRP != 0)continue;  // for PSI2 can only do nharmonic 2,4,6 etc
	  vsum_a[nn][mm]+=cos(nharmonic*(phi[i]-PSI[mm])); 
	  vsumtot_a[nn][mm]+=cos(nharmonic*(phi[i]-PSI[mm])); 
	  double costerm=cos(nharmonic*(phi[i]-PSI[mm]));
	  //	  std::cout<<" cos term="<<costerm<<std::endl;
	  if(nharmonic==2 & mRP==2)costermh->Fill(costerm);
	}
      }
      nvsum_a++;
      nvsumtot_a++;
    } // end particle loop
    for(int nn=0; nn<6 ; nn++){
      int nharmonic=nn+1;
      for(int mm=0; mm<6 ; mm++){
	int mRP=mm+1;
	vsum_a[nn][mm]/=nvsum_a;
      }
    }
    v1_RP1h->Fill(vsum_a[1-1][1-1]);
    v2_RP1h->Fill(vsum_a[2-1][1-1]);
    v3_RP1h->Fill(vsum_a[3-1][1-1]);
    v4_RP1h->Fill(vsum_a[4-1][1-1]);
    v5_RP1h->Fill(vsum_a[5-1][1-1]);
    v6_RP1h->Fill(vsum_a[6-1][1-1]);
    v2_RP2h->Fill(vsum_a[2-1][2-1]);
    v4_RP2h->Fill(vsum_a[4-1][2-1]);
    v6_RP2h->Fill(vsum_a[6-1][2-1]);
    v3_RP3h->Fill(vsum_a[3-1][3-1]);
    v6_RP3h->Fill(vsum_a[6-1][3-1]);
    
  } // end event loop 3
  
  for(int nn=0; nn<6 ; nn++){
    int nharmonic=nn+1;
    for(int mm=0; mm<6 ; mm++){
      int mRP=mm+1;
      vsumtot_a[nn][mm]/=nvsumtot_a;
    }
  }
  if(i_inject>0)std::cout<<std::endl<<" injected flow  v1="<<v_injected[0]<<" v2="<<v_injected[1]<<" v3="<<v_injected[2]<<" v4="<<v_injected[3]<<" v5="<<v_injected[4]<<" v6="<<v_injected[5]<<std::endl<<std::endl;
  std::cout<<" v1_rp1="<<vsumtot_a[1-1][1-1]<<std::endl;
  std::cout<<" v2_RP1="<<vsumtot_a[2-1][1-1]<<std::endl;
  std::cout<<" v3_RP1="<<vsumtot_a[3-1][1-1]<<std::endl;
  std::cout<<" v4_RP1="<<vsumtot_a[4-1][1-1]<<std::endl;
  std::cout<<" v5_RP1="<<vsumtot_a[5-1][1-1]<<std::endl;
  std::cout<<" v6_RP1="<<vsumtot_a[6-1][1-1]<<std::endl;
  std::cout<<" v2_RP2="<<vsumtot_a[2-1][2-1]<<std::endl;
  std::cout<<" v4_RP2="<<vsumtot_a[4-1][2-1]<<std::endl;
  std::cout<<" v6_RP2="<<vsumtot_a[6-1][2-1]<<std::endl;
  std::cout<<" v3_RP3="<<vsumtot_a[3-1][3-1]<<std::endl;
  std::cout<<" v6_RP3="<<vsumtot_a[6-1][3-1]<<std::endl;
  

  //===============================4==============================
  // event loop 4 Find vn's using cumulants
  double bb2bb[6]={0};
  double W2sum=0;
  double c_2[6]={0};  // _2 refers to cumulant order [nn] refers to v1,v2,v3,..
  double v_2[6]={0};
  double bb4bb[6]={0};
  double W4sum=0;
  double c_4[6]={0};  // _4 refers to cumulant order [nn] refers to v1,v2,v3,..
  double v_4[6]={0};
  for(int iev=0; iev<nev; iev++){
    if(iev%ncount==0)std::cout<<" 4 eventnum="<<iev<<std::endl;
    if(ntodo>0 && iev > ntodo)break;
    
    eventnum=eventnum_v[iev];
    number_particles=number_particles_v[iev];
    refmultfxt=refmultfxt_v[iev];
    b_impact=b_impact_v[iev];
    npart_proj=npart_proj_v[iev];
    npart_targ=npart_targ_v[iev];
    nelas_proj=nelas_proj_v[iev];
    ninelas_proj=ninelas_proj_v[iev];
    nelas_targ=nelas_targ_v[iev];
    ninelas_targ=ninelas_targ_v[iev];
    part_id=part_id_v[iev];
    px=px_v[iev];
    py=py_v[iev];
    pz=pz_v[iev];
    mass=mass_v[iev];
    x=x_v[iev];
    y=y_v[iev];
    z=z_v[iev];
    t=t_v[iev];
    energy=energy_v[iev];
    pt=pt_v[iev];
    rapidity=rapidity_v[iev];
    phi=phi_v[iev];
    eta=eta_v[iev];
    
    if(b_impact<bmin || b_impact>bmax)continue; 
    
    TComplex Qnharmonic[12]={0};       

    for (int i=0; i< number_particles; i++){
      //      if(t[i]>tmax)continue; // AMPT only follows stuff 400*0.2 fm =80 fm - leavve out cut allow decays
      if(abs(eta[i])>5)continue;
      // calculate Qvectors
      for(int j=0; j<12; j++){
	Qnharmonic[j]+=TComplex(cos((j+1)*phi[i]),sin((j+1)*phi[i]),false); 
      }
    } // end particle loop

    //    std::cout<<" number_particles="<<number_particles<<std::endl;
    // now figure out up to nn=6
    // 2nd order cumulant
    double q2[6]={0};
    double b2b[6]={0};
    double W2=number_particles*(number_particles-1);
    W2sum+=W2;  
    for(int nn=0; nn<6; nn++){
      int nharmonic=nn+1;
      q2[nn]=Qnharmonic[nn].Rho2();// Q_nharmonic^2
      b2b[nn]=(q2[nn]-number_particles)/number_particles/(number_particles-1);
      //      std::cout<<"1 nn="<<nn<<" q2="<<q2[nn]<<" b2b="<<b2b[nn]<<" W2="<<W2<<std::endl;
      bb2bb[nn]+=W2*b2b[nn];
    }

    // 4th order cumulant
    double q4[6]={0};
    double q2n2[6]={0};
    double q2nnstarnstar[6]={0};    
    double b4b[6]={0};
    //   double W4=number_particles*(number_particles-1)*(number_particles-2)*(number_particles-3);
    //   std::cout<<" 1 W4="<<W4;     // VERY WEIRD PROBLEM 1 W4 and 2 W4 are not equal
    double W4=0;
    W4=number_particles*(number_particles-1); //TEST
    W4*=(number_particles-2)*(number_particles-3); // TEST
    //    std::cout<<" 2 W4="<<W4<<std::endl;
    W4sum+=W4;  
    for(int nn=0; nn<6; nn++){
      int nharmonic=nn+1;
      q4[nn]=q2[nn]*q2[nn];// Q_nharmonic^4
      q2n2[nn]=Qnharmonic[2*nharmonic-1].Rho2();// Q_2*nharmonic^2
      TComplex tmpc=Qnharmonic[2*nharmonic-1]*TComplex::Conjugate(Qnharmonic[nharmonic-1])*TComplex::Conjugate(Qnharmonic[nharmonic-1]);
      q2nnstarnstar[nn]=tmpc.Re();
      b4b[nn] = q4[nn]+q2n2[nn]-2*q2nnstarnstar[nn]-2*2*(number_particles-2)*q2[nn]+2*number_particles*(number_particles-3);
      b4b[nn]=b4b[nn]/(number_particles*(number_particles-1)*(number_particles-2)*(number_particles-3));
      bb4bb[nn]+=W4*b4b[nn];
    }
  } // end event loop 4
  std::cout<<" W2sum="<<W2sum<<std::endl;
  for(int nn=0; nn<6; nn++){
    int nharmonic=nn+1;
    std::cout<<"2 nn="<<nn<<" bb2bb="<<bb2bb[nn]<<std::endl;
    bb2bb[nn]=bb2bb[nn]/W2sum;
    c_2[nn]=bb2bb[nn];
    v_2[nn]=sqrt(c_2[nn]);
    std::cout<<" using 2nd order cumulant nharmonic="<<nn+1<<" c_nharmonic="<<c_2[nn]<<" v_nharmonic="<<v_2[nn]<<std::endl;		 
  }   
  std::cout<<" W4sum="<<W4sum<<std::endl;
  for(int nn=0; nn<6; nn++){
    int nharmonic=nn+1;
    std::cout<<"4 nn="<<nn<<" bb4bb="<<bb4bb[nn]<<std::endl;
    bb4bb[nn]=bb4bb[nn]/W4sum;
    c_4[nn]=bb4bb[nn]-2*bb2bb[nn]*bb2bb[nn];
    v_4[nn]=sqrt(sqrt(-c_4[nn]));
    std::cout<<" using 4th order cumulant nharmonic="<<nn+1<<" c_nharmonic="<<c_4[nn]<<" v_nharmonic="<<v_4[nn]<<std::endl;		 
  }   

  std::cout<<" ybeam="<<ybeam<<std::endl;
  std::cout<<" bmin centrality="<<float(nevbmin)/float(nev)*100<<"%  bmax frac="<<float(nevbmax)/float(nev)*100<<"%"<<std::endl;
  
  //===============================================================
  
  std::cout<<" finished analysis, now write out histograms"<<std::endl;
  // Open a ROOT file and save the histograms
  //
  TFile *myfile = new TFile("jam1ana5.root","RECREATE");
  infoh->Write();
  number_particlesh->Write();
  b_impactallh->Write();
  b_impacth->Write();
  bsq_impacth->Write();
  n_participantsh->Write();
  part_idh->Write();
  part_id2h->Write();

  fxtmulth->Write();        
  number_particles_bsqh->Write();
  fxtmult_bsqh->Write();    
  nparticipants_bsqh->Write(); 
    
  //********************
  rapidityptrawh->Write();
  rapiditypth->Write();
  rapiditypt_piph->Write();
  rapiditypt_pimh->Write();
  rapiditypt_kph->Write(); 
  rapiditypt_kmh->Write(); 
  rapiditypt_prh->Write(); 
  
  xh->Write();
  yh->Write();
  zh->Write();
  th->Write();

  ncollparth->Write();
  formtimeh->Write();
  lastcolltimeh->Write();
  
  ptpzh->Write();    
  etarapidityh->Write();     
  xvsz_allnh->Write();
  xvsrapidity_allnucleonsh->Write();
  xvsy_allnh->Write();     
  xvst_allnh->Write();     
  ptvst_allnh->Write(); 
  etavst_allnh->Write(); 
  rapidityvst_allnh->Write(); 

  pxvst_allnh->Write();
  pyvst_allnh->Write();
  pzvst_allnh->Write();     
  pxvstposrap_allnh->Write();
  pyvstposrap_allnh->Write();
  pzvstposrap_allnh->Write();     
  pxvstybeam_allnh->Write();
  pyvstybeam_allnh->Write();
  pzvstybeam_allnh->Write();     
  
  cth1vst_allnh->Write();
  cth1vstposrap_allnh->Write();
  cth1vstybeam_allnh->Write();

  cth2vst_allnh->Write();
  cth2vstposrap_allnh->Write();
  cth2vstybeam_allnh->Write();

  cth3vst_allnh->Write();
  cth3vstposrap_allnh->Write();
  cth3vstybeam_allnh->Write();

  xvsyh->Write();
  zvsth->Write();
  xvsth->Write();
  yvsth->Write();
  rvsth->Write();

  xvsy_p_cuth->Write();
  xvsy_p_cut2h->Write();

  xvsy_form_time_h->Write();
  xvsy_t_last_coll_h->Write();
  xvsy_ncollpart_h->Write();

  xvsy_ptsmall_ybeamh->Write();      
  xvsy_not_ptsmall_ybeamh->Write();  

  xvsyybeamh->Write();  
  xvszh->Write();
  xvszybeamh->Write();
  part_idybeamh->Write();
  phi_p_cuth->Write();
  phi_p_cut2h->Write();

  phi_p_posraph->Write();
  phi_p_negraph->Write();
  xvsy_p_posraph->Write();
  xvsy_p_negraph->Write();

  xvsy_p_binm4->Write();
  xvsy_p_binm3->Write();
  xvsy_p_binm2->Write();
  xvsy_p_binm1->Write();
  xvsy_p_bin0->Write();
  xvsy_p_binp1->Write();
  xvsy_p_binp2->Write();
  xvsy_p_binp3->Write();
  xvsy_p_binp4->Write();



  xvsy_p_pt0_binm4->Write();
  xvsy_p_pt0_binm3->Write();
  xvsy_p_pt0_binm2->Write();
  xvsy_p_pt0_binm1->Write();
  xvsy_p_pt0_bin0->Write();
  xvsy_p_pt0_binp1->Write();
  xvsy_p_pt0_binp2->Write();
  xvsy_p_pt0_binp3->Write();
  xvsy_p_pt0_binp4->Write();

  xvsy_p_pt1_binm4->Write();
  xvsy_p_pt1_binm3->Write();
  xvsy_p_pt1_binm2->Write();
  xvsy_p_pt1_binm1->Write();
  xvsy_p_pt1_bin0->Write();
  xvsy_p_pt1_binp1->Write();
  xvsy_p_pt1_binp2->Write();
  xvsy_p_pt1_binp3->Write();
  xvsy_p_pt1_binp4->Write();
  
  xvsy_p_pt2_binm4->Write();
  xvsy_p_pt2_binm3->Write();
  xvsy_p_pt2_binm2->Write();
  xvsy_p_pt2_binm1->Write();
  xvsy_p_pt2_bin0->Write();
  xvsy_p_pt2_binp1->Write();
  xvsy_p_pt2_binp2->Write();
  xvsy_p_pt2_binp3->Write();
  xvsy_p_pt2_binp4->Write();

  xvsy_p_pt3_binm4->Write();
  xvsy_p_pt3_binm3->Write();
  xvsy_p_pt3_binm2->Write();
  xvsy_p_pt3_binm1->Write();
  xvsy_p_pt3_bin0->Write();
  xvsy_p_pt3_binp1->Write();
  xvsy_p_pt3_binp2->Write();
  xvsy_p_pt3_binp3->Write();
  xvsy_p_pt3_binp4->Write();

  
  xvsyallh->Write();
  zvstallh->Write();
  xvsyspech->Write();
  xvszspech->Write();
  zvstspech->Write();
  ncollpartspech->Write();
    
  pth->Write();
  rapidityh->Write();
  phih->Write();
  etah->Write();

  pxh->Write();
  pyh->Write();
  pzh->Write();
  
  reactionplanerotationh->Write();  
  PSI1h->Write();
  PSI2h->Write();
  PSI3h->Write();
  PSI4h->Write();

  x0y0h->Write();  
  x0y0posraph->Write();
  x0y0negraph->Write();

  phixyh->Write(); 
  phixyposraph->Write();
  phixynegraph->Write();

  rsqh->Write(); 
  rsqposraph->Write(); 
  rsqnegraph->Write(); 

  rsqavh->Write();
  rsqavposraph->Write(); 
  rsqavnegraph->Write(); 

  eps1h->Write();
  eps2h->Write();
  eps3h->Write();
  eps4h->Write();
  eps5h->Write();
  eps6h->Write();

  eps1posraph->Write(); 
  eps2posraph->Write(); 
  eps3posraph->Write(); 
  eps4posraph->Write(); 
  eps5posraph->Write(); 
  eps6posraph->Write(); 

  eps1negraph->Write(); 
  eps2negraph->Write(); 
  eps3negraph->Write(); 
  eps4negraph->Write(); 
  eps5negraph->Write(); 
  eps6negraph->Write(); 
 
  PSI1posraph->Write();
  PSI1negraph->Write();
  PSI2posraph->Write();
  PSI2negraph->Write();
  PSI3posraph->Write();
  PSI3negraph->Write();

  PSI1_roth->Write();
  PSI2_roth->Write();
  PSI3_roth->Write();
  PSI4_roth->Write();  
  
  v1_RP1h->Write();
  v2_RP1h->Write();
  v3_RP1h->Write();
  v4_RP1h->Write();
  v5_RP1h->Write();
  v6_RP1h->Write();
  v2_RP2h->Write();
  v4_RP2h->Write();
  v6_RP2h->Write();
  v3_RP3h->Write();
  v6_RP3h->Write();

  cth_1_rapidity_all_h->Write();
  cth_1_pt_all_h->Write();
  cth_1_b_all_h->Write();
  cth_1_rapidity_pip_h->Write();
  cth_1_pt_pip_h->Write();
  cth_1_b_pip_h->Write();
  cth_1_rapidity_pim_h->Write();
  cth_1_pt_pim_h->Write();
  cth_1_b_pim_h->Write();
  cth_1_rapidity_kp_h->Write();
  cth_1_pt_kp_h->Write();
  cth_1_b_kp_h->Write();
  cth_1_rapidity_km_h->Write();
  cth_1_pt_km_h->Write();
  cth_1_b_km_h->Write();
  cth_1_rapidity_pr_h->Write();
  cth_1_pt_pr_h->Write();
  cth_1_b_pr_h->Write();
  cth_1_rapidity_deut_h->Write();
  cth_1_pt_deut_h->Write();
  cth_1_b_deut_h->Write();

  cth_2_rapidity_all_h->Write();
  cth_2_pt_all_h->Write();
  cth_2_b_all_h->Write();
  cth_2_rapidity_pip_h->Write();
  cth_2_pt_pip_h->Write();
  cth_2_b_pip_h->Write();
  cth_2_rapidity_pim_h->Write();
  cth_2_pt_pim_h->Write();
  cth_2_b_pim_h->Write();
  cth_2_rapidity_kp_h->Write();
  cth_2_pt_kp_h->Write();
  cth_2_b_kp_h->Write();
  cth_2_rapidity_km_h->Write();
  cth_2_pt_km_h->Write();
  cth_2_b_km_h->Write();
  cth_2_rapidity_pr_h->Write();
  cth_2_pt_pr_h->Write();
  cth_2_b_pr_h->Write();
  cth_2_rapidity_deut_h->Write();
  cth_2_pt_deut_h->Write();
  cth_2_b_deut_h->Write();

  cth_3_rapidity_all_h->Write();
  cth_3_pt_all_h->Write();
  cth_3_b_all_h->Write();
  cth_3_rapidity_pip_h->Write();
  cth_3_pt_pip_h->Write();
  cth_3_b_pip_h->Write();
  cth_3_rapidity_pim_h->Write();
  cth_3_pt_pim_h->Write();
  cth_3_b_pim_h->Write();
  cth_3_rapidity_kp_h->Write();
  cth_3_pt_kp_h->Write();
  cth_3_b_kp_h->Write();
  cth_3_rapidity_km_h->Write();
  cth_3_pt_km_h->Write();
  cth_3_b_km_h->Write();
  cth_3_rapidity_pr_h->Write();
  cth_3_pt_pr_h->Write();
  cth_3_b_pr_h->Write();
  cth_3_rapidity_deut_h->Write();
  cth_3_pt_deut_h->Write();
  cth_3_b_deut_h->Write();

  //---------------------------------------------

  cth_3_rapidity_ptcut1_cent_all_h->Write();
  cth_3_rapidity_ptcut1_cent_pip_h->Write();
  cth_3_rapidity_ptcut1_cent_pim_h->Write();
  cth_3_rapidity_ptcut1_cent_kp_h->Write();
  cth_3_rapidity_ptcut1_cent_km_h->Write();
  cth_3_rapidity_ptcut1_cent_pr_h->Write();
  cth_3_rapidity_ptcut2_cent_pr_h->Write();
  cth_3_rapidity_ptcut1_cent_deut_h->Write();
  cth_3_rapidity_ptcut2_cent_deut_h->Write();
  
  cth_3_rapidity_ptcut1_mid_all_h->Write();
  cth_3_rapidity_ptcut1_mid_pip_h->Write();
  cth_3_rapidity_ptcut1_mid_pim_h->Write();
  cth_3_rapidity_ptcut1_mid_kp_h->Write();
  cth_3_rapidity_ptcut1_mid_km_h->Write();
  cth_3_rapidity_ptcut1_mid_pr_h->Write();
  cth_3_rapidity_ptcut2_mid_pr_h->Write();
  cth_3_rapidity_ptcut1_mid_deut_h->Write();
  cth_3_rapidity_ptcut2_mid_deut_h->Write();
  
  cth_3_rapidity_ptcut1_periph_all_h->Write();
  cth_3_rapidity_ptcut1_periph_pip_h->Write();
  cth_3_rapidity_ptcut1_periph_pim_h->Write();
  cth_3_rapidity_ptcut1_periph_kp_h->Write();
  cth_3_rapidity_ptcut1_periph_km_h->Write();
  cth_3_rapidity_ptcut1_periph_pr_h->Write();
  cth_3_rapidity_ptcut2_periph_pr_h->Write();
  cth_3_rapidity_ptcut1_periph_deut_h->Write();
  cth_3_rapidity_ptcut2_periph_deut_h->Write();
  
  
  cth_3_pt_ycut1_cent_all_h->Write();
  cth_3_pt_ycut1_cent_pip_h->Write();
  cth_3_pt_ycut1_cent_pim_h->Write();
  cth_3_pt_ycut1_cent_kp_h->Write();
  cth_3_pt_ycut1_cent_km_h->Write();
  cth_3_pt_ycut1_cent_pr_h->Write();
  cth_3_pt_ycut1_cent_deut_h->Write();
  
  cth_3_pt_ycut1_mid_all_h->Write();
  cth_3_pt_ycut1_mid_pip_h->Write();
  cth_3_pt_ycut1_mid_pim_h->Write();
  cth_3_pt_ycut1_mid_kp_h->Write();
  cth_3_pt_ycut1_mid_km_h->Write();
  cth_3_pt_ycut1_mid_pr_h->Write();
  cth_3_pt_ycut1_mid_deut_h->Write();
  
  cth_3_pt_ycut1_periph_all_h->Write();
  cth_3_pt_ycut1_periph_pip_h->Write();
  cth_3_pt_ycut1_periph_pim_h->Write();
  cth_3_pt_ycut1_periph_kp_h->Write();
  cth_3_pt_ycut1_periph_km_h->Write();
  cth_3_pt_ycut1_periph_pr_h->Write();
  cth_3_pt_ycut1_periph_deut_h->Write();    
  
  cth_3_b_ycut1_ptcut1_all_h->Write();
  cth_3_b_ycut1_ptcut1_pip_h->Write();
  cth_3_b_ycut1_ptcut1_pim_h->Write();
  cth_3_b_ycut1_ptcut1_kp_h->Write();
  cth_3_b_ycut1_ptcut1_km_h->Write();
  cth_3_b_ycut1_ptcut1_pr_h->Write();
  cth_3_b_ycut1_ptcut2_pr_h->Write();
  cth_3_b_ycut1_ptcut1_deut_h->Write();
  cth_3_b_ycut1_ptcut2_deut_h->Write();


    eps_x_rapidity_ptcut1_cent_pr_h->Write();
  eps_x_rapidity_ptcut2_cent_pr_h->Write();
  eps_x_rapidity_ptcut1_mid_pr_h->Write();
  eps_x_rapidity_ptcut2_mid_pr_h->Write();
  eps_x_rapidity_ptcut1_periph_pr_h->Write();
  eps_x_rapidity_ptcut2_periph_pr_h->Write();
  eps_x_pt_ycut1_cent_pr_h->Write(); 
  eps_x_pt_ycut1_mid_pr_h->Write();  
  eps_x_pt_ycut1_periph_pr_h->Write();
  eps_x_b_ycut1_ptcut1_pr_h->Write(); 
  eps_x_b_ycut1_ptcut2_pr_h->Write(); 
  eps_y_rapidity_ptcut1_cent_pr_h->Write();
  eps_y_rapidity_ptcut2_cent_pr_h->Write();
  eps_y_rapidity_ptcut1_mid_pr_h->Write(); 
  eps_y_rapidity_ptcut2_mid_pr_h->Write(); 
  eps_y_rapidity_ptcut1_periph_pr_h->Write(); 
  eps_y_rapidity_ptcut2_periph_pr_h->Write(); 
  eps_y_pt_ycut1_cent_pr_h->Write(); 
  eps_y_pt_ycut1_mid_pr_h->Write();  
  eps_y_pt_ycut1_periph_pr_h->Write();
  eps_y_b_ycut1_ptcut1_pr_h->Write(); 
  eps_y_b_ycut1_ptcut2_pr_h->Write(); 
  rsq_rapidity_ptcut1_cent_pr_h->Write();
  rsq_rapidity_ptcut2_cent_pr_h->Write();
  rsq_rapidity_ptcut1_mid_pr_h->Write(); 
  rsq_rapidity_ptcut2_mid_pr_h->Write(); 
  rsq_rapidity_ptcut1_periph_pr_h->Write();
  rsq_rapidity_ptcut2_periph_pr_h->Write();
  rsq_pt_ycut1_cent_pr_h->Write(); 
  rsq_pt_ycut1_mid_pr_h->Write();   
  rsq_pt_ycut1_periph_pr_h->Write();
  rsq_b_ycut1_ptcut1_pr_h->Write(); 
  rsq_b_ycut1_ptcut2_pr_h->Write(); 


  
  //------------------------------------------------
  
  cth_1_eta_all_h->Write();
  cth_1_eta_pip_h->Write();
  cth_1_eta_pim_h->Write();
  cth_1_eta_kp_h->Write();
  cth_1_eta_km_h->Write();
  cth_1_eta_pr_h->Write();
  cth_1_eta_deut_h->Write();  

  cth_2_eta_all_h->Write();
  cth_2_eta_pip_h->Write();
  cth_2_eta_pim_h->Write();
  cth_2_eta_kp_h->Write();
  cth_2_eta_km_h->Write();
  cth_2_eta_pr_h->Write();
  cth_2_eta_deut_h->Write();
  
  cth_3_eta_all_h->Write();
  cth_3_eta_pip_h->Write();
  cth_3_eta_pim_h->Write();
  cth_3_eta_kp_h->Write();
  cth_3_eta_km_h->Write();
  cth_3_eta_pr_h->Write();
  cth_3_eta_deut_h->Write();

  rapidityext_all_h->Write();
  etaext_all_h->Write();
  cth_1_rapidityext_all_h->Write();
  cth_1_etaext_all_h->Write();
  rapidityext_pn_h->Write();
  etaext_pn_h->Write();
  cth_1_rapidityext_pn_h->Write();
  cth_1_etaext_pn_h->Write();
  rapidityext_pnbar_h->Write();
  etaext_pnbar_h->Write();
  cth_1_rapidityext_pnbar_h->Write();
  cth_1_etaext_pnbar_h->Write();
  
  costermh->Write();

  myfile->Close();
}

//====================================================================================
// reading stuff in
//====================================================================================
//bool readit(int itype=0, int ido=1, int itrack=1){
bool readit(int itype, int ido, int itrack){
  // itype =1: ampt  2: rqmd 3: my generator  4:smash 5:jam1
  // ido =0:open file, 1:read event 2:read particle -1:close file  3: for SMASH read b_impact after all particles read
  //  std::cout<<" in readit itype="<<itype<<" ido="<<ido<<std::endl;

  if(itype==5)return true; // for itype=5 (jam1) do read in main loop
  
  if(ido==3 && itype !=4) return true;  // only for smash, use ido =3 to get b_impact
  
  static int thrownevents=0; // for itype=3
  static int thrownparticle=0; // for itype=3

  // ---> Define event variables to be read from file
  int evnr=0, ntracks=0, aProj=0, zProj=0, aTarg=0, zTarg=0;
  float b = 0., ekin = 0.;
  
  int ityp=0, i3=0, ichg=0, pid=0;
  float ppx=0., ppy=0., ppz=0., m=0.,theta=0.0, eta=0.0;
  float pp0=0., rr0=0., rrx=0., rry=0., rrz=0.; 
  
  int nKaon=0, nLambda=0, tot=0, tot2=0;
  double px1[10], py1[10], px2[10], py2[10], pt1[10], pt2[10];

	 
  // AMPT
  if(itype==1){
    if(ido==0){
      file.open("amptdata.dat", std::fstream::in);
      if(!file){
      	std::cout<<" amptdata.dat not found"<<std::endl;
      	return false;
      }
    }
    
    if(ido==1){
      if(!(file>>eventnum>>testdummy>>number_particles>>b_impact>>npart_proj>>npart_targ>>nelas_proj>>ninelas_proj>>nelas_targ>>ninelas_targ>>dummy2)){
	std::cout<<" reached end of file "<<std::endl;
	return false;
      }
      //      std::cout<<" in readit "<<number_particles<<std::endl;

    }
    
    if(ido==2){
      if(!(file>>part_id_in>>px_in>>py_in>>pz_in>>mass_in>>x_in>>y_in>>z_in>>t_in)){
	std::cout<<" reached end of file "<<std::endl;
	return false;
      }
    }
    
    if(ido==-1){
      file.close();
      std::cout<<" closed input file"<<std::endl;
    }
  }


  //==============================================================
    // SMASH

  char line[200];
  char line1[200];
  char line2[200];
  char line3[200];
  char line4[200];
  char line5[200];
  char so_scprojtarg[200];
  string string1;

  string tmp;
  
  // so - smash oscar
  double so_t,so_x,so_y,so_z,so_mass,so_p0,so_px,so_py,so_pz,so_bimpact;
  int so_pdg,so_id,so_charge;
  int so_evttmp, so_end;

  if(itype==4){
    if(ido==0){
      file.open("smashdata.oscar", std::fstream::in);
      if(!file){
      	std::cout<<" smashdata.oscar not found"<<std::endl;
      	return false;
      }

      std::getline(file,string1);
      //      std::cout<<string1<<std::endl;
      std::getline(file,string1);
      //      std::cout<<string1<<std::endl;
      std::getline(file,string1);
      //      std::cout<<string1<<std::endl;
  
    }
    
    if(ido==1){

      if(!std::getline(file,string1)){
	std::cout<<" 1c-reached end of file "<<std::endl;
	return false;
      }

      //      std::cout<<" 5 output is string1: "<<string1<<std::endl;
      //	    std::getline(file,string1);
      int nn=string1.length();
      char char_array[nn+1];
      strcpy(char_array,string1.c_str());
      //      std::cout<<" char_array="<<char_array<<std::endl;
      //      for (int i=0; i<nn;i++)std::cout<<char_array[i];
      //      std::cout<<std::endl;
      std::sscanf(char_array,"%s %s %i %s %i",line,line1,&eventnum,line2,&number_particles);
      //      std::cout<<" converted string1 eventnum="<<eventnum<<" number_particles="<<number_particles<<std::endl;
      //      std::cout<<" converted string1  line="<<line<<" line1="<<line1<<" line2="<<line2<<std::endl;
      //      std::cout<<string1<<std::endl;

      //      std::cout<<" line="<<line<<" "<<line1<<std::endl;;
      if(!strncmp(line,"#!OSCAR2013",8)){
	//	std::cout<<" reached end of record "<<std::endl;
	std::getline(file,string1);
	//	std::cout<<string1<<std::endl;
	std::getline(file,string1);
	//	std::cout<<string1<<std::endl;
	if(!(file>>line>>line1>>eventnum>>line2>>number_particles)){
	  std::cout<<" 1a-reached end of file "<<std::endl;
	  return false;
	}
      }

      //      std::cout<<" in readit event number="<<eventnum<<" number particles="<<number_particles<<std::endl;
    } // end ido=1


    //    double p0_in;
    //    int so_id, charge_in;

    //   double form_time_in,xsecfac_in,t_last_coll_in;
    //   int ncoll_in,proc_id_origin_in,proc_type_origin_in,pdg_mother1_in,pdg_mother2_in;

   //ncoll_in,form_time_in,xsecfac_in,proc_id_origin_in,proc_type_origin_in,t_last_coll_in,pdg_mother1_in,pdg_mother2_in

      //ncoll_in<<" "<<form_time_in<<" "<<xsecfac_in<<" "<<proc_id_origin_in<<" "<<proc_type_origin_in<<" "<<t_last_coll_in<<" "<<pdg_mother1_in<<" "<<pdg_mother2_in


      //ncoll_in>>form_time_in>>xsecfac_in>>proc_id_origin_in>>proc_type_origin_in>>t_last_coll_in>>pdg_mother1_in>>pdg_mother2_in

   
    if(ido==2){
      //      if(!(file>>t_in>>x_in>>y_in>>z_in>>mass_in>>p0_in>>px_in>>py_in>>pz_in>>part_id_in>>so_id>>charge_in)){
      if(!(file>>t_in>>x_in>>y_in>>z_in>>mass_in>>p0_in>>px_in>>py_in>>pz_in>>part_id_in>>so_id>>charge_in>>ncoll_in>>form_time_in>>xsecfac_in>>proc_id_origin_in>>proc_type_origin_in>>t_last_coll_in>>pdg_mother1_in>>pdg_mother2_in)){
	std::cout<<" 2-reached end of file "<<std::endl;
	return false;
      } 
      //      std::cout<<" "<<t_in<<" "<<x_in<<" "<<y_in<<" "<<z_in<<" "<<mass_in<<" "<<p0_in<<" "<<px_in<<" "<<py_in<<" "<<pz_in<<" "<<part_id_in<<" "<<so_id<<" "<<charge_in<<std::endl;
      //      std::cout<<" "<<t_in<<" "<<x_in<<" "<<y_in<<" "<<z_in<<" "<<mass_in<<" "<<p0_in<<" "<<px_in<<" "<<py_in<<" "<<pz_in<<" "<<part_id_in<<" "<<so_id<<" "<<charge_in<<" new="<<ncoll_in<<" "<<form_time_in<<" "<<xsecfac_in<<" "<<proc_id_origin_in<<" "<<proc_type_origin_in<<" "<<t_last_coll_in<<" "<<pdg_mother1_in<<" "<<pdg_mother2_in<<std::endl;

      //#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge
      // 200 -38.335 77.2283 92.7543 0.938 1.2793081 -0.265927134 0.535650097 0.631780772 2212 935 1 
      // ncoll form_time xsecfac proc_id_origin proc_type_origin t_last_coll pdg_mother1 pdg_mother2
      // ncoll  form_time   xsecfrac    proc_id_origin   proc_type_origin    t_last_coll  pdg_mother1  pdg_mother2
      // 4      8.65802     1           898               1                  14.293        2114         0

    } // end ido=2
    
    if(ido==-1){
      file.close();
      std::cout<<" closed input file"<<std::endl;
    } // end ido=-1


    int so_evttmp,so_end;
    char so_scprojtarg[10];
    if(ido==3){
       if(!(file>>line>>line1>>so_evttmp>>line2>>so_end>>line3>>b_impact>>line4>>so_scprojtarg)){
	std::cout<<" 3-reached end of file "<<std::endl;
	return true;
      }
       //      std::cout<<" end of event b_impact="<<b_impact<<std::endl;
      //      std::cout<<line<<" "<<line1<<" "<<line2<<" "<<line3<<" "<<line4<<" "<<so_scprojtarg<<std::endl;

      std::getline(file,string1);
      //      std::cout<<" test "<<string1<<std::endl;
    } // end ido=3


  }
  //============================================================================================================


  
  // RQMD  
  if(itype==2){
    if(ido==0){
      std::cout<<" opening input file"<<std::endl;
      fInputFile = fopen("urqmddata.f14","r");
      if(!fInputFile){
	std::cout<<"-E- Can not open input file."<<std::endl;
      	return false;	
	//	exit(1);
      }
      std::cout<<" input file open"<<std::endl;
      Int_t myEvents=0;
    }
    
    if(ido==1){
      //      std::cout<<" reading header"<<std::endl;
    nKaon=0;nLambda=0;
    
    for(int i=0; i<10; i++){
      px1[i]=0.;py1[i]=0.;px2[i]=0.;py2[i]=0.;pt1[i]=0.;pt2[i]=0.;
    }

    // ---> Read and check first event header line from input file
    char read[200];
    //     std::cout<<" 0.1 "<<std::endl;   
    fgets(read, 200, fInputFile);
    //    std::cout<<" read="<<read<<std::endl;
    //    std::cout<<" 0.15 "<<std::endl;   
    if ( feof(fInputFile) ) {
      std::cout << "-I- End of input file reached." << std::endl;
      fclose(fInputFile);
      fInputFile = NULL;
      //   exit(1);
      //      break;
      return false;
    }
    //     std::cout<<" 0.2 "<<std::endl;       
    if ( read[0] != 'U' ) {
      std::cout << "-E- Wrong event header" << std::endl;
      //exit(1);
      return false;
      //      break;
    }

    // ---> Read rest of event header
    //    std::cout<<" 1 "<<std::endl;
    fgets(read, 26, fInputFile);
    fscanf(fInputFile, "%d", &aProj);
    fscanf(fInputFile, "%d", &zProj);
    fgets(read, 25, fInputFile);
    fscanf(fInputFile, "%d", &aTarg);
    fscanf(fInputFile, "%d", &zTarg);
    fgets(read, 200, fInputFile);
    fgets(read, 200, fInputFile);
    fgets(read, 36, fInputFile);
    fscanf(fInputFile, "%f", &b);
    fgets(read, 200, fInputFile);
    fgets(read, 39, fInputFile);
    fscanf(fInputFile, "%e", &ekin);
    fgets(read, 200, fInputFile);
    fgets(read, 7, fInputFile);
    fscanf(fInputFile, "%d", &evnr);
    fgets(read, 200, fInputFile);
    // update # lines from 8 to 11 RKS March 17, 2022
    for (int iline=0; iline<11; iline++)  { fgets(read, 200,fInputFile); }
    fscanf(fInputFile, "%d", &ntracks);
    fgets(read, 200, fInputFile);
    //    std::cout<<" read 2 "<<read<<std::endl;
    fgets(read, 200, fInputFile);
    //    std::cout<<" read 3 "<<read<<std::endl;
    //    std::cout<<aProj<<" "<<zProj<<" "<<aTarg<<" "<<zTarg<<" "<<b<<" "<<ekin<<" "<<evnr<<" "<<ntracks<<" "<<b<<std::endl;
    number_particles=ntracks;
    b_impact=b;
    }

    if(ido==2){
      //      std::cout<<" in readit itype=2 ido=2 itrack="<<itrack<<std::endl;
      // event cuts
      if(number_particles <= 0 || number_particles > 25000){
	//	std::cout<<" number of particles <= 0 or > 25000 -skipping event, number_particles="<<number_particles<<std::endl;
	return false;
      }    // event cut on track numbers to prevent segmentation violation

      //        if(centrality <= 0) continue;    // centrality starts at 1
      
      fscanf(fInputFile, "%e", &rr0);
      fscanf(fInputFile, "%e", &rrx);
      fscanf(fInputFile, "%e", &rry);
      fscanf(fInputFile, "%e", &rrz);
      fscanf(fInputFile, "%e", &pp0);
      
      fscanf(fInputFile, "%e", &ppx);
      fscanf(fInputFile, "%e", &ppy);
      fscanf(fInputFile, "%e", &ppz);
      fscanf(fInputFile, "%e", &m);
      fscanf(fInputFile, "%d", &ityp);
      fscanf(fInputFile, "%d", &i3);
      fscanf(fInputFile, "%d", &ichg);
      char read[200];
      fgets(read, 200, fInputFile);
      
      //      std::cout<<" trk="<<itrack<<" "<<ppx<<" "<<ppy<<" "<<ppz<<" "<<m<<" "<<ityp<<" "<<i3<<" "<<ichg<<std::endl;
      //      std::cout<<" cont "<<rr0<<" "<<rrx<<" "<<rry<<" "<<rrz<<" "<<pp0<<std::endl;

      part_id_in = ityp;
      px_in = ppx;
      py_in = ppy;
      pz_in = ppz;
      x_in = rrx;
      y_in = rry;
      z_in = rrz;
      t_in = rr0;
      Int_t charge = ichg;
      Double_t     pt = sqrt(ppx*ppx+ppy*ppy);
      Double_t pmag = sqrt(ppx*ppx+ppy*ppy+ppz*ppz);

      theta=TMath::ACos(ppz/pmag);
      double eta=- TMath::Log(tan(0.5*theta));
            
      Double_t  phi = TMath::ATan2( ppy, ppx );      
      Double_t     En1 = TMath::Sqrt( ppx*ppx + ppy*ppy + ppz*ppz + m*m );
      double rap = -9999;
      if( (En1 - ppz)!= 0.  ){
	//    Double_t   rap = 0.5*( TMath::Log( (En1+ppz)/(En1-ppz) ) );
	rap = 0.5*( TMath::Log( (En1+ppz)/(En1-ppz) ) );
      }
      Double_t y = rap;
      Double_t      e = pp0;
      mass_in=sqrt(e*e-pmag*pmag);
      //      std::cout<<" mass="<<mass_in<<" m from rqmd="<<m<<std::endl;
      
		// now change PID to standard PID
      //      std::cout<<" 2 part_id_in="<<part_id_in<<std::endl;       
	int pidtmp=part_id_in;
	//	std::cout<<" pidtmp="<<pidtmp<<" "<<" mass="<<mass_in<<" charge="<<charge<<std::endl;
	// protons
	if(pidtmp==1 && charge == 1){part_id_in= 2212;}
	else if(pidtmp==-1 && charge==-1){part_id_in= -2212;}
	// neutron
	else if(pidtmp==1 && charge == 0){part_id_in= 2112;}
	else if(pidtmp==-1 && charge==0){ part_id_in= -2112;}
	// pions
	else if(pidtmp==101 && charge == 1){part_id_in= 211;}
	else if(pidtmp==101 && charge == -1){part_id_in= -211;}
	else if(pidtmp==101 && charge == 0){part_id_in= 111;}
	//K
	else if(pidtmp==106 && charge == 1){part_id_in= 321;}
	else if(pidtmp==106 && charge == -1){part_id_in= -321;}
	else if(pidtmp==106 && charge == 0){part_id_in= -311;}
	else if(pidtmp==-106 && charge==0){part_id_in= 311;}
	else if(pidtmp==-106 && charge==-1){
	  //std::cout<<" pid mass charge="<<pidtmp<<" "<<mass_in<<" "<<charge<<std::endl;
	  part_id_in= -321;} // this makes no sense - do some checks on this
	//eta
	else if(pidtmp==102 && charge == 0){part_id_in= 2211;}
	//lambda
	else if(pidtmp==27 && charge == 0){part_id_in= 3122;}
	else if(pidtmp==-27 && charge == 0){part_id_in= -3122;}
	//sigma
	else if(pidtmp==40 && charge == 1){part_id_in= 3222;}
	else if(pidtmp==40 && charge == -1){part_id_in= 3112;}
	else if(pidtmp==40 && charge == 0){part_id_in= 3212;}
	else if(pidtmp==-40 && charge == -1){part_id_in= -3222;}
	else if(pidtmp==-40 && charge ==  1){part_id_in= -3112;}
	else if(pidtmp==-40 && charge == 0){part_id_in= -3212;}
	//Xsi
	else if(pidtmp==49 && charge == 0){part_id_in= 3322;}
	else if(pidtmp==49 && charge == -1){part_id_in= 3312;}
	else if(pidtmp==-49 && charge == 0){part_id_in= -3322;}
	else if(pidtmp==-49 && charge == 1){part_id_in= -3312;}
	// Omega
	else if(pidtmp==55 && charge == -1){part_id_in= 3334;}
	else if(pidtmp==-55 && charge == 1){part_id_in= -3334;}
	// D	
	else if(pidtmp==133 && charge == 0){part_id_in= 421;}
	else if(pidtmp==-133 && charge == 0){part_id_in= -421;}
	else if(pidtmp==133 && charge == 1){part_id_in= 411;}
	else if(pidtmp==-133 && charge == -1){part_id_in= -411;}

	// ityp from looking a fortran code
	// photon
	else if(pidtmp==100 && charge == 0){part_id_in= 22;}
	// D star for some reason
	else if(pidtmp==134 && charge == 1){part_id_in= 413;}
	else if(pidtmp==-134 && charge == -1){part_id_in= -413;}
	
	else if(pidtmp>1400){
	  //	  std::cout<<" pid mass charge="<<pidtmp<<" "<<mass_in<<" "<<charge<<std::endl;
	  part_id_in=pidtmp-1000;}
	else if(pidtmp<-1400){
	  //	  std::cout<<" pid mass charge="<<pidtmp<<" "<<mass_in<<" "<<charge<<std::endl;
	  part_id_in=pidtmp+1000;}
	
	else{std::cout<<" rqmd input particle id mass charge="<<pidtmp<<" "<<mass_in<<" "<<charge<<std::endl;}
       
    }
    
    if(ido==-1){
      file.close();
      std::cout<<" closed input file"<<std::endl;
    }
  }

  
  // THROW MY OWN EVENT
  if(itype==3){
    
    if(ido==0){
      throwevent(0);
    }
    
    if(ido==1){
      throwevent(1);
      thrownevents++;
      thrownparticle=0;
      eventnum=thrownevents;
      testdummy=0;
      number_particles=nparticles;
      b_impact=0;
      npart_proj=0;
      npart_targ=0;
      nelas_proj=0;
      ninelas_proj=0;
      nelas_targ=0;
      ninelas_targ=0;
      dummy2=0;
      return true;
    }
    //      std::cout<<" in readit "<<number_particles<<std::endl;
    
    
    if(ido==2){
      //    file>>part_id_in>>px_in>>py_in>>pz_in>>mass_in>>x_in>>y_in>>z_in>>t_in
      part_id_in=0;
      px_in=pttrack[thrownparticle]*cos(phitrack[thrownparticle]);
      py_in=pttrack[thrownparticle]*sin(phitrack[thrownparticle]);
      pz_in=pttrack[thrownparticle]*sinh(ytrack[thrownparticle]);
      mass_in=0;
      x_in=1;
      y_in=1;
      z_in=1;
      t_in=1;
      thrownparticle++;
      return true;
    }
    
    if(ido==-1){
      throwevent(-1);
      return true;
    }
  }  
  return true;
}

//====================================================================================
// Throwing my own event
//====================================================================================

int throwevent(int ido){

  // zero for init
  // 1  for throw event

  static bool iprintthrowevent=true;  
  static TF1 *flowf = new TF1("flowf",fitf,0,TWOPI,5);
  static  TRandom3 *ran0 = new TRandom3();
  static TRandom3 *ran1 = new TRandom3();   
  
  static int nparticlesav=200; //100
  static int ptdepv=0;
  /*  
  static double v1=0.;//0.10; //0.02   // 0.05; 
  static double v2=0.;//0.20; //0.05   //0.15;  
  static double v3=0.;//0.15; //0.01  //0.07;  
  static double v4=0.;//0.05; //0.02  //0.03;  
  */
  /*
  static double v1=0.02;   // 0.05; 
  static double v2=0.05;   //0.15;  
  static double v3=0.01;  //0.07;  
  static double v4=0.02;  //0.03;  
  */
  
  static double v1=0.10; //0.02   // 0.05; 
  static double v2=0.20; //0.05   //0.15;  
  static double v3=0.15; //0.01  //0.07;  
  static double v4=0.05; //0.02  //0.03;  
  
  /*
  static double v1=0.40; //0.02   // 0.05; 
  static double v2=0.60; //0.05   //0.15;  
  static double v3=0.30; //0.01  //0.07;  
  static double v4=0.20; //0.02  //0.03;  
  */
  static double PSI2Event=0.0;
  
  if(ido==0){
    std::cout<<" inputtype: throwevents INIT"<<std::endl;
    std::cout<<"Initial values: v1="<<v1<<"  v2="<<v2<<"  v3="<<v3<<"  v4="<<v4
	<<std::endl<<std::endl;
    
    std::cout<<" nparticles/event average(100) max is 10000, for pt dependence need at least 500 for T=300 MeV: ";
    cin>>nparticlesav;
    if(nparticlesav>10000)nparticlesav=10000;
    std::cout<<" nparticles set to "<<nparticlesav<<std::endl;
    
    std::cout<<" pt dependent v's 0-no 1-yes: ";
    cin>>ptdepv;
    std::cout<<" ptdepv set to "<<ptdepv<<std::endl;
    double PSI2=0.;
    
    flowf->SetParNames("PSI2","v1","v2","v3","v4");   
    flowf->SetParameters(PSI2,v1,v2,v3,v4);

    TCanvas *c1 = new TCanvas("c1","The FillRandom example",200,10,700,900);
    c1->SetFillColor(18);
    
    TPad *pad1 = new TPad("pad1","The pad with the function",0.05,0.50,0.95,0.95,21);
    TPad *pad2 = new TPad("pad2","The pad with the histogram",0.05,0.05,0.95,0.45,21);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    
    pad1->SetGridx();
    pad1->SetGridy();
    pad1->GetFrame()->SetFillColor(42);
    pad1->GetFrame()->SetBorderMode(-1);
    pad1->GetFrame()->SetBorderSize(5);
    flowf->SetLineColor(4);
    flowf->SetLineWidth(6);
    flowf->Draw();
    TPaveLabel *lfunction = new TPaveLabel(5,39,9.8,46,"v2 function");
    lfunction->SetFillColor(41);
    lfunction->Draw();
    c1->Update();
    //
    // Create a one dimensional histogram (one float per bin)
    // and fill it following the distribution in function flowf.
    //
    TH1D *h1f = new TH1D("h1f","Test of an Observed Distribution",200,0,TWOPI);
    
    pad2->cd();
    pad2->GetFrame()->SetFillColor(42);
    pad2->GetFrame()->SetBorderMode(-1);
    pad2->GetFrame()->SetBorderSize(5);
    h1f->SetFillColor(45);
    h1f->FillRandom("flowf",1000000);
    h1f->Draw();
    c1->Update();
    
    // the previous stuff was just for illustration
  }  //init


  if(ido==1){

    if(iprintthrowevent)std::cout<<" inputtype: throwevents THROWEVENT"<<std::endl;       
    iprintthrowevent=false;
    //   nparticles=nparticlesav;
   nparticles=ran0->Poisson(nparticlesav);
   if(nparticles<nparticlesmin)std::cout<<" nparticles small "<<nparticles<<std::endl;

   //double PSI2Event=0.0;
        PSI2Event=-0.2; // for testing
     //     PSI2Event=ran0->Rndm()*TWOPI-PI;

     //     PSI2h->Fill(PSI2Event);
     //     PSI2PSI2foundh->Fill(PSI2Event);
       

     flowf->SetParameters(PSI2Event,v1,v2,v3,v4); // for testing &&&
       
     for(int i=0; i<nparticles;i++){
       pttrack[i]=gRandom->Exp(0.300);  // .250 GeV tau
       //       pttrack[i]=gRandom->Rndm()*2.;  // flat in pt
       pttrack[i]+=0.180;
       double vv1=v1;
       double vv2=v2;
       double vv3=v3;
       double vv4=v4;
       
       double ppt=pttrack[i];
       if(ppt>1.)ppt=1;
       vv1=ppt*v1;
       vv2=ppt*v2;
       vv3=ppt*v3;
       vv4=ppt*v4;
       
       if(ptdepv)flowf->SetParameters(PSI2Event,vv1,vv2,vv3,vv4);      
       //       if(!ptdepv)flowf->SetParameters(PSI2Event,v1,v2,v3,v4); // for testing &&&
       
       phitrack[i]=flowf->GetRandom();
       ytrack[i]=(gRandom->Rndm()-.5)*5; // flat over 5 units of rapidity
       //      ytrack[i]=1.; // set all y to 1 now
        //       if(phitrack[i]>1. && phitrack[i]<1.9)i--; // mess up acceptance
      }

  return nparticles;       
  }
  if(ido==-1){
    std::cout<<" For THROWEVENT input values for v1="<<v1<<" v2="<<v2<<" v3="<<v3<<" v4="<<v4<<std::endl;
  }
  return nparticles;       
}

Double_t fitf(Double_t *xx,Double_t *par) {
  double v1=par[1];
  double v2=par[2];
  double v3=par[3];
  double v4=par[4];
  double PSI=par[0];
  
  double x=xx[0];  

  double  fitval = 1./TWOPI*(
			     1.+2*v1*cos(1*(x-PSI))+2*v2*cos(2*(x-PSI))+2*v3*cos(3*(x-PSI))+2*v4*cos(4*(x-PSI))
			     );  
  return fitval;
}

