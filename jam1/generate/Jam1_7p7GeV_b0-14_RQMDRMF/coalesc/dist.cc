#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>
#include <set>
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TFile.h"

#define B_MAXMUL 20000
#include "Cluster_Scanner.h" // Class Scanner
//#define __MNT_NO_TH1D_HIST
#include "inc/Scanner.h"
#define MNT_FAST_FILL
#include "inc/particles.h"
#include "inc/PDGData.h"

using namespace std;

void dist() {
	 ifstream input("filelist");
	 string line;
	 vector<string> InputList;
	 while (input >> line) {
		  InputList.push_back(line);
	 }
    Scanner _sc("filelist", "jam");
    Cluster_Scanner sc(_sc);

    PDGData pdg;

    ///////////////////////////////////////////////////////////////////////////
    // set parameter for coalecsce. 
    sc.mstc106 = 14;
    sc.R0_d = 4.0;
    sc.P0_d = 0.3;
    sc.R0_t = 4.0;
    sc.P0_t = 0.3;
    sc.R0_he4 = 4;
    sc.P0_he4 = 0.3;
    ///////////////////////////////////////////////////////////////////////////
		//
		//pt cut for flow calculation
		float PptCut = 1.0;
		float LptCut = 2.0;
		float DptCut = 2.0;
		float TptCut = 3.0;
		float He4ptCut = 4.0;

    ///////////////////////////////////////////////////////////////////////////
    // Histograms
    TH1D* hrefMult = new TH1D("hrefMult", "", 2000, 0, 2000);
   	TH1D* hrefMult1  = new TH1D("hrefMult1", "", 2000, 0, 2000);
	  TH1D* hb = new TH1D("hb", "b dist.", 150, 0, 15);
	  TH1D* hnpart = new TH1D("hnpart", "", 600, 0, 600);
		TH1D* hcent = new TH1D("hcent", "", 4, -0.5, 3.5);

		// pt and y
		TH2D* p_ptandy = new TH2D("p_ptandy", "p_ptandy", 2000, -5, 5, 2000, 0, 5);
		TH2D* l_ptandy = new TH2D("l_ptandy", "l_ptandy", 2000, -5, 5, 2000, 0, 5);
		TH2D* d_ptandy = new TH2D("d_ptandy", "d_ptandy", 2000, -5, 5, 2000, 0, 5);
		TH2D* t_ptandy = new TH2D("t_ptandy", "t_ptandy", 2000, -5, 5, 2000, 0, 5);
		TH2D* he3_ptandy = new TH2D("he3_ptandy", "he3_ptandy", 2000, -5, 5, 2000, 0, 5);
		TH2D* he4_ptandy = new TH2D("he4_ptandy", "he4_ptandy", 2000, -5, 5, 2000, 0, 5);
		
		Int_t Cent(-1);
	  Double_t ybins[21];
		for(int ii=0; ii<21;ii++) ybins[ii] = ii*0.1-1;

    ///////////////////////////////////////////////////////////////////////////
	 //root file loop
	 TFile* file;
	 TTree* tree;
	 int Entries = 0, refmult = 0, mul = 0;
	 const int MAXMUL=20000;
	 float x[MAXMUL];
	 float ft[MAXMUL];
	 float px[MAXMUL];
	 float py[MAXMUL];
	 float pz[MAXMUL];
	 int pid[MAXMUL];
	 float b;

	 // the initial trees
	 for(int i=0; i<InputList.size(); i++){
		  
		  file = TFile::Open(InputList[i].c_str());
		  file->GetObject("jam", tree);
		  tree->SetBranchAddress("t",ft);
		  tree->SetBranchAddress("x",x);
		  tree->SetBranchAddress("px",px);
		  tree->SetBranchAddress("py",py);
		  tree->SetBranchAddress("pz",pz);
		  tree->SetBranchAddress("mul",&mul);
		  tree->SetBranchAddress("pid",pid);
		  tree->SetBranchAddress("b",&b);
			
		  //Event loop
		  Entries = tree->GetEntries();
		  for(int j=0; j<Entries; j++){
			   tree->GetEntry(j);
					int cent = -1;
					if(b < 2.34) cent = 1;
					else if(b>2.34 && b<8.5) cent = 2; // 5-40%
					else if(b>8.5 && b<12.1) cent = 3; // 40-80%
					if(cent<0) continue;
			
			   refmult = 0;
			   //particle loop
			   for(int k=0; k<=mul; k++){
				
					double pcm  = sqrt(px[k]*px[k]+py[k]*py[k]+pz[k]*pz[k]);
					double pt   = sqrt(px[k]*px[k]+py[k]*py[k]);
					if (pcm > 10.) continue;
					double apid = abs(pid[k]);
					double mass = pdg[apid]->m0;
					double ecm   = sqrt(pcm*pcm + mass*mass); 
					double ycm   = 0.5*log((ecm + pz[k])/(ecm - pz[k]));

					float charge = (pdg[apid]->charge) / 3.0;
					if (charge > 0) ++refmult;
					
					if (pid[k] == 2212 ) { 
						p_ptandy->Fill(ycm, pt);	
			    }//end of proton loop

					if (pid[k] == 3122 ) { 
						l_ptandy->Fill(ycm, pt);	
					} // end of lambda loop
				} // end of track loop
			  hrefMult->Fill(refmult);
		  } //end of  events loop
	 }//end of file loop



   //light nucleus
	// Event loop
	for (int et = 1; sc.hasNextEntry(); ++et) {
		if(et % 10000 == 1) printf("#%d events...\n", et); 
		int mul = sc.mul;
		int refmult1 = 0; // total charged multiplicity

		// Particle Loop
		for (int pr = 0; pr < mul; ++pr) {
			int pid   = sc.pid[pr];
			int apid  = abs(pid);

			float charge = pdg[apid]->charge / 3.0;
			if (charge > 0) ++refmult1;

			hb->Fill(sc.b);
			hnpart->Fill(sc.Npart);
		} // end of particle loop
	
		hrefMult1->Fill(refmult1);

		if(sc.b < 2.34) Cent = 1;
		else if(sc.b>2.34 && sc.b<8.5) Cent = 2; // 5-40%
		else if(sc.b>8.50 && sc.b<11.6) Cent = 3;

		if(Cent<0) continue;

		hcent->Fill(Cent);

		for (int pr = 0; pr < sc.cmul; ++pr) {
			double px    = sc.cpx[pr];
			double py    = sc.cpy[pr];
			double pz    = sc.cpz[pr];
			double mass  = sc.cmass[pr];
		
			int     ncl    = sc.npid[pr];
			int*    cpid   = sc.cpid[pr];
			float* cpx_p   = sc.cpx_p[pr];
			float* cpy_p   = sc.cpy_p[pr];
			float* cpz_p   = sc.cpz_p[pr];
			float* cmass_p = sc.cmass_p[pr];
			float* ct_p    = sc.ct_p[pr];

			double pt2   = px*px + py*py;
			if (pt2 == 0) continue;
			double pt    = sqrt(pt2);
			double pcm2  = pt2 + pz*pz;
			double pcm   = sqrt(pcm2);
			//if (pcm > 10.) continue;
			double ecm   = sqrt(pcm2 + mass*mass);
			double ycm   = 0.5*log((ecm + pz)/(ecm - pz));
			double aycm  = fabs(ycm);
			double eta   = 0.5*log((pcm + pz)/(pcm - pz));
			double aeta  = fabs(eta);

			int clusterType = sc.ClusterType(pr);
			if (clusterType == Cluster_Scanner::ksUnknown) // unwanted type
					continue;

			if (clusterType == Cluster_Scanner::ksDeuteron) {
			  d_ptandy->Fill(ycm, pt);	
			}  // end of deuteron


			if (clusterType == Cluster_Scanner::ksTriton) {
			  t_ptandy->Fill(ycm, pt);	
			} // end of triton


			if (clusterType == Cluster_Scanner::ksHelium3) {
			  he3_ptandy->Fill(ycm, pt);	
			}


			if (clusterType == Cluster_Scanner::ksAlpha) {
			  he4_ptandy->Fill(ycm, pt);	
			}
			
		} // end of cluster loop
 } //end of  Events loop
    ///////////////////////////////////////////////////////////////////////////

  TFile fdist("fdist.root", "recreate");
  fdist.cd();

  hrefMult->Write();
	hrefMult1->Write();
	hb->Write();
	hnpart->Write();
	hcent->Write();

  p_ptandy->Write();
  l_ptandy->Write();
  d_ptandy->Write();
  t_ptandy->Write();
  he3_ptandy->Write();
  he4_ptandy->Write();
  fdist.Close();

} // void dist

int main() {
    dist();
    return 0;
}
