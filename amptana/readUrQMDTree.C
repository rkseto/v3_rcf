// load c++ and c headers
#include <iostream>
#include <stdio.h>

// load ROOT headers
#include "TChain.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TLeaf.h"

// main function
Int_t readUrQMDTree () {
    
	// load events
    char infile [ 200 ]; 
	sprintf ( infile , "/home/seto/Desktop/work/programs/urqmd-3.4/data/urqmd3.4_sim_3gev_100_test_00.root" );
    Int_t nfile = 0; TChain * chain = new TChain ( "Autree" ); nfile += chain -> Add( infile );
    Int_t nentries = chain -> GetEntries (); std::cout << std::endl << "Added " << nfile << " files, " << "# of events is " << nentries << std::endl << std::endl;
    
    // output filename
	TFile * f_out = NULL;
    f_out = new TFile ( "test_histograms.root" , "recreate" );
	
    // histograms
	
	// reference multiplicity
	TH1D * hrefmult = new TH1D ( "hrefmult" , "Reference multiplicity distribution" , 2000 , 0.5 , 2000 + 0.5 );
    
	// real multiplicity
	TH1D * hrealmult = new TH1D ( "hrealmult" , "Real multiplicity distribution" , 2000 , 0.5 , 2000 + 0.5 );
    
	// all tracks
    TH1D * hpx = new TH1D ( "hpx" , "p_{X}" , 1200 , -6.0 , 6.0 );
    hpx -> GetXaxis () -> SetTitle ( "p_{X} (GeV/c)" );
    //hpx -> GetYaxis () -> SetTitle ("(#frac{1}{p_{X}})#frac{dN}{dp_{X}} (GeV/c)^{-2}" );
    hpx -> GetXaxis () -> SetTitle ( "# of Tracks" );
    
    TH1D * hpy = new TH1D ( "hpy" , "p_{Y}" , 1200 , -6.0 , 6.0 );
    hpy -> GetXaxis () -> SetTitle ( "p_{Y} (GeV/c)" );
    //hpy -> GetYaxis () -> SetTitle ("(#frac{1}{p_{Y}})#frac{dN}{dp_{Y}} (GeV/c)^{-2}" );
    hpy -> GetXaxis () -> SetTitle ( "# of Tracks" );
    
    TH1D * hpz = new TH1D ( "hpz" , "p_{Z}" , 1200 , -6.0 , 6.0 );
    hpz -> GetXaxis () -> SetTitle ( "p_{Z} (GeV/c)" );
    //hpz -> GetYaxis () -> SetTitle ("(#frac{1}{p_{Z}})#frac{dN}{dp_{Z}} (GeV/c)^{-2}" );
    hpz -> GetXaxis () -> SetTitle ( "# of Tracks" );
    
    TH1D * hpt = new TH1D ( "hpt" , "p_{T}" , 1000 , 0.0 , 10.0 );
    hpt -> GetXaxis () -> SetTitle ( "p_{T} (GeV/c)" );
    //hpt -> GetYaxis () -> SetTitle ("(#frac{1}{p_{T}})#frac{dN}{dp_{T}} (GeV/c)^{-2}" );
    hpt -> GetXaxis () -> SetTitle ( "# of Tracks" );
    
    TH1D * heta = new TH1D ( "heta" , "#eta" , 1200 , -6.0 , 6.0 );
    heta -> GetXaxis () -> SetTitle ( "Pseudorapidity #eta" );
    heta -> GetYaxis () -> SetTitle ( "# of Tracks" );
    
    TH1D * hy = new TH1D ( "hy" , "y" , 1200 , -6.0 , 6.0 );
    hy -> GetXaxis () -> SetTitle ( "Rapidity y" );
    hy -> GetYaxis () -> SetTitle ( "# of Tracks" );
    
    TH1D * hphi = new TH1D ( "hphi" , "#phi" , 200 , -0.5 * TMath::Pi () , 2.5 * TMath::Pi () );
    hphi -> GetXaxis () -> SetTitle ( "#phi" );
    hphi -> GetYaxis () -> SetTitle ( "# of Tracks" );
	    
	// loop over events
    for ( Long64_t ievent = 0 ; ievent < nentries ; ievent ++ ) {
        
		if ( ievent + 1 > nentries ) break; 
        //if ( ievent + 1 > 1e4 ) break; 
		
		if ( ( ievent + 1 ) % 1000 == 0 ) std::cout << "Analyzed event " << ievent + 1 << std::endl << std::endl; chain -> GetEntry ( ievent );
        
        // read track multiplicity
        TLeaf * leaf_refmult = chain -> GetLeaf ( "tracknumber" );
        Int_t numberOfInputTracks = leaf_refmult -> GetValue ( 0 );
        
        // get event centrality bin
        TLeaf * leaf_centrality = chain -> GetLeaf ( "centrality" );
        Int_t centrality = leaf_centrality -> GetValue ( 0 );
		
		//  0 -  5% centrality = 1
		//  5 - 10% centrality = 2
		// 10 - 20% centrality = 3
		// 20 - 30% centrality = 4
		// 30 - 40% centrality = 5
		// 40 - 50% centrality = 6
		// 50 - 60% centrality = 7
		// 60 - 70% centrality = 8
		// 70 - 80% centrality = 9
		
		Double_t reaction_plane = 0.0;    // Reaction plane angle is always 0 in UrQMD
		
        // event cuts
        if(numberOfInputTracks <= 0 || numberOfInputTracks > 2500) continue;    // event cut on track numbers to prevent segmentation violation
        if(centrality <= 0) continue;    // centrality starts at 1
        
        // read TLeaves that store track parameters
        TLeaf * leaf_PID      = chain -> GetLeaf ( "PID" );
        TLeaf * leaf_Charge   = chain -> GetLeaf ( "Charge" );
        TLeaf * leaf_Px       = chain -> GetLeaf ( "Px" );
        TLeaf * leaf_Py       = chain -> GetLeaf ( "Py" );
        TLeaf * leaf_Pz       = chain -> GetLeaf ( "Pz" );
        TLeaf * leaf_Pt       = chain -> GetLeaf ( "Pt" );
        TLeaf * leaf_Pmag     = chain -> GetLeaf ( "Pmag" );
        TLeaf * leaf_Eta      = chain -> GetLeaf ( "Eta" );
        TLeaf * leaf_Phi      = chain -> GetLeaf ( "Phi" );
        TLeaf * leaf_E        = chain -> GetLeaf ( "E" );
        TLeaf * leaf_Rapidity = chain -> GetLeaf ( "Rapidity" );
        //TLeaf * leaf_M      = chain -> GetLeaf ( "M" );
        
		// define and initialize event/track-wise variables
        // track # and QC RFP, POI counters
        Int_t Ntrack = 0;

		// loop through tracks
        for ( Int_t itrack = 0 ; itrack < numberOfInputTracks ; itrack ++ ) 
		{
            
            Int_t       pid = leaf_PID -> GetValue ( itrack );
            Int_t    charge = leaf_Charge -> GetValue ( itrack );
            Double_t     px = leaf_Px -> GetValue ( itrack );
            Double_t     py = leaf_Py -> GetValue ( itrack );
            Double_t     pz = leaf_Pz -> GetValue ( itrack );
            Double_t     pt = leaf_Pt -> GetValue ( itrack );
            Double_t pmag = leaf_Pmag -> GetValue ( itrack );
            Double_t    eta = leaf_Eta -> GetValue ( itrack );
            Double_t    phi = leaf_Phi -> GetValue ( itrack );
            Double_t      e = leaf_E -> GetValue ( itrack );
            Double_t      y = leaf_Rapidity -> GetValue ( itrack );
            //Double_t    m = leaf_M -> GetValue ( itrack );
            
            if ( phi < 0.0                ) phi += 2.0 * TMath::Pi ();
            if ( phi > 2.0 * TMath::Pi () ) phi -= 2.0 * TMath::Pi ();
            
            // track cuts
            if(e <= 0) continue;
			
			// count track #
            Ntrack ++;
            
			// filling track-wise histograms
            // all tracks
            hpx -> Fill ( px ); hpy -> Fill ( py ); hpz -> Fill ( pz ); hpt -> Fill ( pt );
            heta -> Fill ( eta ); hy -> Fill ( y ); hphi -> Fill ( phi );

	   	    
            // (proton)
            if ( pid == 1   /*&& charge > 0*/ ) {
                
	      //   hpx_proton -> Fill ( px ); hpy_proton -> Fill ( py ); hpz_proton -> Fill ( pz ); hpt_proton -> Fill ( pt );
              //  heta_proton -> Fill ( eta ); hy_proton -> Fill ( y ); hphi_proton -> Fill ( phi );
                
            }
            
            // (pion)
            if ( pid == 101 /*&& charge > 0*/ ) {
                
	      // hpx_pion -> Fill ( px ); hpy_pion -> Fill ( py ); hpz_pion -> Fill ( pz ); hpt_pion -> Fill ( pt );
              //  heta_pion -> Fill ( eta ); hy_pion -> Fill ( y ); hphi_pion -> Fill ( phi );
                
            }
            
            // (kaon)
            if ( pid == 106 /*&& charge > 0*/ ) {
                
              //  hpx_kaon -> Fill ( px ); hpy_kaon -> Fill ( py ); hpz_kaon -> Fill ( pz ); hpt_kaon -> Fill ( pt );
              //  heta_kaon -> Fill ( eta ); hy_kaon -> Fill ( y ); hphi_kaon -> Fill ( phi );
                
            }
            
	  
    }    // track loop ends
        
		// fill event-wise histograms
        hrefmult -> Fill ( numberOfInputTracks );
        
	}    // Event loop ends
    
    f_out -> Write ();
    return 0;
    
}
