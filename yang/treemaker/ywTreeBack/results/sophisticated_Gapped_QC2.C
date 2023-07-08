// load c++ and c headers
#include <iostream>
#include <stdio.h>

// load ROOT headers
#include "TChain.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1D.h"
//#include "TH2D.h"
//#include "TH3D.h"
//#include "TProfile2D.h"
#include "TLeaf.h"
#include "TComplex.h"
#include "TVector3.h"
//#include "TRandom3.h"

// load STAR headers
//#include "StParticleTypes.hh"

Double_t GetCovariance ( Double_t w1 , Double_t w1v1 , Double_t w2 , Double_t w2v2 , Double_t w1w2 , Double_t w1w2v1v2 ) {
    
    // test parameter values to prevent Inf or NaN
    if( w1 == 0.0 || w2 == 0.0 || w1w2 == 0.0 ) return -999.0;
    if( w1 * w2 / w1w2 == 1.0 ) return -999.0;
    
    // compute covariance term
    Double_t cov = ( ( w1w2v1v2 / w1w2 ) - ( w1v1 / w1 ) * ( w2v2 / w2 ) ) / ( w1 * w2 / w1w2 - 1.0 );
    //if ( !TMath::Finite ( cov ) ) cov = -999.0;
    
    // return covariance value
    return cov;
    
}

Double_t GetError ( Double_t w1 , Double_t w1square , Double_t w1v1 , Double_t w1v1square ) {
    
    // test parameter values to prevent Inf or NaN
    if( w1 == 0.0 || w1square == 0.0 ) return -999.0;
    if( w1 * w1 / w1square == 1.0 ) return -999.0;
    
    // compute error term
    Double_t err = ( ( w1v1square / w1 ) - ( w1v1 / w1 ) * ( w1v1 / w1 ) ) / ( w1 * w1 / w1square - 1.0 );
    if ( /*!TMath::Finite ( err ) || w1 * w1 / w1square <= 1.0 ||*/ err <= 0 ) err = -999.0;
    
    // return error term
    return err;
    
}

// main function
Int_t sophisticated_Gapped_QC2 ( Int_t InputFileListIndex = 0 ) {
    
    // load events
    char infile [ 200 ], outfile [ 200 ]; sprintf ( infile , "/star/data01/pwg/ywu27/v3_fxt_qa/results/77A91B590300FCEA974693F69DAF73DA_*%d.root", InputFileListIndex);
    Int_t nfile = 0; TChain * chain = new TChain ( "Autree" ); nfile += chain -> Add( infile );
    Int_t nentries = chain -> GetEntries (); std::cout << "Added " << nfile << " files, " << "# of events is " << nentries << std::endl << std::endl;
    
    // define global variables
    const Int_t Ncentralities =  16;
    const Int_t 	  Ntracks = 196;
    const Int_t         Nbins =  30;
    
    // define flow analysis up down phase space limits
    Double_t rapidityLow = -1.5, rapidityHigh = 1.5;
    Double_t       ptLow =  0.2,       ptHigh = 5.0;
    Double_t midrapidity = -1.0450222;
    Double_t Mass_Proton = 0.938272, Mass_Pion = 0.139568, Mass_Kaon = 0.493646;
    
	// define and initialize global QC variables
    // track and bin multiplicity counters
	Int_t               multCounter [ Ntracks ][ Ncentralities ],
                     binMultCounter [ Nbins ][ Ntracks ][ Ncentralities ];
    
    // variables for averaging over all events
    // reference flow variables
	Double_t                 sumS11 [ Ntracks ][ Ncentralities ],
                            sumbS11 [ Ntracks ][ Ncentralities ],
                             sumM11 [ Ntracks ][ Ncentralities ],
    
                           sumCorr2 [ Ntracks ][ Ncentralities ],
    
                          sumAddon1 [ Ntracks ][ Ncentralities ],
                          sumAddon2 [ Ntracks ][ Ncentralities ],
                          sumAddon3 [ Ntracks ][ Ncentralities ],
                          sumAddon4 [ Ntracks ][ Ncentralities ];
    
    // differential flow variables
    Double_t                  sumMp [ Nbins ][ Ntracks ][ Ncentralities ],
                             sumM01 [ Nbins ][ Ntracks ][ Ncentralities ],
    
                    sumCorr2Reduced [ Nbins ][ Ntracks ][ Ncentralities ],
    
                          sumAddon5 [ Nbins ][ Ntracks ][ Ncentralities ],
                          sumAddon6 [ Nbins ][ Ntracks ][ Ncentralities ];
    
    // squared weights and variables for error estimations
    // reference flow variables
    Double_t          sumS11Squared [ Ntracks ][ Ncentralities ],
                     sumbS11Squared [ Ntracks ][ Ncentralities ],
                      sumM11Squared [ Ntracks ][ Ncentralities ],
    
                    sumCorr2Squared [ Ntracks ][ Ncentralities ],
    
                   sumAddon1Squared [ Ntracks ][ Ncentralities ],
                   sumAddon2Squared [ Ntracks ][ Ncentralities ],
                   sumAddon3Squared [ Ntracks ][ Ncentralities ],
                   sumAddon4Squared [ Ntracks ][ Ncentralities ];
    
    // differential flow variables
    Double_t           sumMpSquared [ Nbins ][ Ntracks ][ Ncentralities ],
                      sumM01Squared [ Nbins ][ Ntracks ][ Ncentralities ],
    
             sumCorr2ReducedSquared [ Nbins ][ Ntracks ][ Ncentralities ],
    
                   sumAddon5Squared [ Nbins ][ Ntracks ][ Ncentralities ],
                   sumAddon6Squared [ Nbins ][ Ntracks ][ Ncentralities ];
    
    // covariance terms
    // reference flow variables
    // multiplications of weights
    Double_t                   sumM11S11 [ Ntracks ][ Ncentralities ],
                              sumM11bS11 [ Ntracks ][ Ncentralities ],
                              sumS11bS11 [ Ntracks ][ Ncentralities ];
    
    // multiplications of variables
    Double_t 				  sumCorr2C1 [ Ntracks ][ Ncentralities ],
             				  sumCorr2C2 [ Ntracks ][ Ncentralities ],
             				  sumCorr2C3 [ Ntracks ][ Ncentralities ],
             				  sumCorr2C4 [ Ntracks ][ Ncentralities ],
    
             					 sumC1C2 [ Ntracks ][ Ncentralities ],
             					 sumC1C3 [ Ntracks ][ Ncentralities ],
             					 sumC1C4 [ Ntracks ][ Ncentralities ],
    
             					 sumC2C3 [ Ntracks ][ Ncentralities ],
             					 sumC2C4 [ Ntracks ][ Ncentralities ],
                                 sumC3C4 [ Ntracks ][ Ncentralities ];
    
    // differential flow variables
    // multiplications of weights
    Double_t                    sumS11mp [ Nbins ][ Ntracks ][ Ncentralities ],
                               sumS11M01 [ Nbins ][ Ntracks ][ Ncentralities ],
    
                               sumbS11mp [ Nbins ][ Ntracks ][ Ncentralities ],
                              sumbS11M01 [ Nbins ][ Ntracks ][ Ncentralities ],
    
                                sumM11mp [ Nbins ][ Ntracks ][ Ncentralities ],
                               sumM11M01 [ Nbins ][ Ntracks ][ Ncentralities ],
    
                                sumM01mp [ Nbins ][ Ntracks ][ Ncentralities ];
    
    // multiplications of variables
    Double_t        sumCorr2Corr2Reduced [ Nbins ][ Ntracks ][ Ncentralities ],
    
             				  sumCorr2C5 [ Nbins ][ Ntracks ][ Ncentralities ],
             				  sumCorr2C6 [ Nbins ][ Ntracks ][ Ncentralities ],
    
             		   sumCorr2ReducedC1 [ Nbins ][ Ntracks ][ Ncentralities ],
             		   sumCorr2ReducedC2 [ Nbins ][ Ntracks ][ Ncentralities ],
                       sumCorr2ReducedC3 [ Nbins ][ Ntracks ][ Ncentralities ],
                       sumCorr2ReducedC4 [ Nbins ][ Ntracks ][ Ncentralities ],
             		   sumCorr2ReducedC5 [ Nbins ][ Ntracks ][ Ncentralities ],
             		   sumCorr2ReducedC6 [ Nbins ][ Ntracks ][ Ncentralities ],
    
                                 sumC1C5 [ Nbins ][ Ntracks ][ Ncentralities ],
             					 sumC1C6 [ Nbins ][ Ntracks ][ Ncentralities ],
    
             					 sumC2C5 [ Nbins ][ Ntracks ][ Ncentralities ],
             					 sumC2C6 [ Nbins ][ Ntracks ][ Ncentralities ],
    
             					 sumC3C5 [ Nbins ][ Ntracks ][ Ncentralities ],
             					 sumC3C6 [ Nbins ][ Ntracks ][ Ncentralities ],
    
             					 sumC4C5 [ Nbins ][ Ntracks ][ Ncentralities ],
             					 sumC4C6 [ Nbins ][ Ntracks ][ Ncentralities ],
    
                                 sumC5C6 [ Nbins ][ Ntracks ][ Ncentralities ];
    
    // differential flows for each centrality and track multiplicity
    // QC{2}
	Double_t      v1_2 [ Nbins ][ Ntracks ][ Ncentralities ],
             v1_2Error [ Nbins ][ Ntracks ][ Ncentralities ];
    // QC{4}
    //              v1_4 [ Nbins ][ Ntracks ][ Ncentralities ],
    //         v1_4Error [ Nbins ][ Ntracks ][ Ncentralities ];
    
    // initialize QC global variables
    for ( Int_t icent = 0 ; icent < Ncentralities ; icent ++ ) {
        
        for ( Int_t itrack = 0 ; itrack < Ntracks ; itrack ++ ) {
            
            multCounter [ itrack ][ icent ] = 0;
            
            sumS11 [ itrack ][ icent ] = 0.0;
            sumbS11 [ itrack ][ icent ] = 0.0;
            sumM11 [ itrack ][ icent ] = 0.0;
            
            sumCorr2 [ itrack ][ icent ] = 0.0;
            
            sumAddon1 [ itrack ][ icent ] = 0.0;
            sumAddon2 [ itrack ][ icent ] = 0.0;
            sumAddon3 [ itrack ][ icent ] = 0.0;
            sumAddon4 [ itrack ][ icent ] = 0.0;
            
            sumS11Squared [ itrack ][ icent ] = 0.0;
            sumbS11Squared [ itrack ][ icent ] = 0.0;
            sumM11Squared [ itrack ][ icent ] = 0.0;
            
            sumCorr2Squared [ itrack ][ icent ] = 0.0;
            
            sumAddon1Squared [ itrack ][ icent ] = 0.0;
            sumAddon2Squared [ itrack ][ icent ] = 0.0;
            sumAddon3Squared [ itrack ][ icent ] = 0.0;
            sumAddon4Squared [ itrack ][ icent ] = 0.0;
            
            sumM11S11 [ itrack ][ icent ] = 0.0;
            sumM11bS11 [ itrack ][ icent ] = 0.0;
            
            sumS11bS11 [ itrack ][ icent ] = 0.0;
            
            sumCorr2C1 [ itrack ][ icent ] = 0.0;
            sumCorr2C2 [ itrack ][ icent ] = 0.0;
            sumCorr2C3 [ itrack ][ icent ] = 0.0;
            sumCorr2C4 [ itrack ][ icent ] = 0.0;
            
            sumC1C2 [ itrack ][ icent ] = 0.0;
            sumC1C3 [ itrack ][ icent ] = 0.0;
            sumC1C4 [ itrack ][ icent ] = 0.0;
            
            sumC2C3 [ itrack ][ icent ] = 0.0;
            sumC2C4 [ itrack ][ icent ] = 0.0;
            
            sumC3C4 [ itrack ][ icent ] = 0.0;
            
			for ( Int_t ibin = 0 ; ibin < Nbins ; ibin ++ ) {
                
                binMultCounter [ ibin ][ itrack ][ icent ] = 0;
                
                sumMp [ ibin ][ itrack ][ icent ] = 0.0;
                sumM01 [ ibin ][ itrack ][ icent ] = 0.0;
                
                sumCorr2Reduced [ ibin ][ itrack ][ icent ] = 0.0;
                
                sumAddon5 [ ibin ][ itrack ][ icent ] = 0.0;
                sumAddon6 [ ibin ][ itrack ][ icent ] = 0.0;
                
                sumMpSquared [ ibin ][ itrack ][ icent ] = 0.0;
                sumM01Squared [ ibin ][ itrack ][ icent ] = 0.0;
                
                sumCorr2ReducedSquared [ ibin ][ itrack ][ icent ] = 0.0;
                
                sumAddon5Squared [ ibin ][ itrack ][ icent ] = 0.0;
                sumAddon6Squared [ ibin ][ itrack ][ icent ] = 0.0;

                sumS11mp [ ibin ][ itrack ][ icent ] = 0.0;
                sumS11M01 [ ibin ][ itrack ][ icent ] = 0.0;

                sumbS11mp [ ibin ][ itrack ][ icent ] = 0.0;
                sumbS11M01 [ ibin ][ itrack ][ icent ] = 0.0;

                sumM11mp [ ibin ][ itrack ][ icent ] = 0.0;
                sumM11M01 [ ibin ][ itrack ][ icent ] = 0.0;
                
                sumM01mp [ ibin ][ itrack ][ icent ] = 0.0;
                
                sumCorr2Corr2Reduced [ ibin ][ itrack ][ icent ] = 0.0;
                
                sumCorr2C5 [ ibin ][ itrack ][ icent ] = 0.0;
                sumCorr2C6 [ ibin ][ itrack ][ icent ] = 0.0;
                
                sumCorr2ReducedC1 [ ibin ][ itrack ][ icent ] = 0.0;
                sumCorr2ReducedC2 [ ibin ][ itrack ][ icent ] = 0.0;
                sumCorr2ReducedC3 [ ibin ][ itrack ][ icent ] = 0.0;
                sumCorr2ReducedC4 [ ibin ][ itrack ][ icent ] = 0.0;
                sumCorr2ReducedC5 [ ibin ][ itrack ][ icent ] = 0.0;
                sumCorr2ReducedC6 [ ibin ][ itrack ][ icent ] = 0.0;
                
                sumC1C5 [ ibin ][ itrack ][ icent ] = 0.0;
                sumC1C6 [ ibin ][ itrack ][ icent ] = 0.0;
                
                sumC2C5 [ ibin ][ itrack ][ icent ] = 0.0;
                sumC2C6 [ ibin ][ itrack ][ icent ] = 0.0;
                
                sumC3C5 [ ibin ][ itrack ][ icent ] = 0.0;
                sumC3C6 [ ibin ][ itrack ][ icent ] = 0.0;
                
                sumC4C5 [ ibin ][ itrack ][ icent ] = 0.0;
                sumC4C6 [ ibin ][ itrack ][ icent ] = 0.0;
                
                sumC5C6 [ ibin ][ itrack ][ icent ] = 0.0;
                
			}
            
		}
        
    }    // QC variables initialization ends
    
    // read in Phi weight for QC from input file
    /*Double_t weight_entries = 0.0, weightPhi [ 134 ]; for ( Int_t ibin = 0 ; ibin < 134 ; ibin ++ ) { weightPhi [ ibin ] = 0.0; }
    
    TFile * f_in = new TFile ( "Gapped_QC2_v1_v2_test_10000_v1_result.root" , "read" );
    
    if ( !f_in -> IsOpen () ) std::cout << "No Phi weight profile!" << std::endl << std::endl;
    
    if ( f_in -> IsOpen () ) {
        
        std::cout << "Phi weight profile loaded!" << std::endl << std::endl;
        
        TH1D * hphi_temp = ( TH1D * ) f_in -> Get ( "hphi" );
        
        for ( Int_t ibin = 34 ; ibin <= 167 ; ibin ++ ) {
            
            weight_entries += hphi_temp -> GetBinContent ( ibin );
            weightPhi [ ibin - 34 ] = ( hphi_temp -> GetBinContent ( ibin ) > 0.0 )? TMath::Power ( hphi_temp -> GetBinContent ( ibin ) , -1.0 ) : 0.0;
            
        }
        
        for ( Int_t ibin = 0 ; ibin < 134 ; ibin ++ ) {
            
            weightPhi [ ibin ] *= ( weight_entries > 0 )? weight_entries / 134.0 : 0.0;
            
        }
        
        f_in -> Close ();
        
    }*/
    
    // output filename
    //Int_t ratio_integer = 100 * ratio;
	sprintf ( outfile , "fxt_3gev_ppiK_v1y_qc24_result_1_000%d.root" , InputFileListIndex );
    TFile * f_out = new TFile ( outfile , "recreate" );
    
    // histograms
    char name [ 200 ], description [ 200 ];
    
    // QC variables stored in these histograms
    /*TH2D * QCvariables2D [ 30 ]; TH3D * QCvariables3D [ 36 ];
    for ( Int_t line = 0 ; line < 36 ; line ++ ) {
        
        if ( line < 30 ) {
            
            sprintf ( name , "qc2d_line%d" , line ); sprintf ( description , "Gapped QC2 variables 2D line %d" , line );
            QCvariables2D [ line ] = new TH2D( name , description , Ntracks , 0.5 , Ntracks + 0.5 , Ncentralities , 0.5 , Ncentralities + 0.5 );
            
        }
        
        sprintf ( name , "qc3d_line%d" , line ); sprintf ( description , "Gapped QC2 variables 3D line %d" , line );
        QCvariables3D [ line ] = new TH3D( name , description , Nbins , 0.5 , Nbins + 0.5 , Ntracks , 0.5 , Ntracks + 0.5 , Ncentralities , 0.5 , Ncentralities + 0.5 );
        
    }*/
    
    // reference multiplicity
    /*TH1D * hrefmult = new TH1D ( "hrefmult" , "Reference multiplicity distribution" , 2000 , 0.5 , 2000 + 0.5 );
    
    // real multiplicity
    TH1D * hrealmult = new TH1D ( "hrealmult" , "Real multiplicity distribution" , 2000 , 0.5 , 2000 + 0.5 );
    
    // impact parameter
    //TH1D * himpact = new TH1D ( "himpact" , "Impact parameter b" , 320 , 0 , 16 );
	
	// reaction plane
	TH1D * hpsi = new TH1D ( "hpsi" , "#psi" , 200 , -0.5 * TMath::Pi () , 2.5 * TMath::Pi () );
    hpsi -> GetXaxis () -> SetTitle ( "#psi" );
    hpsi -> GetYaxis () -> SetTitle ( "# of Events" );
    
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
    
    // proton
    TH1D * hpx_proton = new TH1D ( "hpx_proton" , "Proton p_{X}" , 1200 , -6.0 , 6.0 );
    hpx_proton -> GetXaxis () -> SetTitle ( "p_{X} (GeV/c)" );
    //hpx_proton -> GetYaxis () -> SetTitle ("(#frac{1}{p_{X}})#frac{dN}{dp_{X}} (GeV/c)^{-2}" );
    hpx_proton -> GetXaxis () -> SetTitle ( "# of Tracks" );
    
    TH1D * hpy_proton = new TH1D ( "hpy_proton" , "Proton p_{Y}" , 1200 , -6.0 , 6.0 );
    hpy_proton -> GetXaxis () -> SetTitle ( "p_{Y} (GeV/c)" );
    //hpy_proton -> GetYaxis () -> SetTitle ("(#frac{1}{p_{Y}})#frac{dN}{dp_{Y}} (GeV/c)^{-2}" );
    hpy_proton -> GetXaxis () -> SetTitle ( "p_{Y} (GeV/c)" );
    
    TH1D * hpz_proton = new TH1D ( "hpz_proton" , "Proton p_{Z}" , 1200 , -6.0 , 6.0 );
    hpz_proton -> GetXaxis () -> SetTitle ( "p_{Z} (GeV/c)" );
    //hpz_proton -> GetYaxis () -> SetTitle ("(#frac{1}{p_{Z}})#frac{dN}{dp_{Z}} (GeV/c)^{-2}" );
    hpz_proton -> GetXaxis () -> SetTitle ( "# of Tracks" );
    
    TH1D * hpt_proton = new TH1D ( "hpt_proton" , "Proton p_{T}" , 1000 , 0.0 , 10.0 );
    hpt_proton -> GetXaxis () -> SetTitle ( "p_{T} (GeV/c)" );
    //hpt_proton -> GetYaxis () -> SetTitle ("(#frac{1}{p_{T}})#frac{dN}{dp_{T}} (GeV/c)^{-2}" );
    hpt_proton -> GetXaxis () -> SetTitle ( "# of Tracks" );
    
    TH1D * heta_proton = new TH1D ( "heta_proton" , "Proton #eta" , 1200 , -6.0 , 6.0 );
    heta_proton -> GetXaxis () -> SetTitle ( "Pseudorapidity #eta" );
    heta_proton -> GetYaxis () -> SetTitle ( "# of Tracks" );
    
    TH1D * hy_proton = new TH1D ( "hy_proton" , "Proton y" , 1200 , -6.0 , 6.0 );
    hy_proton -> GetXaxis () -> SetTitle ( "Rapidity y" );
    hy_proton -> GetYaxis () -> SetTitle ( "# of Tracks" );
    
    TH1D * hphi_proton = new TH1D ( "hphi_proton" , "Proton #phi" , 200 , -0.5 * TMath::Pi () , 2.5 * TMath::Pi () );
    hphi_proton -> GetXaxis () -> SetTitle ( "#phi" );
    hphi_proton -> GetYaxis () -> SetTitle ( "# of Tracks" );
    
    // pion
    TH1D * hpx_pion = new TH1D ( "hpx_pion" , "Pion p_{X}" , 1200 , -6.0 , 6.0 );
    hpx_pion -> GetXaxis () -> SetTitle ( "p_{X} (GeV/c)" );
    //hpx_pion -> GetYaxis () -> SetTitle ("(#frac{1}{p_{X}})#frac{dN}{dp_{X}} (GeV/c)^{-2}" );
    hpx_pion -> GetXaxis () -> SetTitle ( "# of Tracks" );
    
    TH1D * hpy_pion = new TH1D ( "hpy_pion" , "Pion p_{Y}" , 1200 , -6.0 , 6.0 );
    hpy_pion -> GetXaxis () -> SetTitle ( "p_{Y} (GeV/c)" );
    //hpy_pion -> GetYaxis () -> SetTitle ("(#frac{1}{p_{Y}})#frac{dN}{dp_{Y}} (GeV/c)^{-2}" );
    hpy_pion -> GetXaxis () -> SetTitle ( "p_{Y} (GeV/c)" );
    
    TH1D * hpz_pion = new TH1D ( "hpz_pion" , "Pion p_{Z}" , 1200 , -6.0 , 6.0 );
    hpz_pion -> GetXaxis () -> SetTitle ( "p_{Z} (GeV/c)" );
    //hpz_pion -> GetYaxis () -> SetTitle ("(#frac{1}{p_{Z}})#frac{dN}{dp_{Z}} (GeV/c)^{-2}" );
    hpz_pion -> GetXaxis () -> SetTitle ( "# of Tracks" );
    
    TH1D * hpt_pion = new TH1D ( "hpt_pion" , "Pion p_{T}" , 1000 , 0.0 , 10.0 );
    hpt_pion -> GetXaxis () -> SetTitle ( "p_{T} (GeV/c)" );
    //hpt_pion -> GetYaxis () -> SetTitle ("(#frac{1}{p_{T}})#frac{dN}{dp_{T}} (GeV/c)^{-2}" );
    hpt_pion -> GetXaxis () -> SetTitle ( "# of Tracks" );
    
    TH1D * heta_pion = new TH1D ( "heta_pion" , "Pion #eta" , 1200 , -6.0 , 6.0 );
    heta_pion -> GetXaxis () -> SetTitle ( "Pseudorapidity #eta" );
    heta_pion -> GetYaxis () -> SetTitle ( "# of Tracks" );
    
    TH1D * hy_pion = new TH1D ( "hy_pion" , "Pion y" , 1200 , -6.0 , 6.0 );
    hy_pion -> GetXaxis () -> SetTitle ( "Rapidity y" );
    hy_pion -> GetYaxis () -> SetTitle ( "# of Tracks" );
    
    TH1D * hphi_pion = new TH1D ( "hphi_pion" , "Pion #phi" , 200 , -0.5 * TMath::Pi () , 2.5 * TMath::Pi () );
    hphi_pion -> GetXaxis () -> SetTitle ( "#phi" );
    hphi_pion -> GetYaxis () -> SetTitle ( "# of Tracks" );
    
    // kaon
    TH1D * hpx_kaon = new TH1D ( "hpx_kaon" , "Kaon p_{X}" , 1200 , -6.0 , 6.0 );
    hpx_kaon -> GetXaxis () -> SetTitle ( "p_{X} (GeV/c)" );
    //hpx_kaon -> GetYaxis () -> SetTitle ("(#frac{1}{p_{X}})#frac{dN}{dp_{X}} (GeV/c)^{-2}" );
    hpx_kaon -> GetXaxis () -> SetTitle ( "# of Tracks" );
    
    TH1D * hpy_kaon = new TH1D ( "hpy_kaon" , "Kaon p_{Y}" , 1200 , -6.0 , 6.0 );
    hpy_kaon -> GetXaxis () -> SetTitle ( "p_{Y} (GeV/c)" );
    //hpy_kaon -> GetYaxis () -> SetTitle ("(#frac{1}{p_{Y}})#frac{dN}{dp_{Y}} (GeV/c)^{-2}" );
    hpy_kaon -> GetXaxis () -> SetTitle ( "p_{Y} (GeV/c)" );
    
    TH1D * hpz_kaon = new TH1D ( "hpz_kaon" , "Kaon p_{Z}" , 1200 , -6.0 , 6.0 );
    hpz_kaon -> GetXaxis () -> SetTitle ( "p_{Z} (GeV/c)" );
    //hpz_kaon -> GetYaxis () -> SetTitle ("(#frac{1}{p_{Z}})#frac{dN}{dp_{Z}} (GeV/c)^{-2}" );
    hpz_kaon -> GetXaxis () -> SetTitle ( "# of Tracks" );
    
    TH1D * hpt_kaon = new TH1D ( "hpt_kaon" , "Kaon p_{T}" , 1000 , 0.0 , 10.0);
    hpt_kaon -> GetXaxis () -> SetTitle ( "p_{T} (GeV/c)" );
    //hpt_kaon -> GetYaxis () -> SetTitle ("(#frac{1}{p_{T}})#frac{dN}{dp_{T}} (GeV/c)^{-2}" );
    hpt_kaon -> GetXaxis () -> SetTitle ( "# of Tracks" );
    
    TH1D * heta_kaon = new TH1D ( "heta_kaon" , "Kaon #eta" , 1200 , -6.0 , 6.0 );
    heta_kaon -> GetXaxis () -> SetTitle ( "Pseudorapidity #eta" );
    heta_kaon -> GetYaxis () -> SetTitle ( "# of Tracks" );
    
    TH1D * hy_kaon = new TH1D ( "hy_kaon" , "Kaon y" , 1200 , -6.0 , 6.0 );
    hy_kaon -> GetXaxis () -> SetTitle ( "Rapidity y" );
    hy_kaon -> GetYaxis () -> SetTitle ( "# of Tracks" );
    
    TH1D * hphi_kaon = new TH1D ( "hphi_kaon" , "Kaon #phi" , 200 , -0.5 * TMath::Pi () , 2.5 * TMath::Pi () );
    hphi_kaon -> GetXaxis () -> SetTitle ( "#phi" );
    hphi_kaon -> GetYaxis () -> SetTitle ( "# of Tracks" );*/
    
    // reference flow
    TH1D * hr2 [ Ncentralities ], * hr4 [ Ncentralities ];
    for ( Int_t icent = 0 ; icent < Ncentralities ; icent ++ ) {
        
        sprintf ( name , "hr2_%d" , icent + 1 ); sprintf ( description , "Centrality bin %d QC2 reference flow" , icent + 1 );
        hr2 [ icent ] = new TH1D ( name , description , Ntracks , 0.5 , Ntracks + 0.5 );
        hr2 [ icent ] -> GetXaxis () -> SetTitle ( "Track multiplicity" );
        hr2 [ icent ] -> GetYaxis () -> SetTitle ( "Reference flow" );
        hr2 [ icent ] -> SetMarkerStyle ( 24 ); hr2[ icent ] -> SetMarkerColor ( 2 ); hr2[ icent ] -> SetLineColor ( 2 );
        
        sprintf ( name , "hr4_%d" , icent + 1 ); sprintf ( description , "Centrality bin %d QC4 reference flow" , icent + 1 );
        hr4 [ icent ] = new TH1D ( name , description , Ntracks , 0.5 , Ntracks + 0.5 );
        hr4 [ icent ] -> GetXaxis () -> SetTitle ( "Track multiplicity" );
        hr4 [ icent ] -> GetYaxis () -> SetTitle ( "Reference flow" );
        hr4 [ icent ] -> SetMarkerStyle ( 24 ); hr4[ icent ] -> SetMarkerColor ( 4 ); hr4[ icent ] -> SetLineColor ( 4 );
        
    }
    
    // differential flow
    TH1D * hv1_2 = new TH1D ( "v1_2" , "v_{1}(#eta)" , Nbins , rapidityLow , rapidityHigh );
    hv1_2 -> GetXaxis () -> SetTitle ( "Pseudorapidity #eta" );
    hv1_2 -> GetYaxis () -> SetTitle ( "v_{1}" );
    hv1_2 -> SetMarkerStyle ( 24 ); hv1_2 -> SetMarkerColor ( 2 ); hv1_2 -> SetLineColor ( 2 );
    
    TH1D * hv1_4 = new TH1D ( "v1_4" , "v_{1}(#eta)" , Nbins , rapidityLow , rapidityHigh );
    hv1_4 -> GetXaxis () -> SetTitle ( "Pseudorapidity #eta" );
    hv1_4 -> GetYaxis () -> SetTitle ( "v_{1}" );
    hv1_4 -> SetMarkerStyle ( 24 ); hv1_4 -> SetMarkerColor ( 4 ); hv1_4 -> SetLineColor ( 4 );
    
    TH1D * hv2_2 = new TH1D ( "v2_2" , "v_{2}(p_{T})" , Nbins , rapidityLow , rapidityHigh );
    hv2_2 -> GetXaxis () -> SetTitle ( "p_{T} (GeV/c)" );
    hv2_2 -> GetYaxis () -> SetTitle ( "v_{2}" );
    hv2_2 -> SetMarkerStyle ( 24 ); hv2_2 -> SetMarkerColor ( 2 ); hv2_2 -> SetLineColor ( 2 );
    
    TH1D * hv2_4 = new TH1D ( "v2_4" , "v_{2}(p_{T})" , Nbins , rapidityLow , rapidityHigh );
    hv2_4 -> GetXaxis () -> SetTitle ( "p_{T} (GeV/c)" );
    hv2_4 -> GetYaxis () -> SetTitle ( "v_{2}" );
    hv2_4 -> SetMarkerStyle ( 24 ); hv2_4 -> SetMarkerColor ( 4 ); hv2_4 -> SetLineColor ( 4 );
	
	TH1D * hv3_2 = new TH1D ( "v3_2" , "v_{3}(p_{T})" , Nbins , rapidityLow , rapidityHigh );
    hv3_2 -> GetXaxis () -> SetTitle ( "p_{T} (GeV/c)" );
    hv3_2 -> GetYaxis () -> SetTitle ( "v_{3}" );
    hv3_2 -> SetMarkerStyle ( 24 ); hv3_2 -> SetMarkerColor ( 2 ); hv3_2 -> SetLineColor ( 2 );
    
    TH1D * hv3_4 = new TH1D ( "v3_4" , "v_{3}(p_{T})" , Nbins , rapidityLow , rapidityHigh );
    hv3_4 -> GetXaxis () -> SetTitle ( "p_{T} (GeV/c)" );
    hv3_4 -> GetYaxis () -> SetTitle ( "v_{3}" );
    hv3_4 -> SetMarkerStyle ( 24 ); hv3_4 -> SetMarkerColor ( 4 ); hv3_4 -> SetLineColor ( 4 );
    
    // flow to compare with
    TH1D * hv1_0 = new TH1D ( "v1_0" , "v_{1}(#eta)" , 60 , -3.0 , 3.0 );
    hv1_0 -> GetXaxis () -> SetTitle ( "Pseudorapidity #eta" );
    hv1_0 -> GetYaxis () -> SetTitle ( "v_{1}" );
    hv1_0 -> SetMarkerStyle ( 20 ); hv1_0 -> SetMarkerColor ( 1 ); hv1_0 -> SetLineColor ( 1 );
    
    TH1D * hv1_1 = new TH1D ( "v1_1" , "v_{1}(#eta)" , 60 , -3.0 , 3.0 );
    hv1_1 -> GetXaxis () -> SetTitle ( "Pseudorapidity #eta" );
    hv1_1 -> GetYaxis () -> SetTitle ( "v_{1}" );
    hv1_1 -> SetMarkerStyle ( 20 ); hv1_1 -> SetMarkerColor ( 3 ); hv1_1 -> SetLineColor ( 3 );
    
    TH1D * hv2_0 = new TH1D ( "v2_0" , "v_{2}(p_{T})" , 60 , 0.0 , 6.0 );
    hv2_0 -> GetXaxis () -> SetTitle ( "p_{T} (GeV/c)" );
    hv2_0 -> GetYaxis () -> SetTitle ( "v_{2}" );
    hv2_0 -> SetMarkerStyle ( 20 ); hv2_0 -> SetMarkerColor ( 1 ); hv2_0 -> SetLineColor ( 1 );
    
    TH1D * hv2_1 = new TH1D ( "v2_1" , "v_{2}(p_{T})" , 60 , 0.0 , 6.0 );
    hv2_1 -> GetXaxis () -> SetTitle ( "p_{T} (GeV/c)" );
    hv2_1 -> GetYaxis () -> SetTitle ( "v_{2}" );
    hv2_1 -> SetMarkerStyle ( 20 ); hv2_1 -> SetMarkerColor ( 3 ); hv2_1 -> SetLineColor ( 3 );
	
	TH1D * hv3_0 = new TH1D ( "v3_0" , "v_{3}(p_{T})" , 60 , 0.0 , 6.0 );
    hv3_0 -> GetXaxis () -> SetTitle ( "p_{T} (GeV/c)" );
    hv3_0 -> GetYaxis () -> SetTitle ( "v_{3}" );
    hv3_0 -> SetMarkerStyle ( 20 ); hv3_0 -> SetMarkerColor ( 1 ); hv3_0 -> SetLineColor ( 1 );
    
    // get particle masses
    /*StProton * Aproton = StProton::instance (); StPionMinus * Apion = StPionMinus::instance (); StKaonMinus * Akaon = StKaonMinus::instance ();
     Double_t mass_proton = Aproton -> mass (), mass_pion = Apion -> mass (), mass_kaon = Akaon -> mass ();*/
    
	// random number generator for phi asymmetry simulation
	//TRandom3 *r0 = new TRandom3(0);
	
    // loop over events
    for ( Long64_t ievent = 0 ; ievent < nentries ; ievent ++ ) {
        
        //if ( ievent + 1 > 10 ) break;
        if ( ( ievent + 1 ) % 1000 == 0 ) std::cout << "Analyzed event " << ievent + 1 << std::endl << std::endl; chain -> GetEntry ( ievent );
        
        TLeaf* leaf_runId = (TLeaf*)chain->GetLeaf("runId");
        Int_t runId = (Int_t)leaf_runId->GetValue(0);
        
        TLeaf* leaf_eventId = (TLeaf*)chain->GetLeaf("eventId");
        Int_t eventId = (Int_t)leaf_eventId->GetValue(0);
        
        TLeaf* leaf_bField = (TLeaf*)chain->GetLeaf("bField");
        Float_t bField = (Float_t)leaf_bField->GetValue(0);
        
        TLeaf* leaf_Vx = (TLeaf*)chain->GetLeaf("Vx");
        TLeaf* leaf_Vy = (TLeaf*)chain->GetLeaf("Vy");
        TLeaf* leaf_Vz = (TLeaf*)chain->GetLeaf("Vz");
        //hist_Vz->Fill((Float_t)leaf_Vz->GetValue(0));
        TVector3 V((Float_t)leaf_Vx->GetValue(0),(Float_t)leaf_Vy->GetValue(0),(Float_t)leaf_Vz->GetValue(0));
        
        TLeaf* leaf_centrality = (TLeaf*)chain->GetLeaf("centrality");
        UShort_t centrality = (UShort_t)leaf_centrality->GetValue(0);
        //hist_centralities->Fill(centrality);
        
        TLeaf* leaf_tracknumber = (TLeaf*)chain->GetLeaf("tracknumber");
        UShort_t tracknumber = (UShort_t)leaf_tracknumber->GetValue(0);
        //hist_trackmult->Fill(tracknumber);
        
        TLeaf* leaf_PID = (TLeaf*)chain->GetLeaf("PID");
        TLeaf* leaf_Charge = (TLeaf*)chain->GetLeaf("Charge");
        TLeaf* leaf_Px = (TLeaf*)chain->GetLeaf("Px");
        TLeaf* leaf_Py = (TLeaf*)chain->GetLeaf("Py");
        TLeaf* leaf_Pz = (TLeaf*)chain->GetLeaf("Pz");
        
        TLeaf* leaf_BBCadc = (TLeaf*)chain->GetLeaf("BBCadc");
        
        TLeaf* leaf_nEPDhits = (TLeaf*)chain->GetLeaf("nEPDhits");
        UShort_t nEPDhits = (UShort_t)leaf_nEPDhits->GetValue(0);
        TLeaf* leaf_EPDid = (TLeaf*)chain->GetLeaf("EPDid");
        TLeaf* leaf_EPDnMip = (TLeaf*)chain->GetLeaf("EPDnMip");
        TLeaf* leaf_EPDadc = (TLeaf*)chain->GetLeaf("EPDadc");
        
        // define and initialize event/track-wise variables
        // track # and QC RFPa, RFPb, POI counters
        Int_t Ntrack = 0, N_RFPa = 0, N_RFPb = 0, N_POI = 0;
        
        // define v_n
        Double_t n = 1.0;
        
        // reference variables
        Double_t  Q_n_1_r = 0.0,  Q_n_1_i = 0.0,  S11 = 0.0,    // track-wise variables
                 bQ_n_1_r = 0.0, bQ_n_1_i = 0.0, bS11 = 0.0;
        
        // differential variables
        Double_t p_n_0_r [ Nbins ], p_n_0_i [ Nbins ], mp [ Nbins ];
        
        // initialize QC variables
        for ( Int_t ibin = 0 ; ibin < Nbins ; ibin ++ ) {
            
            p_n_0_r [ ibin ] = 0.0; p_n_0_i [ ibin ] = 0.0; mp [ ibin ] = 0.0;
            
        }
		
        // loop through tracks
        for(UShort_t i = 0;i < tracknumber;i++) {
            
            UShort_t pid = (UShort_t)leaf_PID->GetValue(i);
            Short_t charge = (Short_t)leaf_Charge->GetValue(i);
            Float_t px = (Float_t)leaf_Px->GetValue(i);
            Float_t py = (Float_t)leaf_Py->GetValue(i);
            Float_t pz = (Float_t)leaf_Pz->GetValue(i);
            Double_t pt = TMath::Sqrt(px*px + py*py);
            Double_t p  = TMath::Sqrt(px*px + py*py + pz*pz);
            TVector3 pMom(px,py,pz);
            Double_t eta = pMom.Eta();
            Double_t phi = pMom.Phi();
            // filling histograms
            //hpx->Fill(px/*,1.0/px*/);
            //hpy->Fill(py/*,1.0/py*/);
            //hpz->Fill(pz/*,1.0/pz*/);
            //hpt->Fill(pt/*,1.0/pt*/);
            //heta->Fill(eta);
            //hphi->Fill(phi);
            Double_t energy = 0.0;
            if(pid == 0) energy = TMath::Sqrt(p*p+Mass_Proton*Mass_Proton);
            if(pid == 1) energy = TMath::Sqrt(p*p+Mass_Pion*Mass_Pion);
            if(pid == 2) energy = TMath::Sqrt(p*p+Mass_Kaon*Mass_Kaon);
            Double_t y = 0.5*TMath::Log((energy+pz)/(energy-pz));
            
            // track cuts
            /*if(pt <= 0) continue;
             if(pmag <= 0) continue;*/
            //if(e <= 0) continue;
			
			// phi asymmetry simulation
			//if(phi > 0.75*TMath::Pi() && phi < TMath::Pi()){
			//	if(r0->Rndm() < ratio) continue;
			//}
            
            // count track #
            Ntrack ++;
            
            // define weight for QC
            Double_t w = 1.0;
            //Double_t w = y;
			//Double_t w = pt;
            //Double_t w = ( pt < 2.0 )? pt : 2.0;
            // find Phi bin number
            /*Int_t Phi_bin = -1; Phi_bin = hphi -> FindBin ( phi ) - 34;
            Double_t w = ( Phi_bin >= 0 && Phi_bin < 134 )? weightPhi [ Phi_bin ] : 0.0;*/
            //w *= ( pt < 2.0 )? pt : 2.0;
            
            Bool_t etaZoneLow = kFALSE, etaZoneMiddle = kFALSE, etaZoneHigh = kFALSE;
            
            if(eta < -1.1) etaZoneLow = kTRUE;
            if(eta > -1.0 && eta < -0.6) etaZoneMiddle = kTRUE;
            if(eta > -0.5) etaZoneHigh = kTRUE;
            
			// RFPa (all) for QC
            if ( pid == 1 && charge > 0 && pt < 1.6 && etaZoneLow ) {
                
                N_RFPa ++;
                w = (y-midrapidity > 0.0)? 1.0 : -1.0;
				// define weight for QC
                //Double_t w = 1.0;
                //Double_t w = y;
                //Double_t w = pt;
                
				// collect QC RFP Q vector parameters
                Q_n_1_r += w * TMath::Cos ( n * phi ); Q_n_1_i += w * TMath::Sin ( n * phi ); S11 += TMath::Abs(w);
                
			}
            
            // RFPb (all) for QC to correlate with POI
            if ( pid == 1 && charge > 0 && pt < 1.6 && etaZoneHigh ) {
                
                N_RFPb++;
                w = (y-midrapidity > 0.0)? 1.0 : -1.0;
                // define weight for QC
                //Double_t w = 1.0;
                //Double_t w = y;
                //Double_t w = pt;
                // collect QC RFP Q vector parameters
                
                bQ_n_1_r += w * TMath::Cos ( n * phi ); bQ_n_1_i += w * TMath::Sin ( n * phi ); bS11 += TMath::Abs(w);
                
            }
            
            // POI (all) for QC
            if( pid == 1 && charge > 0 && pt < 1.6 && etaZoneLow ) {
                
                // find bin number
                //Int_t bin = -1; bin = hv1_4 -> FindBin ( eta ) - 1;
				//Int_t bin = -1; bin = hv2_4 -> FindBin (  pt ) - 1;
				Int_t bin = -1; bin = hv3_4 -> FindBin (  y-midrapidity ) - 1;
                // bin cut to prevent segmentation violation
                if ( bin >= 0 && bin < Nbins ) {
                    
                    N_POI ++;
                    
                    // collect QC POI p vector parameters
                    p_n_0_r [ bin ] += TMath::Cos ( n * phi ); p_n_0_i [ bin ] += TMath::Sin ( n * phi ); mp [ bin ] ++;
                    
				}
                
			}
            
		}    // track loop ends
        
		// get real track multiplicity
        //Int_t tracknumber = Ntrack - 1;
        //Int_t tracknumber = N_RFPb - 1;
        //Int_t tracknumber = 0; centrality = 0;
        
		// event cuts continued
        //if ( tracknumber < 0 || tracknumber >= Ntracks ) continue;    // event cut to prevent segmentation violation
        if ( N_RFPa >= 5 && N_RFPb >= 5 && N_POI >= 1 ) {
            // calculate Reference Flows and corrections
            TComplex  Q_n_1 = TComplex (  Q_n_1_r ,  Q_n_1_i );
            TComplex bQ_n_1 = TComplex ( bQ_n_1_r , bQ_n_1_i );
            
            //Double_t    M11 = S21 - S12;
            Double_t    M11 = bS11 * S11;
            
            //Double_t  corr2 = TMath::Power ( TComplex::Abs ( Q_n_1 ) , 2.0 ) - S12;
            Double_t  corr2 = ( Q_n_1 * TComplex::Conjugate ( bQ_n_1 ) ) . Re ();
            
            Double_t addon1 =  Q_n_1 . Re (); //if ( TMath::Abs ( S11 ) <= 1.0e-6 ) addon1 = 0.0;
            Double_t addon2 =  Q_n_1 . Im (); //if ( TMath::Abs ( S11 ) <= 1.0e-6 ) addon2 = 0.0;
            
            Double_t addon3 = bQ_n_1 . Re (); //if ( TMath::Abs ( bS11 ) <= 1.0e-6 ) addon3 = 0.0;
            Double_t addon4 = bQ_n_1 . Im (); //if ( TMath::Abs ( bS11 ) <= 1.0e-6 ) addon4 = 0.0;
            
            // Collect variables for averaging Reference Flows over all events
            /*if( Q_n_1_r != 0.0 && Q_n_1_i != 0.0 && S11 != 0.0
               && bQ_n_1_r != 0.0 && bQ_n_1_i != 0.0 && bS11 != 0.0
               
               && M11 != 0.0
               
               && corr2 != 0.0
               
               && addon1 != 0.0 && addon2 != 0.0
               
               && addon3 != 0.0 && addon4 != 0.0
               )*/ {
                
                // count # of events for this centrality this track multiplicity
                multCounter [ tracknumber ][ centrality ] ++;
                
                // collect weights and variables
                if (  S11 != 0.0 )   sumS11 [ tracknumber ][ centrality ] += (  S11 );
                if ( bS11 != 0.0 )  sumbS11 [ tracknumber ][ centrality ] += ( bS11 );
                if (  M11 != 0.0 )   sumM11 [ tracknumber ][ centrality ] += (  M11 );
                
                if (  M11 != 0.0 ) sumCorr2 [ tracknumber ][ centrality ] += corr2 * M11 / ( M11 );
                
                // collect variables for statistical errors
                if (  S11 != 0.0 )   sumS11Squared [ tracknumber ][ centrality ] += TMath::Power (  S11 , 2.0 );
                if ( bS11 != 0.0 )  sumbS11Squared [ tracknumber ][ centrality ] += TMath::Power ( bS11 , 2.0 );
                if (  M11 != 0.0 )   sumM11Squared [ tracknumber ][ centrality ] += TMath::Power (  M11 , 2.0 );
                
                if (  M11 != 0.0 ) sumCorr2Squared [ tracknumber ][ centrality ] += TMath::Power ( corr2 , 2.0 ) / ( M11 );
                
                // detector inefficiency correction
                if (  S11 != 0.0 ) sumAddon1 [ tracknumber ][ centrality ] += addon1 *  S11 / (  S11 );
                if (  S11 != 0.0 ) sumAddon2 [ tracknumber ][ centrality ] += addon2 *  S11 / (  S11 );
                if ( bS11 != 0.0 ) sumAddon3 [ tracknumber ][ centrality ] += addon3 * bS11 / ( bS11 );
                if ( bS11 != 0.0 ) sumAddon4 [ tracknumber ][ centrality ] += addon4 * bS11 / ( bS11 );
                
                if (  S11 != 0.0 ) sumAddon1Squared [ tracknumber ][ centrality ] += TMath::Power ( addon1 , 2.0 ) / (  S11 );
                if (  S11 != 0.0 ) sumAddon2Squared [ tracknumber ][ centrality ] += TMath::Power ( addon2 , 2.0 ) / (  S11 );
                if ( bS11 != 0.0 ) sumAddon3Squared [ tracknumber ][ centrality ] += TMath::Power ( addon3 , 2.0 ) / ( bS11 );
                if ( bS11 != 0.0 ) sumAddon4Squared [ tracknumber ][ centrality ] += TMath::Power ( addon4 , 2.0 ) / ( bS11 );
                
                // covariances
                if ( S11 != 0.0 && bS11 != 0.0 ) sumS11bS11 [ tracknumber ][ centrality ] += ( S11 * bS11 );
                if ( M11 != 0.0 &&  S11 != 0.0 )  sumM11S11 [ tracknumber ][ centrality ] += ( M11 *  S11 );
                if ( M11 != 0.0 && bS11 != 0.0 ) sumM11bS11 [ tracknumber ][ centrality ] += ( M11 * bS11 );
                
                if ( M11 != 0.0 &&  S11 != 0.0 ) sumCorr2C1 [ tracknumber ][ centrality ] += ( corr2 * M11 / ( M11 ) ) * ( addon1 *  S11 / (  S11 ) );
                if ( M11 != 0.0 &&  S11 != 0.0 ) sumCorr2C2 [ tracknumber ][ centrality ] += ( corr2 * M11 / ( M11 ) ) * ( addon2 *  S11 / (  S11 ) );
                if ( M11 != 0.0 && bS11 != 0.0 ) sumCorr2C3 [ tracknumber ][ centrality ] += ( corr2 * M11 / ( M11 ) ) * ( addon3 * bS11 / ( bS11 ) );
                if ( M11 != 0.0 && bS11 != 0.0 ) sumCorr2C4 [ tracknumber ][ centrality ] += ( corr2 * M11 / ( M11 ) ) * ( addon4 * bS11 / ( bS11 ) );
                
                if ( S11 != 0.0 )                   sumC1C2 [ tracknumber ][ centrality ] += ( addon1 *  S11 / (  S11 ) ) * ( addon2 *  S11 / (  S11 ) );
                if ( S11 != 0.0 && bS11 != 0.0 )    sumC1C3 [ tracknumber ][ centrality ] += ( addon1 *  S11 / (  S11 ) ) * ( addon3 * bS11 / ( bS11 ) );
                if ( S11 != 0.0 && bS11 != 0.0 )    sumC1C4 [ tracknumber ][ centrality ] += ( addon1 *  S11 / (  S11 ) ) * ( addon4 * bS11 / ( bS11 ) );
                
                if ( S11 != 0.0 && bS11 != 0.0 )    sumC2C3 [ tracknumber ][ centrality ] += ( addon2 *  S11 / (  S11 ) ) * ( addon3 * bS11 / ( bS11 ) );
                if ( S11 != 0.0 && bS11 != 0.0 )    sumC2C4 [ tracknumber ][ centrality ] += ( addon2 *  S11 / (  S11 ) ) * ( addon4 * bS11 / ( bS11 ) );
                
                if ( bS11 != 0.0 )                  sumC3C4 [ tracknumber ][ centrality ] += ( addon3 * bS11 / ( bS11 ) ) * ( addon4 * bS11 / ( bS11 ) );
            }
            
            // compute Differential Flows and corrections
            // loop through y/pT bins
            TComplex p_n_0;
            for ( Int_t ibin = 0 ; ibin < Nbins ; ibin ++ ) {
                
                p_n_0 = TComplex ( p_n_0_r [ ibin ] , p_n_0_i [ ibin ] );
                
                //Double_t    M01 = mp [ ibin ] * S11 - s11 [ ibin ];
                Double_t    M01 = mp [ ibin ] * bS11;
                
                //Double_t corr2Reduced = ( p_n_0 * TComplex::Conjugate ( Q_n_1 ) - s11 [ ibin ] ) . Re ();
                Double_t corr2Reduced = ( p_n_0 * TComplex::Conjugate ( bQ_n_1 ) ) . Re ();
                
                Double_t  addon5 = p_n_0 . Re ();
                Double_t  addon6 = p_n_0 . Im ();
                
                // Collect variables for averaging Differential Flows over all events
                /*if ( Q_n_1_r != 0.0 && Q_n_1_i != 0.0 && S11 != 0.0
                   
                   && bQ_n_1_r != 0.0 && bQ_n_1_i != 0.0 && bS11 != 0.0
                   
                   && M11 != 0.0
                   
                   && corr2 != 0.0
                   
                   && addon1 != 0.0 && addon2 != 0.0
                   && addon3 != 0.0 && addon4 != 0.0
                   
                   && p_n_0_r [ ibin ] != 0.0 && p_n_0_i [ ibin ] != 0.0 && mp [ ibin ] != 0.0
                   
                   && M01 != 0.0
                   
                   && corr2Reduced != 0.0
                   
                   && addon5 != 0.0 && addon6 != 0.0
                   )*/ {
                    
                    // count # of tracks for this centrality this track multiplicity this bin #
                    binMultCounter [ ibin ][ tracknumber ][ centrality ]++;
                    
                    // collect weights and variables
                    if ( mp [ ibin ] != 0.0 )      sumMp [ ibin ][ tracknumber ][ centrality ] += ( mp [ ibin ] );
                    if ( M01 != 0.0 )             sumM01 [ ibin ][ tracknumber ][ centrality ] += ( M01         );
                    
                    if ( M01 != 0.0 )    sumCorr2Reduced [ ibin ][ tracknumber ][ centrality ] += corr2Reduced * M01 / ( M01 );
                    
                    // Collect variables for statistical errors
                    if ( mp [ ibin ] != 0.0 )      sumMpSquared [ ibin ][ tracknumber ][ centrality ] += TMath::Power ( mp [ ibin ] , 2.0 );
                    if ( M01 != 0.0 )             sumM01Squared [ ibin ][ tracknumber ][ centrality ] += TMath::Power ( M01         , 2.0 );
                    
                    if ( M01 != 0.0 ) sumCorr2ReducedSquared [ ibin ][ tracknumber ][ centrality ] += TMath::Power ( corr2Reduced , 2.0 ) / ( M01 );
                    
                    // detector inefficiency correction
                    if ( mp [ ibin ] != 0.0 ) sumAddon5 [ ibin ][ tracknumber ][ centrality ] += addon5 * mp [ ibin ] / ( mp [ ibin ] );
                    if ( mp [ ibin ] != 0.0 ) sumAddon6 [ ibin ][ tracknumber ][ centrality ] += addon6 * mp [ ibin ] / ( mp [ ibin ] );
                    
                    if ( mp [ ibin ] != 0.0 ) sumAddon5Squared [ ibin ][ tracknumber ][ centrality ] += TMath::Power ( addon5 , 2.0 ) / ( mp [ ibin ] );
                    if ( mp [ ibin ] != 0.0 ) sumAddon6Squared [ ibin ][ tracknumber ][ centrality ] += TMath::Power ( addon6 , 2.0 ) / ( mp [ ibin ] );
                    
                    // Covariances
                    if (  M11 != 0.0 && M01 != 0.0         )  sumM11M01 [ ibin ][ tracknumber ][ centrality ] += (  M11 * M01         );
                    if (  M11 != 0.0 && mp [ ibin ] != 0.0 )   sumM11mp [ ibin ][ tracknumber ][ centrality ] += (  M11 * mp [ ibin ] );
                    if (  S11 != 0.0 && M01 != 0.0         )  sumS11M01 [ ibin ][ tracknumber ][ centrality ] += (  S11 * M01         );
                    if (  S11 != 0.0 && mp [ ibin ] != 0.0 )   sumS11mp [ ibin ][ tracknumber ][ centrality ] += (  S11 * mp [ ibin ] );
                    if ( bS11 != 0.0 && M01 != 0.0         ) sumbS11M01 [ ibin ][ tracknumber ][ centrality ] += ( bS11 * M01         );
                    if ( bS11 != 0.0 && mp [ ibin ] != 0.0 )  sumbS11mp [ ibin ][ tracknumber ][ centrality ] += ( bS11 * mp [ ibin ] );
                    if (  M01 != 0.0 && mp [ ibin ] != 0.0 )   sumM01mp [ ibin ][ tracknumber ][ centrality ] += (  M01 * mp [ ibin ] );
                    
                    if (  M11 != 0.0 && M01 != 0.0         ) sumCorr2Corr2Reduced [ ibin ][ tracknumber ][ centrality ] += ( corr2 * M11 / ( M11 ) )
                                                                                                                           * ( corr2Reduced * M01 / ( M01 ) );
                    if (  M11 != 0.0 && mp [ ibin ] != 0.0 )           sumCorr2C5 [ ibin ][ tracknumber ][ centrality ] += ( corr2 * M11 / ( M11 ) )
                                                                                                                           * ( addon5 * mp [ ibin ] / ( mp [ ibin ] ) );
                    if (  M11 != 0.0 && mp [ ibin ] != 0.0 )           sumCorr2C6 [ ibin ][ tracknumber ][ centrality ] += ( corr2 * M11 / ( M11 ) )
                                                                                                                           * ( addon6 * mp [ ibin ] / ( mp [ ibin ] ) );
                    
                    if (  M01 != 0.0 &&      S11 != 0.0 )    sumCorr2ReducedC1 [ ibin ][ tracknumber ][ centrality ] += ( corr2Reduced * M01 / ( M01 ) )
                                                                                                                        * ( addon1 *  S11 / (  S11 ) );
                    if (  M01 != 0.0 &&      S11 != 0.0 )    sumCorr2ReducedC2 [ ibin ][ tracknumber ][ centrality ] += ( corr2Reduced * M01 / ( M01 ) )
                                                                                                                        * ( addon2 *  S11 / (  S11 ) );
                    if (  M01 != 0.0 &&     bS11 != 0.0 )    sumCorr2ReducedC3 [ ibin ][ tracknumber ][ centrality ] += ( corr2Reduced * M01 / ( M01 ) )
                                                                                                                        * ( addon3 * bS11 / ( bS11 ) );
                    if (  M01 != 0.0 &&     bS11 != 0.0 )    sumCorr2ReducedC4 [ ibin ][ tracknumber ][ centrality ] += ( corr2Reduced * M01 / ( M01 ) )
                                                                                                                        * ( addon4 * bS11 / ( bS11 ) );
                    if (  M01 != 0.0 && mp [ ibin ] != 0.0 )    sumCorr2ReducedC5 [ ibin ][ tracknumber ][ centrality ] += ( corr2Reduced * M01 / ( M01 ) )
                                                                                                                           * ( addon5 * mp [ ibin ] / ( mp [ ibin ] ) );
                    if (  M01 != 0.0 && mp [ ibin ] != 0.0 )    sumCorr2ReducedC6 [ ibin ][ tracknumber ][ centrality ] += ( corr2Reduced * M01 / ( M01 ) )
                                                                                                                           * ( addon6 * mp [ ibin ] / ( mp [ ibin ] ) );
                    
                    if (  S11 != 0.0 && mp [ ibin ] != 0.0 )              sumC1C5 [ ibin ][ tracknumber ][ centrality ] += ( addon1 *  S11 / (  S11 ) )
                                                                                                                           * ( addon5 * mp [ ibin ] / ( mp [ ibin ] ) );
                    if (  S11 != 0.0 && mp [ ibin ] != 0.0 )              sumC1C6 [ ibin ][ tracknumber ][ centrality ] += ( addon1 *  S11 / (  S11 ) )
                                                                                                                           * ( addon6 * mp [ ibin ] / ( mp [ ibin ] ) );
                    
                    if (  S11 != 0.0 && mp [ ibin ] != 0.0 )              sumC2C5 [ ibin ][ tracknumber ][ centrality ] += ( addon2 *  S11 / (  S11 ) )
                                                                                                                           * ( addon5 * mp [ ibin ] / ( mp [ ibin ] ) );
                    if (  S11 != 0.0 && mp [ ibin ] != 0.0 )              sumC2C6 [ ibin ][ tracknumber ][ centrality ] += ( addon2 *  S11 / (  S11 ) )
                                                                                                                           * ( addon6 * mp [ ibin ] / ( mp [ ibin ] ) );
                    
                    if ( bS11 != 0.0 && mp [ ibin ] != 0.0 )              sumC3C5 [ ibin ][ tracknumber ][ centrality ] += ( addon3 * bS11 / ( bS11 ) )
                                                                                                                           * ( addon5 * mp [ ibin ] / ( mp [ ibin ] ) );
                    if ( bS11 != 0.0 && mp [ ibin ] != 0.0 )              sumC3C6 [ ibin ][ tracknumber ][ centrality ] += ( addon3 * bS11 / ( bS11 ) )
                                                                                                                           * ( addon6 * mp [ ibin ] / ( mp [ ibin ] ) );
                    
                    if ( bS11 != 0.0 && mp [ ibin ] != 0.0 )              sumC4C5 [ ibin ][ tracknumber ][ centrality ] += ( addon4 * bS11 / ( bS11 ) )
                                                                                                                           * ( addon5 * mp [ ibin ] / ( mp [ ibin ] ) );
                    if ( bS11 != 0.0 && mp [ ibin ] != 0.0 )              sumC4C6 [ ibin ][ tracknumber ][ centrality ] += ( addon4 * bS11 / ( bS11 ) )
                                                                                                                           * ( addon6 * mp [ ibin ] / ( mp [ ibin ] ) );
                    
                    if ( mp [ ibin ] != 0.0 && mp [ ibin ] != 0.0 )          sumC5C6 [ ibin ][ tracknumber ][ centrality ] += ( addon5 * mp [ ibin ] / ( mp [ ibin ] ) )
                                                                                                                              * ( addon6 * mp [ ibin ] / ( mp [ ibin ] ) );
                    
                }
                
            }    // y/pT bins loop ends
        }
        
        // fill event-wise histograms
        //hrefmult -> Fill ( numberOfInputTracks ); hrealmult -> Fill ( N_RFPb ); /*himpact -> Fill ( b );*/ hpsi -> Fill ( reaction_plane );
        
        // count number of events for current centrality and track multiplicity
        //multCounter [ tracknumber ][ centrality ] ++;
        
    }    // event loop ends
    
    // save QC variables
    /*Int_t line;
    for ( Int_t icent = 0 ; icent < Ncentralities ; icent ++ ) {
        
        for ( Int_t itrack = 0 ; itrack < Ntracks ; itrack ++ ) {
            
            line = 0;
            
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , multCounter [ itrack ][ icent] ); line ++;
            
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumS11 [ itrack ][ icent] ); line ++;
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumbS11 [ itrack ][ icent] ); line ++;
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumM11 [ itrack ][ icent] ); line ++;
            
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumCorr2 [ itrack ][ icent] ); line ++;
            
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumAddon1 [ itrack ][ icent] ); line ++;
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumAddon2 [ itrack ][ icent] ); line ++;
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumAddon3 [ itrack ][ icent] ); line ++;
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumAddon4 [ itrack ][ icent] ); line ++;
            
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumS11Squared [ itrack ][ icent] ); line ++;
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumbS11Squared [ itrack ][ icent] ); line ++;
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumM11Squared [ itrack ][ icent] ); line ++;
            
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumCorr2Squared [ itrack ][ icent] ); line ++;
            
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumAddon1Squared [ itrack ][ icent] ); line ++;
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumAddon2Squared [ itrack ][ icent] ); line ++;
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumAddon3Squared [ itrack ][ icent] ); line ++;
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumAddon4Squared [ itrack ][ icent] ); line ++;
            
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumM11S11 [ itrack ][ icent] ); line ++;
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumM11bS11 [ itrack ][ icent] ); line ++;
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumS11bS11 [ itrack ][ icent] ); line ++;
            
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumCorr2C1 [ itrack ][ icent] ); line ++;
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumCorr2C2 [ itrack ][ icent] ); line ++;
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumCorr2C3 [ itrack ][ icent] ); line ++;
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumCorr2C4 [ itrack ][ icent] ); line ++;
            
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumC1C2 [ itrack ][ icent] ); line ++;
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumC1C3 [ itrack ][ icent] ); line ++;
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumC1C4 [ itrack ][ icent] ); line ++;
            
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumC2C3 [ itrack ][ icent] ); line ++;
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumC2C4 [ itrack ][ icent] ); line ++;
            
            QCvariables2D [ line ] -> SetBinContent ( itrack + 1 , icent + 1 , sumC3C4 [ itrack ][ icent] ); line ++;
            
            for ( Int_t ibin = 0 ; ibin < Nbins ; ibin ++ ) {
                
                line = 0;
                
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , binMultCounter [ ibin ][ itrack ][ icent] ); line ++;
                
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumMp [ ibin ][ itrack ][ icent] ); line ++;
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumM01 [ ibin ][ itrack ][ icent] ); line ++;
                
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumCorr2Reduced [ ibin ][ itrack ][ icent] ); line ++;
                
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumAddon5 [ ibin ][ itrack ][ icent] ); line ++;
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumAddon6 [ ibin ][ itrack ][ icent] ); line ++;
                
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumMpSquared [ ibin ][ itrack ][ icent] ); line ++;
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumM01Squared [ ibin ][ itrack ][ icent] ); line ++;
                
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumCorr2ReducedSquared [ ibin ][ itrack ][ icent] ); line ++;
                
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumAddon5Squared [ ibin ][ itrack ][ icent] ); line ++;
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumAddon6Squared [ ibin ][ itrack ][ icent] ); line ++;
                
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumM11M01 [ ibin ][ itrack ][ icent] ); line ++;
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumM11mp [ ibin ][ itrack ][ icent] ); line ++;
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumS11M01 [ ibin ][ itrack ][ icent] ); line ++;
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumS11mp [ ibin ][ itrack ][ icent] ); line ++;
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumbS11M01 [ ibin ][ itrack ][ icent] ); line ++;
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumbS11mp [ ibin ][ itrack ][ icent] ); line ++;
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumM01mp [ ibin ][ itrack ][ icent] ); line ++;
                
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumCorr2Corr2Reduced [ ibin ][ itrack ][ icent] ); line ++;
                
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumCorr2C5 [ ibin ][ itrack ][ icent] ); line ++;
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumCorr2C6 [ ibin ][ itrack ][ icent] ); line ++;
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumCorr2ReducedC1 [ ibin ][ itrack ][ icent] ); line ++;
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumCorr2ReducedC2 [ ibin ][ itrack ][ icent] ); line ++;
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumCorr2ReducedC3 [ ibin ][ itrack ][ icent] ); line ++;
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumCorr2ReducedC4 [ ibin ][ itrack ][ icent] ); line ++;
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumCorr2ReducedC5 [ ibin ][ itrack ][ icent] ); line ++;
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumCorr2ReducedC6 [ ibin ][ itrack ][ icent] ); line ++;
                
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumC1C5 [ ibin ][ itrack ][ icent] ); line ++;
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumC1C6 [ ibin ][ itrack ][ icent] ); line ++;
                
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumC2C5 [ ibin ][ itrack ][ icent] ); line ++;
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumC2C6 [ ibin ][ itrack ][ icent] ); line ++;
                
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumC3C5 [ ibin ][ itrack ][ icent] ); line ++;
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumC3C6 [ ibin ][ itrack ][ icent] ); line ++;
                
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumC4C5 [ ibin ][ itrack ][ icent] ); line ++;
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumC4C6 [ ibin ][ itrack ][ icent] ); line ++;
                
                QCvariables3D [ line ] -> SetBinContent ( ibin + 1 , itrack + 1 , icent + 1 , sumC5C6 [ ibin ][ itrack ][ icent] ); line ++;
                
            }    // bin loop ends
            
        }    // track loop ends
        
    }*/    // centrality loop ends
    
    std::cout << "Computing flows" << std::endl << std::endl;
    
    for ( Int_t icent = 0 ; icent < Ncentralities ; icent ++ ) {
        
        for ( Int_t itrack = 0 ; itrack < Ntracks ; itrack ++ ) {
            
            // get mean
			Double_t reference2_1 = (  sumM11 [ itrack ][ icent ] > 0.0 &&  sumCorr2 [ itrack ][ icent ] != 0.0 )?  sumCorr2 [ itrack ][ icent ] /  sumM11 [ itrack ][ icent ] : -999.0;
            //if(!TMath::Finite(reference2_1)) reference2_1 = -999.0;
            
			Double_t  correction1 = (  sumS11 [ itrack ][ icent ] > 0.0 && sumAddon1 [ itrack ][ icent ] != 0.0 )? sumAddon1 [ itrack ][ icent ] /  sumS11 [ itrack ][ icent ] : -999.0;
            //if(!TMath::Finite( correction1))  correction1 = -999.0;
			Double_t  correction2 = (  sumS11 [ itrack ][ icent ] > 0.0 && sumAddon2 [ itrack ][ icent ] != 0.0 )? sumAddon2 [ itrack ][ icent ] /  sumS11 [ itrack ][ icent ] : -999.0;
            //if(!TMath::Finite( correction2))  correction2 = -999.0;
            Double_t  correction3 = ( sumbS11 [ itrack ][ icent ] > 0.0 && sumAddon3 [ itrack ][ icent ] != 0.0 )? sumAddon3 [ itrack ][ icent ] / sumbS11 [ itrack ][ icent ] : -999.0;
            //if(!TMath::Finite( correction3))  correction3 = -999.0;
            Double_t  correction4 = ( sumbS11 [ itrack ][ icent ] > 0.0 && sumAddon4 [ itrack ][ icent ] != 0.0 )? sumAddon4 [ itrack ][ icent ] / sumbS11 [ itrack ][ icent ] : -999.0;
            //if(!TMath::Finite( correction4))  correction4 = -999.0;
            
			// get error on mean
            Double_t reference2Error1 = GetError (  sumM11 [ itrack ][ icent ] ,  sumM11Squared [ itrack ][ icent ] ,  sumCorr2 [ itrack ][ icent ] ,  sumCorr2Squared [ itrack ][ icent ] );
            
            Double_t correction1Error = GetError (  sumS11 [ itrack ][ icent ] ,  sumS11Squared [ itrack ][ icent ] , sumAddon1 [ itrack ][ icent ] , sumAddon1Squared [ itrack ][ icent ] );
            Double_t correction2Error = GetError (  sumS11 [ itrack ][ icent ] ,  sumS11Squared [ itrack ][ icent ] , sumAddon2 [ itrack ][ icent ] , sumAddon2Squared [ itrack ][ icent ] );
            Double_t correction3Error = GetError ( sumbS11 [ itrack ][ icent ] , sumbS11Squared [ itrack ][ icent ] , sumAddon3 [ itrack ][ icent ] , sumAddon3Squared [ itrack ][ icent ] );
            Double_t correction4Error = GetError ( sumbS11 [ itrack ][ icent ] , sumbS11Squared [ itrack ][ icent ] , sumAddon4 [ itrack ][ icent ] , sumAddon4Squared [ itrack ][ icent ] );
            
			// zero bad averages and errors
            if ( reference2_1 == -999.0 || reference2Error1  == -999.0 || reference2_1 * reference2Error1 == 0.0 ) { reference2_1 = -999.0; reference2Error1 = -999.0; }
            
            if (  correction1 == -999.0 || correction1Error  == -999.0 ||  correction1 * correction1Error == 0.0 ) {  correction1 = -999.0; correction1Error = -999.0; }
            if (  correction2 == -999.0 || correction2Error  == -999.0 ||  correction2 * correction2Error == 0.0 ) {  correction2 = -999.0; correction2Error = -999.0; }
            if (  correction3 == -999.0 || correction3Error  == -999.0 ||  correction3 * correction3Error == 0.0 ) {  correction3 = -999.0; correction3Error = -999.0; }
            if (  correction4 == -999.0 || correction4Error  == -999.0 ||  correction4 * correction4Error == 0.0 ) {  correction4 = -999.0; correction4Error = -999.0; }
            // Covariance terms
            Double_t covRef2C1 = GetCovariance (   sumM11 [ itrack ][ icent ] ,   sumCorr2 [ itrack ][ icent ] ,
                                                   sumS11 [ itrack ][ icent ] ,  sumAddon1 [ itrack ][ icent ] ,
                                                sumM11S11 [ itrack ][ icent ] , sumCorr2C1 [ itrack ][ icent ] );
            Double_t covRef2C2 = GetCovariance (   sumM11 [ itrack ][ icent ] ,   sumCorr2 [ itrack ][ icent ] ,
                                                   sumS11 [ itrack ][ icent ] ,  sumAddon2 [ itrack ][ icent ] ,
                                                sumM11S11 [ itrack ][ icent ] , sumCorr2C2 [ itrack ][ icent ] );
            Double_t covRef2C3 = GetCovariance (    sumM11 [ itrack ][ icent ] ,   sumCorr2 [ itrack ][ icent ] ,
                                                   sumbS11 [ itrack ][ icent ] ,  sumAddon3 [ itrack ][ icent ] ,
                                                sumM11bS11 [ itrack ][ icent ] , sumCorr2C3 [ itrack ][ icent ] );
            Double_t covRef2C4 = GetCovariance (    sumM11 [ itrack ][ icent ] ,   sumCorr2 [ itrack ][ icent ] ,
                                                   sumbS11 [ itrack ][ icent ] ,  sumAddon4 [ itrack ][ icent ] ,
                                                sumM11bS11 [ itrack ][ icent ] , sumCorr2C4 [ itrack ][ icent ] );

            Double_t covC1C2 = GetCovariance (       sumS11 [ itrack ][ icent ] , sumAddon1 [ itrack ][ icent ] ,
                                                     sumS11 [ itrack ][ icent ] , sumAddon2 [ itrack ][ icent ] ,
                                              sumS11Squared [ itrack ][ icent ] ,   sumC1C2 [ itrack ][ icent ] );
            Double_t covC1C3 = GetCovariance (    sumS11 [ itrack ][ icent ] , sumAddon1 [ itrack ][ icent ] ,
                                                 sumbS11 [ itrack ][ icent ] , sumAddon3 [ itrack ][ icent ] ,
                                              sumS11bS11 [ itrack ][ icent ] ,   sumC1C3 [ itrack ][ icent ] );
            Double_t covC1C4 = GetCovariance (    sumS11 [ itrack ][ icent ] , sumAddon1 [ itrack ][ icent ] ,
                                                 sumbS11 [ itrack ][ icent ] , sumAddon4 [ itrack ][ icent ] ,
                                              sumS11bS11 [ itrack ][ icent ] ,   sumC1C4 [ itrack ][ icent ] );
            
            Double_t covC2C3 = GetCovariance (    sumS11 [ itrack ][ icent ] , sumAddon2 [ itrack ][ icent ] ,
                                                 sumbS11 [ itrack ][ icent ] , sumAddon3 [ itrack ][ icent ] ,
                                              sumS11bS11 [ itrack ][ icent ] ,   sumC2C3 [ itrack ][ icent ] );
            Double_t covC2C4 = GetCovariance (    sumS11 [ itrack ][ icent ] , sumAddon2 [ itrack ][ icent ] ,
                                                 sumbS11 [ itrack ][ icent ] , sumAddon4 [ itrack ][ icent ] ,
                                              sumS11bS11 [ itrack ][ icent ] ,   sumC2C4 [ itrack ][ icent ] );
            
            Double_t covC3C4 = GetCovariance (       sumbS11 [ itrack ][ icent ] , sumAddon3 [ itrack ][ icent ] ,
                                                     sumbS11 [ itrack ][ icent ] , sumAddon4 [ itrack ][ icent ] ,
                                              sumbS11Squared [ itrack ][ icent ] ,   sumC3C4 [ itrack ][ icent ] );
            
            // Corrected Reference Flows
            Double_t reference2_2 = -999.0;
            if( reference2_1 != -999.0 && correction1 != -999.0 && correction2 != -999.0 && correction3 != -999.0 && correction4 != -999.0 ) {
                reference2_2 = reference2_1 - correction1 * correction3 - correction2 * correction4;
            }
            //Double_t reference2_2 = ( reference2_1 != -999.0 && correction1 != -999.0 && correction2 != -999.0 && correction3 != -999.0 && correction4 != -999.0 )?
            //                        reference2_1 - correction1 * correction3 - correction2 * correction4 : -999.0;
            /*if(!TMath::Finite(reference2_2) || reference2_1 == -999.0
               || correction1 == -999.0 || correction2 == -999.0 || correction3 == -999.0 || correction4 == -999.0) reference2_2 = -999.0;*/
            
			// prepare reference flow to be divided
            //Double_t sign2 = ( reference2_2 > 0.0 )? 1.0 : -1.0;
            Double_t sign2 = 1.0;
            
            Double_t reference2_3 = -999.0;
            if( reference2_2 != -999.0 && sign2 * reference2_2 > 0.0 ) {
                reference2_3 = TMath::Sqrt ( sign2 * reference2_2 );
            }
            //if ( /*!TMath::Finite ( reference2_3 ) || reference2_2 == -999.0 ||*/ TMath::Abs ( reference2_3 ) < 1.0e-44 ) reference2_3 = -999.0;
            
			// differentiate QC2 reference flow
            Double_t dr23r2 = -999.0;
            if( reference2_3 != -999.0 ) {
                dr23r2 = 0.5 * TMath::Power ( sign2 * reference2_2 , -0.5 ) * sign2 * 1.0;
            }
            //if(!TMath::Finite(dr23r2)) dr23r2 = -999.0;
            Double_t dr23c1 = -999.0;
            if( reference2_3 != -999.0 ) {
                dr23c1 = 0.5 * TMath::Power ( sign2 * reference2_2 , -0.5 ) * sign2 * (-correction3);
            }
            //if(!TMath::Finite(dr23c1)) dr23c1 = -999.0;
            Double_t dr23c2 = -999.0;
            if( reference2_3 != -999.0 ) {
                dr23c2 = 0.5 * TMath::Power ( sign2 * reference2_2 , -0.5 ) * sign2 * (-correction4);
            }
            //if(!TMath::Finite(dr23c2)) dr23c2 = -999.0;
            Double_t dr23c3 = -999.0;
            if( reference2_3 != -999.0 ) {
                dr23c3 = 0.5 * TMath::Power ( sign2 * reference2_2 , -0.5 ) * sign2 * (-correction1);
            }
            //if(!TMath::Finite(dr23c3)) dr23c3 = -999.0;
            Double_t dr23c4 = -999.0;
            if( reference2_3 != -999.0 ) {
                dr23c4 = 0.5 * TMath::Power ( sign2 * reference2_2 , -0.5 ) * sign2 * (-correction2);
            }
            //if(!TMath::Finite(dr23c4)) dr23c4 = -999.0;
            
			// error propagation
			Double_t reference2Error3 = -999.0;
			if ( dr23r2 != -999.0 && reference2Error1 != -999.0
                
                && dr23c1 != -999.0 && correction1Error != -999.0
                && dr23c2 != -999.0 && correction2Error != -999.0
                && dr23c3 != -999.0 && correction3Error != -999.0
                && dr23c4 != -999.0 && correction4Error != -999.0
                
                && covRef2C1 != -999.0
                && covRef2C2 != -999.0
                && covRef2C3 != -999.0
                && covRef2C4 != -999.0
                
                && covC1C2 != -999.0
                && covC1C3 != -999.0
                && covC1C4 != -999.0
                
                && covC2C3 != -999.0
                && covC2C4 != -999.0
                
                && covC3C4 != -999.0
                ) {
                
				reference2Error3 = TMath::Sqrt( TMath::Power ( dr23r2 , 2.0 ) * reference2Error1
                                               
                                               + TMath::Power ( dr23c1 , 2.0 ) * correction1Error
                                               + TMath::Power ( dr23c2 , 2.0 ) * correction2Error
                                               + TMath::Power ( dr23c3 , 2.0 ) * correction3Error
                                               + TMath::Power ( dr23c4 , 2.0 ) * correction4Error
                                               
                                               + 2.0 * dr23r2 * dr23c1 * covRef2C1
                                               + 2.0 * dr23r2 * dr23c2 * covRef2C2
                                               + 2.0 * dr23r2 * dr23c3 * covRef2C3
                                               + 2.0 * dr23r2 * dr23c4 * covRef2C4
                                               
                                               + 2.0 * dr23c1 * dr23c2 * covC1C2
                                               + 2.0 * dr23c1 * dr23c3 * covC1C3
                                               + 2.0 * dr23c1 * dr23c4 * covC1C4
                                               
                                               + 2.0 * dr23c2 * dr23c3 * covC2C3
                                               + 2.0 * dr23c2 * dr23c4 * covC2C4
                                               
                                               + 2.0 * dr23c3 * dr23c4 * covC3C4
                                               );
                
                if ( !TMath::Finite ( reference2Error3 ) ) reference2Error3 = -999.0;

			}
            
            if ( reference2_3 == -999.0 || reference2Error3 == -999.0 || reference2_3 * reference2Error3 == 0.0
                //|| TMath::Abs ( reference2_3 ) <= reference2Error3 
				) { reference2_3 = -999.0; reference2Error3 = -999.0; }
            
			// QC2 reference flow
			if ( reference2_3 != -999.0 && reference2Error3 != -999.0
                //TMath::Finite(reference2_3) && TMath::Finite(reference2Error3)
                ) {
                hr2 [ icent ] -> SetBinContent ( itrack + 1 , reference2_3 ); hr2 [ icent ] -> SetBinError ( itrack + 1 , reference2Error3 );
            }
            
            // compute differential flows
            for ( Int_t ibin = 0 ; ibin < Nbins ; ibin ++ ) {
                
                // Differential Flows and their correction terms
				Double_t differential2_1 = ( sumM01 [ ibin ][ itrack ][ icent ] > 0.0 && sumCorr2Reduced [ ibin ][ itrack ][ icent ] != 0.0 )?
                                           sumCorr2Reduced [ ibin ][ itrack ][ icent ] / sumM01 [ ibin ][ itrack ][ icent ] : -999.0;
                
                //if(!TMath::Finite(differential2_1)) differential2_1 = -999.0;
				Double_t     correction5 = (  sumMp [ ibin ][ itrack ][ icent ] > 0.0 &&       sumAddon5 [ ibin ][ itrack ][ icent ] != 0.0 )?
                                                 sumAddon5 [ ibin ][ itrack ][ icent ] /  sumMp [ ibin ][ itrack ][ icent ] : -999.0;
                //if(!TMath::Finite(    correction5))     correction5 = -999.0;
				Double_t	 correction6 = (  sumMp [ ibin ][ itrack ][ icent ] > 0.0 &&       sumAddon6 [ ibin ][ itrack ][ icent ] != 0.0 )?
                                                 sumAddon6 [ ibin ][ itrack ][ icent ] /  sumMp [ ibin ][ itrack ][ icent ] : -999.0;
                //if(!TMath::Finite(    correction6))     correction6 = -999.0;
                
				// error on mean
                Double_t differential2Error1 = GetError (         sumM01 [ ibin ][ itrack ][ icent ] ,          sumM01Squared [ ibin ][ itrack ][ icent ] ,
                                                         sumCorr2Reduced [ ibin ][ itrack ][ icent ] , sumCorr2ReducedSquared [ ibin ][ itrack ][ icent ] );
                
                Double_t    correction5Error = GetError (    sumMp [ ibin ][ itrack ][ icent ] ,     sumMpSquared [ ibin ][ itrack ][ icent ] ,
                                                         sumAddon5 [ ibin ][ itrack ][ icent ] , sumAddon5Squared [ ibin ][ itrack ][ icent ] );
                Double_t    correction6Error = GetError (    sumMp [ ibin ][ itrack ][ icent ] ,     sumMpSquared [ ibin ][ itrack ][ icent ] ,
                                                         sumAddon6 [ ibin ][ itrack ][ icent ] , sumAddon6Squared [ ibin ][ itrack ][ icent ] );
                
                if ( differential2_1 == -999.0 || differential2Error1  == -999.0 || differential2_1 * differential2Error1 == 0.0 ) { differential2_1 = -999.0; differential2Error1 = -999.0; }
                if (     correction5 == -999.0 ||    correction5Error  == -999.0 ||     correction5 *    correction5Error == 0.0 ) {     correction5 = -999.0;    correction5Error = -999.0; }
                if (     correction6 == -999.0 ||    correction6Error  == -999.0 ||     correction6 *    correction6Error == 0.0 ) {     correction6 = -999.0;    correction6Error = -999.0; }
                
                // Covariance terms
                Double_t covRef2Def2 = GetCovariance (   sumM11 [ itrack ][ icent ]         ,             sumCorr2 [ itrack ][ icent ]         ,
                                                         sumM01 [ ibin ][ itrack ][ icent ] ,      sumCorr2Reduced [ ibin ][ itrack ][ icent ] ,
                                                      sumM11M01 [ ibin ][ itrack ][ icent ] , sumCorr2Corr2Reduced [ ibin ][ itrack ][ icent ] );
                //if ( TMath::Abs ( 1.0 - sumM11M01 [ ibin ][ itrack ][ icent ] / ( sumM11 [ itrack ][ icent ] * sumM01 [ ibin ][ itrack ][ icent ] ) ) <= 1.0e-6 ) covRef2Def2 = -999.0;
                
                Double_t covRef2C5 = GetCovariance (  sumM11 [ itrack ][ icent ]         ,   sumCorr2 [ itrack ][ icent ]         ,
                                                       sumMp [ ibin ][ itrack ][ icent ] ,  sumAddon5 [ ibin ][ itrack ][ icent ] ,
                                                    sumM11mp [ ibin ][ itrack ][ icent ] , sumCorr2C5 [ ibin ][ itrack ][ icent ] );
                Double_t covRef2C6 = GetCovariance (  sumM11 [ itrack ][ icent ]         ,   sumCorr2 [ itrack ][ icent ]         ,
                                                       sumMp [ ibin ][ itrack ][ icent ] ,  sumAddon6 [ ibin ][ itrack ][ icent ] ,
                                                    sumM11mp [ ibin ][ itrack ][ icent ] , sumCorr2C6 [ ibin ][ itrack ][ icent ] );
                
                Double_t covDef2C1 = GetCovariance (   sumM01 [ ibin ][ itrack ][ icent ] ,   sumCorr2Reduced [ ibin ][ itrack ][ icent ] ,
                                                       sumS11 [ itrack ][ icent ]         ,         sumAddon1 [ itrack ][ icent ]         ,
                                                    sumS11M01 [ ibin ][ itrack ][ icent ] , sumCorr2ReducedC1 [ ibin ][ itrack ][ icent ] );
                Double_t covDef2C2 = GetCovariance (   sumM01 [ ibin ][ itrack ][ icent ] ,   sumCorr2Reduced [ ibin ][ itrack ][ icent ] ,
                                                       sumS11 [ itrack ][ icent ]         ,         sumAddon2 [ itrack ][ icent ]         ,
                                                    sumS11M01 [ ibin ][ itrack ][ icent ] , sumCorr2ReducedC2 [ ibin ][ itrack ][ icent ] );
                Double_t covDef2C3 = GetCovariance (    sumM01 [ ibin ][ itrack ][ icent ] ,   sumCorr2Reduced [ ibin ][ itrack ][ icent ] ,
                                                       sumbS11 [ itrack ][ icent ]         ,         sumAddon3 [ itrack ][ icent ]         ,
                                                    sumbS11M01 [ ibin ][ itrack ][ icent ] , sumCorr2ReducedC3 [ ibin ][ itrack ][ icent ] );
                Double_t covDef2C4 = GetCovariance (    sumM01 [ ibin ][ itrack ][ icent ] ,   sumCorr2Reduced [ ibin ][ itrack ][ icent ] ,
                                                       sumbS11 [ itrack ][ icent ]         ,         sumAddon4 [ itrack ][ icent ]         ,
                                                    sumbS11M01 [ ibin ][ itrack ][ icent ] , sumCorr2ReducedC4 [ ibin ][ itrack ][ icent ] );
                Double_t covDef2C5 = GetCovariance (  sumM01 [ ibin ][ itrack ][ icent ] ,   sumCorr2Reduced [ ibin ][ itrack ][ icent ] ,
                                                       sumMp [ ibin ][ itrack ][ icent ] ,         sumAddon5 [ ibin ][ itrack ][ icent ] ,
                                                    sumM01mp [ ibin ][ itrack ][ icent ] , sumCorr2ReducedC5 [ ibin ][ itrack ][ icent ] );
                Double_t covDef2C6 = GetCovariance (  sumM01 [ ibin ][ itrack ][ icent ] ,   sumCorr2Reduced [ ibin ][ itrack ][ icent ] ,
                                                       sumMp [ ibin ][ itrack ][ icent ] ,         sumAddon6 [ ibin ][ itrack ][ icent ] ,
                                                    sumM01mp [ ibin ][ itrack ][ icent ] , sumCorr2ReducedC6 [ ibin ][ itrack ][ icent ] );
                
                Double_t covC1C5 = GetCovariance (  sumS11 [ itrack ][ icent ]         , sumAddon1 [ itrack ][ icent ]         ,
                                                     sumMp [ ibin ][ itrack ][ icent ] , sumAddon5 [ ibin ][ itrack ][ icent ] ,
                                                  sumS11mp [ ibin ][ itrack ][ icent ] ,   sumC1C5 [ ibin ][ itrack ][ icent ] );
                Double_t covC1C6 = GetCovariance (  sumS11 [ itrack ][ icent ]         , sumAddon1 [ itrack ][ icent ]         ,
                                                     sumMp [ ibin ][ itrack ][ icent ] , sumAddon6 [ ibin ][ itrack ][ icent ] ,
                                                  sumS11mp [ ibin ][ itrack ][ icent ] ,   sumC1C6 [ ibin ][ itrack ][ icent ] );
                
                Double_t covC2C5 = GetCovariance (  sumS11 [ itrack ][ icent ]         , sumAddon2 [ itrack ][ icent ]         ,
                                                     sumMp [ ibin ][ itrack ][ icent ] , sumAddon5 [ ibin ][ itrack ][ icent ] ,
                                                  sumS11mp [ ibin ][ itrack ][ icent ] ,   sumC2C5 [ ibin ][ itrack ][ icent ] );
                Double_t covC2C6 = GetCovariance (  sumS11 [ itrack ][ icent ]         , sumAddon2 [ itrack ][ icent ]         ,
                                                     sumMp [ ibin ][ itrack ][ icent ] , sumAddon6 [ ibin ][ itrack ][ icent ] ,
                                                  sumS11mp [ ibin ][ itrack ][ icent ] ,   sumC2C6 [ ibin ][ itrack ][ icent ] );
                
                Double_t covC3C5 = GetCovariance (  sumbS11 [ itrack ][ icent ]         , sumAddon3 [ itrack ][ icent ]         ,
                                                      sumMp [ ibin ][ itrack ][ icent ] , sumAddon5 [ ibin ][ itrack ][ icent ] ,
                                                  sumbS11mp [ ibin ][ itrack ][ icent ] ,   sumC3C5 [ ibin ][ itrack ][ icent ] );
                Double_t covC3C6 = GetCovariance (  sumbS11 [ itrack ][ icent ]         , sumAddon3 [ itrack ][ icent ]         ,
                                                      sumMp [ ibin ][ itrack ][ icent ] , sumAddon6 [ ibin ][ itrack ][ icent ] ,
                                                  sumbS11mp [ ibin ][ itrack ][ icent ] ,   sumC3C6 [ ibin ][ itrack ][ icent ] );
                
                Double_t covC4C5 = GetCovariance (  sumbS11 [ itrack ][ icent ]         , sumAddon4 [ itrack ][ icent ]         ,
                                                      sumMp [ ibin ][ itrack ][ icent ] , sumAddon5 [ ibin ][ itrack ][ icent ] ,
                                                  sumbS11mp [ ibin ][ itrack ][ icent ] ,   sumC4C5 [ ibin ][ itrack ][ icent ] );
                Double_t covC4C6 = GetCovariance (  sumbS11 [ itrack ][ icent ]         , sumAddon4 [ itrack ][ icent ]         ,
                                                      sumMp [ ibin ][ itrack ][ icent ] , sumAddon6 [ ibin ][ itrack ][ icent ] ,
                                                  sumbS11mp [ ibin ][ itrack ][ icent ] ,   sumC4C6 [ ibin ][ itrack ][ icent ] );
                
                Double_t covC5C6 = GetCovariance (       sumMp [ ibin ][ itrack ][ icent ] , sumAddon5 [ ibin ][ itrack ][ icent ] ,
                                                         sumMp [ ibin ][ itrack ][ icent ] , sumAddon6 [ ibin ][ itrack ][ icent ] ,
                                                  sumMpSquared [ ibin ][ itrack ][ icent ] ,   sumC5C6 [ ibin ][ itrack ][ icent ] );
                
                // Corrected Differential Flows
                Double_t differential2_2 = -999.0;
                if( differential2_1 != -999.0 && correction3 != -999.0 && correction4 != -999.0 && correction5 != -999.0 && correction6 != -999.0 ) {
                    differential2_2 = differential2_1 - correction5*correction3 - correction6*correction4;
                }
                /*if(!TMath::Finite(differential2_2) || differential2_1 == -999.0
                   || correction3 == -999.0 || correction4 == -999.0 || correction5 == -999.0 || correction6 == -999.0) differential2_2 = -999.0;*/
                
				// Final flows
				v1_2 [ ibin ][ itrack ][ icent ] = ( differential2_2 != -999.0 && reference2_3 != -999.0 )? sign2 * differential2_2 / reference2_3 : -999.0;
                //if(!TMath::Finite(v1_2 [ ibin ][ itrack ][ icent ]) || differential2_2 == -999.0 || reference2_3 == -999.0) v1_2 [ ibin ][ itrack ][ icent ] = -999.0;
                
				// differentiate QC2 differential flow
                Double_t df2dr2 = -999.0;
                if( v1_2 [ ibin ][ itrack ][ icent ] != -999.0 ) {
                    df2dr2 = TMath::Power ( sign2 * reference2_2 , -1.0 ) * ( - differential2_2 * 0.5 * TMath::Power ( sign2 * reference2_2 , -0.5 ) * 1.0 );
                }
                //if(!TMath::Finite(df2dr2)) df2dr2 = -999.0;
                Double_t df2dd2 = -999.0;
                if( v1_2 [ ibin ][ itrack ][ icent ] != -999.0 ) {
                    df2dd2 = TMath::Power ( sign2 * reference2_2 , -1.0 ) * ( TMath::Power ( sign2 * reference2_2 , 0.5 ) * sign2 * 1.0 );
                }
                //if(!TMath::Finite(df2dd2)) df2dd2 = -999.0;
                Double_t df2dc1 = -999.0;
                if( v1_2 [ ibin ][ itrack ][ icent ] != -999.0 ) {
                    df2dc1 = TMath::Power ( sign2 * reference2_2 , -1.0 ) * ( - differential2_2 * 0.5 * TMath::Power ( sign2 * reference2_2 , -0.5 ) * ( -correction3 ) );
                }
                //if(!TMath::Finite(df2dc1)) df2dc1 = -999.0;
                Double_t df2dc2 = -999.0;
                if( v1_2 [ ibin ][ itrack ][ icent ] != -999.0 ) {
                    df2dc2 = TMath::Power ( sign2 * reference2_2 , -1.0 ) * ( - differential2_2 * 0.5 * TMath::Power ( sign2 * reference2_2 , -0.5 ) * ( -correction4 ) );
                }
                //if(!TMath::Finite(df2dc2)) df2dc2 = -999.0;
                Double_t df2dc3 = -999.0;
                if( v1_2 [ ibin ][ itrack ][ icent ] != -999.0 ) {
                    df2dc3 = TMath::Power ( sign2 * reference2_2 , -1.0 ) * ( TMath::Power ( sign2 * reference2_2 , 0.5 ) * sign2 * ( -correction5 )
                                                                             - differential2_2 * 0.5 * TMath::Power ( sign2 * reference2_2 , -0.5 ) * ( -correction1 ) );
                }
                //if(!TMath::Finite(df2dc3)) df2dc3 = -999.0;
                Double_t df2dc4 = -999.0;
                if( v1_2 [ ibin ][ itrack ][ icent ] != -999.0 ) {
                    df2dc4 = TMath::Power ( sign2 * reference2_2 , -1.0 ) * ( TMath::Power ( sign2 * reference2_2 , 0.5 ) * sign2 * ( -correction6 )
                                                                             - differential2_2 * 0.5 * TMath::Power ( sign2 * reference2_2 , -0.5 ) * ( -correction2 ) );
                }
                //if(!TMath::Finite(df2dc4)) df2dc4 = -999.0;
                Double_t df2dc5 = -999.0;
                if( v1_2 [ ibin ][ itrack ][ icent ] != -999.0 ) {
                    df2dc5 = TMath::Power ( sign2 * reference2_2 , -1.0 ) * ( TMath::Power ( sign2 * reference2_2 , 0.5 ) * sign2 * ( -correction3 ) );
                }
                //if(!TMath::Finite(df2dc5)) df2dc5 = -999.0;
                Double_t df2dc6 = -999.0;
                if( v1_2 [ ibin ][ itrack ][ icent ] != -999.0 ) {
                    df2dc6 = TMath::Power ( sign2 * reference2_2 , -1.0 ) * ( TMath::Power ( sign2 * reference2_2 , 0.5 ) * sign2 * ( -correction4 ) );
                }
                //if(!TMath::Finite(df2dc6)) df2dc6 = -999.0;
                
				// error propagation
				v1_2Error [ ibin ][ itrack ][ icent ] = -999.0;
				if ( df2dr2 != -999.0 && reference2Error1 != -999.0
                    && df2dd2 != -999.0 && differential2Error1 != -999.0
                    
                    && df2dc1 != -999.0 && correction1Error != -999.0
                    && df2dc2 != -999.0 && correction2Error != -999.0
                    && df2dc3 != -999.0 && correction3Error != -999.0
                    && df2dc4 != -999.0 && correction4Error != -999.0
                    && df2dc5 != -999.0 && correction5Error != -999.0
                    && df2dc6 != -999.0 && correction6Error != -999.0
                    
                    && covRef2Def2 != -999.0
                    && covRef2C1 != -999.0
                    && covRef2C2 != -999.0
                    && covRef2C3 != -999.0
                    && covRef2C4 != -999.0
                    && covRef2C5 != -999.0
                    && covRef2C6 != -999.0
                    
                    && covDef2C1 != -999.0
                    && covDef2C2 != -999.0
                    && covDef2C3 != -999.0
                    && covDef2C4 != -999.0
                    && covDef2C5 != -999.0
                    && covDef2C6 != -999.0
                    
                    && covC1C2 != -999.0
                    && covC1C3 != -999.0
                    && covC1C4 != -999.0
                    && covC1C5 != -999.0
                    && covC1C6 != -999.0
                    
                    && covC2C3 != -999.0
                    && covC2C4 != -999.0
                    && covC2C5 != -999.0
                    && covC2C6 != -999.0
                    
                    && covC3C4 != -999.0
                    && covC3C5 != -999.0
                    && covC3C6 != -999.0
                    
                    && covC4C5 != -999.0
                    && covC4C6 != -999.0
                    
                    && covC5C6 != -999.0
                    ) {
                    
                    v1_2Error [ ibin ][ itrack ][ icent ] = TMath::Sqrt ( TMath::Power ( df2dr2 , 2.0 ) * reference2Error1
                                                                         + TMath::Power ( df2dd2 , 2.0 ) * differential2Error1
                                                                         
                                                                         + TMath::Power ( df2dc1 , 2.0 ) * correction1Error
                                                                         + TMath::Power ( df2dc2 , 2.0 ) * correction2Error
                                                                         + TMath::Power ( df2dc3 , 2.0 ) * correction3Error
                                                                         + TMath::Power ( df2dc4 , 2.0 ) * correction4Error
                                                                         + TMath::Power ( df2dc5 , 2.0 ) * correction5Error
                                                                         + TMath::Power ( df2dc6 , 2.0 ) * correction6Error
                                                                         
                                                                         + 2.0 * df2dr2 * df2dd2 * covRef2Def2
                                                                         + 2.0 * df2dr2 * df2dc1 * covRef2C1
                                                                         + 2.0 * df2dr2 * df2dc2 * covRef2C2
                                                                         + 2.0 * df2dr2 * df2dc3 * covRef2C3
                                                                         + 2.0 * df2dr2 * df2dc4 * covRef2C4
                                                                         + 2.0 * df2dr2 * df2dc5 * covRef2C5
                                                                         + 2.0 * df2dr2 * df2dc6 * covRef2C6
                                                                         
                                                                         + 2.0 * df2dd2 * df2dc1 * covDef2C1
                                                                         + 2.0 * df2dd2 * df2dc2 * covDef2C2
                                                                         + 2.0 * df2dd2 * df2dc3 * covDef2C3
                                                                         + 2.0 * df2dd2 * df2dc4 * covDef2C4
                                                                         + 2.0 * df2dd2 * df2dc5 * covDef2C5
                                                                         + 2.0 * df2dd2 * df2dc6 * covDef2C6
                                                                         
                                                                         + 2.0 * df2dc1 * df2dc2 * covC1C2
                                                                         + 2.0 * df2dc1 * df2dc3 * covC1C3
                                                                         + 2.0 * df2dc1 * df2dc4 * covC1C4
                                                                         + 2.0 * df2dc1 * df2dc5 * covC1C5
                                                                         + 2.0 * df2dc1 * df2dc6 * covC1C6
                                                                         
                                                                         + 2.0 * df2dc2 * df2dc3 * covC2C3
                                                                         + 2.0 * df2dc2 * df2dc4 * covC2C4
                                                                         + 2.0 * df2dc2 * df2dc5 * covC2C5
                                                                         + 2.0 * df2dc2 * df2dc6 * covC2C6
                                                                         
                                                                         + 2.0 * df2dc3 * df2dc4 * covC3C4
                                                                         + 2.0 * df2dc3 * df2dc5 * covC3C5
                                                                         + 2.0 * df2dc3 * df2dc6 * covC3C6
                                                                         
                                                                         + 2.0 * df2dc4 * df2dc5 * covC4C5
                                                                         + 2.0 * df2dc4 * df2dc6 * covC4C6
                                                                         
                                                                         + 2.0 * df2dc5 * df2dc6 * covC5C6
                                                                         );
                    
                    if ( !TMath::Finite ( v1_2Error [ ibin ][ itrack ][ icent ] ) ) v1_2Error [ ibin ][ itrack ][ icent ] = -999.0;
                    
				}
                
                if ( v1_2 [ ibin ][ itrack ][ icent ] == -999.0 || v1_2Error [ ibin ][ itrack ][ icent ] == -999.0
                    || v1_2 [ ibin ][ itrack ][ icent ] * v1_2Error [ ibin ][ itrack ][ icent ] == 0.0
                    || multCounter [ itrack ][ icent ] < 5 || binMultCounter [ ibin ][ itrack ][ icent ] < 5 ) {
                    v1_2 [ ibin ][ itrack ][ icent ] = -999.0; v1_2Error [ ibin ][ itrack ][ icent ] = -999.0;
                }
				
				// QC2 differential flow
                /*if ( //reference2_3 != -999.0 && reference2Error3 != -999.0
                    !TMath::Finite(v1_2 [ ibin ][ itrack ][ icent ]) || !TMath::Finite(v1_2Error [ ibin ][ itrack ][ icent ])
                    ) {
                    v1_2 [ ibin ][ itrack ][ icent ] = -999.0; v1_2Error [ ibin ][ itrack ][ icent ] = -999.0;
                }*/
                
                /*std::cout << icent << "  " << itrack << "  " << ibin << "  " << reference2_1 << "  " << sign2*reference2_2 << "  " << reference2_3 << "  "
                            << differential2_1 << "  " << sign2 * differential2_2 << "  " << v1_2 [ ibin ][ itrack ][ icent ] << "  " << v1_2Error [ ibin ][ itrack ][ icent ] 
                            << std::endl << std::endl;*/
                
			}    // bin loop ends
            
		}    // track loop ends
        
		std::cout << "Finished computing centrality bin " << icent + 1 << std::endl << std::endl;
        
    }    // centrality loop ends
    
    // average over centrality bins
    Double_t /*content1 [ Nbins ], error1 [ Nbins ],*/ content2 [ Nbins ], error2 [ Nbins ], content3 [ Nbins ], error3 [ Nbins ], content4 [ Nbins ], error4 [ Nbins ];
    for ( Int_t ibin = 0; ibin < Nbins; ibin ++ ) {
        
        /*content1 [ ibin ] = 0.0; error1 [ ibin ] = 0.0;*/ content2 [ ibin ] = 0.0; error2 [ ibin ] = 0.0; content4 [ ibin ] = 0.0; error4 [ ibin ] = 0.0;
        
        for ( Int_t icent = 0; icent < Ncentralities; icent ++ ) {   // (0-80)%
        //for ( Int_t icent = 1; icent < 5; icent++ ) {   // (10-30)%

			for ( Int_t itrack = 0; itrack < Ntracks; itrack ++ ) {
                
                if( v1_2 [ ibin ][ itrack ][ icent ] != -999.0 && v1_2Error [ ibin ][ itrack ][ icent ] != -999.0 ) {
                    content2 [ ibin ] += v1_2 [ ibin ][ itrack ][ icent ] / TMath::Power ( v1_2Error [ ibin ][ itrack ][ icent ] , 2.0 );
                    error2 [ ibin ] += 1.0 / TMath::Power ( v1_2Error [ ibin ][ itrack ][ icent ] , 2.0 );
                }
                
                
                /*content4 [ ibin ] += ( v1_4 [ ibin ][ itrack ][ icent ] != -999.0 && v1_4Error [ ibin ][ itrack ][ icent ] != -999.0 )?
                                     v1_4 [ ibin ][ itrack ][ icent ] / TMath::Power ( v1_4Error [ ibin ][ itrack ][ icent ] , 2.0 ) : 0.0;
                  error4 [ ibin ] += ( v1_4 [ ibin ][ itrack ][ icent ] != -999.0 && v1_4Error [ ibin ][ itrack ][ icent ] != -999.0 )?
                                     1.0 / TMath::Power ( v1_4Error [ ibin ][ itrack ][ icent ] , 2.0 ) : 0.0;*/
                
            }
            
        }
        if( error2 [ ibin ] > 0.0 ) {
            content2 [ ibin ] = content2 [ ibin ] / error2 [ ibin ];
            //if ( !TMath::Finite ( content2 [ ibin ] ) ) content2 [ ibin ] = -999.0;
              error2 [ ibin ] = TMath::Sqrt ( 1.0 / error2 [ ibin ] );
        }
        
        //if ( !TMath::Finite ( error2 [ ibin ] ) ) error2 [ ibin ] = -999.0;
        if ( content2 [ ibin ] == -999.0 || error2 [ ibin ] == -999.0 || content2 [ ibin ] * error2 [ ibin ] == 0.0 ) { content2 [ ibin ] = 0.0; error2 [ ibin ] = 0.0; }
        
        /*content4 [ ibin ] = ( error4 [ ibin ] > 0.0 )? content4 [ ibin ] / error4 [ ibin ] : -999.0;
        //if ( !TMath::Finite ( content4 [ ibin ] ) ) content4 [ ibin ] = -999.0;
          error4 [ ibin ] = ( error4 [ ibin ] > 0.0 )? TMath::Sqrt ( 1.0 / error4 [ ibin ] ) : -999.0;
        //if ( !TMath::Finite ( error4 [ ibin ] ) ) error4 [ ibin ] = -999.0;
        if ( content4 [ ibin ] == -999.0 || error4 [ ibin ] == -999.0 || content4 [ ibin ] * error4 [ ibin ] == 0.0 ) { content4 [ ibin ] = 0.0; error4 [ ibin ] = 0.0; }*/
        
		// filling v1 plots
        hv1_2 -> SetBinContent ( ibin + 1 , content2 [ ibin ] ); hv1_2 -> SetBinError ( ibin + 1 , error2 [ ibin ] );
        //hv1_4 -> SetBinContent ( ibin + 1 , content4 [ ibin ] ); hv1_4 -> SetBinError ( ibin + 1 , error4 [ ibin ] );
        
		// filling v2 plots
        hv2_2 -> SetBinContent ( ibin + 1 , content2 [ ibin ] ); hv2_2 -> SetBinError ( ibin + 1 , error2 [ ibin ] );
        //hv2_4 -> SetBinContent ( ibin + 1 , content4 [ ibin ] ); hv2_4 -> SetBinError ( ibin + 1 , error4 [ ibin ] );
		
		// filling v3 plots
        hv3_2 -> SetBinContent ( ibin + 1 , content2 [ ibin ] ); hv3_2 -> SetBinError ( ibin + 1 , error2 [ ibin ] );
        //hv3_4 -> SetBinContent ( ibin + 1 , content4 [ ibin ] ); hv3_4 -> SetBinError ( ibin + 1 , error4 [ ibin ] );
        
      }
    
	/*Double_t resolution[Ncentralities], resolutionError[Ncentralities], resolution2[Ncentralities], resolution2Error[Ncentralities];
	TH2D *temp1 = new TH2D("temp1","temp1",Ncentralities,0.5,Ncentralities+0.5,60,-3.,3.);
	TH2D *temp2 = new TH2D("temp2","temp2",Ncentralities,0.5,Ncentralities+0.5,60,0.,6.);
    for ( Int_t icent = 0; icent < Ncentralities; icent++ )    // (0-80)%
    //for ( Int_t icent = 4; icent < 7; icent++ )    // (10-40)%
      {
        // get Event Plane resolution <cos1*1[psi_1^a - psi_1^b]> k=1 m=1
        resolution[icent] = 0.0; resolutionError[icent] = 0.0;
        resolution[icent] = hres->GetBinContent(icent+1,5); resolutionError[icent] = hres->GetBinError(icent+1,5);
        // get sub Event Plane resolution
        resolutionError[icent] = TMath::Abs(1.*0.5*TMath::Power ( 1.*resolution[icent],-0.5)*resolutionError[icent]); resolution[icent] = TMath::Sqrt(1.*resolution[icent]);
        // get full Event Plane resolution
        Double_t sub_res_chi = chi(1,resolution[icent]);
        resolution2[icent] = resEventPlane(1,TMath::Sqrt(2.0)*sub_res_chi);
        Double_t sub_res_delta = 0.005, full_res_delta = 0.0;
        Double_t sub_chi_delta = chi(1,(resolution[icent]+sub_res_delta));
        full_res_delta = resEventPlane(1,TMath::Sqrt(2.0)*sub_chi_delta);
        resolution2Error[icent] = resolutionError[icent] * TMath::Abs(resolution2[icent] - full_res_delta) / sub_res_delta;
        // replace resolutions
        resolution[icent] = resolution2[icent];
        resolutionError[icent] = resolution2Error[icent];
        for ( Int_t ibin = 0; ibin < 60; ibin++ )
          {
            content1[ibin] = 0.0; error1[ibin] = 0.0;
			content1[ibin] = hprof1->GetBinContent(icent+1,ibin+1); error1[ibin] = hprof1->GetBinError(icent+1,ibin+1);
			error1[ibin] = TMath::Sqrt( TMath::Power ( error1[ibin]/resolution[icent],2.)
                                        + TMath::Power (  -content1[ibin]*resolutionError[icent]/TMath::Power ( resolution[icent],2.) , 2.) );
            content1[ibin] = content1[ibin]/resolution[icent];
			temp1->SetBinContent(icent+1,ibin+1,content1[ibin]); temp1->SetBinError(icent+1,ibin+1,error1[ibin]);
		  }
        // get Event Plane resolution <cos2*1[psi_1^a - psi_1^b]> k=2 m=1
        resolution[icent] = 0.0; resolutionError[icent] = 0.0;
        resolution[icent] = hres->GetBinContent(icent+1,6); resolutionError[icent] = hres->GetBinError(icent+1,6);
        // get sub Event Plane resolution
        resolutionError[icent] = TMath::Abs(1.*0.5*TMath::Power ( 1.*resolution[icent],-0.5)*resolutionError[icent]); resolution[icent] = TMath::Sqrt(1.*resolution[icent]);
        // get full Event Plane resolution
        sub_res_chi = chi(2,resolution[icent]);
        resolution2[icent] = resEventPlane(2,TMath::Sqrt(2.0)*sub_res_chi);
        sub_res_delta = 0.005; full_res_delta = 0.0;
        sub_chi_delta = chi(2,(resolution[icent]+sub_res_delta));
        full_res_delta = resEventPlane(2,TMath::Sqrt(2.0)*sub_chi_delta);
        resolution2Error[icent] = resolutionError[icent] * TMath::Abs(resolution2[icent] - full_res_delta) / sub_res_delta;
        // replace resolutions
        resolution[icent] = resolution2[icent];
        resolutionError[icent] = resolution2Error[icent];
        for ( Int_t ibin = 0; ibin < 60; ibin++ )
          {
            content2[ibin] = 0.0; error2[ibin] = 0.0;
            content2[ibin] = hprof2->GetBinContent(icent+1,ibin+1); error2[ibin] = hprof2->GetBinError(icent+1,ibin+1);
            error2[ibin] = TMath::Sqrt( TMath::Power ( error2[ibin]/resolution[icent],2.)
                                       + TMath::Power (  -content2[ibin]*resolutionError[icent]/TMath::Power ( resolution[icent],2.) , 2.) );
            content2[ibin] = content2[ibin]/resolution[icent];
            temp2->SetBinContent(icent+1,ibin+1,content2[ibin]); temp2->SetBinError(icent+1,ibin+1,error2[ibin]);
          }
	  }*/
    
    f_out -> Write ();
    return 0;
    
}
