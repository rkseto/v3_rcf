#include <cstdio>
#include <iostream>
#include "TFile.h"
#include "TList.h"
#include "TH1.h"
#include "TKey.h"
#include "TObject.h"
#include "TObjArray.h"

Int_t combine (/*Int_t run,*/TString jobId,Int_t fileNumber, Int_t keyJumpIndex = 6) {
    // output file
	char name[1500];
	//sprintf(name,"run%d_combined.root",run);
    TString outputFilename = jobId;
    outputFilename.Append("_total.root");
	TFile *outputFile = new TFile(outputFilename,"recreate");
    // input file
    TFile *tmpFile[1500]; TObjArray *ObjArray = new TObjArray(300); Int_t Nkeys = 0;
    TString tmpFilename;
	for(Int_t i=0;i<fileNumber;i++) {
        sprintf(name,"_%d.root",i);
        tmpFilename = jobId;
        tmpFilename.Append(name);
        tmpFile[i] = new TFile(tmpFilename,"read");
        if(i == 0) Nkeys = tmpFile[i]->GetNkeys();
        //if(Nkeys == 0) return -1;
        if(i == 0) std::cout<<"There're "<<Nkeys-keyJumpIndex+1<<" keys"<<std::endl;
        TList *list = tmpFile[i]->GetListOfKeys();
        TIter next((TList*)list);
        TKey  *key;
        Int_t keyCounter = 0;
        Int_t keyJumper  = 0;
        while ( ( key = (TKey*)next() ) ) {
            keyJumper++;
            if(keyJumper < keyJumpIndex) continue;
            TObject *obj = tmpFile[i]->Get(key->GetName());
            //if(!obj) return -1;
            if((obj->IsA() == TDirectory::Class())) continue;
            if(obj->InheritsFrom(TH1::Class())) {
                if(i == 0) ObjArray->Add((TH1*)obj);
                if(i >= 1) ((TH1*)ObjArray->At(keyCounter))->Add((TH1*)obj);
                keyCounter++;
            }
        }
        // clean up
        if(i >= 1) {
            tmpFile[i]->Close();
            delete tmpFile[i];
        }
        std::cout<<"File "<<i+1<<" read in"<<std::endl;
    }
    // write merged histograms
    for(Int_t k=0;k<=Nkeys-keyJumpIndex;k++) {
        outputFile->cd();
        ((TH1*)ObjArray->At(k))->Write();
        std::cout<<"Histogram "<<k+1<<" merged"<<std::endl;
    }
	return 0;
}
