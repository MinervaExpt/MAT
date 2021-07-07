#ifndef POTCounter_cxx
#define POTCounter_cxx 1

#include "POTCounter.h"

#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TCollection.h>
#include <TLeaf.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TBranch.h>

using namespace PlotUtils;

double POTCounter::getPOTfromPlaylist(std::string playlist, bool getTotal, int batch_number, bool useFast)
{
    // Create TChain form Playlist and Initialize
    TChain* fChain = new TChain("Meta");
    Init(playlist, fChain);

    // Check TChain
    if (!fChain || fChain == 0){
        std::cout<<"PlotUtils::POTCounter -- Can NOT find TChain from input playlist! -- Returning -1"<<std::endl;
        return -1;
    }else{
        std::cout<<"PlotUtils::POTCounter -- Initialized fChain using the Input playlist!"<<std::endl;
    }

    // Get POT_Used and return
    double sumPOTUsed;
    if (useFast) sumPOTUsed = getPOTfromTChain_Fast(fChain, getTotal, batch_number);
    else sumPOTUsed = getPOTfromTChain(fChain, getTotal,batch_number);

    return sumPOTUsed;
}

double POTCounter::getPOTfromTChain(TChain* ch, bool getTotal, int batch_number)
{
    double sumPOTUsed = 0;
    char leafName[50];

    if( batch_number == -1 ) {
        std::cout<<"PlotUtils::POTCounter -- Counting POT"<<std::endl;
        if(getTotal) sprintf(leafName,"POT_Total");
        else sprintf(leafName,"POT_Used");
    }
    else if( 1 <= batch_number && batch_number <= 6 ) {
        std::cout<<"PlotUtils::POTCounter -- Counting POT from batch "<< batch_number<<std::endl;
        if(getTotal){
          std::cout <<"PlotUtils::POTCounter -- for getPOTfromTChain, there is no POT_Total for a batch"<<std::endl;
          std::cout <<"PlotUtils::POTCounter -- Please use batch number -1 "<<std::endl;
          return sumPOTUsed;
         } 
        else sprintf(leafName,"POT_Used_batch%d",batch_number);
    }
    else{
        std::cout <<"PlotUtils::POTCounter -- for getPOTfromTChain, batch number is 1 through 6"<<std::endl;
        std::cout <<"PlotUtils::POTCounter -- inputted batch number: "<<batch_number<<std::endl;
        return sumPOTUsed;
    }

    TObjArray* fileElements = ch->GetListOfFiles();
    TIter next(fileElements);
    TChainElement* chEl=0;

    std::cout<<"PlotUtils::POTCounter -- It can take some time depending on the size of the TChain"<<std::endl;

    while (( chEl=(TChainElement*)next() )) {
        TFile f(chEl->GetTitle());
        TTree* t = (TTree*)f.Get("Meta");
        if (!t){
            std::cout<<"PlotUtils::POTCounter -- No Meta tree in file "<<chEl->GetTitle()<<std::endl;
            continue;
        }

        // Loop Over all Entries 
        int n_entries = t->GetEntries();
        for (int i=0;i<n_entries;i++){
            t->GetEntry(i);
            TLeaf* POT_Used = t->GetLeaf(leafName);

            if (POT_Used){
                sumPOTUsed = sumPOTUsed + POT_Used->GetValue();
            }
        }
    }

    return sumPOTUsed;
}

double POTCounter::getPOTfromTChain_Fast(TChain* ch, bool getTotal, int batch_number)
{
    double sumPOTUsed = 0;

    std::cout<<"PlotUtils::POTCounter -- Called getPOTfromTChain_Fast()"<<std::endl;
    std::cout<<"PlotUtils::POTCounter -- This method works only for non-merged files, if your playlist includes merged ROOT files,use getPOTfromPlaylist(playlist, false)"<<std::endl;
    std::cout<<"PlotUtils::POTCounter -- This can't get batch POT";
    std::cout<<"PlotUtils::POTCounter -- Counting POT Fast"<<std::endl;
    std::cout<<"PlotUtils::POTCounter -- It can take some time depending on the size of the TChain"<<std::endl;

    Long64_t nentries = ch->GetEntriesFast();

    for (Long64_t jentry=0; jentry < nentries; jentry++) {

        Long64_t ientry = ch->GetEntry(jentry);

        if (ientry == 0) {
            std::cout<<"\tGetEntry failure "<<jentry<<std::endl;
            break;
        }

        if( batch_number == -1 ) {
            if(getTotal) sumPOTUsed = sumPOTUsed + POT_Total;
            else sumPOTUsed = sumPOTUsed + POT_Used; 
            
        }
        else if( 1 <= batch_number && batch_number <= 6 ) {
            if(getTotal){
              std::cout <<"PlotUtils::POTCounter -- for getPOTfromTChain, there is no POT_Total for a batch"<<std::endl;
              std::cout <<"PlotUtils::POTCounter -- Please use batch number -1 "<<std::endl;
              return sumPOTUsed;
             } 
             if(batch_number == 1) sumPOTUsed = sumPOTUsed + POT_Used_batch1;
             if(batch_number == 2) sumPOTUsed = sumPOTUsed + POT_Used_batch2;
             if(batch_number == 3) sumPOTUsed = sumPOTUsed + POT_Used_batch3;
             if(batch_number == 4) sumPOTUsed = sumPOTUsed + POT_Used_batch4;
             if(batch_number == 5) sumPOTUsed = sumPOTUsed + POT_Used_batch5;
             if(batch_number == 6) sumPOTUsed = sumPOTUsed + POT_Used_batch6;
        }
        else{
            std::cout <<"PlotUtils::POTCounter -- for getPOTfromTChain, batch number is 1 through 6"<<std::endl;
            std::cout <<"PlotUtils::POTCounter -- inputted batch number: "<<batch_number<<std::endl;
            return sumPOTUsed;
        }
        
        // Progress Message on Terminal
        int msg_entry = 10000;
        if (jentry%msg_entry == 0) std::cout<<"\tCurrent Total POT = "<<sumPOTUsed<<std::endl;
    }
    return sumPOTUsed;
}

void POTCounter::Init(std::string playlist, TChain* fChain)
{
  std::ifstream input_pl(playlist.c_str());
    std::string filename;

    if( !input_pl.is_open() ){
        std::cerr<<"PlotUtils::POTCounter -- Cannot open Playlist File!"<<std::endl;
        exit(1);
    }else{
        std::cout<<"PlotUtils::POTCounter -- Reading Playlist: "<<playlist.c_str()<<std::endl;
    }

    /* 
     * Loop input playlist and get file names 
     *     Assumption: file names start with '/' character 
     */
    while (!input_pl.eof()) {
        getline(input_pl,filename);

        if (filename[0] != '/') continue;

        fChain->Add( filename.c_str() );
        //std::cout<<"PlotUtils::POTCounter -- Added "<<filename.c_str()<<std::endl;
    }

    fChain->SetMakeClass(1);
    fChain->SetBranchAddress("POT_Total", &POT_Total, &b_POT_Total);
    fChain->SetBranchAddress("POT_Used", &POT_Used, &b_POT_Used);
    fChain->SetBranchAddress("POT_Used_batch1", &POT_Used_batch1, &b_POT_Used_batch1);
    fChain->SetBranchAddress("POT_Used_batch2", &POT_Used_batch2, &b_POT_Used_batch2);
    fChain->SetBranchAddress("POT_Used_batch3", &POT_Used_batch3, &b_POT_Used_batch3);
    fChain->SetBranchAddress("POT_Used_batch4", &POT_Used_batch4, &b_POT_Used_batch4);
    fChain->SetBranchAddress("POT_Used_batch5", &POT_Used_batch5, &b_POT_Used_batch5);
    fChain->SetBranchAddress("POT_Used_batch6", &POT_Used_batch6, &b_POT_Used_batch6);
    fChain->SetBranchAddress("nEntries_Header", &nEntries_Header, &b_nEntries_Header);
    fChain->SetBranchAddress("nEntries_Truth", &nEntries_Truth, &b_nEntries_Truth);

    input_pl.close();
}

POTCounter::POTCounter()
{
    // Do Nothing!
}

#endif

