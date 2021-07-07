/*
================================================================================
Class: PlotUtils::POTCounter 
    A Class designed to Count POT_Used for a given 
        1) TChain
        2) Playlist - A list of root files

Example Usage:
    std::string playlist = "Input/Playlists/pl_MC_All.dat";
    
    POTCounter pot_counter;
    double totalPOT = pot_counter.getPOTfromPlaylist(playlist);
    
    std::cout<<"Total POT = "<<totalPOT<<std::endl;

Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/

#ifndef POTCounter_h
#define POTCounter_h 1

// C++ Libraries
#include <string>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <cstdlib>

// ROOT Libraries
#include <TROOT.h>

class TChain;
class TBranch;

namespace PlotUtils{

    class POTCounter {

        public:
            // Default Constructor
            POTCounter();
            double getPOTfromTChain(TChain* ch,bool getTotal = false,int batch_number = -1);
            double getPOTfromTChain_Fast(TChain* ch,bool getTotal=false, int batch_number = -1);
            double getPOTfromPlaylist(std::string playlist, bool getTotal=false,int batch_number = -1, bool useFast=false);            

        private:
            void Init(std::string playlist, TChain* fChain);
             // List of Variables in Meta Tree 
            Double_t        POT_Total;
            Double_t        POT_Used;
            Double_t        POT_Used_batch1;
            Double_t        POT_Used_batch2;
            Double_t        POT_Used_batch3;
            Double_t        POT_Used_batch4;
            Double_t        POT_Used_batch5;
            Double_t        POT_Used_batch6;
            Int_t           nEntries_Header;
            Int_t           nEntries_Truth;

            // List of branches
            TBranch        *b_POT_Total;   //!
            TBranch        *b_POT_Used;   //!
            TBranch        *b_POT_Used_batch1;
            TBranch        *b_POT_Used_batch2;
            TBranch        *b_POT_Used_batch3;
            TBranch        *b_POT_Used_batch4;
            TBranch        *b_POT_Used_batch5;
            TBranch        *b_POT_Used_batch6;
            TBranch        *b_nEntries_Header;   //!
            TBranch        *b_nEntries_Truth;   //!

    }; // end of class POTCounter
} // end of namespace PlotUtils

#endif

