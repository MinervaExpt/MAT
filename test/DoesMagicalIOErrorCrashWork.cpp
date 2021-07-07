//File: DoesMagicalIOErrorCrashWork.cpp
//Brief: Test simplest use case of ErrorHandler.h to crash a job processing MINERvA AnaTuples.
//Usage: DoesMagicalIOErrorCrashWork.cpp CCQENu eventID /path/to/some/AnaTuple.root || echo "Test of ErrorHandler.h passed."
//
//       OR to test root macros:
//
//       root -l -b -q ../setup/.rootlogon.C 'DoesMagicalIOErrorCrashWork.cpp+g("CCQENu", "eventID", "root://fndca1.fnal.gov:1094///pnfs/fnal.gov/usr/minerva/persistent/users/drut1186/CCQENu_Anatuples/LowEnergy_BugFix/MC/ccqenu_mc_minerva13C_FinalProc/grid/central_value/minerva/ana/v10r8p9/00/01/32/12/SIM_minerva_00013212_Subruns_0112_CCQENu_Ana_Tuple_v10r8p9.root")' || echo "Job crashed."
//
//Author: Andrew Olivier aolivier@ur.rochester.edu

#define USAGE "DoesIOErrorHandlingWork: Test whether ErrorHandler.h works in PlotUtils.\n\n"\
              "Usage: DoesIOErrorHandlingWork <name of AnaTuple TTree> <name of a branch to read> /path/to/some/AnaTuple.root [more .root files...]\n\n"\
              "\tDoesIOErrorHandlingWork demonstrates how to use the ErrorHandler.h header from\n"\
              "\tPlotUtils to react to I/O errors while reading from a TTree.  It accepts exactly 2\n"\
              "\targuments:\n"\
              "\t\t1) The name of an AnaTuple TTree to read\n"\
              "\t\t2) The name of a branch to read\n"\
              "\t\t3) A MINERvA AnaTuple .root file to process\n"\
              "\tIt returns nonzero on failing to read a file, a TTree, or the specified branch.\n"\
              "\tA return code of 2 means that you did something wrong on the command line\n"\
              "\tand saw this message.\n"

#undef NDEBUG //TODO: Remove this.  I'm too lazy to put the framework in debug mode.

//ROOT includes
#include "TFile.h"
#include "TTree.h"

//c++ includes
#include <iostream>
#include <cassert>

//C UNIX include (POSIX?)
#include <unistd.h>
#include <stdio.h>

//PlotUtils includes I'm testing
#include "PlotUtils/CrashOnROOTMessage.h"

//The majority of the work is done in this function that can be used as a ROOT macro.
int DoesMagicalIOErrorCrashWork(const std::string& treeName, const std::string& branchName, const std::vector<std::string>& filesToRead)
{
  for(std::vector<std::string>::const_iterator whichFile = filesToRead.begin(); whichFile != filesToRead.end(); ++whichFile)
  {
    const std::string& fileName = *whichFile;
    #ifndef NDEBUG
      std::cout << "Reading " << fileName << "...\n";
    #endif //NDEBUG

    //See DoesIOErrorHandlingWork.cpp for an example of how to write code that reacts to specific failure modes.
    TFile* input = TFile::Open(fileName.c_str());

    //Normally, I'd check whether input is nullptr here, but ErrorHandler.h should cause an exception
    //to be thrown in all such cases that I know of.  These cases indicate that ErrorHandler.h isn't
    //working, but they're not actually the cases I'm testing for.
    assert(input != NULL && "ErrorHandler.h failed to intercept a file that could not be opened!");

    TTree* tree = (TTree*)input->Get(treeName.c_str());
    if(tree == NULL)
    {
      std::cerr << "Failed to read a TTree named " << treeName << "from file " << fileName << ".  ErrorHandler.h"
                << "isn't actually capable of catching this failure mode, so you'll have to check for it yourself.\n";
      return 1;
    }

    //Read some standard branch that every AnaTuple has
    int eventID = -1;
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus(branchName.c_str(), 1);
    tree->Branch(branchName.c_str(), &eventID);

    for(int whichEntry = 0; tree->GetEntry(whichEntry) > 0; ++whichEntry)
    {
      #ifndef NDEBUG
        if(whichEntry % 1000 == 0) std::cout << "Processing entry " << whichEntry << "...\n";
      #endif //NDEBUG
      //Your analysis code would normally go here
    }

    //Deleting input automatically deletes tree.
    //In anything but Eroica, just say std::unique_ptr<>.  Seriously.
    delete input;
  }

  //If no I/O errors
  return 0;
}

int DoesMagicalIOErrorCrashWork(const std::string& treeName, const std::string& branchName, const std::string& fileToRead)
{
  std::vector<std::string> onlyOneFile(1, fileToRead);
  return DoesMagicalIOErrorCrashWork(treeName, branchName, onlyOneFile);
}

#ifndef __CINT__
int main(const int argc, const char** argv)
{
  const int nOtherArgs = 2 + 1; //2 arguments documented plus the path to this executable

  if(argc < nOtherArgs || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))
  {
    std::cerr << argc - 1 << " command line arguments given, but at least " << nOtherArgs - 1 << " expected.\n\n" << USAGE;
    return 2;
  }

  std::vector<std::string> filesToRead(argv + nOtherArgs, argv + argc);
  //Attempt to read file names from STDIN.
  //This lets me pipe the output of find or
  //sed to this program to select files to read.
  if(!isatty(fileno(stdin)))
  {
    std::string line;
    while(std::getline(std::cin, line)) filesToRead.push_back(line);
  }

  return DoesMagicalIOErrorCrashWork(argv[1], argv[2], filesToRead);
}
#endif //__CINT__
