//File: DoesIOErrorHandlingWork.cpp
//Brief: Test using ErrorHandler.h to react to and crash a job processing MINERvA AnaTuples.
//       If this test works, it will return to the operating system with non-zero exit status.
//       So, if you're using it in a testing framework, you want to require that this program
//       FAILS.
//Usage: DoesIOErrorHandlingWork /path/to/some/AnaTuple.root || echo "Test of ErrorHandler.h passed."
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
              "\tIt will return -1 for the special case of failing to access a .root file via xrootd.\n"\
              "\tA return code of -2 means that you did something wrong on the command line\n"\
              "\tand saw this message.\n"

#undef NDEBUG //TODO: Remove this.  I'm too lazy to put the framework in debug mode.

//ROOT includes
#include "TFile.h"
#include "TTree.h"

//c++ includes
#include <iostream>
#include <memory>
#include <cassert>

//C UNIX include (POSIX?)
#include <unistd.h>
#include <stdio.h>

//PlotUtils includes I'm testing
#include "PlotUtils/CrashOnROOTMessage.h"

int main(const int argc, const char** argv)
{
  const int nOtherArgs = 2 + 1; //2 arguments documented plus the path to this executable

  if(argc < nOtherArgs || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))
  {
    std::cerr << argc - 1 << " command line arguments given, but at least " << nOtherArgs - 1 << " expected.\n\n" << USAGE;
    return -2;
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

  bool nSkippedFiles = 0; //Number of files I skipped because of errors from ROOT

  for(std::vector<std::string>::const_iterator whichFile = filesToRead.begin(); whichFile != filesToRead.end(); ++whichFile)
  {
    const std::string& fileName = *whichFile;
    #ifndef NDEBUG
      std::cout << "Reading " << fileName << "...\n";
    #endif //NDEBUG

    //You could react to bad things happening within ROOT like this try-catch block.
    //I'm going to keep going with a warning message if I have a problem reading a file.
    //However, if there's a problem with xrootd, I'm going to bail immediately with a
    //special return code of -1.
    try
    {
      std::unique_ptr<TFile> input(TFile::Open(fileName.c_str()));
      //Normally, I'd check whether input is nullptr here, but ErrorHandler.h should cause an exception
      //to be thrown in all such cases that I know of.  These cases indicate that ErrorHandler.h isn't
      //working, but they're not actually the cases I'm testing for.
      assert(input != nullptr && "ErrorHandler.h failed to intercept a file that could not be opened!");

      TTree* tree = (TTree*)input->Get(argv[1]);
      //assert(tree != nullptr && "ErrorHandler.h failed to intercept a missing TTree!");

      //Huh, ROOT just happily returns a nullptr from TFile::Get() if the object doesn't
      //exist.  No Error()s or Warning()s at all.
      //This is a good opportunity to point out that you too can throw ROOT::errors.
      if(tree == nullptr)
      {
        throw ROOT::error(1, "TDirectory::Get", (std::string("No such object: ") + argv[1]).c_str());
      }

      //Mispell some standard branch that every AnaTuple has
      const char* someCommonBranch = argv[2];
      int eventID = -1;
      tree->SetBranchStatus("*", 0);
      tree->SetBranchStatus(someCommonBranch, 1);
      tree->Branch(someCommonBranch, &eventID);

      for(int whichEntry = 0; tree->GetEntry(whichEntry) > 0; ++whichEntry)
      {
        #ifndef NDEBUG
          if(whichEntry % 1000 == 0) std::cout << "Processing entry " << whichEntry << "...\n";
        #endif //NDEBUG
        //Your analysis code would normally go here
      }
    }
    catch(const ROOT::warning& warn)
    {
      std::cerr << "Got a warning from ROOT of severity " << warn.level << ": " << warn.what() << "\n";

      //Check err.location for xrootd-related problems
      if(warn.where.find("TNetXNGFile") != std::string::npos)
      {
        std::cerr << "This is an xrootd error which is really bad.  So, I'm bailing on the whole job!\n";
        return -1;
      }

      std::cerr << "Skipping " << fileName << "...\n";
      ++nSkippedFiles;
      continue;
    }
    catch(const ROOT::error& err)
    {
      std::cerr << "Got an error from ROOT of severity " << err.level << ": " << err.what() << "\n";

      //Check err.location for xrootd-related problems
      if(err.where.find("TNetXNGFile") != std::string::npos)
      {
        std::cerr << "This is an xrootd error which is really bad.  So, I'm bailing on the whole job!\n";
        return -1;
      }
                                                                                                        
      std::cerr << "This would normally crash a ROOT job, but I'm just going to skip " << fileName << "\n";
      ++nSkippedFiles;
      continue;
    }

    //Or, you could let the exceptions just crash your job for you by
    //doing nothing besides including ErrorHandler.h.  You could emulate
    //the normal behavior of ROOT by only catching ROOT::warnings and
    //omitting the ROOT::error catch block.
  }

  //Basically, don't use your output files if your analysis returns nonzero.
  //I'm returning the number of files I skipped.
  return nSkippedFiles;
}
