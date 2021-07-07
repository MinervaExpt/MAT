//File: DoesChainWrapperWorkWithXROOTD.cpp
//Brief: Test whether PlotUtils::ChainWrapper works with regular expressions and xrootd.
//       Can either be built into an executable or interpretted as a ROOT macro.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//c++ includes
#ifndef __CINT__
#include <iostream>
#endif

//PlotUtils includes
#ifndef __CINT__
#include "PlotUtils/ChainWrapper.h"
#endif

//This function does all of the work.
bool DoesChainWrapperWorkWithXROOTD(const std::string& regularExpression = "") //TODO: Put a file here that will probably "always" exist
{
  //Necessary for incomplete Reflex support in ROOT version < 6
  //TODO: ifndefs for ROOT version
  gSystem->Load( "libCintex.so" );
  ROOT::Cintex::Cintex::Enable();

  PlotUtils::ChainWrapper chw("Truth"); //The Truth tree should be in any MINERvA AnaTool output file.
  chw.AddFiles(regularExpression.c_str());
  if(chw.GetEntries() <= 0)
  {
    std::cerr << "chw has no entries after trying construct it with a file!\n";
    return false;
  }

  TCanvas allMine("DoesChainWrapperWorkWithXROOTD");
  chw.GetTree()->Draw("mc_vtx[0]:mc_vtx[1]", "", "colz"); //A nice histogram that should work with any MINERvA AnaTuple.
  allMine.Print("DoesChainWrapperWorkWithXROOTD.png");

  return true;
}

//Run the above ROOT macro in an executable instead.  Returns to the OS:
//  0 if ChainWrapper works with xrootd
//  1 if I get to the point of setting up ChainWrapper but can't Add() any files
//  2 if command line arguments are incorrect
//Using return codes should make this easy to incorporate in a test.
int main(const int argc, const char** argv)
{
  if(argc != 2)
  {
    std::cerr << "Used wrong number of arguments to ChainWrapper test!\n"
              << "Usage: DoesChainWrapperWorkWithXROOTD <regular expression to test>\n";
    return 2; //Returning 2 indicates a failure to parse arguments.
  }

  return DoesChainWrapperWorkWithXROOTD(argv[1]);
}
