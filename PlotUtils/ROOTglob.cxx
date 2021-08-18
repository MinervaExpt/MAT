//File: ROOTglob.cxx
//Brief: Use regular expressions, like */*.root, to generate a list of files for
//       any kind of filesystem ROOT supports.  This includes streaming files over
//       xrootd.  Unlike TChain::Add(), this function does a "recursive glob".
//Author: Andrew Olivier aolivier@ur.rochester.edu

//PlotUtils includes
#include "PlotUtils/ROOTglob.h"

//ROOT includes
#include "TRegexp.h"
#include "TSystem.h"
#include "TString.h"

//I'm only going to print things out in debug builds
#ifndef NDEBUG
#include <iostream>
#endif

//c++ includes
#include <algorithm>

namespace PlotUtils
{
  //Returns: A list of files that matched the regular expression regex.  If
  //         os is a TNetXNGSystem or something more exotic that uses special
  //         URLs, these should be appropriate URLs for use with TFile::Open().
  //
  //Parameters:
  //path: Regular expression for files to find.  Could be an absolute or a
  //      relative path.  Could even be just a single file.
  //
  //os:   (Polymorphic) reference to a TSystem object.  You probably want to use
  //      *gSystem here most of the time.  Provides an interface to the operating system
  std::vector<std::string> glob(const std::string& path, TSystem& os)
  {
    std::vector<std::string> result;
    std::vector<std::string> dirNameStack; //Names of directories to traverse
    const size_t beginHostPos = path.find("//") + 2;
    const size_t beginPathPos = path.find("//", beginHostPos) + 2;
    const size_t firstRegexChar = path.find_first_of("^$.[]*+?", beginPathPos);
    dirNameStack.push_back(path.substr(0, path.rfind("/", firstRegexChar)));
                                                                                                                                         
    //Translation of recursive glob into a FIFO-based loop
    while(!dirNameStack.empty())
    {
      std::string currentDir = dirNameStack.back();
      dirNameStack.pop_back();

      //There might have been expanded regular expression matches in path to get to currentDir.
      //I never cross /s in regular expressions, so count /s to figure out where to start looking
      //in path.
      const size_t howManySlashes = std::count(currentDir.begin(), currentDir.end(), '/');
      size_t startPos = 0;
      for(size_t nSlashes = 0; startPos < path.length() && nSlashes <= howManySlashes; ++nSlashes) startPos = path.find("/", startPos) + 1;

      //I don't want to look for a specific subdirectory in each directory.  xrootd directory queries are probably pretty expensive, and
      //that would be BORING besides.  Skip to the next regular expression candidate.  Define the regular expression itself as the
      //thing bracketed by /s at the next regex character.
      const size_t nextRegexChar = path.find_first_of("^$.[]*+?", startPos);
      const size_t lowerSlash = path.rfind("/", nextRegexChar) + 1;
      const size_t upperSlash = path.find("/", nextRegexChar);
      currentDir += path.substr(startPos-1, lowerSlash - startPos); //At the end of the day, THIS is the next directory I really need to do regex matches in.
      const bool recurse = (upperSlash != std::string::npos); //Don't keep going forever when I get a path that ends with a regular expression!

      const TRegexp regex(path.substr(lowerSlash, upperSlash - lowerSlash).c_str(), kTRUE);
      void* dir = os.OpenDirectory(currentDir.c_str());
      if(dir)
      {
        const char* entryName;
        while((entryName = os.GetDirEntry(dir))) //Loop over entries in dir
        {
          if(!strcmp(entryName, ".") || !strcmp(entryName, "..")) continue; //Avoid endless recursion
                                                                                                                                         
          TString entryString(entryName); //entryName seems to be a path relative to dir
          if(entryString.Index(regex) >= 0)
          {
            if(entryString.EndsWith(".root")) result.push_back(currentDir + "/" + entryName);
            else if(recurse) dirNameStack.push_back(currentDir + "/" + entryName);
          } //If matches regular expression
        } //While there's another directory entry to view
                                                                                                                                         
        os.FreeDirectory(dir);
      } //If could open directory
      #ifndef NDEBUG //Only print out messages about missing directories in debug builds
      else
      {
        std::cerr << "PlotUtils::glob(): Failed to open a directory named " << currentDir << " that was there a minute ago!  "
                  << "Was it moved or deleted while this program was running?  This is just a debugging message, so "
                  << "I'll just ignore this directory and move on..." << std::endl;
      }
      #endif
    } //While there are more directories to process

    return result;
  }
}
