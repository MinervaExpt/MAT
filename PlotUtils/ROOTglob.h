//File: ROOTglob.h
//Brief: Use regular expressions, like */*.root, to generate a list of files for
//       any kind of filesystem ROOT supports.  This includes streaming files over
//       xrootd.  Unlike TChain::Add(), this function does a "recursive glob".
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef PLOTUTILS_ROOTGLOB_H
#define PLOTUTILS_ROOTGLOB_H

//c++ includes
#include <vector>
#include <string>

//Forward declarations for ROOT classes
class TSystem;

namespace MAT
{
  //Returns: A list of files that matched the regular expression regex.  If
  //         os is a TNetXNGSystem or something more exotic that uses special
  //         URLs, these should be appropriate URLs for use with TFile::Open().
  //
  //Parameters:
  //regex: Regular expression for files to find.  Could be an absolute or a
  //       relative path.  Could even be just a single file.
  //
  //os:    (Polymorphic) reference to a TSystem object.  You probably want to use
  //       *gSystem here most of the time.  Provides an interface to the operating system
  std::vector<std::string> glob(const std::string& regex, TSystem& os);
}

#endif //PLOTUTILS_ROOTGLOB_H
