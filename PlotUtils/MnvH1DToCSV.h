#ifndef MNV1HDTOCSV_H

#include "MnvH1D.h"

#include <string>

// functions to dump histograms as csv files - from Cheryl Patrick


namespace MAT{
  

  void MnvH1DToCSV(MAT::MnvH1D *hist, std::string name, std::string directory, double scale, bool fullprecision=true, bool errors=true);
}

#endif
