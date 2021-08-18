#ifndef MNV1HDTOCSV_H

#include "PlotUtils/MnvH1D.h"

#include <string>

// functions to dump histograms as csv files - from Cheryl Patrick


namespace PlotUtils{
  

  void MnvH1DToCSV(PlotUtils::MnvH1D *hist, std::string name, std::string directory, double scale, bool fullprecision=true, bool errors=true);
}

#endif
