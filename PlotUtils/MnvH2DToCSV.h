#ifndef MNVH2DTOCSV_H


#include "PlotUtils/MnvH2D.h"
#include <string>

// functions to dump histograms as csv files - from Cheryl Patrick


namespace PlotUtils{
  
  void MnvH2DToCSV(PlotUtils::MnvH2D *hist, std::string name, std::string directory, double scale, bool fullprecision=true, bool syserrors=true);

}

#endif
