#ifndef MNV_ANABINNING_h
#define MNV_ANABINNING_h 1

#include <vector>

namespace PlotUtils{

  struct axis_binning{
    std::vector<double> bin_edges; 
    unsigned int nbins;
    double min; double max; 
    bool uniform;
    double default_width;
  };

  class AnaBinning {

    public:

      //! Default Constructor
      AnaBinning(){};

      //! singleton gettor
      static AnaBinning& Get();

      axis_binning GetProtonEnergyBinsGeV()const;

      axis_binning GetBlobEnergyBinsGeV()const;
      axis_binning GetBlobEnergyBinsMeV()const;

      axis_binning GetMuonEnergyBinsGeV()const;
      axis_binning GetMuonEnergyUniformBinsGeV()const;
      axis_binning GetMuonCosineAngleBins()const;
      axis_binning GetMuonAngleBinsDegrees()const;

      axis_binning GetThetaXBinsDegrees()const;
      axis_binning GetThetaYBinsDegrees()const;
      axis_binning GetPhiBinsDegrees()const;

      axis_binning GetQ2BinsGeV()const;
      axis_binning GetQ2UniformBinsGeV()const;
      axis_binning GetQ2FineBinsGeV()const;
      axis_binning GetQ2FineBinsGeV_05()const;

      axis_binning GetNeutrinoEnergyBinsGeV()const;
      axis_binning GetNeutrinoEnergyUniformBinsGeV()const;

      axis_binning GetVertexR2Bins()const;
      axis_binning GetVertexXYBins()const;
      axis_binning GetVertexZBins()const;
      axis_binning GetVertexZBinsModule()const;
      axis_binning GetVertexTBins()const;
      axis_binning GetVertexTBinsZoom()const;

      axis_binning GetMinosVertexXYBins()const;
      axis_binning GetMinosQPBins()const;

      //! Add to bin_edges n_bins bins each with width bin_width.  Add a bin at 0 if bin_edges is empty.
      void AddBins( std::vector<double>& bin_edges, const int n_bins, const double bin_width ) const;

      //! Add to bin_edges bins of width bin_width, from lowBin to high_bin (inclusive)
      void AddBins( std::vector<double>& bin_edges, const double bin_width, const double low_bin, const double high_bin ) const;

  }; //end of AnaBinning

} //end of PlotUtils

#endif
