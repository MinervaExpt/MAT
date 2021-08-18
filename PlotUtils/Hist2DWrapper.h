//=============================================================================
/*!
  @brief A wrapper around MnvH1Ds which allows one to interface more easily,
         intuitively, and directly with error bands and their constituent 
         universes.

  This class allows one to more easily and intuitively interface with an
  MnvH1D's systematic error bands and their constituent universes. 

  It was built to eliminate the distinction between lateral and vertical error
  bands, eliminate FillVert/LatErrorBand functions, and standardize as many
  MINERvA errors as possible, and linearize the process of looping events and
  filling error bands.

  e.g. To create an MnvH1D and fill it with the full set of standard genie
  systematics:

      std::map<std::string, std::vector<CVUniverse*>> error_bands = GetGenieSystematics<CVUniverse>(my_event_chain);
      PlotUtils::HistWrapper<CVUniverse> my_hw("E_{#nu}", nbins, xmin, xmax, error_bands);
      loop my_event_chain 
        loop band in error_bands 
          loop universe in band
            if(PassesCuts(universe))
              my_hw.univHist(universe)->Fill(universe->GetEnu(), universe->GetWeight);

  Pass a container of universe objects, and the constructor will make an
  MnvH1D and populate it with MnvVertErrorBands.
  Then, in your loop over events and universes, HistWrapper::univHist accesses
  the correct TH1 so you can fill it directly.

  The user must write her own CVUniverse class, which can be very simple,
  containing only these members:

  - virtual std::string ShortName() const = 0;
  - virtual std::string LatexName() const = 0;
  - virtual double GetWeight() const;
  - PlotUtils::ChainWrapper& m_chw;
  - double m_nsigma;
  - Long64_t m_entry;

  See a complete example [here](http://cdcvs0.fnal.gov/cgi-bin/public-cvs/cvsweb-public.cgi/Personal/bmesserl/SystematicsFramework/?cvsroot=mnvsoft)

  @author Ben Messerly
*/
//=============================================================================

#ifndef HIST2DWRAPPER_H
#define HIST2DWRAPPER_H

#include "PlotUtils/MnvH2D.h"
#include <vector>
#include <map>

namespace PlotUtils{
  template <typename T>
  struct Hist2DWrapper
  {
    // Default Constructor
    Hist2DWrapper();

    //! Constructor from a vector of universes
    Hist2DWrapper(const char* hist_name, const char* title, int nBinsX, double xmin, double xmax, int nBinsY, double ymin, double ymax, std::vector<T*>& univs);

    //! Constructor from a pre-existing MnvH2D and a vector of universes
    Hist2DWrapper(MnvH2D* h2d, std::vector<T*>& univs, bool clear_error_bands = false );

    //! Constructor from a map of errors/universes
    Hist2DWrapper(const char* hist_name, const char* title, int nBins, double xmin, double xmax, int nBinsY, double ymin, double ymax, std::map< std::string, std::vector<T*> >& bands);

    //! Constructor from a pre-existing MnvH2D and a map of errors/universes
    Hist2DWrapper(MnvH2D* h2d, std::map< std::string, std::vector<T*> >& bands, bool clear_error_bands = false);

    // The next 2 constructors are hidden from our python bindings because they use a c++11 feature: delegating constructors
    #ifndef __GCCXML__
    //! Constructor from a vector of universes and variable bin widths
    Hist2DWrapper(const char* name, const char* title, const std::vector<double> xBins, const std::vector<double> yBins, std::vector<T*>& univs);

    //! Constructor from a map of universes and variable bin widths
    Hist2DWrapper(const char* name, const char* title, const std::vector<double> xBins, const std::vector<double> yBins, std::map< std::string, std::vector<T*> >& bands);
    #endif //__GCCXML__

    // Data members
    PlotUtils::MnvH2D* hist;                        /*!< The MnvH2D that is created and filled */
    std::map<const T*, TH2D*> univToHistMap;        /*!< The map between univs and TH2 */ 
    std::map<std::string, int> nhistsAssignedSoFar; /*!< Counter for nunivs in each band.*/ 

    //! Synchronize the MnvH2D's CV with each of its error band's CVs.
    void SyncCVHistos();

    //! Access universe TH2 given a universe object
    TH2D* univHist(const T* univ) const;

    //! Fill universe TH2
    void FillUniverse(const T& univ, double valueX, double valueY, double weight);

    //! Fill universe TH2
    void FillUniverse(const T* univ, double valueX, double valueY, double weight);

    private:
      // Add vertical error bands to the MnvH2D and associate them with universes
      void FillErrorBandsWithSysUni( std::vector<T*> univs, int nhists, std::string name );  
  };
}

// Template classes must have function definitions available in header.
#include "Hist2DWrapper.cxx"

#endif // HISTWRAPPER_H
