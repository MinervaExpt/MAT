//=============================================================================
/*! @brief A wrapper around MnvH1Ds which allows one to interface more easily,
           intuitively, and directly with error bands and their constituent
           universes.

This class allows one to more easily and intuitively interface with an
MnvH1D's systematic error bands and their constituent universes.

It was built to eliminate the distinction between lateral and vertical error
bands, eliminate FillVert/LatErrorBand functions, and standardize as many
MINERvA errors as possible, and linearize the process of looping events and
filling error bands.

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

See a complete example
[here](http://cdcvs0.fnal.gov/cgi-bin/public-cvs/cvsweb-public.cgi/Personal/bmesserl/MasterAnaMacro/?cvsroot=mnvsoft)

@author Ben Messerly
*/
//=============================================================================

#ifndef HISTWRAPPER_H
#define HISTWRAPPER_H

#include <map>
#include <vector>

#include "MnvH1D.h"

namespace PlotUtils {

template <typename T>
struct HistWrapper {
  // Default Constructor
  HistWrapper();

  //! Constructor from a map of errors/universes, uniform bins
  HistWrapper(const char* hist_name, const char* title, int nBins, double xmin,
              double xmax, std::map<std::string, std::vector<T*> >& bands);

  //! Constructor from a map of errors/universes, variable bins
  HistWrapper(const char* hist_name, const char* title, int nBins,
              std::vector<double>,
              std::map<std::string, std::vector<T*> >& bands);

  //! Constructor from a vector of universes, uniform binning
  HistWrapper(const char* hist_name, const char* title, int nBins, double xmin,
              double xmax, std::vector<T*>& univs);

  //! Constructor from a vector of universes, variable binning
  HistWrapper(const char* hist_name, const char* title, int nbins,
              std::vector<double> bins, std::vector<T*>& univs);

  //! Constructor from a pre-existing MnvH1D and a vector of universes
  HistWrapper(MnvH1D* h1d, std::vector<T*>& univs,
              bool clear_error_bands = false);

  //! Constructor from a pre-existing MnvH1D and a map of errors/universes
  HistWrapper(MnvH1D* h1d, std::map<std::string, std::vector<T*> >& bands,
              bool clear_error_bands = false);

// The next 2 constructors are hidden from our python bindings because they use
// a c++11 feature for code reuse: delegating constructors.
#ifndef __GCCXML__
  //! Constructor from a map of errors/universes and variable-sized bins
  HistWrapper(const char* name, const char* title,
              const std::vector<double> binLowEdges,
              std::map<std::string, std::vector<T*> >& bands);

  //! Constructor from a vector of universes and variable-sized bins
  HistWrapper(const char* name, const char* title,
              const std::vector<double> binLowEdges, std::vector<T*>& univs);
#endif  //__gccxml__

  //! Constructor from a histogram to add universes later, uniform bins. HMS
  HistWrapper(const char* hist_name, const char* title, int nBins, double xmin,
              double xmax);

  //! Constructor from a histogram to add universes later, variable bins. HMS
  HistWrapper(const char* hist_name, const char* title, int nBins,
              std::vector<double> bins);

  //! Synchronize the MnvH1D's CV with each of its error band's CVs.
  void SyncCVHistos();

  //! Access universe TH1 given a universe object
  TH1D* univHist(const T* univ) const;

  //! Fill universe TH1
  void FillUniverse(const T& univ, double value, double weight);

  //! Fill universe TH1
  void FillUniverse(const T* univ, double value, double weight);

  //! add a universe to the list
  void AddUniverses(const std::string name, const T* univ, const int nhists,
                    const int i);

  // Data members
  PlotUtils::MnvH1D* hist; /*!< The MnvH1D that is created and filled */
  std::map<const T*, TH1D*> univToHistMap; /*!< The map between univs and TH1 */
  std::map<std::string, int>
      nhistsAssignedSoFar; /*!< Counter for nunivs in each band.*/

 private:
  // Add vertical error bands to the MnvH1D and associate them with universes
  void FillErrorBandsWithSysUni(std::vector<T*> univs, int nhists,
                                std::string name);
};

}  // namespace PlotUtils

// Template classes must have function definitions available in header.
#include "HistWrapper.cxx"

#endif  // HISTWRAPPER_H
