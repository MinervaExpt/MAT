#ifndef HISTWRAPPER_CXX
#define HISTWRAPPER_CXX

#include "PlotUtils/HistWrapper.h"

#include "PlotUtils/FluxSystematics.cxx"  // PlotUtils::flux_reweighter

using namespace PlotUtils;

// Default Constructor
template <typename T>
HistWrapper<T>::HistWrapper()
    : hist(nullptr), univToHistMap(), nhistsAssignedSoFar() {}

//! Constructor from a map of errors/universes, uniform bins
template <typename T>
HistWrapper<T>::HistWrapper(const char* hist_name, const char* title, int nBins,
                            double xmin, double xmax,
                            std::map<std::string, std::vector<T*> >& bands) {
  hist = new PlotUtils::MnvH1D(hist_name, title, nBins, xmin, xmax);
  hist->SetDirectory(0);

  for (typename std::map<std::string, std::vector<T*> >::const_iterator band =
           bands.begin();
       band != bands.end(); ++band) {
    std::vector<T*> univs = band->second;
    const std::string name(univs.front()->ShortName());

    FillErrorBandsWithSysUni(univs, univs.size(), name);
  }

  assert(hist);
}

//! Constructor from a map of errors/universes, variable bins
template <typename T>
HistWrapper<T>::HistWrapper(const char* hist_name, const char* title, int nBins,
                            std::vector<double> bins,
                            std::map<std::string, std::vector<T*> >& bands) {
  hist = new PlotUtils::MnvH1D(hist_name, title, nBins, &bins[0]);
  hist->SetDirectory(0);

  for (typename std::map<std::string, std::vector<T*> >::const_iterator band =
           bands.begin();
       band != bands.end(); ++band) {
    std::vector<T*> univs = band->second;
    const std::string name(univs.front()->ShortName());

    FillErrorBandsWithSysUni(univs, univs.size(), name);
  }

  assert(hist);
}

// Constructor from a pre-existing MnvH1D & map< string, vector<universes> >
template <typename T>
HistWrapper<T>::HistWrapper(MnvH1D* h1d,
                            std::map<std::string, std::vector<T*> >& bands,
                            bool clear_error_bands) {
  hist = new MnvH1D(*h1d);
  if (clear_error_bands) {
    if (hist->GetNErrorSources() > 0)
      printf(
          "HistWrapper: Clearing error bands from MnvH1D %s with existing "
          "error bands\n",
          hist->GetName());
    hist->ClearAllErrorBands();
  }
  hist->SetDirectory(0);

  for (typename std::map<std::string, std::vector<T*> >::const_iterator band =
           bands.begin();
       band != bands.end(); ++band) {
    std::vector<T*> univs = band->second;
    const std::string name(univs.front()->ShortName());

    FillErrorBandsWithSysUni(univs, univs.size(), name);
  }

  assert(hist);
}

// Constructor from unorganized vector<universe>, uniform binning
// Default 2 universes for any and all bands included.
template <typename T>
HistWrapper<T>::HistWrapper(const char* hist_name, const char* title, int nBins,
                            double xmin, double xmax, std::vector<T*>& univs) {
  hist = new PlotUtils::MnvH1D(hist_name, title, nBins, xmin, xmax);
  hist->SetDirectory(0);
  std::vector<T*> tmpunivs;

  for (unsigned int iUniv = 0; iUniv < univs.size(); ++iUniv) {
    tmpunivs.clear();
    T* univ = univs[iUniv];
    tmpunivs.push_back(univ);
    const std::string name(univ->ShortName());

    FillErrorBandsWithSysUni(tmpunivs, 2, name);
  }

  assert(hist);
}

// Constructor from unorganized vector<universe>, variable binning
// Default 2 universes for any and all bands included.
template <typename T>
HistWrapper<T>::HistWrapper(const char* hist_name, const char* title, int nBins,
                            std::vector<double> bins, std::vector<T*>& univs) {
  double* dbins = &bins[0];

  hist = new PlotUtils::MnvH1D(hist_name, title, nBins, dbins);
  hist->SetDirectory(0);
  std::vector<T*> tmpunivs;

  for (unsigned int iUniv = 0; iUniv < univs.size(); ++iUniv) {
    tmpunivs.clear();
    T* univ = univs[iUniv];
    tmpunivs.push_back(univ);
    const std::string name(univ->ShortName());

    FillErrorBandsWithSysUni(tmpunivs, 2, name);
  }

  assert(hist);
}

// Constructor from an MnvH1D & vector of universes
template <typename T>
HistWrapper<T>::HistWrapper(MnvH1D* h1d, std::vector<T*>& univs,
                            bool clear_error_bands) {
  hist = new MnvH1D(*h1d);
  if (clear_error_bands) {
    if (hist->GetNErrorSources() > 0)
      printf(
          "HistWrapper: Clearing error bands from MnvH1D %s with existing "
          "error bands\n",
          hist->GetName());
    hist->ClearAllErrorBands();
  }
  hist->SetDirectory(0);
  std::vector<T*> tmpunivs;

  for (unsigned int iUniv = 0; iUniv < univs.size(); ++iUniv) {
    tmpunivs.clear();
    T* univ = univs[iUniv];
    tmpunivs.push_back(univ);
    const std::string name(univ->ShortName());

    FillErrorBandsWithSysUni(tmpunivs, 2, name);
  }

  assert(hist);
}

//! Constructor from a histogram to add universes later, uniform bins. HMS
template <typename T>
HistWrapper<T>::HistWrapper(const char* hist_name, const char* title, int nBins,
                            std::vector<double> bins) {
  double* dbins = &bins[0];
  hist = new PlotUtils::MnvH1D(hist_name, title, nBins, dbins);
  hist->SetDirectory(0);
  assert(hist);
}

//! Constructor from a histogram to add universes later, variable bins. HMS
template <typename T>
HistWrapper<T>::HistWrapper(const char* hist_name, const char* title, int nBins,
                            double xmin, double xmax) {
  hist = new PlotUtils::MnvH1D(hist_name, title, nBins, xmin, xmax);
  hist->SetDirectory(0);
  assert(hist);
}

// The next 2 constructors are hidden from our python bindings because they use
// a c++11 feature: delegating constructors
#ifndef __GCCXML__
// Construct from a map of universes and variable-sized bins.
// TODO: If we could stand to make a breaking change by changing the order of
// HistWrapper<>
//       constructor arguments, I think we could get rid of all but 1 of the
//       above constructors.
template <class T>
HistWrapper<T>::HistWrapper(const char* name, const char* title,
                            const std::vector<double> binLowEdges,
                            std::map<std::string, std::vector<T*> >& bands)
    : HistWrapper(
          new MnvH1D(name, title, binLowEdges.size() - 1, binLowEdges.data()),
          bands) {}

// Construct from a vector of universes and variable-sized bins.
template <class T>
HistWrapper<T>::HistWrapper(const char* name, const char* title,
                            const std::vector<double> binLowEdges,
                            std::vector<T*>& univs)
    : HistWrapper(
          new MnvH1D(name, title, binLowEdges.size() - 1, binLowEdges.data()),
          univs) {}
#endif  //__GCCXML__

// Copy an MnvH1D's CV histo to each of its vertical error band CV histos
template <typename T>
void HistWrapper<T>::SyncCVHistos() {
  TH1D* theCVHisto = new TH1D(*hist);
  theCVHisto->SetDirectory(0);
  std::vector<std::string> bandnames = hist->GetErrorBandNames();
  // std::cout << "Synching Error Band CV's with MnvH1D's CV" << std::endl;
  // loop this MnvH1D's bands
  for (std::vector<std::string>::const_iterator bandname = bandnames.begin();
       bandname != bandnames.end(); ++bandname) {
    // std::cout << *bandname << std::endl;
    // band is a reference to the error band's CV histo
    PlotUtils::MnvVertErrorBand& band =
        *(hist->GetVertErrorBand((*bandname).c_str()));
    // set the BAND's CV histo to the MnvH1D's CV histo
    band.TH1D::operator=(*theCVHisto);
  }
}

// Fill
template <typename T>
void HistWrapper<T>::FillUniverse(const T& univ, const double value,
                                  const double weight) {
  try {
    univHist(&univ)->Fill(value, weight);
  } catch (std::out_of_range) {
    std::cerr << "From HistWrapper::FillUniverse\n";
    throw;
  }
}

template <typename T>
void HistWrapper<T>::FillUniverse(const T* univ, const double value,
                                  const double weight) {
  FillUniverse(*univ, value, weight);
}

// Access to the TH1 corresponding to the given universe
template <typename T>
TH1D* HistWrapper<T>::univHist(const T* univ) const {
  try {
    return univToHistMap.at(univ);
  } catch (std::out_of_range) {
    assert(hist);
    std::cerr << "HistWrapper::univHist out_of_range error.\n"
              << "The universe passed here to univHist does not exist in this "
              << "HistWrapper (" << hist->GetName() << ").\n"
              << "The universe in question is " << univ->ShortName()
              << " with m_nsigma " << univ->GetSigma() << ".\nMake sure this "
              << "exact universe was used to construct this HistWrapper.\n";
    throw;
  }
}

//! add a universe to the list
template <typename T>
void HistWrapper<T>::AddUniverses(const std::string name, const T* univ,
                                  const int nhists, const int histno) {
  assert(histno < nhists);
  // Special-case for central value
  if (name == "cv") {
    if (nhists != 1) {
      std::cerr << "HistWrapper constructor ERROR: CV should not have "
                   "more than one universe!"
                << std::endl;
      std::exit(2);
    }

    univToHistMap[univ] = hist;
    return;
  }
  std::cout << "try to add a " << univ->ShortName() << " as " << name
            << std::endl;

  // add error band if not already there
  if (!hist->HasVertErrorBand(name)) {
    std::cerr << "adding vert error band with " << name << " " << nhists
              << std::endl;
    hist->AddVertErrorBand(name, nhists);
  }
  // Connect each band's universe to the corresponding MnvH1D's TH1
  univToHistMap[univ] = hist->GetVertErrorBand(name)->GetHist(histno);

  assert(hist);
}

// Private, used by CTORs
template <typename T>
void HistWrapper<T>::FillErrorBandsWithSysUni(std::vector<T*> univs, int nhists,
                                              std::string name) {
  // Special-case for central value
  if (name == "cv") {
    if (univs.size() != 1) {
      std::cerr << "HistWrapper constructor ERROR: CV should not have "
                   "more than one universe!"
                << std::endl;
      std::exit(2);
    }
    T* cv = univs.front();
    univToHistMap[cv] = hist;
    return;
  }

  // Add this band to the MnvH1D
  if (!hist->HasVertErrorBand(name)) {
    if (univs.front()->ShortName() == "Flux" &&
        univs.front()->UseNuEConstraint()) {
      const int nflux_universes = univs.front()->GetNFluxUniverses();
      if(nhists != nflux_universes) {
        std::cout << "WARNING from HistWrapper::FillErrorBandsWithSysUni\n" 
            << "  You're attempting to make a HistWrapper with " << nhists 
            << " flux universes\n  but this is at odds with the number of flux "
            << "universes you have\n  associated with your CVUniverse "
            << "(GetNFluxUniverses) " << nflux_universes  << ".\n  Your HW will "
            << "be constructed with " << nflux_universes << " universes.\n  "
            << "Change this with MinervaUniverse::SetNFluxUniverses.\n";
      }
      PlotUtils::flux_reweighter(univs.front()->GetPlaylist(),
                                 univs.front()->GetAnalysisNuPDG(), true,
                                 nflux_universes)
          .AddFluxErrorBand(hist);
    } else
      hist->AddVertErrorBand(name, nhists);
  }

  // Connect each band's universe to the corresponding MnvH1D's TH1
  for (unsigned int iUniv = 0; iUniv < univs.size(); ++iUniv) {
    T* univ = univs[iUniv];
    if (univ->ShortName() != name) {
      std::cerr << "HistWrapper constructor ERROR: faulty construction "
                   "of systematic bands map. All of a bands' universes "
                   "should have same name."
                << std::endl;
      std::exit(2);
    }
    int& nAssigned = nhistsAssignedSoFar[name];
    univToHistMap[univ] = hist->GetVertErrorBand(name)->GetHist(nAssigned);
    nAssigned++;
  }
}

#endif  // HISTWRAPPER_CXX
