//==============================================================================
// A HistFolio is a smart histogram container
//
// Component histograms can be accessed by name (string) or by a any type you
// choose, for example an enum corresponding to signal and background
// categories.
//
// This is a template class which at minimum needs to know what kind of
// underlying histogram object you would like to use, e.g. MnvH2D, TH1,
// HistWrapper, etc.
//
// TODO Capacity to select subsections of the folio for plotting
// TODO More accessors, probably.
//
// Example script showing how to use this can be found at:
// Personal/bmesserl/MasterAnaMacro/runEventLoop_HistFolio.C
//==============================================================================
#ifndef HISTFOLIO_H
#define HISTFOLIO_H

#include "MnvColors.h"
#include "MnvH1D.h"
#include "NamedCategory.h"
#include "TFile.h"
#include "TObjArray.h"

namespace PlotUtils {

template <class HIST = PlotUtils::MnvH1D, class CATEGORY = int>
class HistFolio {
 public:
  //============================================================================
  // CTORs
  //============================================================================
  // CTOR -- only binning
  template <class... HISTARGS>
  HistFolio(std::string folio_name, HISTARGS... hist_args);

  // CTOR -- source/template hist
  HistFolio(std::string folio_name, HIST* hist);

  // CTOR -- copy/clone
  HistFolio(const HistFolio& rhs, std::string new_name);

  // CTOR -- binning and categories (strings)
  template <class... HISTARGS>
  HistFolio(const std::vector<std::string>& names, std::string folio_name,
            HISTARGS... hist_args);
 
  // CTOR -- binning and categories (strings) vector 
  template <class... HISTARGS>  
  HistFolio(const std::vector<std::string>& names, std::string folio_name,
            std::vector<double> Bins, HISTARGS... hist_args);

  // CTOR -- binning and categories (NamedCategories)
  template <class... HISTARGS>
  HistFolio(const std::vector<NamedCategory<CATEGORY>>& categories,
            std::string folio_name, const HISTARGS... hist_args);
  
  // 1D hist CTOR -- Vector binning and categories (NamedCategories)
  template <class... HISTARGS>
  HistFolio(const std::vector<NamedCategory<CATEGORY>>& categories,
            std::string folio_name, std::vector<double> Bins, const HISTARGS... hist_args);

  // 2D hist CTOR -- Vector binning and categories (NamedCategories)
  template <class... HISTARGS>
   HistFolio(const std::vector<NamedCategory<CATEGORY>>& categories,
            std::string folio_name, std::vector<double> Bins1, std::vector<double> Bins2, const HISTARGS... hist_args);


  // CTOR -- default.
  HistFolio();

  //============================================================================
  // Add hists to folio ...
  //============================================================================
  // ... from a string
  void AddComponentHist(const std::string category_name);

  // ... from category
  void AddComponentHist(const NamedCategory<CATEGORY> category,
                        int n_hists_so_far = -99);

  // ... from a hist
  void AddComponentHist(const std::string category_name, HIST* component_hist);

  //============================================================================
  // Setters and Getters
  //============================================================================
  // Number of component histograms
  int GetSize() const;

  // Name of folio
  std::string GetName() const;
  // std::string GetTitle() const;
  // void SetTitle(const std::string);
  // void SetName(const std::string);

  // Access component histogram. T = string or CATEGORY.
  template <class T>
  HIST* GetComponentHist(const T key, bool be_quiet = true);

  // Useful for MnvPlotter::DrawDataStackedMC
  TObjArray GetHistArray() const;

  std::map<NamedCategory<CATEGORY>, HIST*> GetHistMap() const;

  //============================================================================
  // Colors
  //============================================================================
  void ApplyColorPalette(
      MnvColors::EColorPalettes palette = MnvColors::kOkabeItoPalette);
  void ApplyColorPalette(std::vector<int> palette =
                             MnvColors::GetColors(MnvColors::kOkabeItoPalette));

  //============================================================================
  // Write and save to file
  //============================================================================
  void WriteToFile(TFile& f) const;

 private:
  //============================================================================
  // Data Members
  //============================================================================
  mutable std::string m_folio_name;
  HIST* m_source_hist;
  std::map<NamedCategory<CATEGORY>, HIST*> m_hist_map;
  TObjArray m_hist_array;

  //============================================================================
  // Component Histogram Initialization
  //============================================================================
  // Set axis labels, titles, colors, and add to member containers
  void InitializeComponentHist(NamedCategory<CATEGORY> category,
                               HIST* component_hist, const int n_hists_so_far);

  //============================================================================
  // Histogram access helpers
  //============================================================================
  // Make sure new component additions aren't duplicates
  bool CheckForDuplicates(const NamedCategory<CATEGORY> category);

  // Search for histo in map by string
  typename std::map<NamedCategory<CATEGORY>, HIST*>::const_iterator
  SearchMapForHist(const std::string name) const;

  // Search for histo in map by a CATEGORY
  typename std::map<NamedCategory<CATEGORY>, HIST*>::const_iterator
  SearchMapForHist(const CATEGORY value) const;

  // Helper predicates to search the hist map
  struct by_name {
   public:
    by_name(std::string name);
    bool operator()(
        const std::pair<NamedCategory<CATEGORY>, HIST*>& element) const;

   private:
    const std::string m_name;
  };

  struct by_value {
   public:
    by_value(CATEGORY value);
    bool operator()(
        const std::pair<NamedCategory<CATEGORY>, HIST*>& element) const;

   private:
    const CATEGORY m_value;
  };
};

//==============================================================================
// Load HistFolio from file (not a member function)
//==============================================================================
template <class HIST>
PlotUtils::HistFolio<HIST, int> LoadHistFolioFromFile(TFile& f,
                                                      std::string folio_name);

// Helper to LoadHistFolioFromFile
template <class HIST>
HIST* FindSourceHist(TFile& f, std::string folio_name);

}  // namespace PlotUtils

#include "HistFolio.cxx"

#endif  // HISTFOLIO_H
