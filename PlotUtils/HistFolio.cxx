#ifndef HISTFOLIO_CXX
#define HISTFOLIO_CXX

#include "HistFolio.h"

#include <cassert>
#include <exception>

#include "TKey.h"  // For loading HistFolio from file

using namespace MAT;

//============================================================================
//
// CTORs
//
//============================================================================
// CTOR -- only binning
template <class HIST, class CATEGORY>
template <class... HISTARGS>
HistFolio<HIST, CATEGORY>::HistFolio(std::string folio_name,
                                     HISTARGS... hist_args)
    : m_folio_name(folio_name),
      m_source_hist(new HIST(folio_name.c_str(), hist_args...)) {
  assert(
      ("ERROR: If you do not specify a NamedCategory list in this CTOR"
       "then you cannot specify a template argument. Do "
       "HistFolio<MAT::MnvH1D>, or provide a category list.",
       std::is_same<int, CATEGORY>::value));
}



// CTOR -- source/template histogram
template <class HIST, class CATEGORY>
HistFolio<HIST, CATEGORY>::HistFolio(std::string folio_name, HIST* hist)
    : m_folio_name(folio_name), m_source_hist(hist->Clone(folio_name.c_str())) {
  m_source_hist->SetTitle(m_folio_name.c_str());
}

// CTOR -- copy/clone
template <class HIST, class CATEGORY>
HistFolio<HIST, CATEGORY>::HistFolio(const HistFolio& rhs,
                                     std::string new_name) {
  HistFolio<HIST, CATEGORY> lhs(rhs);  // Default
  lhs.m_folio_name = new_name;
  for (auto h : lhs.m_hist_map) {
    std::string component_hist_name = lhs.m_folio_name + "_" + h.first.m_name;
    h.second->SetName(component_hist_name.c_str());
  }
  m_folio_name = lhs.m_folio_name;
  m_source_hist = lhs.m_source_hist;
  m_hist_map = lhs.m_hist_map;
  m_hist_array = lhs.m_hist_array;
}

// CTOR -- binning and categories (strings)
template <class HIST, class CATEGORY>
template <class... HISTARGS>
HistFolio<HIST, CATEGORY>::HistFolio(const std::vector<std::string>& names,
                                     std::string folio_name,
                                     HISTARGS... hist_args)
    : m_folio_name(folio_name),
      m_source_hist(new HIST(folio_name.c_str(), "source_hist", hist_args...)) {
  for (const auto& name : names) AddComponentHist(name);
}

//CTOR binning and categorier intialized with vector 
template <class HIST, class CATEGORY>
template <class... HISTARGS>
HistFolio<HIST, CATEGORY>::HistFolio(const std::vector<std::string>& names,
                                     std::string folio_name,std::vector<double> Bins, HISTARGS... hist_args)
    : m_folio_name(folio_name),
      m_source_hist(new HIST(folio_name.c_str(), hist_args...  ,Bins.size()-1, Bins.data() )) {
  for (const auto& name : names) AddComponentHist(name);
}



// CTOR -- binning and categories (NamedCategories)
template <class HIST, class CATEGORY>
template <class... HISTARGS>
HistFolio<HIST, CATEGORY>::HistFolio(
    const std::vector<NamedCategory<CATEGORY>>& categories,
    std::string folio_name, const HISTARGS... hist_args)
    : m_folio_name(folio_name),
      m_source_hist(new HIST(folio_name.c_str(), "source_hist", hist_args...)) {
  for (const auto& cat : categories) AddComponentHist(cat);
}

//CTOR -- same as above but intialized with  vector 
template <class HIST, class CATEGORY>
template <class... HISTARGS>
HistFolio<HIST, CATEGORY>::HistFolio(  const std::vector<NamedCategory<CATEGORY>>& categories,
                                       std::string folio_name, 
                                       std::vector<double> Bins, 
                                      const HISTARGS... hist_args)
    : m_folio_name(folio_name),
      m_source_hist(new HIST(folio_name.c_str(), hist_args..., Bins.size()-1, Bins.data() )) {
  for (const auto& cat : categories) AddComponentHist(cat);
}


//for 2D His with vectorBins
template <class HIST, class CATEGORY>
template <class... HISTARGS>
HistFolio<HIST, CATEGORY>::HistFolio(  const std::vector<NamedCategory<CATEGORY>>& categories,
                                       std::string folio_name,
                                       std::vector<double> Bins1,
                                       std::vector<double> Bins2,
                                      const HISTARGS... hist_args)
    : m_folio_name(folio_name),
      m_source_hist(new HIST(folio_name.c_str(), hist_args..., Bins1.size()-1, Bins1.data(),Bins2.size()-1,Bins2.data() )) {
  for (const auto& cat : categories) AddComponentHist(cat);
}



// CTOR -- Default
template <class HIST, class CATEGORY>
HistFolio<HIST, CATEGORY>::HistFolio(){};

//============================================================================
//
// ADD HISTS TO FOLIO (PUBLIC)
//
//============================================================================
// ... from a string
template <class HIST, class CATEGORY>
void HistFolio<HIST, CATEGORY>::AddComponentHist(const std::string name) {
  assert(
      ("ERROR: this specific function is not supported for HistFolio if the "
       "CATEGORY is not int (default).",
       std::is_same<int, CATEGORY>::value));
  int n_hists_so_far = GetSize();
  NamedCategory<CATEGORY> category(n_hists_so_far, name.c_str());
  AddComponentHist(category, n_hists_so_far);
}

// ... from category
template <class HIST, class CATEGORY>
void HistFolio<HIST, CATEGORY>::AddComponentHist(
    const NamedCategory<CATEGORY> category, int n_hists_so_far) {
  CheckForDuplicates(category);
  if (n_hists_so_far == -99) n_hists_so_far = GetSize();
  std::string component_hist_name = m_folio_name + "_" + category.m_name;
  // std::cout << "Creating the " << n_hists_so_far << "th histo in folio for "
  //          << component_hist_name << "\n";
  HIST* component_hist =
      (HIST*)m_source_hist->Clone(component_hist_name.c_str());
  InitializeComponentHist(category, component_hist, n_hists_so_far);
}

// ... from a hist
template <class HIST, class CATEGORY>
void HistFolio<HIST, CATEGORY>::AddComponentHist(const std::string name,
                                                 HIST* component_hist) {
  const int n_hists_so_far = GetSize();
  NamedCategory<CATEGORY> category(n_hists_so_far, name.c_str());
  CheckForDuplicates(category);

  // std::string component_hist_name = m_folio_name + "_" + category.m_name;
  // std::cout << "Adding the " << n_hists_so_far << "th histo in folio for "
  //          << component_hist_name << "\n";
  InitializeComponentHist(category, component_hist, n_hists_so_far);
}

//============================================================================
//
// SET AND GET (PUBLIC)
//
//============================================================================
template <class HIST, class CATEGORY>
int HistFolio<HIST, CATEGORY>::GetSize() const {
  int size = static_cast<int>(m_hist_map.size());
  assert(("ERROR: hist map and array size mismatch",
          size == m_hist_array.GetEntries()));
  return size;
}

template <class HIST, class CATEGORY>
std::string HistFolio<HIST, CATEGORY>::GetName() const {
  return m_folio_name;
}

// Get component histogram either by name or CATEGORY
template <class HIST, class CATEGORY>
template <class T>
HIST* HistFolio<HIST, CATEGORY>::GetComponentHist(const T key, bool be_quiet) {
  auto it = SearchMapForHist(key);
  if (it != m_hist_map.end()) {
    return (*it).second;
  } else {
    if (!be_quiet) {
      std::cerr << "GetComponentHist Error: " << key << " histogram was not "
                << "found in hist folio " << m_folio_name << "!\n"
                << "Returning a nullptr.\n";
    }
    return nullptr;
  }
}

template <class HIST, class CATEGORY>
TObjArray HistFolio<HIST, CATEGORY>::GetHistArray() const {
  return m_hist_array;
}

template <class HIST, class CATEGORY>
std::map<NamedCategory<CATEGORY>, HIST*> HistFolio<HIST, CATEGORY>::GetHistMap()
    const {
  return m_hist_map;
}

//============================================================================
//
// WRITE TO FILE (PUBLIC)
//
//============================================================================
template <class HIST, class CATEGORY>
void HistFolio<HIST, CATEGORY>::WriteToFile(TFile& f) const {
  f.cd();
  for (auto h : GetHistMap()) ((HIST*)h.second)->Write();
}

//==============================================================================
//
// COLORS (PUBLIC)
//
//==============================================================================
template <class HIST, class CATEGORY>
void HistFolio<HIST, CATEGORY>::ApplyColorPalette(
    MnvColors::EColorPalettes mnv_palette) {
  std::vector<int> palette = MnvColors::GetColors(mnv_palette);
  ApplyColorPalette(palette);
}

template <class HIST, class CATEGORY>
void HistFolio<HIST, CATEGORY>::ApplyColorPalette(std::vector<int> palette) {
  for (int i = 0; i < m_hist_array.GetEntries(); ++i) {
    ((HIST*)m_hist_array.At(i))->SetFillColor(palette[i]);
    ((HIST*)m_hist_array.At(i))->SetLineColor(palette[i]);
  }
}

//============================================================================
//
// PRIVATE
//
//============================================================================
// Component Histogram Initialization
// Set axis labels, titles, colors, and add to member containers
// TODO check that component_hist and m_source_hist binning are compatible.
template <class HIST, class CATEGORY>
void HistFolio<HIST, CATEGORY>::InitializeComponentHist(
    NamedCategory<CATEGORY> category, HIST* component_hist,
    const int n_hists_so_far) {
  std::string component_hist_name = m_folio_name + "_" + category.m_name;
  component_hist->GetXaxis()->SetTitle(m_source_hist->GetXaxis()->GetTitle());
  component_hist->GetYaxis()->SetTitle(m_source_hist->GetYaxis()->GetTitle());
  component_hist->SetTitle(category.m_name.c_str());
  component_hist->SetName(component_hist_name.c_str());

  // Default Colors
  std::vector<int> fill_colors =
                   MnvColors::GetColors(MnvColors::kOkabeItoPalette);
  std::vector<int> line_colors =
                    MnvColors::GetColors(MnvColors::kOkabeItoDarkPalette);
  component_hist->SetLineColor(line_colors.at(n_hists_so_far));
  component_hist->SetFillColor(fill_colors.at(n_hists_so_far));

  // Add histograms to containers
  m_hist_map[category] = component_hist;
  m_hist_array.Add(component_hist);
}

// Histogram access helpers
template <class HIST, class CATEGORY>
bool HistFolio<HIST, CATEGORY>::CheckForDuplicates(
    const NamedCategory<CATEGORY> category) {
  bool already_contains = GetComponentHist(category.m_value) != nullptr ||
                          GetComponentHist(category.m_name) != nullptr;
  if (already_contains)
    std::cerr << "WARNING: HistFolio " << GetName()
              << " already contains component hist with name \""
              << category.m_name << "\"!\n  Did you mean to add a duplicate?\n";

  return already_contains;
}

// Search for histo in map by string
template <class HIST, class CATEGORY>
typename std::map<NamedCategory<CATEGORY>, HIST*>::const_iterator
HistFolio<HIST, CATEGORY>::SearchMapForHist(const std::string name) const {
  return find_if(m_hist_map.begin(), m_hist_map.end(), by_name(name));
}

// Search for histo in map by a CATEGORY
template <class HIST, class CATEGORY>
typename std::map<NamedCategory<CATEGORY>, HIST*>::const_iterator
HistFolio<HIST, CATEGORY>::SearchMapForHist(const CATEGORY value) const {
  return find_if(m_hist_map.begin(), m_hist_map.end(), by_value(value));
}

// Helper predicate CTOR - string
template <class HIST, class CATEGORY>
HistFolio<HIST, CATEGORY>::by_name::by_name(std::string name) : m_name(name) {}

// Helper predicate comparison - string
template <class HIST, class CATEGORY>
bool HistFolio<HIST, CATEGORY>::by_name::operator()(
    const std::pair<NamedCategory<CATEGORY>, HIST*>& element) const {
  return element.first.m_name == m_name;
}

// Helper predicate CTOR - value
template <class HIST, class CATEGORY>
HistFolio<HIST, CATEGORY>::by_value::by_value(CATEGORY value)
    : m_value(value){};

// Helper predicate comparison - value
template <class HIST, class CATEGORY>
bool HistFolio<HIST, CATEGORY>::by_value::operator()(
    const std::pair<NamedCategory<CATEGORY>, HIST*>& element) const {
  return element.first.m_value == m_value;
}

//==============================================================================
//
// LOAD HISTFOLIO FROM FILE (NOT A MEMBER FUNCTION)
//
//==============================================================================
// LoadHistFolio (not a member function)
#ifndef __CINT__
template <class HIST>
HistFolio<HIST, int> MAT::LoadHistFolioFromFile(TFile& f,
                                                      std::string folio_name) {
  // 1. First find a source histogram with name == folio_name
  HIST* source_hist = FindSourceHist<HIST>(f, folio_name);

  if (!source_hist) {
    HIST* dummy = new HIST();
    std::string input_hist_type = std::string(typeid(*dummy).name());
    std::cerr << "  ERROR LoadHistFolioFromFile: No objects resembling folio \""
              << folio_name << "\" of type \"" << input_hist_type
              << "\" were found in file " << f.GetName() << ".\n";
    std::cerr << "  Returning empty HistFolio.\n";
    HistFolio<HIST> ret(folio_name, dummy);
    return ret;
  }

  // 2. Create return HistFolio from source_hist
  HistFolio<HIST, int> ret(folio_name, source_hist);

  // 3. Then loop the file and add components to folio
  TObject* obj;
  TKey* key;
  TIter next(f.GetListOfKeys());
  while ((key = (TKey*)next())) {
    obj = f.Get(key->GetName());
    std::string full_name = std::string(key->GetName());
    if (full_name.find(folio_name) != std::string::npos) {  // found a match
      std::string category_name =
          full_name.substr(full_name.find(folio_name) + folio_name.size() + 1);
      HIST* component_hist = (HIST*)obj;
      ret.AddComponentHist(category_name, component_hist);
    }
  }  // end file objects loop

  std::cout << "Successfully loaded HistFolio " << folio_name
            << " from file "
               " with "
            << ret.GetSize() << " component Hists.\n";

  return ret;
}
#endif  // __CINT__

// Helper to LoadHistFolioFromFile (not a member function)
template <class HIST>
HIST* MAT::FindSourceHist(TFile& f, std::string folio_name) {
  TObject* obj;
  TKey* key;
  TIter next(f.GetListOfKeys());
  HIST* source_hist = nullptr;

  // Find a source_hist in the file
  while ((key = (TKey*)next())) {
    obj = f.Get(key->GetName());
    std::string full_name = std::string(key->GetName());
    if (full_name.find(folio_name) != std::string::npos) {  // found a match
      source_hist = dynamic_cast<HIST*>(obj);

      // Make sure that source hist is of correct type HIST
      try {
        typeid(*source_hist).name();  // badcast if type mismatch
      } catch (std::exception& e) {
        HIST* dummy = new HIST();
        std::string input_hist_type = std::string(typeid(*dummy).name());
        delete dummy;
        std::cerr << "  ERROR FindSourceHist: Histogram template parameter "
                  << "type " << input_hist_type << " differs from type found "
                  << "in file!\n";
      }
      break;
    }
  }
  return source_hist;
}

#endif  // HISTFOLIO
