#ifndef HIST2DWRAPPER_CXX
#define HIST2DWRAPPER_CXX

#include "Hist2DWrapper.h"
#include "FluxSystematics.h" //flux_reweighter

using namespace PlotUtils;

// Default Constructor
template <typename T>
Hist2DWrapper<T>::Hist2DWrapper() :
  hist(nullptr),
  univToHistMap(),
  nhistsAssignedSoFar()
{}

// Constructor from map< string, vector<universe> >
template <typename T>
Hist2DWrapper<T>::Hist2DWrapper(const char* hist_name, const char* title, 
                             int nBinsX, double xmin, double xmax, 
                             int nBinsY, double ymin, double ymax, 
                             std::map< std::string, std::vector<T*> >& bands) {
  hist=new PlotUtils::MnvH2D(hist_name, title, nBinsX, xmin, xmax, nBinsY, ymin, ymax );
  hist->SetDirectory(0);
  
  for( typename std::map< std::string, std::vector<T*> >::const_iterator band = bands.begin();
         band != bands.end(); ++band){
    std::vector<T*> univs = band->second;
    const std::string name( univs.front()->ShortName() );
    
    FillErrorBandsWithSysUni( univs, univs.size() , name );  
  }
}


// Constructor from an MnvH2D & map< string, vector<universes> >
template <typename T> 
Hist2DWrapper<T>::Hist2DWrapper(MnvH2D* h2d, 
                                std::map< std::string, std::vector<T*> >& bands,
                                bool clear_error_bands) {
  hist = new MnvH2D(*h2d);
  if(clear_error_bands){
    if( hist->GetNErrorSources() > 0 ) 
      printf("Hist2DWrapper: Clearing error bands from MnvH2D %s with existing "
             "error bands",hist->GetName());
    hist->ClearAllErrorBands();
  }
  hist->SetDirectory(0);
  
  for( typename std::map< std::string, std::vector<T*> >::const_iterator 
           band = bands.begin();
           band != bands.end(); ++band){
    std::vector<T*> univs = band->second;
    const std::string name( univs.front()->ShortName() );
    
    FillErrorBandsWithSysUni( univs, univs.size() , name );  
  }
}


// Constructor from unorganized vector<universe>
// Default 2 universes for any and all bands included.
template <typename T>
Hist2DWrapper<T>::Hist2DWrapper(const char* hist_name, const char* title, 
                            int nBinsX, double xmin, double xmax, 
                            int nBinsY, double ymin, double ymax, 
                            std::vector<T*>& univs) {
  hist=new PlotUtils::MnvH2D(hist_name, title, nBinsX, xmin, xmax, nBinsY, ymin, ymax);
  hist->SetDirectory(0);
  std::vector<T*> tmpunivs;

  for(unsigned int iUniv=0; iUniv<univs.size(); ++iUniv){
    tmpunivs.clear();
    T* univ=univs[iUniv];
    tmpunivs.push_back(univ);
    const std::string name(univ->ShortName());

    FillErrorBandsWithSysUni( tmpunivs, 2, name );  
  }
}


// Constructor from an MnvH2D & vector of universes
template <typename T>
Hist2DWrapper<T>::Hist2DWrapper(MnvH2D* h2d, std::vector<T*>& univs,
                                bool clear_error_bands) {
  hist = new MnvH2D(*h2d);
  if(clear_error_bands){
    if( hist->GetNErrorSources() > 0 ) 
      printf("Hist2DWrapper: Clearing error bands from MnvH2D %s with existing "
             "error bands",hist->GetName());
    hist->ClearAllErrorBands();
  }
  hist->SetDirectory(0);
  std::vector<T*> tmpunivs;

  for(unsigned int iUniv=0; iUniv<univs.size(); ++iUniv){
    tmpunivs.clear();
    T* univ=univs[iUniv];
    tmpunivs.push_back(univ);
    const std::string name(univ->ShortName());

    FillErrorBandsWithSysUni( tmpunivs, 2, name );  
  }
}

// The next 2 constructors are hidden from our python bindings because they use a c++11 feature: delegating constructors
#ifndef __GCCXML__
// Constructor from a vector of universes and variable bin widths
template <typename T>
Hist2DWrapper<T>::Hist2DWrapper(const char* name, const char* title, const std::vector<double> xBins,
                                const std::vector<double> yBins,
                                std::vector<T*>& univs): Hist2DWrapper(new MnvH2D(name, title, xBins.size() - 1, xBins.data(),
                                                                                  yBins.size() - 1, yBins.data()),
                                                                       univs)
{
}


//! Constructor from a map of universes and variable bin widths
template <typename T>
Hist2DWrapper<T>::Hist2DWrapper(const char* name, const char* title, const std::vector<double> xBins,
                                const std::vector<double> yBins,
                                std::map< std::string, std::vector<T*> >& bands): Hist2DWrapper(new MnvH2D(name, title, xBins.size() - 1, xBins.data(),
                                                                                                           yBins.size() - 1, yBins.data()),
                                                                                                bands)
{
}
#endif //__GCCXML__


// Copy an MnvH2D's CV histo to each of its vertical error band CV histos
template <typename T>
void Hist2DWrapper<T>::SyncCVHistos() {
  TH2D* theCVHisto = new TH2D( *hist ); 
  std::vector<std::string> bandnames = hist->GetErrorBandNames();
  //std::cout << "Synching Error Band CV's with MnvH2D's CV" << std::endl;
  // loop this MnvH2D's bands
  for(std::vector<std::string>::const_iterator bandname = bandnames.begin(); 
      bandname != bandnames.end(); ++bandname) {
    //std::cout << bandname << std::endl;
    // band is a reference to the error band's CV histo
    PlotUtils::MnvVertErrorBand2D& band = *(hist->GetVertErrorBand((*bandname).c_str()));
    // set the BAND's CV histo to the MnvH2D's CV histo
    band.TH2D::operator=( *theCVHisto );
  }
}


// Fill
template <typename T>
void Hist2DWrapper<T>::FillUniverse(const T& univ, double valueX, double valueY, double weight){
  univHist(&univ)->Fill(valueX, valueY, weight);
}


template <typename T>
void Hist2DWrapper<T>::FillUniverse(const T* univ, double valueX, double valueY, double weight){
  univHist(univ)->Fill(valueX, valueY, weight);
}


// Access to the TH2 corresponding to the given universe
template <typename T> 
TH2D* Hist2DWrapper<T>::univHist(const T* univ) const { 
  try{
    return univToHistMap.at(univ);
  }
  catch(std::out_of_range){
    std::cerr << "Hist2DWrapper::univHist out_of_range error: " 
              << "Universe " << univ->ShortName() << " with m_nsigma " 
              << univ->GetSigma() << " does not exist in this Hist2DWrapper.\n";
    std::exit(2);
  }
}


// Private, used by CTORs
template <typename T>
void Hist2DWrapper<T>::FillErrorBandsWithSysUni(std::vector<T*> univs, 
                                              int nhists, 
                                              std::string name) {
  // Special-case for central value
  if(name=="cv"){
    if( univs.size() != 1 ){ 
      std::cerr << "Hist2DWrapper constructor ERROR: CV should not have "
                   "more than one universe!" << std::endl; 
      std::exit(2); 
    }
    T* cv = univs.front();
    univToHistMap[cv]=hist;
    return;
  }

  // Add this band to the MnvH2D 
  if(!hist->HasVertErrorBand(name)){
    if( univs.front()->ShortName() == "Flux" && univs.front()->UseNuEConstraint() ){
      PlotUtils::flux_reweighter(univs.front()->GetPlaylist(), univs.front()->GetAnalysisNuPDG(), true, univs.front()->GetNFluxUniverses()).AddFluxErrorBand( hist );
    }
    else hist->AddVertErrorBand( name, nhists );
  }
  
  // Connect each band's universe to the corresponding MnvH2D's TH2
  for(unsigned int iUniv=0; iUniv<univs.size(); ++iUniv){
    T* univ=univs[iUniv];
    if( univ->ShortName() != name){ 
      std::cerr << "Hist2DWrapper constructor ERROR: faulty construction "
                   "of systematic bands map. All of a bands' universes "
                   "should have same name."  << std::endl;
      std::exit(2);
    }
    int& nAssigned=nhistsAssignedSoFar[name];
    univToHistMap[univ]=hist->GetVertErrorBand(name)->GetHist(nAssigned);
    nAssigned++;
  }
}

#endif // HIST2DWRAPPER_CXX

