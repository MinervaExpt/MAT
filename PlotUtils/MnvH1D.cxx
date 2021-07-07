#ifndef MNV_MnvH1D_cxx
#define MNV_MnvH1D_cxx 1

#include "MnvH1D.h"
#include "MnvH2D.h"
#include "HistogramUtils.h"

#include <TMath.h>
#include <cassert>
#include <fstream>
using namespace PlotUtils;

//#define MNV1DBG  Turn this on to debug native covariance matrices


//==================================================================================
// CONSTRUCTORS
//==================================================================================

//--------------------------------------------------------
// Copy constructors with default normalized bin width
//--------------------------------------------------------
MnvH1D::MnvH1D() : 
  TH1D(),
  fNormBinWidth(1.)
{ 
  SilentSumw2(); 
}

MnvH1D::MnvH1D( const TVectorD& v ) :
  TH1D( v ),
  fNormBinWidth(1.)
{ 
  SilentSumw2(); 
}

MnvH1D::MnvH1D( const TH1D& h1d ) :
  TH1D( h1d ) 
{ 
  fNormBinWidth = h1d.GetBinWidth(1);
  SilentSumw2();
}


//--------------------------------------------------------
// Copy constructors with specified normalized bin width
//--------------------------------------------------------
MnvH1D::MnvH1D(Double_t normBinWidth) :
  TH1D(),
  fNormBinWidth(normBinWidth)
{
  SilentSumw2();
}

MnvH1D::MnvH1D( const TVectorD& v, Double_t normBinWidth ) :
  TH1D( v ),
  fNormBinWidth(normBinWidth)
{
  SilentSumw2(); 
}

MnvH1D::MnvH1D( const TH1D& h1d, Double_t normBinWidth ) :
  TH1D( h1d ),
  fNormBinWidth(normBinWidth)
{ 
  SilentSumw2();
}

MnvH1D::MnvH1D( const MnvH1D& h ) :
  TH1D( h )
{
  // Deep copy the variables
  DeepCopy( h );
}

MnvH1D::~MnvH1D()
{
  while (fVertErrorBandMap.size() > 0)
  {
    delete fVertErrorBandMap.begin()->second;
    fVertErrorBandMap.erase(fVertErrorBandMap.begin(), ++fVertErrorBandMap.begin());
  }

  while (fLatErrorBandMap.size() > 0)
  {
    delete fLatErrorBandMap.begin()->second;
    fLatErrorBandMap.erase(fLatErrorBandMap.begin(), ++fLatErrorBandMap.begin());
  }
  
  while (fSysErrorMatrix.size() > 0)
  {
    delete fSysErrorMatrix.begin()->second;
    fSysErrorMatrix.erase(fSysErrorMatrix.begin(), ++fSysErrorMatrix.begin());
  }

  while (fRemovedSysErrorMatrix.size() > 0)
  {
    delete fRemovedSysErrorMatrix.begin()->second;
    fRemovedSysErrorMatrix.erase(fRemovedSysErrorMatrix.begin(), ++fRemovedSysErrorMatrix.begin());
  }

  while (fUncorrErrorMap.size() > 0)
  {
    delete fUncorrErrorMap.begin()->second;
    fUncorrErrorMap.erase(fUncorrErrorMap.begin(), ++fUncorrErrorMap.begin());
  }
}

MnvH1D* MnvH1D::Clone(const char* name) const
{
  MnvH1D* a = new MnvH1D(*this);
  if (TString(name) != TString("")){
    a->SetName(name);
  }
  return a;
}

MnvH1D& MnvH1D::operator=( const MnvH1D& h )
{
  // If this is me, then no copy is necessary
  if( this == &h )
    return *this;

  //call the base class assignment operator
  this->TH1D::operator=(h);

  // Delete and clear all vert and lat error bands
  for( std::map<std::string, MnvVertErrorBand*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
    delete it->second;
  fVertErrorBandMap.clear();

  for( std::map<std::string, MnvLatErrorBand*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
    delete it->second;
  fLatErrorBandMap.clear();

  //delete and clear all uncorr errors
  for( std::map<std::string, TH1D*>::iterator it = fUncorrErrorMap.begin(); it != fUncorrErrorMap.end(); ++it )
    delete it->second;
  fUncorrErrorMap.clear();

  for( std::map<std::string, TMatrixD*>::iterator it = fSysErrorMatrix.begin(); it != fSysErrorMatrix.end(); ++it )
    delete it->second;
  fSysErrorMatrix.clear();

  // Then deeop copy the variables
  DeepCopy(h);

  return *this;
}

void MnvH1D::DeepCopy( const MnvH1D& h )
{
  // Set bin norm width
  fNormBinWidth = h.GetNormBinWidth();

  // Copy the vert and lat error bands
  std::vector<std::string> vertNames = h.GetVertErrorBandNames();
  for( std::vector<std::string>::iterator name = vertNames.begin(); name != vertNames.end(); ++name )
    fVertErrorBandMap[*name] = new MnvVertErrorBand( *h.GetVertErrorBand(*name) );

  std::vector<std::string> latNames = h.GetLatErrorBandNames();
  for( std::vector<std::string>::iterator name = latNames.begin(); name != latNames.end(); ++name )
    fLatErrorBandMap[*name] = new MnvLatErrorBand( *h.GetLatErrorBand(*name) );

  //copy all uncorr errors
  std::vector<std::string> uncorrNames = h.GetUncorrErrorNames();
  for( std::vector<std::string>::iterator name = uncorrNames.begin(); name != uncorrNames.end(); ++name )
    fUncorrErrorMap[*name] = new TH1D( *h.GetUncorrError(*name) );
  
  std::vector<std::string> sysErrorMatrixNames = h.GetCovMatricesNames();
  for( std::vector<std::string>::iterator name = sysErrorMatrixNames.begin(); name != sysErrorMatrixNames.end(); ++name ){
    fSysErrorMatrix[*name] = new TMatrixD( h.GetSysErrorMatrix(*name) );
  }
}

//--------------------------------------------------------
// Constructors with default normalization bin width
//--------------------------------------------------------
MnvH1D::MnvH1D( const char* name, const char* title, Int_t nbinsx, const Float_t* xbins) :
  TH1D( name, title, nbinsx, xbins) 
{ 
  // Default normlization bin width is the first bin width
  fNormBinWidth = xbins[1] - xbins[0];
  SilentSumw2();
}

MnvH1D::MnvH1D( const char* name, const char* title, Int_t nbinsx, const Double_t* xbins) :
  TH1D( name, title, nbinsx, xbins) 
{ 
  // Default normlization bin width is the first bin width
  fNormBinWidth = xbins[1] - xbins[0];
  SilentSumw2();
}

MnvH1D::MnvH1D( const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup):
  TH1D( name, title, nbinsx, xlow, xup )
{ 
  // Default normalization bin width is the constant width of the bins
  fNormBinWidth = (xup - xlow) / float(nbinsx);
  SilentSumw2();
}

//--------------------------------------------------------
// Constructors with specified normalization bin width
//--------------------------------------------------------
MnvH1D::MnvH1D( const char* name, const char* title, Int_t nbinsx, const Float_t* xbins, Double_t normBinWidth) :
  TH1D( name, title, nbinsx, xbins),
  fNormBinWidth(normBinWidth)
{ 
  SilentSumw2();
}

MnvH1D::MnvH1D( const char* name, const char* title, Int_t nbinsx, const Double_t* xbins, Double_t normBinWidth) :
  TH1D( name, title, nbinsx, xbins),
  fNormBinWidth(normBinWidth)
{ 
  SilentSumw2();
}

MnvH1D::MnvH1D( const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Double_t normBinWidth):
  TH1D( name, title, nbinsx, xlow, xup ),
  fNormBinWidth(normBinWidth)
{
  SilentSumw2();
}


//--------------------------------------------------------
// trivial helper functions
//--------------------------------------------------------
bool MnvH1D::HasEnding (std::string const &fullString, std::string const &ending) const
{
  TString a(fullString);
  return a.EndsWith( ending.c_str() );
}


void MnvH1D::SilentSumw2()
{
  if( 0 == GetSumw2N() )
    Sumw2();
}


//! Rename all histograms inside MnvH1D + Error Bands
void MnvH1D::RenameHistosAndErrorBands( const std::string& name )
{

  this->SetName( name.c_str() );

  std::vector<std::string> vert_errBandNames = this->GetVertErrorBandNames();
  std::vector<std::string> lat_errBandNames  = this->GetLatErrorBandNames();
  std::vector<std::string> uncorr_errBandNames  = this->GetUncorrErrorNames();

  for (std::vector<std::string>::iterator itName = vert_errBandNames.begin(); itName != vert_errBandNames.end(); ++itName) 
  {
    MnvVertErrorBand* tmp_band = this->GetVertErrorBand(*itName);
    std::string band_name = std::string(name + "_" + *itName);
    tmp_band->SetName( band_name.c_str() );
    for (unsigned int i=0; i < tmp_band->GetNHists(); ++i)
      tmp_band->GetHist(i)->SetName( Form("%s_universe%d",band_name.c_str(),i) );
  }

  for (std::vector<std::string>::iterator itName = lat_errBandNames.begin(); itName != lat_errBandNames.end(); ++itName) 
  {
    MnvLatErrorBand* tmp_band = this->GetLatErrorBand(*itName);
    std::string band_name = std::string(name + "_" + *itName);
    tmp_band->SetName( band_name.c_str() );
    for (unsigned int i=0; i < tmp_band->GetNHists(); ++i)
      tmp_band->GetHist(i)->SetName( Form("%s_universe%d",band_name.c_str(),i) );
  }

  for (std::vector<std::string>::iterator itName = uncorr_errBandNames.begin(); itName != uncorr_errBandNames.end(); ++itName)
  {
    TH1D* tmp_band = this->GetUncorrError(*itName);
    std::string band_name = std::string(name + "_" + *itName);
    tmp_band->SetName( band_name.c_str() );
  }
  
  // UPDATE can't rename matrices
}

//This removes SumW2 from the error band hists.  This is saves about half the space, but
// if you care about the stat error of your error bands, calling this will change that error
void MnvH1D::UnSumw2Universes()
{
  std::vector<std::string> vert_errBandNames = this->GetVertErrorBandNames();

  for (std::vector<std::string>::iterator itName = vert_errBandNames.begin(); itName != vert_errBandNames.end(); ++itName) 
  {
    MnvVertErrorBand* band = this->GetVertErrorBand(*itName);
    for (unsigned int i=0; i < band->GetNHists(); ++i){
      band->GetHist(i)->Sumw2(false);
    }
  }

  std::vector<std::string> lat_errBandNames  = this->GetLatErrorBandNames();

  for (std::vector<std::string>::iterator itName = lat_errBandNames.begin(); itName != lat_errBandNames.end(); ++itName) 
  {
    MnvLatErrorBand* band = this->GetLatErrorBand(*itName);
    for (unsigned int i=0; i < band->GetNHists(); ++i){
      band->GetHist(i)->Sumw2(false);
    }
  }

}

//! Delete all Error Bands 
void MnvH1D::ClearAllErrorBands()
{

  ClearSysErrorMatrices( );

  // Delete and clear all vert and lat error bands
  for( std::map<std::string, MnvVertErrorBand*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
    delete it->second;
  fVertErrorBandMap.clear();

  for( std::map<std::string, MnvLatErrorBand*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
    delete it->second;
  fLatErrorBandMap.clear();

  for( std::map<std::string, TH1D*>::iterator it = fUncorrErrorMap.begin(); it != fUncorrErrorMap.end(); ++it )
    delete it->second;
  fUncorrErrorMap.clear();
  
  for( std::map<std::string, TMatrixD*>::iterator it = fSysErrorMatrix.begin(); it != fSysErrorMatrix.end(); ++it )
    delete it->second;
  fSysErrorMatrix.clear();

}

bool MnvH1D::AddLatErrorBand( const std::string& name, const int nhists /* = -1 */ )
{
  // Make sure there are no ErrorBands with this name already
  if( HasErrorBand( name ) )
  {
    Warning("MnvH1D::AddLatErrorBand", Form("There is already an error band with name \"%s\".  Doing nothing.", name.c_str()) );

    return false;
  }

  // Error bands we own have this MnvH1D's name as a prefix
  const std::string errName( std::string(GetName()) + "_" + name );

  // non-positive nhists means you want to use the LatErrorBandDefault
  if( nhists > 0 )
    fLatErrorBandMap[name] = new MnvLatErrorBand( errName, (TH1D*)this, nhists );
  else
    fLatErrorBandMap[name] = new MnvLatErrorBand( errName, (TH1D*)this );

  return true;
}

bool MnvH1D::AddLatErrorBand( const std::string& name, const std::vector<TH1D*>& base )
{
  // Make sure there are no ErrorBands with this name already
  if( HasErrorBand( name ) )
  {
    Warning("MnvH1D::AddLatErrorBand", Form("There is already an error band with name \"%s\".  Doing nothing.", name.c_str()) );
    return false;
  }

  // Error bands we own have this MnvH1D's name as a prefix
  const std::string errName( std::string(GetName()) + "_" + name );

  // Set the ErrorBand
  fLatErrorBandMap[name] = new MnvLatErrorBand( errName, (TH1D*)this, base );

  return true;
}

bool MnvH1D::AddLatErrorBandAndFillWithCV( const std::string& name, const int nhists )
{
  // Make sure there are no ErrorBands with this name already
  if( HasErrorBand( name ) )
  {
    Warning("MnvH1D::AddLatErrorBandAndFillWithCV", Form("There is already an error band with name \"%s\".  Doing nothing.", name.c_str()) );
    return false;
  }

  // Make sure number of requested number of histos is not negative 
  if( nhists < 0 )
  {
    Warning("MnvH1D::AddLatErrorBandAndFillWithCV", Form("Passing a negative number of universes to create error band with name \"%s\".  Doing nothing.", name.c_str()) );
    return false;
  }

  // Make a vector of histos with the CV 
  std::vector<TH1D*> histos;
  for( int universe=0; universe<nhists; ++universe )
  {
    TH1D* histo = new TH1D( *this );
    histo->SetName( Form( "tmp_universe_%i", universe ) );
    histos.push_back( histo );
  }

  // Add the error band and fill it with the vector of histos
  bool ok = this->AddLatErrorBand( name, histos );

  // Clean vector of histos
  for( std::vector<TH1D*>::iterator it=histos.begin(); it!=histos.end(); ++it )
  {
    if( *it ) delete *it;
  }

  return ok;
}

bool MnvH1D::AddVertErrorBand( const std::string& name, const int nhists /* = -1 */ )
{
  // Make sure there are no ErrorBands with this name already
  if( HasErrorBand( name ) )
  {
    Warning("MnvH1D::AddVertErrorBand", Form("There is already an error band with name \"%s\".  Doing nothing.", name.c_str()) );
  
    return false;
  }

  // Error bands we own have this MnvH1D's name as a prefix
  const std::string errName( std::string(GetName()) + "_" + name );

  // non-positive nhists means you want to use the VertErrorBand's default
  if( nhists > 0 )
    fVertErrorBandMap[name] = new MnvVertErrorBand( errName, (TH1D*)this, nhists );
  else
    fVertErrorBandMap[name] = new MnvVertErrorBand( errName, (TH1D*)this );

  return true;
}

bool MnvH1D::AddVertErrorBand( const std::string& name, const std::vector<TH1D*>& base )
{
  // Make sure there are no ErrorBands with this name already
  if( HasErrorBand( name ) )
  {
    Warning("MnvH1D::AddVertErrorBand", Form("There is already an error band with name \"%s\".  Doing nothing.", name.c_str()) );
    return false;
  }

  // Error bands we own have this MnvH1D's name as a prefix
  const std::string errName( std::string(GetName()) + "_" + name );

  // Set the ErrorBand
  fVertErrorBandMap[name] = new MnvVertErrorBand( errName, (TH1D*)this, base );

  return true;
}


bool MnvH1D::AddUncorrError( const std::string& name )
{
  // Make sure there are no ErrorBands with this name already
  if( HasErrorBand( name ) )
  {
    Warning("MnvH1D::AddUncorrError", Form("There is already an error band with name \"%s\".  Doing nothing.", name.c_str()) );
    return false;
  }

  // Histogram will have this MnvH1D's name as a prefix
  const std::string errName( std::string(GetName()) + "_" + name );

  //create the uncorrelated error histogram by making a copy of self
  //reset the histo because we intend to fill it ourselves
  TH1D *uncorrHist = dynamic_cast<TH1D*>( this->TH1D::Clone( errName.c_str() ) );
  uncorrHist->Reset();

  // Set the ErrorBand in the map
  fUncorrErrorMap[name] = uncorrHist;

  return true;
}

bool MnvH1D::AddCovMatrix( const std::string& name )
{
  return AddCovMatrixAndFillWithCV(name);  // CV is zero by the way
}


bool MnvH1D::AddUncorrError( const std::string& name, const TH1D *hist, bool errInContent /*=false*/ )
{
  // Make sure there are no ErrorBands with this name already
  if( HasErrorBand( name ) )
  {
    Warning("MnvH1D::AddUncorrError", Form("There is already an error band with name \"%s\".  Doing nothing.", name.c_str()) );
    return false;
  }

  // Histogram will have this MnvH1D's name as a prefix
  const std::string errName( std::string(GetName()) + "_" + name );

  //create the uncorrelated error histogram by making a copy of self
  //keep content intact since we assume no more filling will happen
  TH1D *uncorrHist = dynamic_cast<TH1D*>( this->TH1D::Clone( errName.c_str() ) );

  //replace the errors with those from provided hist
  for( int i = 0; i <= hist->GetNbinsX()+1; ++i )
  {
    const double binErr = errInContent ?
      hist->GetBinContent(i):
      hist->GetBinError(i);
    uncorrHist->SetBinError(i, binErr);
  }

  // Set the ErrorBand in the map
  fUncorrErrorMap[name] = uncorrHist;

  return true;
}

bool MnvH1D::AddCovMatrix( const std::string& name, const TMatrixD * m, bool errInContent /*=false*/ )
{
  // Make sure there are no ErrorBands with this name already
  if( HasErrorBand( name ) )
  {
    Warning("MnvH1D::AddCovMatrix", Form("There is already an cov matrix with name \"%s\".  Doing nothing.", name.c_str()) );
    return false;
  }

  return FillSysErrorMatrix(name,*m);  // this actually does what Add should have done. Keep the fill method for backwards compatibility

  return true;
}

bool MnvH1D::AddUncorrErrorAndFillWithCV( const std::string& name )
{
  // Make sure there are no Errors with this name already
  if( HasErrorBand( name ) )
  {
    Warning("MnvH1D::AddUncorrErrorAndFillWithCV", Form("There is already an uncorr error with name \"%s\".  Doing nothing.", name.c_str()) );
    return false;
  }

  // Histogram will have this MnvH1D's name as a prefix
  const std::string errName( std::string(GetName()) + "_" + name );

  //create the uncorrelated error histogram by making a copy of self
  //keep content intact since we assume no more filling will happen
  TH1D *uncorrHist = dynamic_cast<TH1D*>( this->TH1D::Clone( errName.c_str() ) );

  //replace the errors with 0
  for( int i = 0; i <= uncorrHist->GetNbinsX()+1; ++i )
    uncorrHist->SetBinError(i, 0.);

  // Set the ErrorBand in the map
  fUncorrErrorMap[name] = uncorrHist;

  return true;
}




bool MnvH1D::AddVertErrorBandAndFillWithCV( const std::string& name, const int nhists )
{
  // Make sure there are no ErrorBands with this name already
  if( HasErrorBand( name ) )
  {
    Warning("MnvH1D::AddVertErrorBandAndFillWithCV", Form("There is already an error band with name \"%s\".  Doing nothing.", name.c_str()) );
    return false;
  }

  // Make sure requested number of histos is not negative 
  if( nhists < 0 )
  {
    Warning("MnvH1D::AddVertErrorBandAndFillWithCV", Form("Passing a negative number of universes to create error band with name \"%s\".  Doing nothing.", name.c_str()) );
    return false;
  }

  // Make a vector of histos with the CV 
  std::vector<TH1D*> histos;
  for( int universe=0; universe<nhists; ++universe )
  {
    TH1D* histo = new TH1D( *this );
    histo->SetName( Form( "tmp_universe_%i", universe ) );
    histos.push_back( histo );
  }

  // Add the error band and fill it with the vector of histos
  bool ok = this->AddVertErrorBand( name, histos );

  // Clean vector of histos
  for( std::vector<TH1D*>::iterator it=histos.begin(); it!=histos.end(); ++it )
  {
    if( *it ) delete *it;
  }

  return ok;
}


bool MnvH1D::AddCovMatrixAndFillWithCV( const std::string& name )
{
  // Make sure there are no ErrorBands with this name already
  if( HasErrorBand( name ) )
  {
    Warning("MnvH1D::AddCovMatrixAndFillWithCV", Form("There is already an error band with name \"%s\".  Doing nothing.", name.c_str()) );
    return false;
  }

  int nrows = this->GetNbinsX() + 2;
  TMatrixD tmp(nrows,nrows);

  // Add the error band and fill it with the vector of histos
  bool ok = this->FillSysErrorMatrix(name, tmp );

  return ok;
}




bool MnvH1D::AddMissingErrorBandsAndFillWithCV( const MnvH1D& ref )
{
  // Declare container for error band names
  std::vector<std::string> names = ref.GetVertErrorBandNames();

  // Add vertical error bands found in reference MnvH1D
  for( std::vector<std::string>::iterator name=names.begin(); name!=names.end(); ++name )
  {
    // Skip already added vertical error bands
    if( HasVertErrorBand( *name ) )
      continue;
    unsigned int nunis = ref.GetVertErrorBand( *name )->GetNHists();
    if( ! this->AddVertErrorBandAndFillWithCV( *name, nunis ) )
      return false;
    if ( ref.GetVertErrorBand(*name)->GetUnivWgts() ) GetVertErrorBand(*name )->SetUnivWgts( *(ref.GetVertErrorBand(*name)->GetUnivWgts()) );
  }

  // Add lateral error bands found in reference MnvH1D
  names = ref.GetLatErrorBandNames();
  for( std::vector<std::string>::iterator name=names.begin(); name!=names.end(); ++name )
  {
    // Skip already added lateral error bands
    if( HasLatErrorBand( *name ) )
      continue;
    unsigned int nunis = ref.GetLatErrorBand( *name )->GetNHists();
    if( ! this->AddLatErrorBandAndFillWithCV( *name, nunis ) )
      return false;
    if ( ref.GetLatErrorBand(*name)->GetUnivWgts() ) GetLatErrorBand(*name )->SetUnivWgts( *(ref.GetLatErrorBand(*name)->GetUnivWgts()) );
  }

  // Add uncorrlated error bands found in reference MnvH1D
  names = ref.GetUncorrErrorNames();
  for( std::vector<std::string>::iterator name=names.begin(); name!=names.end(); ++name )
  {
    // Skip already added vertical and lateral error bands
    if( HasUncorrError( *name ) )
      continue;
    if( ! this->AddUncorrErrorAndFillWithCV( *name ) )
      return false;
  }
  
  names = ref.GetCovMatricesNames();
  for( std::vector<std::string>::iterator name=names.begin(); name!=names.end(); ++name )
  {
    if( HasErrorMatrix( *name ) )
      continue;
    if (! this->AddCovMatrixAndFillWithCV(*name))
      return false;
  }
  
  


  return true;

}

bool MnvH1D::AddMissingErrorBandsAndFillWithCV( const MnvH2D& ref )
{
  // Declare container for error band names
  std::vector<std::string> names = ref.GetVertErrorBandNames();

  // Add vertical error bands found in reference MnvH1D
  for( std::vector<std::string>::iterator name=names.begin(); name!=names.end(); ++name )
  {
    // Skip already added vertical error bands
    if( HasVertErrorBand( *name ) )
      continue;
    unsigned int nunis = ref.GetVertErrorBand( *name )->GetNHists();
    if( ! this->AddVertErrorBandAndFillWithCV( *name, nunis ) )
      return false;
    if ( ref.GetVertErrorBand(*name)->GetUnivWgts() ) GetVertErrorBand(*name )->SetUnivWgts( *(ref.GetVertErrorBand(*name)->GetUnivWgts()) );
  }

  // Add lateral error bands found in reference MnvH1D
  names = ref.GetLatErrorBandNames();
  for( std::vector<std::string>::iterator name=names.begin(); name!=names.end(); ++name )
  {
    // Skip already added lateral error bands
    if( HasLatErrorBand( *name ) )
      continue;
    unsigned int nunis = ref.GetLatErrorBand( *name )->GetNHists();
    if( ! this->AddLatErrorBandAndFillWithCV( *name, nunis ) )
      return false;
    if ( ref.GetLatErrorBand(*name)->GetUnivWgts() ) GetLatErrorBand(*name )->SetUnivWgts( *(ref.GetLatErrorBand(*name)->GetUnivWgts()) );
  }
  return true;

}



bool MnvH1D::HasLatErrorBand( const std::string& name ) const
{
  // Check the MnvLatErrorBands
  if( fLatErrorBandMap.find( name ) != fLatErrorBandMap.end() )
    return true;

  return false;
}

bool MnvH1D::HasVertErrorBand( const std::string& name ) const
{
  // Check the MnvVertErrorBands
  if( fVertErrorBandMap.find( name ) != fVertErrorBandMap.end() )
    return true;

  return false;
}

bool MnvH1D::HasUncorrError( const std::string& name ) const
{
  // Check the uncorr errors
  if( fUncorrErrorMap.find( name ) != fUncorrErrorMap.end() )
    return true;

  return false;
}


bool MnvH1D::HasErrorBand( const std::string& name ) const
{
  // Check the MnvLatErrorBands
  if( HasLatErrorBand( name ) )
    return true;

  // Check the MnvVertErrorBands
  if( HasVertErrorBand( name ) )
    return true;

  // Check the uncorrelated errors
  if( HasUncorrError( name ) )
    return true;

  return false;
}

bool MnvH1D::HasErrorMatrix( const std::string& name ) const
{
  // Check the fSysErrorMatrix
  if( fSysErrorMatrix.find( name ) != fSysErrorMatrix.end() )
    return true;

  return false;
}

bool MnvH1D::HasRemovedErrorMatrix( const std::string& name ) const
{
  // Check the fSysErrorMatrix
  if( fRemovedSysErrorMatrix.find( name ) != fRemovedSysErrorMatrix.end() )
    return true;

  return false;
}

MnvLatErrorBand* MnvH1D::GetLatErrorBand( const std::string& name )
{
  std::map<std::string, MnvLatErrorBand*>::iterator i = fLatErrorBandMap.find( name );
  if( i == fLatErrorBandMap.end() )
  {
    Warning( "MnvH1D::GetLatErrorBand", Form("There is no lateral error band with name \"%s\".  Returning NULL.", name.c_str() ));

    return NULL;
  }

  return i->second;
}

MnvVertErrorBand* MnvH1D::GetVertErrorBand( const std::string& name )
{
  std::map<std::string, MnvVertErrorBand*>::iterator i = fVertErrorBandMap.find( name );
  if( i == fVertErrorBandMap.end() )
  {
    Warning( "MnvH1D::GetVertErrorBand", Form("There is no vertical error band with name \"%s\".  Returning NULL.", name.c_str() ));

    return NULL;
  }

  return i->second;
}


const MnvLatErrorBand* MnvH1D::GetLatErrorBand( const std::string& name ) const
{
  std::map<std::string, MnvLatErrorBand*>::const_iterator i = fLatErrorBandMap.find( name );
  if( i == fLatErrorBandMap.end() )
  {
    Warning( "MnvH1D::GetLatErrorBand", Form("There is no lateral error band with name \"%s\".  Returning NULL.", name.c_str() ));
 
    return NULL;
  }

  return i->second;
}

const MnvVertErrorBand* MnvH1D::GetVertErrorBand( const std::string& name ) const
{
  std::map<std::string, MnvVertErrorBand*>::const_iterator i = fVertErrorBandMap.find( name );
  if( i == fVertErrorBandMap.end() )
  {
    Warning( "MnvH1D::GetVertErrorBand", Form("There is no vertical error band with name \"%s\".  Returning NULL.", name.c_str() ));
 
    return NULL;
  }

  return i->second;
}

TH1D* MnvH1D::GetUncorrError( const std::string& name )
{
  std::map<std::string, TH1D*>::const_iterator i = fUncorrErrorMap.find( name );
  if( i == fUncorrErrorMap.end() )
  {
    Warning( "MnvH1D::GetUncorrError", Form("There is no uncorrelated error with name \"%s\".  Returning NULL.", name.c_str() ));
 
    return NULL;
  }

  return i->second;
}

const TH1D* MnvH1D::GetUncorrError( const std::string& name ) const
{
  std::map<std::string, TH1D*>::const_iterator i = fUncorrErrorMap.find( name );
  if( i == fUncorrErrorMap.end() )
  {
    Warning( "MnvH1D::GetUncorrError", Form("There is no uncorrelated error with name \"%s\".  Returning NULL.", name.c_str() ));
 
    return NULL;
  }

  return i->second;
}

TH1D MnvH1D::GetUncorrErrorAsHist( const std::string& name, bool asFrac ) const
{
  //copy self as a blank TH1D
  TH1D rval( *dynamic_cast<const TH1D*>(this) );
  rval.Reset();

  //if the uncorrelated error exists copy the error from its bins to the content of rval
  const TH1D *uncorr = GetUncorrError(name);
  if(uncorr)
  {
    int lowBin = 0;
    int highBin = uncorr->GetNbinsX()+1;
    for( int i = lowBin; i <= highBin; ++i )
    {
      const double err = uncorr->GetBinError(i);

      if(asFrac)
      {
        //not sure what to do if err!=0 but content=0
        //use 0 for now
        const double content = this->GetBinContent(i);
        const double fracErr = ( 0.==content ) ? 0. : fabs(err/content); 
        rval.SetBinContent(i, fracErr );
      }
      else
        rval.SetBinContent( i, err );
    }
  }

  return rval;
}

//==================================
// Transfer error bands
//==================================
bool MnvH1D::TransferErrorBands( MnvH1D *hist, bool removeFromMe )
{
  bool allOK = true;

  std::vector<std::string> names = GetVertErrorBandNames();
  for( std::vector<std::string>::iterator i = names.begin(); i != names.end(); ++i )
    allOK = TransferVertErrorBand( hist, *i, removeFromMe ) && allOK;

  names = GetLatErrorBandNames();
  for( std::vector<std::string>::iterator i = names.begin(); i != names.end(); ++i )
    allOK = TransferLatErrorBand( hist, *i, removeFromMe ) && allOK;
    
  names = GetUncorrErrorNames();  
  for( std::vector<std::string>::iterator i = names.begin(); i != names.end(); ++i )
    allOK = TransferUncorrErrorBand( hist, *i, removeFromMe ) && allOK;

  return allOK;
}

bool MnvH1D::TransferVertErrorBand( MnvH1D *hist, const std::string& errName, bool removeFromMe )
{
  MnvVertErrorBand *errBand(0);
  if( removeFromMe )
    errBand = PopVertErrorBand( errName );
  else
  {
    MnvVertErrorBand *myErrBand = GetVertErrorBand( errName );
    if( 0 != myErrBand )
      errBand = new MnvVertErrorBand( *myErrBand );
  }

  if( 0 == errBand )
  {
    Error( "MnvH1D::TransferVertErrorBand", Form("Could not find error band %s to transfer.", errName.c_str() ) );
    return false;
  }

  //deep scale the error band by ratio CV histos
  TH1D ratio( *( dynamic_cast<TH1D*>(this) ) );
  ratio.Divide( hist, this );
  errBand->MultiplySingle( errBand, &ratio);

  bool pushOK = hist->PushErrorBand( errName, errBand );

  return pushOK;
}

bool MnvH1D::TransferLatErrorBand( MnvH1D *hist, const std::string& errName, bool removeFromMe )
{
  MnvLatErrorBand *errBand(0);
  if( removeFromMe )
    errBand = PopLatErrorBand( errName );
  else
  {
    MnvLatErrorBand *myErrBand = GetLatErrorBand( errName );
    if( 0 != myErrBand )
      errBand = new MnvLatErrorBand( *myErrBand );
  }

  if( 0 == errBand )
  {
    Error( "MnvH1D::TransferLatErrorBand", Form("Could not find error band %s to transfer.", errName.c_str() ) );
    return false;
  }

  //deep scale the error band by ratio CV histos
  TH1D ratio( *( dynamic_cast<TH1D*>(this) ) );
  ratio.Divide( hist, this );
  errBand->MultiplySingle( errBand, &ratio);


  bool pushOK = hist->PushErrorBand( errName, errBand );

  return pushOK;
}

bool MnvH1D::TransferUncorrErrorBand( MnvH1D *hist, const std::string& errName, bool removeFromMe )
{
  TH1D *errBand(0);
  if( removeFromMe )
    errBand = PopUncorrError( errName );
  
  else
  {
    TH1D *myErrBand = GetUncorrError( errName );
    if( 0 != myErrBand )
      errBand = new TH1D( *myErrBand );
  }

  if( 0 == errBand )
  {
    Error( "MnvH1D::TransferUncorrErrorBand", Form("Could not find error band %s to transfer.", errName.c_str() ) );
    return false;
  }

  //deep scale the error band by ratio CV histos
  //TH1D ratio( *( dynamic_cast<TH1D*>(this) ) );
  //ratio.Divide( hist, this );
  //errBand->Multiply( errBand, &ratio);


  bool pushOK = hist->PushUncorrError( errName, errBand );

  return pushOK;
}

MnvLatErrorBand* MnvH1D::PopLatErrorBand( const std::string& name )
{
  std::map<std::string, MnvLatErrorBand*>::iterator i = fLatErrorBandMap.find( name );
  if( i == fLatErrorBandMap.end() )
  {
    Warning( "MnvH1D::PopLatErrorBand", Form("There is no lateral error band with name \"%s\".  Returning NULL.", name.c_str() ));

    return NULL;
  }

  //get a pointer to the error band and remove it from the MnvH1D's vector
  MnvLatErrorBand* rval = i->second;
  fLatErrorBandMap.erase(i);

  return rval;
}



MnvVertErrorBand* MnvH1D::PopVertErrorBand( const std::string& name )
{
  std::map<std::string, MnvVertErrorBand*>::iterator i = fVertErrorBandMap.find( name );
  if( i == fVertErrorBandMap.end() )
  {
    Warning( "MnvH1D::PopVertErrorBand", Form("There is no vertical error band with name \"%s\".  Returning NULL.", name.c_str() ));
  
    return NULL;
  }

  //get a pointer to the error band and remove it from the MnvH1D's vector
  MnvVertErrorBand* rval = i->second;
  fVertErrorBandMap.erase(i);

  return rval;
}

TH1D* MnvH1D::PopUncorrError( const std::string& name )
{ 
  std::map<std::string, TH1D*>::iterator i = fUncorrErrorMap.find( name );
  if( i == fUncorrErrorMap.end() )
  {
    Warning( "MnvH1D::PopUncorrError", Form("There is no uncorrelated error with name \"%s\".  Returning NULL.", name.c_str() ));
    return NULL;
  }

  //get a pointer to the uncorr error and remove it from the MnvH1D's vector
  TH1D* rval = i->second;
  fUncorrErrorMap.erase(i);

  return rval;
} 

// UPDATE2019

TMatrixD* MnvH1D::PopSysErrorMatrix( const std::string& name )
{
  std::map<std::string, TMatrixD*>::iterator i = fSysErrorMatrix.find( name );
  if( i == fSysErrorMatrix.end() )
    {
    Warning( "MnvH1D::popSysErrorMatrix", Form("There is no systematic error matrix with name \"%s\".  Returning NULL.", name.c_str() ));
    
    return NULL;
    }
  
  //get a pointer to the error band and remove it from the MnvH1D's vector
  TMatrixD* rval = i->second;
  fSysErrorMatrix.erase(i);
  
  return rval;
}




bool MnvH1D::PushErrorBand( const std::string& name, MnvVertErrorBand* err )
{
  if( HasVertErrorBand(name) )
  {
    Warning( "MnvH1D::PushErrorBand", Form("I already had vert error band %s.  I'm deleting it and adding the new one.", name.c_str() ) );
    delete PopVertErrorBand(name);
  }
  fVertErrorBandMap[name] = err;
  return true;
}

bool MnvH1D::PushErrorBand( const std::string& name, MnvLatErrorBand* err )
{ 
  if( HasLatErrorBand(name) )
  {
    Warning( "MnvH1D::PushErrorBand", Form("I already had lat error band %s.  I'm deleting it and adding the new one.", name.c_str() ) );
    delete PopLatErrorBand(name);
  }
  fLatErrorBandMap[name] = err;
  return true;
}

bool MnvH1D::PushUncorrError( const std::string& name, TH1D* err )
{ 
  if( HasUncorrError(name) )
  {
    Warning( "MnvH1D::PushUncorrError", Form("I already had uncorrelated error %s.  I'm deleting it and adding the new one.", name.c_str() ) );
    delete PopUncorrError(name);
  }
  //make sure it has sumw2 so div/mult/etc work correctly even though the errors probably arent sqrt(w*w)
  if( 0 == err->GetSumw2N() )
    err->Sumw2();

  fUncorrErrorMap[name] = err;
  return true;
}


std::vector<std::string> MnvH1D::GetErrorBandNames() const
{
  std::vector<std::string> rval;
  for( std::map<std::string, MnvVertErrorBand*>::const_iterator i = fVertErrorBandMap.begin(); i != fVertErrorBandMap.end(); ++i )
    rval.push_back( i->first );
  for( std::map<std::string, MnvLatErrorBand*>::const_iterator i = fLatErrorBandMap.begin(); i != fLatErrorBandMap.end(); ++i )
    rval.push_back( i->first );
  return rval;
}

std::vector<std::string> MnvH1D::GetVertErrorBandNames() const
{
  std::vector<std::string> rval;
  for( std::map<std::string, MnvVertErrorBand*>::const_iterator i = fVertErrorBandMap.begin(); i != fVertErrorBandMap.end(); ++i )
    rval.push_back( i->first );
  return rval;
}

std::vector<std::string> MnvH1D::GetLatErrorBandNames() const
{
  std::vector<std::string> rval;
  for( std::map<std::string, MnvLatErrorBand*>::const_iterator i = fLatErrorBandMap.begin(); i != fLatErrorBandMap.end(); ++i )
    rval.push_back( i->first );
  return rval;
}

std::vector<std::string> MnvH1D::GetUncorrErrorNames() const
{
  std::vector<std::string> rval;
  for( std::map<std::string, TH1D*>::const_iterator i = fUncorrErrorMap.begin(); i != fUncorrErrorMap.end(); ++i )
    rval.push_back( i->first );
  return rval;
}


std::vector<std::string> MnvH1D::GetSysErrorMatricesNames() const
{

  std::vector<std::string> rval;
  //Vertical Errors
  for( std::map<std::string, MnvVertErrorBand*>::const_iterator i = fVertErrorBandMap.begin(); i != fVertErrorBandMap.end(); ++i )
    rval.push_back( i->first );

  //Lateral Errors
  for( std::map<std::string, MnvLatErrorBand*>::const_iterator i = fLatErrorBandMap.begin(); i != fLatErrorBandMap.end(); ++i )
    rval.push_back( i->first );

  //Uncorrelated Errors
  for( std::map<std::string, TH1D*>::const_iterator i = fUncorrErrorMap.begin(); i != fUncorrErrorMap.end(); ++i )
    rval.push_back( i->first );

  //Special Errors
  for( std::map<std::string, TMatrixD*>::const_iterator i = fSysErrorMatrix.begin(); i != fSysErrorMatrix.end(); ++i )
  {
    if ( !HasEnding(i->first, "_asShape") )
      rval.push_back( i->first );
  }

  return rval;

}

std::vector<std::string> MnvH1D::GetCovMatricesNames() const
{
  std::vector<std::string> rval;
  //Special Errors
  for( std::map<std::string, TMatrixD*>::const_iterator i = fSysErrorMatrix.begin(); i != fSysErrorMatrix.end(); ++i )
    {
    if ( !HasEnding(i->first, "_asShape") )
      rval.push_back( i->first );
    }
  
  return rval;
  
}

std::vector<std::string> MnvH1D::GetRemovedSysErrorMatricesNames() const
{
  std::vector<std::string> rval;
  for( std::map<std::string, TMatrixD*>::const_iterator i = fRemovedSysErrorMatrix.begin(); i != fRemovedSysErrorMatrix.end(); ++i )
  {
    if ( !HasEnding(i->first, "_asShape") )
      rval.push_back( i->first );
  }

  return rval;
}

bool MnvH1D::FillLatErrorBand( const std::string& name, const double val, const std::vector<double>& shifts, const double cvweight /* = 1.0 */, const bool fillcv /* = true */, const double *weights /* = 0 */ )
{
  return FillLatErrorBand( name, val, &(shifts[0]), cvweight, fillcv, weights );
}

bool MnvH1D::FillLatErrorBand( const std::string& name, const double val, const double * shifts, const double cvweight /* = 1.0 */, const bool fillcv /* = true */, const double * weights /* = 0 */ )
{
  // Try to fill a lateral error band
  MnvLatErrorBand *lat = GetLatErrorBand( name );
  if( lat )
    return lat->Fill( val, shifts, cvweight, fillcv, weights );

  Warning( "MnvH1D::FillLatErrorBand", Form("Could not find a lateral error band to fill with name = %s", name.c_str() ) );
  return false;
}

bool MnvH1D::FillLatErrorBand( const std::string& name, const double val, const double shiftDown, const double shiftUp, const double cvweight  /*= 1.0*/, const bool fillcv /* = true */ )
{
  // Try to fill a lateral error band
  MnvLatErrorBand *lat = GetLatErrorBand( name );
  if( lat )
    return lat->Fill( val, shiftDown, shiftUp, cvweight, fillcv );

  Warning( "MnvH1D::FillLatErrorBand", Form("Could not find a lateral error band to fill with name = %s", name.c_str() ) );
  return false;
}


bool MnvH1D::FillVertErrorBand( const std::string& name, const double val, const std::vector<double>& weights, const double cvweight /* = 1.0 */, double cvWeightFromMe /*= 1.*/)
{
  return FillVertErrorBand( name, val, &(weights[0]), cvweight, cvWeightFromMe );
}

bool MnvH1D::FillVertErrorBand( const std::string& name, const double val, const double * weights, const double cvweight /* = 1.0 */, double cvWeightFromMe /*= 1.*/)
{
  // Try to fill a vertical error band
  MnvVertErrorBand* vert = GetVertErrorBand( name );
  if( vert )
    return vert->Fill( val, weights, cvweight, cvWeightFromMe );

  Warning( "MnvH1D::FillVertErrorBand", Form("Could not find a vertical error band to fill with name = %s", name.c_str() ) );
  return false;
}

bool MnvH1D::FillVertErrorBand( const std::string& name, const double val, const double weightDown, const double weightUp, const double cvweight  /*= 1.0*/, double cvWeightFromMe /*= 1.*/ )
{
  // Try to fill a vertical error band
  MnvVertErrorBand *vert = GetVertErrorBand( name );
  if( vert )
    return vert->Fill( val, weightDown, weightUp, cvweight, cvWeightFromMe );

  Warning( "MnvH1D::FillVertErrorBand", Form("Could not find a vertical error band to fill with name = %s", name.c_str() ) );
  return false;
}


bool MnvH1D::FillUncorrError( const std::string& name, const double val, const double err, const double cvweight /*= 1.0*/ )
{
  TH1D *hist = GetUncorrError(name);
  if( !hist )
  {
    Warning("MnvH1D::FillUncorrError", Form("Could not find an uncorrelated error to fill with name = %s", name.c_str() ) );
    return false;
  }

  //find the bin
  int bin = hist->FindBin(val);
  //add to bin content, which doesn't change error
  hist->AddBinContent( bin, cvweight );
  //add to bin error so that binErr/binContent = avg( err/cvweight ) always
  hist->SetBinError( bin, err + hist->GetBinError(bin) );

  //may have to copy global stats from CV
  //it can matter when averaging in some cases
  //not sure how long it will take
  double stats[4] = {0.};
  this->GetStats(stats);
  hist->PutStats(stats);

  return true;
}

bool MnvH1D::FillSysErrorMatrix(const std::string& name, const TMatrixD& matrix){
  
  
  if (!HasErrorMatrix(name)){
    std::cout << "MnvH1D::FillSysErrorMatrix: " << GetName() << " creating systematic error matrix "<< name << endl;
    fSysErrorMatrix[name] = new TMatrixD(matrix);
    return true;
  }
  std::cout << "MnvH1D::FillSysErrorMAtrix: " << GetName() << " modifying systematic error matrix " << name << endl;
  delete fSysErrorMatrix[name];
  fSysErrorMatrix[name] = new TMatrixD(matrix);
  return true;
}


void MnvH1D::SetUseSpreadErrorAll( bool use )
{
  for( std::map<std::string, MnvVertErrorBand*>::iterator i = fVertErrorBandMap.begin(); i != fVertErrorBandMap.end(); ++i )
    i->second->SetUseSpreadError( use );

  for( std::map<std::string, MnvLatErrorBand*>::iterator i = fLatErrorBandMap.begin(); i != fLatErrorBandMap.end(); ++i )
    i->second->SetUseSpreadError( use );
}


TH1D MnvH1D::GetCVHistoWithStatError() const
{
  //! @todo I think GetCVHistoWithStatError is unnecessary because it is identical to the TH1D part of self.
  // Get the stat. error band
  TH1D err = GetStatError( false );

  // Create a copy of this histogram and rename it
  TH1D rval( *this );
  std::string tmpName( std::string( GetName() ) + "_CV_WithStatErr" );
  rval.SetName( tmpName.c_str() );

  for( int iBin = 0; iBin <= GetNbinsX(); ++iBin )
    rval.SetBinError( iBin, err.GetBinContent(iBin) );

  return rval;
}

TH1D MnvH1D::GetCVHistoWithError( bool includeStat /* = true */ , bool cov_area_normalize /* = false */) const
{
  // Get the error band
  TH1D err = GetTotalError( includeStat , false, cov_area_normalize);

  // Create a copy of this histogram and rename it
  TH1D rval( *this );
  std::string tmpName( std::string( GetName() ) + "_CV_WithErr" );
  rval.SetName( tmpName.c_str() );

  for( int iBin = 0; iBin <= GetNbinsX(); ++iBin )
    rval.SetBinError( iBin, err.GetBinContent(iBin) );

  return rval;
}

TH1D MnvH1D::GetTotalError(
    bool includeStat /* = true */, 
    bool asFrac /* = false */ ,
    bool cov_area_normalize /*= false */) const
{
  // Make a copy of this histogram as a TH1D and rename it
  TH1D err( *this );
  err.Reset();
  std::string tmpName( std::string(GetName()) + "_TotalError");
  err.SetName( tmpName.c_str() );

  const int highBin = GetNbinsX() + 1;
  const int lowBin = 0;

  // Get the Total Error Matrix
  TMatrixD errMatrix = GetTotalErrorMatrix(includeStat, asFrac, cov_area_normalize);

  for( int iBin = lowBin; iBin <= highBin; ++iBin )
  {
    double derr = errMatrix[iBin][iBin];
    err.SetBinContent( iBin, ( derr > 0 ) ? sqrt(derr): 0. );
  }

  //if you manually set the MnvH1D min/max it exists on err.
  //use default min/max for err
  err.SetMinimum();
  err.SetMaximum();

  return err;
}

//------------------------------------------------------------------------
// Get the Total Covariance Matrix
//------------------------------------------------------------------------
TMatrixD MnvH1D::GetTotalErrorMatrix(
    bool includeStat /*= true*/, 
    bool asFrac /*= false*/, 
    bool cov_area_normalize /*= false*/ ) const
{
  const int highBin = GetNbinsX() + 1;
  const int lowBin = 0;
  TMatrixD covmx(highBin+1,highBin+1);

  std::vector<std::string> names = GetSysErrorMatricesNames();
  for (std::vector<std::string>::const_iterator itName = names.begin() ; itName != names.end() ; ++itName)
  {

    covmx += GetSysErrorMatrix(*itName, false, cov_area_normalize);
  }

  if (includeStat)
    covmx += GetStatErrorMatrix();

  if (asFrac)
  {
    for( int i = lowBin; i <= highBin; ++i )
    {
      for( int k = i; k <= highBin; ++k )
      { 
        //Gettting the the CV value for bin i
        const double cv_i = GetBinContent(i);
        const double cv_k = GetBinContent(k);

        covmx[i][k]= ( (cv_i != 0.) && (cv_k !=0.) )? covmx[i][k]/(cv_i * cv_k) : 0.;
        covmx[k][i]=covmx[i][k];
      }
    }
  }

  return covmx;
}

TMatrixD MnvH1D::GetTotalCorrelationMatrix( bool cov_area_normalize /*= false*/, bool includeStat /*=false*/ ) const
{
  TMatrixD covmx = GetTotalErrorMatrix(includeStat, false, cov_area_normalize);
  const int size = covmx.GetNrows();
  TMatrixD corrmx(size,size);

  for( int i = 0; i < size; ++i )
  {
    for( int k = 0; k < size; ++k )
    {
      corrmx[i][k] = ( covmx[i][i] == 0. || covmx[k][k] == 0. ) ? 0. : covmx[i][k]/sqrt(covmx[i][i]*covmx[k][k]);
    }
  }
  return corrmx;
}

TH2D MnvH1D::GetTotalCorrelationMatrixTH2D( bool cov_area_normalize /*= false*/, bool includeStat /*=false*/ ) const
{
  TMatrixD corrmx = GetTotalCorrelationMatrix(cov_area_normalize,includeStat);
  TH2D convertedCorrMx = TH2D(corrmx);
  convertedCorrMx.GetZaxis()->SetRangeUser(-1.0,1.0);
  return convertedCorrMx;
}

//------------------------------------------------------------------------
// Calculate systematic covariance matrix with another MnvH1D
//------------------------------------------------------------------------
TMatrixD MnvH1D::GetCombinedSysErrorMatrix( std::string errband_name,
                                        MnvH1D *h2, 
					bool   asFrac /*= false*/, 
					bool   cov_area_normalize /*= false*/ ) const
{

  // Combine the two histograms into one histogram with arbitrary bin sizes
  int nbins1 = GetNbinsX();
  int nbins2 = h2->GetNbinsX();
  TH1D *h_temp = new TH1D("h_temp2","h_temp2",nbins1+nbins2,0,nbins1+nbins2); 

  for(int i = 0; i<nbins1; i++) {
    h_temp->SetBinContent(i+1,GetBinContent(i+1));
    h_temp->SetBinError(i+1,GetBinError(i+1));
  }
  for(int i = nbins1; i<nbins1+nbins2; i++) {
    h_temp->SetBinContent(i+1,h2->GetBinContent(i+1-nbins1));
    h_temp->SetBinError(i+1,h2->GetBinError(i+1-nbins1));
  }

  // Create an MnvH1D out of the combined histo
  MnvH1D *mnv_temp = new MnvH1D(*h_temp);

  // Extract the universes from the error bad and for an error band
  // for the combined MnvH1D

  // vert error bands
  if(HasVertErrorBand(errband_name) || h2->HasVertErrorBand(errband_name)) {

    int nhists1 =  0;
    if(HasVertErrorBand(errband_name))
      nhists1 = GetVertErrorBand(errband_name)->GetHists().size();
    int nhists2 = 0;
    if(h2->HasVertErrorBand(errband_name))
       nhists2 = h2->GetVertErrorBand(errband_name)->GetHists().size();
    
    // one of the histograms has got to have an error band for this to work
    if(nhists1!=0 || nhists2!=0) {
    
      // figure out the number of histograms that should be included in the
      // combined error band

      // in most cases, both will have the same number of histograms
      int nhists = nhists1;
    
      // but we have to consider special cases:
      if(nhists1==0) {
	std::cout<<"WARNING: you've asked for the combined covariance matrix of the MnvH1D's, but the first histogram doesn't contain the "<<errband_name<<" error band.  Assuming errors are zero for that histogram and errorband."<<std::endl;
	  nhists = nhists2;
      }
      if(nhists2==0) {
	std::cout<<"WARNING: you've asked for the combined covariance matrix of the MnvH1D's, but the second histogram doesn't contain the "<<errband_name<<" error band.  Assuming errors are zero for that histogram and errorband."<<std::endl;
	  nhists = nhists1;
      }
      if(nhists1!=0 && nhists2!=0 && nhists1<nhists2) {
	std::cout<<"WARNING: you've asked for the combined covariance matrix of the MnvH1D's, but the number of histograms in the "<<errband_name<<" error band don't agree... Only considering the first "<<nhists1<<" histograms."<<std::endl;
	nhists = nhists1;
      }
      if(nhists1!=0 && nhists2!=0 && nhists2<nhists1) {
	std::cout<<"WARNING: you've asked for the combined covariance matrix of the MnvH1D's, but the number of histograms in the "<<errband_name<<" error band don't agree... Only considering the first "<<nhists1<<" histograms."<<std::endl;
	nhists = nhists2;
      }

      std::vector<TH1D*> newhists;

      for(int j = 0; j<nhists; j++) {
	
	std::stringstream ss; ss<<j;
	string temp_name = "combined_"+errband_name+"_"+ss.str();

	TH1D *temp = (TH1D*)h_temp->Clone(temp_name.c_str());
	temp->Clear();

	// Fill in the first half of the error band
	for(int k = 0; k<nbins1; k++) {
	  if(nhists1 !=0) {
	    
	    double area_scale = Integral()/GetVertErrorBand(errband_name)->GetHists()[j]->Integral();
	    if(cov_area_normalize)
	      temp->SetBinContent(k+1,area_scale * GetVertErrorBand(errband_name)->GetHists()[j]->GetBinContent(k+1));
	    else
	      temp->SetBinContent(k+1,GetVertErrorBand(errband_name)->GetHists()[j]->GetBinContent(k+1));
	  }
	  else {
	    // if the error band doesn't exist, just fill in the CV's.
	    temp->SetBinContent(k+1,GetBinContent(k+1));
	  }
	}


	// Fill in second half of error band
	for(int k = nbins1; k<nbins1+nbins2; k++) {
	  if(nhists2 !=0) {
	    
	    double area_scale = h2->Integral()/h2->GetVertErrorBand(errband_name)->GetHists()[j]->Integral();
	    if(cov_area_normalize)
	      temp->SetBinContent(k+1,area_scale * h2->GetVertErrorBand(errband_name)->GetHists()[j]->GetBinContent(k+1-nbins1));
	    else
	      temp->SetBinContent(k+1,h2->GetVertErrorBand(errband_name)->GetHists()[j]->GetBinContent(k+1-nbins1));
	  }
	  else {
	    // if the error band doesn't exist, just fill in the CV's.
	    temp->SetBinContent(k+1,h2->GetBinContent(k+1-nbins1));
	  }
	}
	newhists.push_back(temp);
      }
      mnv_temp->AddVertErrorBand(errband_name,newhists);
      mnv_temp->GetVertErrorBand(errband_name)->SetUseSpreadError(GetVertErrorBand(errband_name)->GetUseSpreadError());
      // Set universe weights; does not work if universes are different for two histograms
      if(nhists1==nhists2) {
	for(int j = 0; j<nhists1; j++) {
	  mnv_temp->GetVertErrorBand(errband_name)->SetUnivWgt(j,GetVertErrorBand(errband_name)->GetUnivWgt(j));
	}
      }
    }
  } 

  // lat error bands
  if(HasLatErrorBand(errband_name) || h2->HasLatErrorBand(errband_name)) {

    int nhists1 =  0;
    if(HasLatErrorBand(errband_name))
      nhists1 = GetLatErrorBand(errband_name)->GetHists().size();
    int nhists2 = 0;
    if(h2->HasLatErrorBand(errband_name))
       nhists2 = h2->GetLatErrorBand(errband_name)->GetHists().size();
    
    // one of the histograms has got to have an error band for this to work
    if(nhists1!=0 || nhists2!=0) {
    
      // figure out the number of histograms that should be included in the
      // combined error band

      // in most cases, both will have the same number of histograms
      int nhists = nhists1;
    
      // but we have to consider special cases:
      if(nhists1==0) {
	std::cout<<"WARNING: you've asked for the combined covariance matrix of the MnvH1D's, but the first histogram doesn't contain the "<<errband_name<<" error band.  Assuming errors are zero for that histogram and errorband."<<std::endl;
	  nhists = nhists2;
      }
      if(nhists2==0) {
	std::cout<<"WARNING: you've asked for the combined covariance matrix of the MnvH1D's, but the second histogram doesn't contain the "<<errband_name<<" error band.  Assuming errors are zero for that histogram and errorband."<<std::endl;
	  nhists = nhists1;
      }
      if(nhists1!=0 && nhists2!=0 && nhists1<nhists2) {
	std::cout<<"WARNING: you've asked for the combined covariance matrix of the MnvH1D's, but the number of histograms in the "<<errband_name<<" error band don't agree... Only considering the first "<<nhists1<<" histograms."<<std::endl;
	nhists = nhists1;
      }
      if(nhists1!=0 && nhists2!=0 && nhists2<nhists1) {
	std::cout<<"WARNING: you've asked for the combined covariance matrix of the MnvH1D's, but the number of histograms in the "<<errband_name<<" error band don't agree... Only considering the first "<<nhists1<<" histograms."<<std::endl;
	nhists = nhists2;
      }

      std::vector<TH1D*> newhists;

      for(int j = 0; j<nhists; j++) {
	
	std::stringstream ss; ss<<j;
	string temp_name = "combined_"+errband_name+"_"+ss.str();

	TH1D *temp = (TH1D*)h_temp->Clone(temp_name.c_str());
	temp->Clear();
	//cout<<"temp_name "<<temp_name<<std::endl;

	// Fill in the first half of the error band
	for(int k = 0; k<nbins1; k++) {
	  if(nhists1 !=0) {
	    
	    double area_scale = Integral()/GetLatErrorBand(errband_name)->GetHists()[j]->Integral();
	    if(cov_area_normalize)
	      temp->SetBinContent(k+1,area_scale * GetLatErrorBand(errband_name)->GetHists()[j]->GetBinContent(k+1));
	    else
	      temp->SetBinContent(k+1,GetLatErrorBand(errband_name)->GetHists()[j]->GetBinContent(k+1));
	  }
	  else {
	    // if the error band doesn't exist, just fill in the CV's.
	    temp->SetBinContent(k+1,GetBinContent(k+1));
	  }
	}


	// Fill in second half of error band
	for(int k = nbins1; k<nbins1+nbins2; k++) {
	  if(nhists2 !=0) {
	    
	    double area_scale = Integral()/h2->GetLatErrorBand(errband_name)->GetHists()[j]->Integral();
	    if(cov_area_normalize)
	      temp->SetBinContent(k+1,area_scale * h2->GetLatErrorBand(errband_name)->GetHists()[j]->GetBinContent(k+1-nbins1));
	    else
	      temp->SetBinContent(k+1,h2->GetLatErrorBand(errband_name)->GetHists()[j]->GetBinContent(k+1-nbins1));
	  }
	  else {
	    // if the error band doesn't exist, just fill in the CV's.
	    temp->SetBinContent(k+1,h2->GetBinContent(k+1-nbins1));
	  }
	}
	newhists.push_back(temp);
      }
      mnv_temp->AddLatErrorBand(errband_name,newhists);
      mnv_temp->GetLatErrorBand(errband_name)->SetUseSpreadError(GetLatErrorBand(errband_name)->GetUseSpreadError());
      // Set universe weights; does not work if universes are different for two histograms
      if(nhists1==nhists2) {
	for(int j = 0; j<nhists1; j++) {
          mnv_temp->GetLatErrorBand(errband_name)->SetUnivWgt(j,GetLatErrorBand(errband_name)->GetUnivWgt(j));
	}
      }
    }
  }   
  
  if(HasUncorrError(errband_name) || h2->HasUncorrError(errband_name)) {
  
    string temp_name = "combined_"+errband_name;

    TH1D *temp = (TH1D*)h_temp->Clone(temp_name.c_str());
    temp->Clear();

    for(int k = 0; k<nbins1; k++) {
      if(HasUncorrError(errband_name)) 
	temp->SetBinContent(k+1,GetUncorrError(errband_name)->GetBinContent(k+1));
      else
	temp->SetBinContent(k+1,0);
    }
       
    for(int k = nbins1; k<nbins1+nbins2; k++) {
      if(h2->HasUncorrError(errband_name)) 
	temp->SetBinContent(k+1,h2->GetUncorrError(errband_name)->GetBinContent(k+1-nbins1));
      else
	temp->SetBinContent(k+1,0);
    }
  
    mnv_temp->AddUncorrError(errband_name,temp);
  }

  // it doesn't make sense to area normalize the final histo; this was done above
  
  TMatrixD covmx = mnv_temp->GetSysErrorMatrix(errband_name, asFrac, false);
  delete h_temp;
  delete mnv_temp;

  return covmx;

}

//------------------------------------------------------------------------
// Get the Total Covariance Matrix
//------------------------------------------------------------------------
TMatrixD MnvH1D::GetCombinedTotalErrorMatrix(
					     MnvH1D *h2,					     
					     bool includeStat /*= true*/, 
					     bool asFrac /*= false*/, 
					     bool cov_area_normalize /*= false*/ ) const
{

  // Combine the two histograms into one histogram with arbitrary bin sizes
  int nbins1 = GetNbinsX();
  int nbins2 = h2->GetNbinsX();
  TH1D *h_temp = new TH1D("h_temp","h_temp",nbins1+nbins2,0,nbins1+nbins2); 

  for(int i = 0; i<nbins1; i++) {
    h_temp->SetBinContent(i+1,GetBinContent(i+1));
    if(includeStat)
      h_temp->SetBinError(i+1,GetBinError(i+1));
  }
  for(int i = nbins1; i<nbins1+nbins2; i++) {
    h_temp->SetBinContent(i+1,h2->GetBinContent(i+1-nbins1));
    if(includeStat)
      h_temp->SetBinError(i+1,h2->GetBinError(i+1-nbins1));
  }

  // Create an MnvH1D out of the combined histo
  MnvH1D *mnv_temp = new MnvH1D(*h_temp);

  TMatrixD covmx = mnv_temp->GetStatErrorMatrix();

  
  std::vector<std::string> names = GetSysErrorMatricesNames();
  for (std::vector<std::string>::const_iterator itName = names.begin() ; itName != names.end() ; ++itName)
  {
    covmx += GetCombinedSysErrorMatrix(*itName, h2, false,cov_area_normalize);
  }
  

  if (asFrac)
  {
    for( int i = 0; i <= covmx.GetNrows(); ++i )
    {
      for( int k = i; k <= covmx.GetNrows(); ++k )
      { 
        //Gettting the the CV value for bin i
        const double cv_i = h_temp->GetBinContent(i);
        const double cv_k = h_temp->GetBinContent(k);

        covmx[i][k]= ( (cv_i != 0.) && (cv_k !=0.) )? covmx[i][k]/(cv_i * cv_k) : 0.;
        covmx[k][i]=covmx[i][k];
      }
    }
  }
  return covmx;
  }


//------------------------------------------------------------------------
// Get the Total Covariance Matrix
//------------------------------------------------------------------------
TMatrixD MnvH1D::GetCombinedStatErrorMatrix( MnvH1D *h2) const			{

  // Combine the two histograms into one histogram with arbitrary bin sizes
  int nbins1 = GetNbinsX();
  int nbins2 = h2->GetNbinsX();
  TH1D *h_temp = new TH1D("h_temp","h_temp",nbins1+nbins2,0,nbins1+nbins2); 

  for(int i = 0; i<nbins1; i++) {
    h_temp->SetBinContent(i+1,GetBinContent(i+1));
    h_temp->SetBinError(i+1,GetBinError(i+1));
  }
  for(int i = nbins1; i<nbins1+nbins2; i++) {
    h_temp->SetBinContent(i+1,h2->GetBinContent(i+1-nbins1));
    h_temp->SetBinError(i+1,h2->GetBinError(i+1-nbins1));
  }

  // Create an MnvH1D out of the combined histo
  MnvH1D *mnv_temp = new MnvH1D(*h_temp);

  TMatrixD covmx = mnv_temp->GetStatErrorMatrix();

  return covmx;
}

//------------------------------------------------------------------------
// Calculate covariance matrix with another MnvH1D
// Assumes stat errors are uncorrelated
//------------------------------------------------------------------------
TMatrixD MnvH1D::GetCombinedTotalErrorMatrix_old(
					MnvH1D *h2, 
					bool   includeStat /*= true*/, 
					bool   asFrac /*= false*/, 
					bool   cov_area_normalize /*= false*/ ) const
{

  // Combine the two histograms into one histogram with arbitrary bin sizes
  int nbins1 = GetNbinsX();
  int nbins2 = h2->GetNbinsX();
  TH1D *h_temp = new TH1D("h_temp","h_temp",nbins1+nbins2,0,nbins1+nbins2); 

  for(int i = 0; i<nbins1; i++) {
    h_temp->SetBinContent(i+1,GetBinContent(i+1));
    h_temp->SetBinError(i+1,GetBinError(i+1));
  }
  for(int i = nbins1; i<nbins1+nbins2; i++) {
    h_temp->SetBinContent(i+1,h2->GetBinContent(i+1-nbins1));
    h_temp->SetBinError(i+1,h2->GetBinError(i+1-nbins1));
  }

  // Create an MnvH1D out of the combined histo
  MnvH1D *mnv_temp = new MnvH1D(*h_temp);

  // Extract the universes from each histo and create error bands
  // for the combined MnvH1D
  
  // vert error bands

  for(unsigned int i = 0; i< GetVertErrorBandNames().size(); i++) {
    
    std::string errband_name = GetVertErrorBandNames()[i];
    int nhists =  GetVertErrorBand(errband_name)->GetHists().size();

    std::vector<TH1D*> newhists;

    for(int j = 0; j<nhists; j++) {
 
      std::stringstream ss; ss<<j;
      string temp_name = "combined_"+errband_name+"_"+ss.str();

      TH1D *temp = (TH1D*)h_temp->Clone(temp_name.c_str());

      temp->Clear();
    
      for(int k = 0; k<nbins1; k++) {
	double area_scale = GetVertErrorBand(errband_name)->GetHists()[j]->Integral()/Integral();
	if(cov_area_normalize)
	  temp->SetBinContent(k+1,area_scale * GetVertErrorBand(errband_name)->GetHists()[j]->GetBinContent(k+1));
	else
	  temp->SetBinContent(k+1,GetVertErrorBand(errband_name)->GetHists()[j]->GetBinContent(k+1));
      }
      for(int k = nbins1; k<nbins1+nbins2; k++) {
	if(h2->HasVertErrorBand(errband_name)) {
	  double area_scale = h2->GetVertErrorBand(errband_name)->GetHists()[j]->Integral()/h2->Integral();
	  if(cov_area_normalize)
	    temp->SetBinContent(k+1,area_scale * h2->GetVertErrorBand(errband_name)->GetHists()[j]->GetBinContent(k+1-nbins1));
	  else
	    temp->SetBinContent(k+1,h2->GetVertErrorBand(errband_name)->GetHists()[j]->GetBinContent(k+1-nbins1));
	}
	else {
	  temp->SetBinContent(k+1,h2->GetBinContent(k+1-nbins1));
	}
      }
      newhists.push_back(temp);
      
    }
    mnv_temp->AddVertErrorBand(errband_name,newhists);
    mnv_temp->GetVertErrorBand(errband_name)->SetUseSpreadError(GetVertErrorBand(errband_name)->GetUseSpreadError());
    // Set universe weights (from nu+e flux constraint, probably)
    for(int j = 0; j<nhists; j++) {
      mnv_temp->GetVertErrorBand(errband_name)->SetUnivWgt(j,GetVertErrorBand(errband_name)->GetUnivWgt(j));
    }
  }

   
  // look for error bands that are in second hist but not the first
  for(unsigned int i = 0; i<h2->GetVertErrorBandNames().size(); i++) {
    std::string errband_name = h2->GetVertErrorBandNames()[i];
    if(!mnv_temp->HasVertErrorBand(errband_name)) {
      int nhists =  GetVertErrorBand(errband_name)->GetHists().size();
      std::vector<TH1D*> newhists;	 
      for(int j=0; j<nhists; j++) {
	std::stringstream ss; ss<<j;
	string temp_name = "combined_"+errband_name+"_"+ss.str();
	
	TH1D *temp = (TH1D*)h_temp->Clone(temp_name.c_str());
	temp->Clear();
	
	for(int k = 0; k<nbins1; k++) {
	  temp->SetBinContent(k+1,GetBinContent(k+1));
	}
	for(int k = nbins1; k<nbins1+nbins2; k++) {
	  double area_scale = h2->GetVertErrorBand(errband_name)->GetHists()[j]->Integral()/h2->Integral();
	  if(cov_area_normalize)
	    temp->SetBinContent(k+1,area_scale * h2->GetVertErrorBand(errband_name)->GetHists()[j]->GetBinContent(k+1-nbins1));
	  else
	    temp->SetBinContent(k+1,h2->GetVertErrorBand(errband_name)->GetHists()[j]->GetBinContent(k+1-nbins1));
	}
	newhists.push_back(temp);
      }
      mnv_temp->AddVertErrorBand(errband_name,newhists);
      mnv_temp->GetVertErrorBand(errband_name)->SetUseSpreadError(GetVertErrorBand(errband_name)->GetUseSpreadError());
      // Set universe weights (from nu+e flux constraint, probably)
      for(int j = 0; j<nhists; j++) {
	mnv_temp->GetVertErrorBand(errband_name)->SetUnivWgt(j,h2->GetVertErrorBand(errband_name)->GetUnivWgt(j));
      }
    }
  }

  // lat error bands
  for(unsigned int i = 0; i< GetLatErrorBandNames().size(); i++) {

    std::string errband_name = GetLatErrorBandNames()[i];
    int nhists =  GetLatErrorBand(errband_name)->GetHists().size();

    std::vector<TH1D*> newhists;


    for(int j = 0; j<nhists; j++) {
 
      std::stringstream ss; ss<<j;
      string temp_name = "combined_"+errband_name+"_"+ss.str();

      TH1D *temp = (TH1D*)h_temp->Clone(temp_name.c_str());
      temp->Clear();

    
      for(int k = 0; k<nbins1; k++) {
	double area_scale = GetLatErrorBand(errband_name)->GetHists()[j]->Integral()/Integral();
	if(cov_area_normalize)
	  temp->SetBinContent(k+1,area_scale * GetLatErrorBand(errband_name)->GetHists()[j]->GetBinContent(k+1));
	else
	  temp->SetBinContent(k+1,GetLatErrorBand(errband_name)->GetHists()[j]->GetBinContent(k+1));
      }
      for(int k = nbins1; k<nbins1+nbins2; k++) {
	if(h2->HasLatErrorBand(errband_name)) {
	  double area_scale = h2->GetLatErrorBand(errband_name)->GetHists()[j]->Integral()/h2->Integral();
	  if(cov_area_normalize)
	    temp->SetBinContent(k+1,area_scale * h2->GetLatErrorBand(errband_name)->GetHists()[j]->GetBinContent(k+1-nbins1));
	  else
	    temp->SetBinContent(k+1,h2->GetLatErrorBand(errband_name)->GetHists()[j]->GetBinContent(k+1-nbins1));
	}
	else {
	  temp->SetBinContent(k+1,h2->GetBinContent(k+1-nbins1));
	}
      }
      newhists.push_back(temp);
    }
    mnv_temp->AddLatErrorBand(errband_name,newhists);
    mnv_temp->GetLatErrorBand(errband_name)->SetUseSpreadError(GetLatErrorBand(errband_name)->GetUseSpreadError());
    // Set universe weights (from nu+e flux constraint, probably)
    for(int j = 0; j<nhists; j++) {
      mnv_temp->GetLatErrorBand(errband_name)->SetUnivWgt(j,GetLatErrorBand(errband_name)->GetUnivWgt(j));
    }
  }
   
  // look for error bands that are in second hist but not the first
  for(unsigned int i = 0; i<h2->GetLatErrorBandNames().size(); i++) {
    std::string errband_name = h2->GetLatErrorBandNames()[i];
    if(!mnv_temp->HasLatErrorBand(errband_name)) {
      int nhists =  GetLatErrorBand(errband_name)->GetHists().size();
      std::vector<TH1D*> newhists;	 
      for(int j=0; j<nhists; j++) {
	std::stringstream ss; ss<<j;
	string temp_name = "combined_"+errband_name+"_"+ss.str();

	TH1D *temp = (TH1D*)h_temp->Clone(temp_name.c_str());
	temp->Clear();

	for(int k = 0; k<nbins1; k++) {
	  temp->SetBinContent(k+1,GetBinContent(k+1));
	}
	for(int k = nbins1; k<nbins1+nbins2; k++) {
	  double area_scale = h2->GetLatErrorBand(errband_name)->GetHists()[j]->Integral()/h2->Integral();
	  if(cov_area_normalize)
	    temp->SetBinContent(k+1,area_scale * h2->GetLatErrorBand(errband_name)->GetHists()[j]->GetBinContent(k+1-nbins1));
	  else
	    temp->SetBinContent(k+1,h2->GetLatErrorBand(errband_name)->GetHists()[j]->GetBinContent(k+1-nbins1));
	}
	newhists.push_back(temp);
      }
      mnv_temp->AddLatErrorBand(errband_name,newhists);
      mnv_temp->GetLatErrorBand(errband_name)->SetUseSpreadError(GetLatErrorBand(errband_name)->GetUseSpreadError());
      // Set universe weights (from nu+e flux constraint, probably)
      for(int j = 0; j<nhists; j++) {
	mnv_temp->GetLatErrorBand(errband_name)->SetUnivWgt(j,h2->GetLatErrorBand(errband_name)->GetUnivWgt(j));
      }
    }   
  }

  // uncorr error bands
  for(unsigned int i = 0; i< GetUncorrErrorNames().size(); i++) {

    std::string errband_name = GetUncorrErrorNames()[i];

    string temp_name = "combined_"+errband_name;

    TH1D *temp = (TH1D*)h_temp->Clone(temp_name.c_str());
    temp->Clear();
    
    for(int k = 0; k<nbins1; k++) {
      temp->SetBinContent(k+1,GetUncorrError(errband_name)->GetBinContent(k+1));
    }
    for(int k = nbins1; k<nbins1+nbins2; k++) {
      if(h2->HasUncorrError(errband_name)) {
	temp->SetBinContent(k+1,h2->GetUncorrError(errband_name)->GetBinContent(k+1-nbins1));
      }
      else {
	temp->SetBinContent(k+1,0);
      }
    }
    mnv_temp->AddUncorrError(errband_name,temp);
  }
   
  // look for error bands that are in second hist but not the first
  for(unsigned int i = 0; i<h2->GetUncorrErrorNames().size(); i++) {
    std::string errband_name = h2->GetUncorrErrorNames()[i];
    if(!mnv_temp->HasUncorrError(errband_name)) {
      string temp_name = "combined_"+errband_name;
      
      TH1D *temp = (TH1D*)h_temp->Clone(temp_name.c_str());
      temp->Clear();

      for(int k = 0; k<nbins1; k++) {
	temp->SetBinContent(k+1,0);
      }
      for(int k = nbins1; k<nbins1+nbins2; k++) {
	temp->SetBinContent(k+1,h2->GetUncorrError(errband_name)->GetBinContent(k+1-nbins1));
      }
      mnv_temp->AddUncorrError(errband_name,temp);
    }
  }

  delete h_temp;
  
  // it doesn't make sense to area normalize the final histo; this was done above
  return mnv_temp->GetTotalErrorMatrix(includeStat,asFrac,false);
}




Int_t MnvH1D::WriteTotalErrorMatrix(const std::string& name, bool includeStat /*= true*/, bool asFrac /*= false*/, bool cov_area_normalize /*= false*/)
{
  TMatrixD covmx = GetTotalErrorMatrix(includeStat, asFrac, cov_area_normalize);
  Int_t state = covmx.Write(name.c_str());
  return state;
}

Int_t MnvH1D::WriteTotalCorrelationMatrix(const std::string& name, bool cov_area_normalize /*= false*/)
{
  TMatrixD corrmx = GetTotalCorrelationMatrix(cov_area_normalize);
  Int_t state = corrmx.Write(name.c_str());
  return state;
}

TMatrixD MnvH1D::GetSysCorrelationMatrix(
    const std::string& name,
    bool cov_area_normalize /*= false*/ ) const
{
  TMatrixD covmx = GetSysErrorMatrix(name, false, cov_area_normalize);
  const int size = covmx.GetNrows();
  TMatrixD corrmx(size,size);

  for( int i = 0; i < size; ++i )
  {
    for( int k = 0; k < size; ++k )
    {
      corrmx[i][k] = ( covmx[i][i] == 0. || covmx[k][k] == 0. ) ? 0. : covmx[i][k]/sqrt(covmx[i][i]*covmx[k][k]);
    }
  }
  return corrmx;
}

//--------------------------------------------------------------------------------
// Push a Covariance Matrix to the Systematic/Statistical fSys/StatErrorMatrix map
//--------------------------------------------------------------------------------
bool MnvH1D::PushCovMatrix(
    const std::string& name,
    TMatrixD covmx, 
    bool cov_area_normalize /*=false*/)
{

  if ( HasEnding(name,"_asShape") )
  {   
    std::cout << "Warning [MnvH1D::pushCovMatrix] : You cannot push a Covariance matrix with prefix: \"_asShape\", at the end. Doing Nothing."<< std::endl;
    std::cout << "\nTo Add an \"only shape\" error Matrix, set the last boolean in this constructor to true:\n        PushCovMatrix(\"errorMatrix_name\", TMatrixD errorMatrix, true).\n  "<< std::endl;
    if(fStrict){
      assert(0);
    }

    return false;
  }

  const std::string name_condition = (cov_area_normalize)? "_asShape" : "";
  const std::string fname = name + name_condition ;

  if ( HasRemovedErrorMatrix( fname ) )
  {
    std::cout << "Warning [MnvH1D::pushCovMatrix] : The systematic you are trying to push is on the Removed Systematics Map List, if you want to re-add this systematic. \nPlease use the function: \nMnvH1D::UnRemoveSysErrorMatrix( std::string name ) "<< std::endl;
    return false;
  }

  const int highBin = GetNbinsX() + 1;

  if (covmx.GetNrows()!=highBin+1 || covmx.GetNcols()!=highBin+1 )
  {
    std::cout << "Warning [MnvH1D::pushCovMatrix] : The pushed covariance matrix dimensions are incorrect (it should be a " << highBin+1 << "x" << highBin+1 << " matrix). Doing nothing." <<std::endl;
    if(fStrict){
      assert(0);
    }
    return false;
  }

  // Make sure there are no Covariance Matrices with this name already
  if( HasErrorMatrix( fname ) )
  {
    std::cout << "Warning [MnvH1D::PushCovMatrix] : There is a Matrix with name \"" << fname << "\" in " << GetName() << "already.  Doing nothing." << std::endl;
    return false;
  }
  else if ( HasErrorBand( name ) )
  {
    std::cout << "Warning [MnvH1D::PushCovMatrix] : " << name << " is already in either a Vertical or Lateral Error" << std::endl;
    return false;
  }

  // Adding matrix to fSysErrorMatrix
  //! @todo can't we just use the copy constructor when creating a new TMatrixD?
  TMatrixD* temp = new TMatrixD(covmx.GetNrows(), covmx.GetNcols() );
  *temp = covmx;
  fSysErrorMatrix[fname] = temp;

  return true;
}

//------------------------------------------------------------------------
// Modify statistical uncertainty based on studies by DGR related to unfolding and finite MC sizes. See docDB 28992,28899
//------------------------------------------------------------------------

void MnvH1D::ModifyStatisticalUnc(double factor, std::string covmatrixname){

  
  TMatrixD cov = *this->PopSysErrorMatrix(covmatrixname);
  if(this->HasErrorMatrix("unfoldingCov") && covmatrixname!="unfoldingCov"){
    throw("Your MnvH object has both the default covariance (unfoldingCov) and a custom matrix you are modifying. This is probably not correct! There should be only one migration covariance matrix and MnvUnfold has been modified to provide you one");
  }

  double uncmod = (1+1/factor);//To apply to the cov
  double sqrtmod = sqrt(uncmod);//To apply to the diagonal error stored in the CV TH

  //covariance first
  cov*=uncmod;
  this->PushCovMatrix(covmatrixname,cov,false);
  
  //TH
  for(int i=0;i<this->GetNbinsX()+2;i++) this->SetBinError(i,this->GetBinError(i)*sqrtmod);
}


//------------------------------------------------------------------------
// Get a Specific Covariance Matrix from the map
//------------------------------------------------------------------------
TMatrixD MnvH1D::GetSysErrorMatrix(const std::string& name, bool asFrac /*= false*/, bool cov_area_normalize /*= false*/) const
{
  if (HasEnding(name,"_asShape") )
    std::cout << "Warning [MnvH1D::GetSysErrorMatrix]: You are calling the error Matrix: " << name <<".\nAssuming the error Band wanted is: " << name.substr(0,name.length()-8) << " with cov_area_normalize = true" << std::endl;

  const std::string name_condition = ( (cov_area_normalize) && !(HasEnding(name,"_asShape")) )?  "_asShape" : "";
  const std::string fname = name + name_condition ; 
  const std::string errName = HasEnding(fname,"_asShape") ? fname.substr(0,fname.length()-8) : fname;

  const int highBin = GetNbinsX() + 1;
  const int lowBin = 0;
  TMatrixD covmx(highBin+1,highBin+1);

  if ( HasErrorMatrix( fname ) )
    covmx = *(fSysErrorMatrix.find(fname)->second);
  else if( fLatErrorBandMap.find( errName ) != fLatErrorBandMap.end() )
  {
    std::map<std::string, MnvLatErrorBand*>::const_iterator it = fLatErrorBandMap.find( errName );
    covmx = it->second->CalcCovMx( ( HasEnding(fname,"_asShape") ) );
  }
  else if( fVertErrorBandMap.find( errName ) != fVertErrorBandMap.end() )
  {
    std::map<std::string, MnvVertErrorBand*>::const_iterator it = fVertErrorBandMap.find( errName );
    covmx = it->second->CalcCovMx( ( HasEnding(fname,"_asShape") ) );
  }
  else if( fUncorrErrorMap.find( errName ) != fUncorrErrorMap.end() )
  {
    std::map<std::string, TH1D*>::const_iterator it = fUncorrErrorMap.find( errName );
    covmx = MnvHist::GetErrorsAsMatrix( it->second );
  }
  else
    std::cout << "Warning [MnvH1D::GetSysErrorMatrix]: There is no Covariance Matrix with name " << fname << " in " << GetName() <<  ".Returning and empty Matrix." << std::endl;

  if (asFrac)
  {
    for( int i = lowBin; i <= highBin; ++i )
    {
      for( int k = i; k <= highBin; ++k )
      { 
        //Gettting the the CV value for bin i
        const double cv_i = GetBinContent(i);
        const double cv_k = GetBinContent(k);

        covmx[i][k]= ( (cv_i != 0.) && (cv_k !=0.) )? covmx[i][k]/(cv_i * cv_k) : 0.;
        covmx[k][i]=covmx[i][k];
      }
    }
  }

  return covmx;
}

bool MnvH1D::RemoveSysErrorMatrix(const std::string& name)
{
  const std::string shapeName = name + "_asShape";
  if ( !HasErrorMatrix( name) && !HasErrorMatrix( shapeName ) )
  {
    std::cout << "Warning [MnvH1D::RemoveSysErrorMatrix] : The Systematic \"" << name << "\" was not found on the systematic map in neither absolute nor shape component.\nNotice Vertical/Lateral Errors cannot be removed. Doing nothing." << std::endl;
    return false;
  }

  if ( HasErrorMatrix( name ) )
  {  
    // Removing systematic from fSysErrorMatrix
    // and Adding it to fRemovedSysErrorMatrix
    fRemovedSysErrorMatrix[name] = fSysErrorMatrix[name];
    fSysErrorMatrix.erase( fSysErrorMatrix.find( name ) );
  }

  if ( HasErrorMatrix(shapeName) )
  {
    // Doing the same for shape component matrix (if any)
    fRemovedSysErrorMatrix[shapeName] = fSysErrorMatrix[shapeName];
    fSysErrorMatrix.erase( fSysErrorMatrix.find( shapeName ) );
  }

  return true;
}

bool MnvH1D::UnRemoveSysErrorMatrix(const std::string& name)
{
  // Moving Removed Systematics (Absolute and Shape) from
  // fRemovedSysErrorMatrix back to fSysErrorMatrix

  const std::string shapeName = name + "_asShape";
  if ( !HasRemovedErrorMatrix(name) && !HasRemovedErrorMatrix(shapeName) )
  {
    std::cout << "Warning [MnvH1D::UnRemoveSysErrorMatrix] : The Systematic \"" << name << "\" was not found on the map of removed systematics neither in absolute nor shape component. Doing nothing." << std::endl;
    return false;
  }

  if( HasRemovedErrorMatrix( name ) )
  {
    fSysErrorMatrix[name] = fRemovedSysErrorMatrix[name];
    fRemovedSysErrorMatrix.erase( fRemovedSysErrorMatrix.find( name ) );
  }

  if( HasRemovedErrorMatrix( shapeName ) )
  {
    fSysErrorMatrix[name] = fRemovedSysErrorMatrix[shapeName];
    fRemovedSysErrorMatrix.erase( fRemovedSysErrorMatrix.find( shapeName ) );
  }
  return true;
}

void MnvH1D::ClearSysErrorMatrices()
{
  // Loop over all SysErrorMatrices and delete them, even the RemovedSysErrorMatrices.
  for( std::map<std::string, TMatrixD*>::const_iterator i = fSysErrorMatrix.begin(); i != fSysErrorMatrix.end(); ++i )
  {
    if( 0 != i->second )
      delete i->second;
  }
  fSysErrorMatrix.clear();

  for( std::map<std::string, TMatrixD*>::const_iterator i = fRemovedSysErrorMatrix.begin(); i != fRemovedSysErrorMatrix.end(); ++i )
  {
    if( 0 != i->second )
      delete i->second;
  }
  fRemovedSysErrorMatrix.clear();
}

TMatrixD MnvH1D::GetStatErrorMatrix( bool asFrac /* =false */ ) const
{
  const int highBin = GetNbinsX() + 1;
  const int lowBin = 0;
  TMatrixD covmx(highBin+1,highBin+1);

  // stat error
  for (int iBin=lowBin; iBin<=highBin; ++iBin )
    covmx[iBin][iBin]= GetBinError(iBin);

  if (asFrac)
  { 
    for( int iBin = lowBin; iBin <= highBin; ++iBin )
    {
      double binCon = GetBinContent(iBin);
      double binErr = ( binCon != 0. ? covmx[iBin][iBin]/binCon : 0. ); 
      covmx[iBin][iBin] = binErr;
    }
  }

  return covmx*covmx;
}

//------------------------------------------------------------------------
// Get a histogram of the statistical errors only.
//------------------------------------------------------------------------
TH1D MnvH1D::GetStatError( bool asFrac /* = false */ ) const
{
  // Make a copy of this histogram as a TH1D and rename it
  TH1D err( *this );
  err.Reset();
  std::string tmpName( std::string(GetName()) + "_StatError");
  err.SetName( tmpName.c_str() );

  const int highBin = GetNbinsX() + 1;
  const int lowBin = 0;

  // stat error
  for( int iBin = lowBin; iBin <= highBin; ++iBin )
    err.SetBinContent( iBin, GetBinError(iBin) );

  if( asFrac )
  {
    for( int iBin = lowBin; iBin <= highBin; ++iBin )
    {
      double binCon = GetBinContent(iBin);
      double binErr = ( binCon != 0. ? err.GetBinContent(iBin) / binCon : 0. );
      err.SetBinContent( iBin, binErr );
    }
  }

  //if you manually set the MnvH1D min/max it exists on err.
  //use default min/max for err
  err.SetMinimum();
  err.SetMaximum();

  return err;
}

//------------------------------------------------------------------------
// Get an MnvH1D which has its bin content and errors normalized to bin width so it looks smooth
//------------------------------------------------------------------------
MnvH1D MnvH1D::GetBinNormalizedCopy( Double_t normBinWidth /* = fNormBinWidth */ ) const
{
  // If normBinWidth is negative, use the default bin width for this MnvH1D
  if( normBinWidth <= 0 )
    normBinWidth = fNormBinWidth;

  MnvH1D rval(*this);

  // If the default bin width for this MnvH1D is negative and the user wanted to use that default, just return self (already normalized or normalization not appropriate).
  // Otherwise, normalize to the this bin width
  if( normBinWidth > 0 )
    rval.Scale( normBinWidth, "width" );

  return rval;
}


Double_t MnvH1D::GetAreaNormFactor(const MnvH1D * h_data) const
{
  //considering overflow
  const int highBin_data = h_data->GetNbinsX()+1;
  const int highBin_mc   = this->GetNbinsX()+1;
  //considering underflow
  const int lowBin = 0;
  Double_t area_scale = 1.;

  if( ! CheckConsistency( this, h_data ) ){
    Warning("GetAreaNormFactor", "Data and MC axes do not match" );
  }
  if (this->Integral( lowBin, highBin_mc )==0 ){
    Warning("GetAreaNormFactor","MC Area Histogram is zero. No Scale Factor calculated");
  }
  else {
    area_scale = h_data->Integral( lowBin, highBin_data ) / this->Integral( lowBin, highBin_mc ) ;
  }
  return area_scale;

}

//=======================================================================
// DrawBinNormalized
//=======================================================================
MnvH1D* MnvH1D::DrawBinNormalized( Option_t* option /* = ""*/, Double_t normBinWidth /* = -1 */ ) const
{
  // Create a (probably) unique name for the clone to avoid potential memory leaks
  static int nNormDraws = 0;
  const char *cloneName = Form( "%s_normClone_%d", GetName(), ++nNormDraws );

  // Create the clone from the bin normalized version of this MnvH1D.
  // Clone is a new object, which is necessary so the drawn histogram stays drawn after the function ends.
  MnvH1D *clone = (MnvH1D*)this->GetBinNormalizedCopy( normBinWidth ).Clone( cloneName );

  // Call draw with the clone
  clone->Draw( option );

  // Return the clone so that the caller can delete the object when it is no longer needed
  return clone;
}


//======================================================================
// Replacements of ROOT versions of functions
//======================================================================

void MnvH1D::Scale( Double_t c1 /*= 1.*/, Option_t* option /*=""*/, Bool_t scaleUniv /*=true*/ )
{
  // Scale yourself using TH1D::Scale
  this->TH1D::Scale( c1, option );

  // ... but only scale the other universes if not disabled
  if (!scaleUniv)
    return;

  // Scale the lateral error bands
  for( std::map<std::string, MnvLatErrorBand*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
    it->second->Scale( c1, option );

  // Scale the vertical error bands
  for( std::map<std::string, MnvVertErrorBand*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
    it->second->Scale( c1, option );

  // Scale the uncorr errors
  for( std::map<std::string, TH1D*>::iterator it = fUncorrErrorMap.begin(); it != fUncorrErrorMap.end(); ++it )
    it->second->Scale( c1, option );


  // Scale special covariance matrices including the removed ones
  for( std::map<std::string, TMatrixD*>::const_iterator it = fSysErrorMatrix.begin(); it != fSysErrorMatrix.end(); ++it ){
#ifdef MNV1DBG
    std::cout << "MnvH1D::Scale covmatrix " << it->first << " scaled by " << c1 << " " << option << std::endl;
#endif
    if( 0 != it->second ){
      if ( !HasEnding(it->first, "_asShape") ){
        *(it->second) *=c1*c1;
        TMatrixD * m = it->second;
        // implement bin width correction!
        if (TString(option).Contains("width")){
          int nx = GetNbinsX()+2;
          for (int ix = 0; ix < nx; ix++){
            double widix = GetXaxis()->GetBinWidth(ix);
            for (int jx = 0; jx < nx; jx++){
              double widjx = GetXaxis()->GetBinWidth(jx);
              if (widix*widjx!= 0.0){
                (*m)[ix][jx] /=(widix*widjx);
              }
              else{
                (*m)[ix][jx] = 0.0;
              }
              
            } //jx
          } //ix
        } //ok
      } // exists
    } // loop over errors
  }
    
    for( std::map<std::string, TMatrixD*>::const_iterator it = fRemovedSysErrorMatrix.begin(); it != fRemovedSysErrorMatrix.end(); ++it ){
      if( 0 != it->second ){
        if ( !HasEnding(it->first, "_asShape") ){
          *(it->second) *=c1*c1;
          TMatrixD * m = it->second;
          // implement bin width correction!
          if (TString(option).Contains("width")){
            int nx = GetNbinsX()+2;
            for (int ix = 0; ix < nx; ix++){
              double widix = GetXaxis()->GetBinWidth(ix);
              for (int jx = 0; jx < nx; jx++){
                double widjx = GetXaxis()->GetBinWidth(jx);
                if (widix*widjx!= 0.0){
                  (*m)[ix][jx] /=(widix*widjx);
                }
                else{
                  (*m)[ix][jx] = 0.0;
                }
                
              } //jx
            } //ix
          } //ok
        } // exists
      } // loop over errors
    }
  }

void MnvH1D::Divide( const MnvH1D* h1, const MnvH1D* h2, Double_t c1 /*= 1*/, Double_t c2 /*= 1*/, Option_t* option /*=""*/ )
{
  //! @todo Would love to return a bool in MnvH1D::Divide, but we want this Divide to override TH1's and that is void
#ifdef MNV1DBG
  std::cout << "MnvH1D::Divide h1,h2,c1,c2,option " << h1->GetName() <<" , "<<h2->GetName() << ", " << c1 << ", " << c2 << ", " << option << std::endl;
#endif
  // Call the TH1D Divide
  bool status = this->TH1D::Divide( (TH1D*)h1, (TH1D*)h2, c1, c2, option );
#ifdef MNVDBG1
  if (!status) {
    std::cout << " MnvH1D::Divide  had problems " << h1->GetName() <<" , "<<h2->GetName() << ", " << c1 << ", " << c2 << ", " << option << std::endl;
    abort();
  }
#endif

  // Call Divide for all lateral error bands
  for( std::map<std::string, MnvLatErrorBand*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
  {
    const MnvLatErrorBand* err1 = h1->GetLatErrorBand( it->first );
    const MnvLatErrorBand* err2 = h2->GetLatErrorBand( it->first );
    if( !err1 || !err2 )
    {
      Error("Divide", Form("Could not divide MnvH1Ds because they all don't have the %s MnvLatErrorBand", it->first.c_str()) );
      return;
    }
    status =it->second->Divide( err1, err2, c1, c2, option );
  #ifdef MNVDBG1
    if (!status) {
      std::cout << " MnvH1D::Divide  had problems " << h1->GetName() <<" , "<<h2->GetName() << ", " << c1 << ", " << c2 << ", " << option << std::endl;
      abort();
    }
  #endif
  }

  // Divide the vertical error bands
  for( std::map<std::string, MnvVertErrorBand*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
  {
    const MnvVertErrorBand* err1 = h1->GetVertErrorBand( it->first );
    const MnvVertErrorBand* err2 = h2->GetVertErrorBand( it->first );
    if( !err1 || !err2 )
    {
      Error("Divide", Form("Could not divide MnvH1Ds because they all don't have the %s MnvVertErrorBand", it->first.c_str()) );
      return;
    }
    status = it->second->Divide( err1, err2, c1, c2, option );
  #ifdef MNVDBG1
    if (!status) {
      std::cout << " MnvH1D::Divide  had problems " << h1->GetName() <<" , "<<h2->GetName() << ", " << c1 << ", " << c2 << ", " << option << std::endl;
      abort();
    }
  #endif
  }

  // Divide the uncorr errors
  for( std::map<std::string, TH1D*>::iterator it = fUncorrErrorMap.begin(); it != fUncorrErrorMap.end(); ++it )
  {
    const TH1D* err1 = h1->GetUncorrError( it->first );
    const TH1D* err2 = h2->GetUncorrError( it->first );
    if( !err1 || !err2 )
    {
      Error("Divide", Form("Could not divide MnvH1Ds because they all don't have the %s uncorrelated error", it->first.c_str()) );
      return;
    }
    status =  it->second->Divide( err1, err2, c1, c2, option );
  #ifdef MNVDBG1
    if (!status) {
      std::cout << " MnvH1D::Divide  had problems " << h1->GetName() <<" , "<<h2->GetName() << ", " << c1 << ", " << c2 << ", " << option << std::endl;
      abort();
    }
  #endif
  }


//  // Do we have special errors in the systemics?
//  if ( !( fSysErrorMatrix.empty() ) || !( fRemovedSysErrorMatrix.empty() ) )
//  {
//    Warning("MnvH1D::Divide", "Customized error matrices were found (errors that come from neither Vertical nor Lateral Error Bands). They will be cleared.");
//    if(fStrict){
//      assert(0);
//    }
//    ClearSysErrorMatrices( );
//  }
  
  
  for( std::map<std::string, TMatrixD*>::const_iterator i = fSysErrorMatrix.begin(); i != fSysErrorMatrix.end(); ++i ){
    
    std::cout << " MnvH1D::Divide " << h1->GetName() << " by " << h2->GetName() << "error matrix from h1 is scaled by " << c1 << "^2/" << c2 << "^2" << std::endl;
    // pull out the old matrix so I can replace it.
   
    
    TMatrixD newmatrix = h1->GetSysErrorMatrix(i->first);
    
#ifdef MNV1DBG
    std::cout << " MnvH1D  Covmatrix before " << std::endl;
    newmatrix.Print();
#endif
    for (int i = 0; i < newmatrix.GetNrows(); i++){
      double hi = h2->GetBinContent(i);
      for (int j = 0; j < newmatrix.GetNcols(); j++){
        double hj = h2->GetBinContent(j);
        if (hj*hi!=0){
          newmatrix[i][j] /=  hi*hj*c2*c2;
          newmatrix[i][j] *=c1*c1;
        }
        else{
          newmatrix[i][j] = 0;
        }
      }
    }
    // if the one you are dividing by also has this matrix do a little error analysis
    if(h2->HasErrorMatrix(i->first)){
      TMatrixD new2 = h2->GetSysErrorMatrix(i->first);
      for (int i = 0; i < newmatrix.GetNrows(); i++){
        double hi = h2->GetBinContent(i);
        for (int j = 0; j < newmatrix.GetNcols(); j++){
          double hj = h2->GetBinContent(j);
          if (hj*hi!=0){
            newmatrix[i][j] -= new2[i][j]*c1*c1*GetBinContent(i)*GetBinContent(j)/(hi*hi*hj*hj*c2*c2);
          }
        }
      }
    }
#ifdef MNV1DBG
    std::cout << " MnvH1D covmatrix after " <<  std::endl;
    newmatrix.Print();
    std::cout << " MnvH1D::Divide replace CovErrorMatrix with new one " << i->first <<  std::endl;
#endif
    this->FillSysErrorMatrix(i->first,newmatrix);
    
  }

  return;
}


Bool_t MnvH1D::Divide( TF1 *f1, Double_t c1 )
{
	bool success = true;

	// cv histo first
	success = success && this->TH1D::Divide(f1, c1);

	// Call Divide for all lateral error bands
	for( std::map<std::string, MnvLatErrorBand*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
		success = success && it->second->Divide( f1, c1 );

	// Call Divide for all vertical error bands
	for( std::map<std::string, MnvVertErrorBand*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
		success = success && it->second->Divide( f1, c1 );

	return success;
}

//=====================================================================
// Divide the the numerator and its error bands by the same denominator
//=====================================================================
void MnvH1D::DivideSingle( const MnvH1D* h1, const TH1* h2, Double_t c1 /*= 1*/, Double_t c2 /*= 1*/, Option_t* option /*=""*/ )
{
  // Call the TH1D Divide
  bool status = this->TH1D::Divide( (TH1D*)h1, h2, c1, c2, option );
  #ifdef MNVDBG1
    if (!status) {
      std::cout << " MnvH1D::Divide  had problems " << h1->GetName() <<" , "<<h2->GetName() << ", " << c1 << ", " << c2 << ", " << option << std::endl;
      abort();
    }
  #endif

  // Call Divide for all lateral error bands
  for( std::map<std::string, MnvLatErrorBand*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
  {
    const MnvLatErrorBand* err1 = h1->GetLatErrorBand( it->first );
    if( !err1  )
    {
      Error("Divide", Form("Could not divide MnvH1Ds because they all don't have the %s MnvLatErrorBand", it->first.c_str()) );
      return;
    }
    status = it->second->DivideSingle( err1, h2, c1, c2, option );
  #ifdef MNVDBG1
    if (!status) {
      std::cout << " MnvH1D::Divide  had problems " << h1->GetName() <<" , "<<h2->GetName() << ", " << c1 << ", " << c2 << ", " << option << std::endl;
      abort();
    }
  #endif
  }

  // Divide the vertical error bands
  for( std::map<std::string, MnvVertErrorBand*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
  {
    const MnvVertErrorBand* err1 = h1->GetVertErrorBand( it->first );
    if( !err1 )
    {
      Error("Divide", Form("Could not divide MnvH1Ds because they all don't have the %s MnvVertErrorBand", it->first.c_str()) );
      return;
    }
    status = it->second->DivideSingle( err1, h2, c1, c2, option );
  #ifdef MNVDBG1
    if (!status) {
      std::cout << " MnvH1D::Divide  had problems " << err1->GetName() <<" , "<<h2->GetName() << ", " << c1 << ", " << c2 << ", " << option << std::endl;
      abort();
    }
  #endif
  }

  // Divide the uncorr errors
  for( std::map<std::string, TH1D*>::iterator it = fUncorrErrorMap.begin(); it != fUncorrErrorMap.end(); ++it )
  {
    const TH1D* err1 = h1->GetUncorrError( it->first );
    if( !err1 )
    {
      Error("Divide", Form("Could not divide MnvH1Ds because they all don't have the %s uncorrelated error", it->first.c_str()) );
      return;
    }
    status = it->second->Divide( err1, h2, c1, c2, option );
  #ifdef MNVDBG1
    if (!status) {
      std::cout << " MnvH1D::Divide  had problems " << err1->GetName() <<" , "<<h2->GetName() << ", " << c1 << ", " << c2 << ", " << option << std::endl;
      abort();
    }
  #endif
  
  }

  for( std::map<std::string, TMatrixD*>::const_iterator i = fSysErrorMatrix.begin(); i != fSysErrorMatrix.end(); ++i ){
    
    std::cout << " Dividing " << h1->GetName() << " by " << h2->GetName() << "error matrix from h1 is scaled by h2^2" << std::endl;
    
    TMatrixD newmatrix = h1->GetSysErrorMatrix(i->first);
    
    for (int i = 0; i < newmatrix.GetNrows(); i++){
      double hi = h2->GetBinContent(i);
      for (int j = 0; j < newmatrix.GetNcols(); j++){
        double hj = h2->GetBinContent(j);
        if (hj*hi!=0){
          newmatrix[i][j] /=  hi*hj*c2*c2;
          newmatrix[i][j] *= c1*c1;
        }
        else{
          newmatrix[i][j] = 0;
        }
      }
    }
    #ifdef MNV1DBG
        std::cout << " MnvH1D::Divide(TH1) replace CovErrorMatrix with new one " << i->first <<  std::endl;
    #endif
    this->PushCovMatrix(i->first,newmatrix);
    
  }

//  // Do we have special errors in the systemics?
//  if ( !( fSysErrorMatrix.empty() ) || !( fRemovedSysErrorMatrix.empty() ) )
//  {
//  Warning( "MnvH1D::DivideSingle", "Customized error matrices were found (errors that come from neither Vertical nor Lateral Error Bands) . They will be cleared.");
//    if(fStrict){
//      assert(0);
//    }
//    ClearSysErrorMatrices( );
//  }

  return;
}


void MnvH1D::Multiply( const MnvH1D* h1, const MnvH1D* h2, Double_t c1 /*= 1*/, Double_t c2 /*= 1*/ )
{
  // Call the TH1D Multiply 
  this->TH1D::Multiply( (TH1D*)h1, (TH1D*)h2, c1, c2 );

  // Call Multiply for all lateral error bands
  for( std::map<std::string, MnvLatErrorBand*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
  {
    const MnvLatErrorBand* err1 = h1->GetLatErrorBand( it->first );
    const MnvLatErrorBand* err2 = h2->GetLatErrorBand( it->first );
    if( !err1 || !err2 )
    {
      Error("Multiply", Form("Could not multiply MnvH1Ds because they all don't have the %s MnvLatErrorBand", it->first.c_str()) );
      return;
    }
    it->second->Multiply( err1, err2, c1, c2 );
  }

  // multiply the vertical error bands
  for( std::map<std::string, MnvVertErrorBand*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
  {
    const MnvVertErrorBand* err1 = h1->GetVertErrorBand( it->first );
    const MnvVertErrorBand* err2 = h2->GetVertErrorBand( it->first );
    if( !err1 || !err2 )
    {
      Error("Multiply", Form("Could not multiply MnvH1Ds because they all don't have the %s MnvVertErrorBand", it->first.c_str()) );
      return;
    }
    it->second->Multiply( err1, err2, c1, c2 );
  }

  // multiply the uncorr errors
  for( std::map<std::string, TH1D*>::iterator it = fUncorrErrorMap.begin(); it != fUncorrErrorMap.end(); ++it )
  {
    const TH1D* err1 = h1->GetUncorrError( it->first );
    const TH1D* err2 = h2->GetUncorrError( it->first );
    if( !err1 || !err2 )
    {
      Error("Multiply", Form("Could not multiply MnvH1Ds because they all don't have the %s uncorrelated error", it->first.c_str()) );
      return;
    }
    it->second->Multiply( err1, err2, c1, c2 );
  }

  for( std::map<std::string, TMatrixD*>::const_iterator i = fSysErrorMatrix.begin(); i != fSysErrorMatrix.end(); ++i ){
    
    std::cout << " Multiply " << h1->GetName() << " by " << h2->GetName() << "error matrix from h1 is scaled by h2^2" << std::endl;
    
    TMatrixD newmatrix = h1->GetSysErrorMatrix(i->first);
    
    for (int i = 0; i < newmatrix.GetNrows(); i++){
      double hi = h2->GetBinContent(i);
      for (int j = 0; j < newmatrix.GetNcols(); j++){
        double hj = h2->GetBinContent(j);
        if (hj!=0){
          newmatrix[i][j] *=  hi*hj*c1*c2*c1*c2;
        }
        else{
          newmatrix[i][j] = 0;
        }
      }
    }
    #ifdef MNV1DBG
        std::cout << " MnvH1D::Multiply replace CovErrorMatrix with new one " << i->first <<  std::endl;
    #endif
    this->PushCovMatrix(i->first,newmatrix);
    
  }
//  // Do we have special errors in the systemics?
//  if ( !( fSysErrorMatrix.empty() ) || !( fRemovedSysErrorMatrix.empty() ) )
//  {
//    Warning("MnvH1D::Multiply", "Customized error matrices were found (errors that come from neither Vertical nor Lateral Error Bands) . They will be cleared.");
//    if(fStrict){
//      assert(0);
//    }
//    ClearSysErrorMatrices( );
//  }

  return;
}



//=====================================================================
// Multiply the the hist and its error bands by the same denominator
//=====================================================================
void MnvH1D::MultiplySingle( const MnvH1D* h1, const TH1* h2, const Double_t c1 /*= 1*/, const Double_t c2 /*= 1*/ )
{
  // Call the TH1D Multiply
  this->TH1D::Multiply( (TH1D*)h1, h2, c1, c2 );

  // Call Multiply for all lateral error bands
  for( std::map<std::string, MnvLatErrorBand*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
  {
    const MnvLatErrorBand* err1 = h1->GetLatErrorBand( it->first );
    if( !err1  )
    {
      Error("MultiplySingle", Form("Could not multiply MnvH1Ds because they all don't have the %s MnvLatErrorBand", it->first.c_str()) );
      return;
    }
    it->second->MultiplySingle( err1, h2, c1, c2 );
  }

  // Call Multiply for all vertical error bands
  for( std::map<std::string, MnvVertErrorBand*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
  {
    const MnvVertErrorBand* err1 = h1->GetVertErrorBand( it->first );
    if( !err1 )
    {
      Error("MultiplySingle", Form("Could not multiply MnvH1Ds because they all don't have the %s MnvVertErrorBand", it->first.c_str()) );
      return;
    }
    it->second->MultiplySingle( err1, h2, c1, c2 );
  }

  // multiply the uncorr errors
  for( std::map<std::string, TH1D*>::iterator it = fUncorrErrorMap.begin(); it != fUncorrErrorMap.end(); ++it )
  {
    const TH1D* err1 = h1->GetUncorrError( it->first );
    if( !err1 )
    {
      Error("MultiplySingle", Form("Could not multiply MnvH1Ds because they all don't have the %s uncorrelated error", it->first.c_str()) );
      return;
    }
    it->second->Multiply( err1, h2, c1, c2 );
  }

  for( std::map<std::string, TMatrixD*>::const_iterator i = fSysErrorMatrix.begin(); i != fSysErrorMatrix.end(); ++i ){
    
    std::cout << " MultiplySingle " << h1->GetName() << " by " << h2->GetName() << "error matrix from h1 is scaled by h2^2" << std::endl;
    
    TMatrixD newmatrix = h1->GetSysErrorMatrix(i->first);
    
    for (int i = 0; i < newmatrix.GetNrows(); i++){
      double hi = h2->GetBinContent(i);
      for (int j = 0; j < newmatrix.GetNcols(); j++){
        double hj = h2->GetBinContent(j);
        if (hj!=0){
          newmatrix[i][j] *=  hi*hj*c1*c2*c1*c2;
        }
        else{
          newmatrix[i][j] = 0;
        }
      }
    }
    #ifdef MNV1DBG
        std::cout << " MnvH1D::MultiplySingle replace CovErrorMatrix with new one " << i->first <<  std::endl;
    #endif
    this->PushCovMatrix(i->first,newmatrix);
    
  }
//  // Do we have special errors in the systemics?
//  if ( !( fSysErrorMatrix.empty() ) || !( fRemovedSysErrorMatrix.empty() ) )
//  {
//    Warning( "MnvH1D::MultiplySingle", "Customized error matrices were found (errors that come from neither Vertical nor Lateral Error Bands) . They will be cleared.");
//    if(fStrict){
//      assert(0);
//    }
//    ClearSysErrorMatrices( );
//  }

  return;
}

// In ROOT 5.34, the signature of TH1::Add changed to return a
// bool. When we no longer need to build against ROOT 5.30, this
// horrible preprocessor jiggerypokery can go away
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,34,0)
Bool_t MnvH1D::Add( const TH1* h1, const Double_t c1 /*= 1.*/ )
#else
void MnvH1D::Add( const TH1* h1, const Double_t c1 /*= 1.*/ )
#endif
{
  // Try to cast the input TH1 to a MnvH1D
  // If cast doesn't work, I don't know what to do. How do we add a MnvH1D and TH1?
  const MnvH1D *mnv1 = dynamic_cast<const MnvH1D*>(h1);

  if( mnv1 )
  {
    // If cast works...
    //! @todo do the consistency checks on axes before adding

    //save orig stat uncertainty in case we need it for averaging uncorrelated errors
    const TH1D origStatErr = this->GetStatError();

    // Add as a TH1D
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,34,0)
    if(!this->TH1D::Add( h1, c1 )) return kFALSE;
#else
    this->TH1D::Add( h1, c1 );
#endif

    // Call Add for all lateral error bands
    for( std::map<std::string, MnvLatErrorBand*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
    {
      Bool_t ok = false;
      const MnvLatErrorBand* err1 = mnv1->GetLatErrorBand( it->first );
      if( !err1  )
      {
           Warning("MnvH1D::Add", Form("Additive MnvH1D '%s' lacks %s MnvLatErrorBand.  Add central value to all universes.", mnv1->GetName(), it->first.c_str()) );
        throw std::runtime_error("missing lat error band in Add");
        ok = it->second->AddSingle( h1, c1 );
      }
      else
        ok = it->second->Add( err1, c1 );

      if( ! ok )
      {
        Error("Add", Form("Could not add MnvH1Ds because histogram add failed for MnvLatErrorBand %s ", it->first.c_str() ) );
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,34,0)
        return kFALSE;
#else
        return;
#endif
      }
    }//done adding Lat errors

    // Call Add for all vertical error bands
    for( std::map<std::string, MnvVertErrorBand*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
    {
      Bool_t ok = false;
      const MnvVertErrorBand* err1 = mnv1->GetVertErrorBand( it->first );
      if( !err1  )
      {
        throw std::runtime_error("missing vert error band in Add");
        Warning("MnvH1D::Add", Form("Additive MnvH1D lacks %s MnvVertErrorBand.  Add central value to all universes.", it->first.c_str()) );
        ok = it->second->AddSingle( h1, c1 );
      }
      else
        ok = it->second->Add( err1, c1 );

      if( ! ok )
      {
        Error("Add", Form("Could not add MnvH1Ds because histogram add failed for MnvVertErrorBand %s ", it->first.c_str() ) );
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,34,0)
        return kFALSE;
#else
        return;
#endif
      }
    }//done adding Vert errors

    //call add for all uncorrelated errors
    for( std::map<std::string, TH1D*>::iterator it = fUncorrErrorMap.begin(); it != fUncorrErrorMap.end(); ++it )
    {
      //note - TH1D::Add is void so we just hope that it works
      const TH1D* err1 = mnv1->GetUncorrError( it->first );
      if( !err1  )
      {
        //todo: a better way would be to add content but keep error equal to the one hist that has error
        Warning("MnvH1D::Add", Form("Additive MnvH1D lacks %s uncorrelated error.  Add central value to all universes.", it->first.c_str()) );
        it->second->TH1D::Add( (const TH1D*)h1, c1 );
      }
      else
      {
        //if this is an average add then we need to force correct behavior
        if( it->second->TestBit( TH1::kIsAverage ) && err1->TestBit( TH1::kIsAverage ) )
        {
          int nbins = it->second->GetNbinsX()+1;
          for( int ibin = 0; ibin <= nbins; ++ibin )
          {
            //apply same operations to error as were applied to CV, propagage the uncorr error
            double uncorrA = it->second->GetBinError(ibin);
            double uncorrB = err1->GetBinError(ibin);
            double errA = origStatErr.GetBinContent(ibin);
            double errB = h1->GetBinError(ibin);
            double wA = ( 0. < errA ) ? 1./(errA*errA) : 1.E200; //use err=sqrt(w) or very large value
            double wB = ( 0. < errB ) ? c1/(errB*errB) : 1.E200; //use err=sqrt(w) or very large value
            double errPieceA = (uncorrA*wA) / (wA+wB);
            double errPieceB = (uncorrB*wB) / (wA+wB);
            double err = sqrt( TMath::Power(errPieceA,2) +  TMath::Power(errPieceB,2) );
            it->second->SetBinError(ibin,err);
            //content is same as CV
            it->second->SetBinContent(ibin, this->GetBinContent(ibin) );
          }
        }
        else
          it->second->TH1D::Add( (const TH1D*)err1, c1 );
      }

    }//done adding uncorrelated errors

   
    // Do we have special errors in the systemics?
  
  // HMS add for error matrices.
      for( std::map<std::string, TMatrixD*>::const_iterator i = fSysErrorMatrix.begin(); i != fSysErrorMatrix.end(); ++i ){
        
        std::cout << " Adding " << mnv1->GetName() << " to " << GetName() << "error matrix from h1 is scaled by c1" << std::endl;
        
          TMatrixD newmatrix = GetSysErrorMatrix(i->first);
          TMatrixD addmatrix = mnv1->GetSysErrorMatrix(i->first);
          addmatrix*=c1;
          newmatrix += addmatrix;
#ifdef MNV1DBG
        std::cout << "MnvH1D::Add modify CovMatrix " << i->first << std::endl;
#endif
        
          FillSysErrorMatrix(i->first,newmatrix);
        
      }
  
//    if ( !( fSysErrorMatrix.empty() ) || !( fRemovedSysErrorMatrix.empty() ) )
//    {
//      Warning( "MnvH1D::Add(TH1)", "Customized error matrices were found (errors that come from neither Vertical nor Lateral Error Bands). They will be cleared.");
//      if (fStrict){
//        assert(0);
//      }
//      ClearSysErrorMatrices( );
//   }


  }// end if cast to MnvH1D worked
  else
  {
    Error( "MnvH1D::Add", Form("Unable to add histogram to MnvH1D '%s' because it could not be cast to an MnvH1D.  Did nothing.", GetName()) );
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,34,0)
    return kFALSE;
#else
    return;
#endif
  }

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,34,0)
  return kTRUE;
#else
  return;
#endif
}

Bool_t MnvH1D::AddSingle( const TH1* h1, const Double_t c1 )
{
	// add to CV histogram
	this->TH1D::Add(h1, c1);

	// now handle error bands
	for( std::map<std::string, MnvVertErrorBand*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
	{
		Bool_t success = it->second->AddSingle(h1, c1);
		if (!success)
			return success;
	}
	for( std::map<std::string, MnvLatErrorBand*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
	{
		Bool_t success = it->second->AddSingle(h1, c1);
		if (!success)
			return success;
	}

  // as you are adding a histogram, there is nothing to be done for the covariance matrices
	return true;
}

TH1* MnvH1D::Rebin(  Int_t ngroup /*= 2*/, const char *newname /*= ""*/, const Double_t *xbins /*= 0*/ )
{
  MnvH1D * clone = NULL;
  if ( (newname && strlen(newname) > 0) || xbins )
    clone = new MnvH1D(*this);

  // Call TH1D's Rebin and store the return value, which we return if everything went well
  TH1 * rval;
  if (clone)
  {
    // we'll be creating a temporary copy generated by the Rebin() method.
    rval = this->TH1D::Rebin( ngroup, newname, xbins );
    if (rval)
    {
      dynamic_cast<TH1D*>(rval)->Copy(*clone);  // for some reason TH1::Copy is protected, but TH1D::Copy isn't
      delete rval;
      rval = clone;
    }
  }
  else
    rval = this->TH1D::Rebin( ngroup );

  //we only need the rebin warning once
  const int oldVerbosity = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kError;

  // Call Rebin for all lateral error bands
  for( std::map<std::string, MnvLatErrorBand*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
  {
    Bool_t ok;
    if (clone)
    {
      MnvLatErrorBand * newEB = it->second->Rebin( ngroup, xbins );
      ok = (newEB != NULL);
      if (ok)
      {
        MnvLatErrorBand * old = clone->PopLatErrorBand(it->first);
        clone->PushErrorBand(it->first, newEB);
        delete old;
      }
    }
    else
      ok = (it->second->Rebin( ngroup ) != NULL);

    if( ! ok )
    {
      Error("Rebin", Form("Could not rebin MnvH1Ds because Rebin failed for MnvLatErrorBand %s ", it->first.c_str() ) );
      gErrorIgnoreLevel = oldVerbosity;
      return NULL;
    }
  }//done adding Lat errors

  // Call Rebin for all vertical error bands
  for( std::map<std::string, MnvVertErrorBand*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
  {
    Bool_t ok;
    if (clone)
    {
      MnvVertErrorBand * newEB = it->second->Rebin( ngroup, xbins );
      ok = (newEB != NULL);
      if (ok)
      {
        MnvVertErrorBand * old = clone->PopVertErrorBand(it->first);
        clone->PushErrorBand(it->first, newEB);
        delete old;
      }
    }
    else
      ok = (it->second->Rebin( ngroup ) != NULL);

    if( ! ok )
    {
      Error("Rebin", Form("Could not rebin MnvH1Ds because Rebin failed for MnvVertErrorBand %s ", it->first.c_str() ) );
      gErrorIgnoreLevel = oldVerbosity;
      return NULL;
    }
  }//done adding Vert errors

  // Call Rebin for all uncorrelated errors
  for( std::map<std::string, TH1D*>::iterator it = fUncorrErrorMap.begin(); it != fUncorrErrorMap.end(); ++it )
  {
    Bool_t ok;
    TH1 * newH = it->second->TH1D::Rebin( ngroup, "tmp", xbins );
    ok = (newH != NULL);
    if( ! ok )
    {
      Error("Rebin", Form("Could not rebin MnvH1Ds because Rebin failed for uncorrelated error %s ", it->first.c_str() ) );
      gErrorIgnoreLevel = oldVerbosity;
      return NULL;
    }

    if (clone)
    {
      dynamic_cast<TH1D*>(newH)->Copy(*clone->GetUncorrError(it->first));
      clone->GetUncorrError(it->first)->SetName(it->second->GetName());
      delete newH;
    }
  }//done adding uncorr errors

  // the systematic error matrices on a rebinned copy
  // will be bunk.  we need to recalculate them.
  if (clone)
    clone->ClearSysErrorMatrices();

  // Do we have special errors in the systemics?
  if ( !( fSysErrorMatrix.empty() ) || !( fRemovedSysErrorMatrix.empty() ) )
  {
    Warning("MnvH1D::Rebin", "Customized error matrices were found (errors that come from neither Vertical nor Lateral Error Bands). They will be cleared.");
  
    if (fStrict){
      assert(0);
    }
    ClearSysErrorMatrices( );
  }

  gErrorIgnoreLevel = oldVerbosity;

  return rval;
}


//=========================================
// Reset
//=========================================
void MnvH1D::Reset(Option_t *option)
{
  //reset the base class
  this->TH1D::Reset(option);

  //reset all vert and lat error bands, but do not remove them
  for( std::map<std::string, MnvLatErrorBand*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
    it->second->Reset(option);
  for( std::map<std::string, MnvVertErrorBand*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
    it->second->Reset(option);
  for( std::map<std::string, TH1D*>::iterator it = fUncorrErrorMap.begin(); it != fUncorrErrorMap.end(); ++it )
    it->second->Reset(option);

  //if there are any error matrices, clear them
  ClearSysErrorMatrices( );
}


//===========================================
// SetBit
//===========================================
void MnvH1D::SetBit(UInt_t f, Bool_t set)
{
  //set the base class bit
  this->TH1D::SetBit(f,set);

  //reset all vert and lat error bands
  for( std::map<std::string, MnvLatErrorBand*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
    it->second->SetBit(f,set);
  for( std::map<std::string, MnvVertErrorBand*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
    it->second->SetBit(f,set);
  for( std::map<std::string, TH1D*>::iterator it = fUncorrErrorMap.begin(); it != fUncorrErrorMap.end(); ++it )
    it->second->SetBit(f,set);

  //maybe we need to clear error matrices.
  //there are cases where you wouldn't, so leave this to the caller
}


void MnvH1D::MnvH1DToCSV(std::string name, std::string directory, double scale, bool fullprecision, bool syserrors, bool percentage, bool binwidth){
  
  //std::cout << "entering 1DToCSV " << name << " scale = " << scale << " percentage = " << percentage << " binwidth = " << binwidth << std::endl;

  std::ofstream *f_values =new std::ofstream();
  std::ofstream *f_err =new std::ofstream();
  std::ofstream *f_staterr =new std::ofstream();
  std::ofstream *f_syserr =new std::ofstream();
  std::ofstream *f_bins =new std::ofstream();
  std::ofstream *f_corr =new std::ofstream();
  std::ofstream *f_cov =new std::ofstream();
  
  
  
  f_values->open((directory+name+"_1d.csv").c_str());
  f_err->open((directory+name+"_errors_1d.csv").c_str());
  f_staterr->open((directory+name+"_staterrors_1d.csv").c_str());
  f_syserr->open((directory+name+"_syserrors_1d.csv").c_str());
  f_bins->open((directory+name+"_bins_1d.csv").c_str());
  f_corr->open((directory+name+"_correlation.csv").c_str());
  f_cov->open((directory+name+"_covariance.csv").c_str());
  
  TH1D stat=GetStatError(); //stat error
  TH1D total=GetCVHistoWithError(); // CV with total error
  TH1D sys=GetTotalError(false); //sys error only
  
  //    *f_bins<<GetXaxis()->GetBinLowEdge(1); //<<std::endl;
  *f_bins<<"Bins  \t" << GetName();
  /*
  *f_values << "Values\t"<< GetName();
  *f_err << "err\t"<< GetName();
  *f_staterr << "staterr\t"<< GetName();
  *f_syserr << "syserr\t"<< GetName();
  
  
  *f_bins << std::endl;
  *f_values << std::endl;
  *f_err << std::endl;
  *f_staterr << std::endl;
  *f_syserr << std::endl;
  */
  *f_bins<<GetXaxis()->GetBinLowEdge(1)<< ","; //<<std::endl;
  if(fullprecision){
    
    
    for (int i=1;i<=GetXaxis()->GetNbins();i++)
      {
      
      if (i>1) {
        *f_values << ",";
        *f_bins << ",";
        *f_err << ",";
        *f_staterr << ",";
        *f_syserr << ",";
      }
      double bincor = 1.0;
      if (binwidth) bincor = (GetXaxis()->GetBinWidth(i));
         
      // Bin width normalize if not enu when we want a total x sec
      *f_bins<<GetXaxis()->GetBinUpEdge(i);//<<std::endl;
      *f_values<<Form("%.17e",total.GetBinContent(i)/bincor*scale);
      *f_err<<Form("%.17e",total.GetBinError(i)/bincor*scale);
      *f_staterr<<Form("%.17e",stat.GetBinContent(i)/bincor*scale);
      *f_syserr<<Form("%.17e",sys.GetBinContent(i)/bincor*scale);
      
      }
  }
  else{
    //    *f_bins<<GetXaxis()->GetBinLowEdge(1); //<<std::endl;
    /**f_bins<<"Bins\t";
     *f_bins<<GetXaxis()->GetBinLowEdge(1)<< "\t"; //<<std::endl;
     *f_values << "Values\t";
     *f_err << "err\t";
     *f_staterr << "staterr\t";
     *f_syserr << "syserr\t";
     */
    
    for (int i=1;i<=GetXaxis()->GetNbins();i++)
      {
      if (i>1) {
        *f_values << ",";
        *f_bins << ",";
        *f_err << ",";
        *f_staterr << ",";
        *f_syserr << ",";
      }
      double bincor = 1.0;
      if (binwidth) bincor = (GetXaxis()->GetBinWidth(i));
         
      *f_bins<<GetXaxis()->GetBinUpEdge(i);//<<std::endl;
      // Bin width normalize if not enu when we want a total x sec
      *f_values<<Form("%.2f",total.GetBinContent(i)/bincor*scale);
      *f_err<<Form("%.2f",total.GetBinError(i)/bincor*scale);
      *f_staterr<<Form("%.2f",stat.GetBinContent(i)/bincor*scale);
      *f_syserr<<Form("%.2f",sys.GetBinContent(i)/bincor*scale);
      
      }
  }
  *f_bins << std::endl;
  *f_values << std::endl;
  *f_err << std::endl;
  *f_staterr << std::endl;
  *f_syserr << std::endl;
  f_values->close();
  f_err->close();
  f_staterr->close();
  f_syserr->close();
  f_bins->close();
  
  //    TMatrixD correlation_matrix= GetTotalCorrelationMatrix();
  TMatrixD correlation_matrix= GetTotalCorrelationMatrix();
  TMatrixD covariance_matrix= GetTotalErrorMatrix();
  correlation_matrix *= (scale*scale); // scale by factor of 10^41
  
  int nbins_x=GetNbinsX();
  
  int totalbins=(nbins_x+2);
  
  
  //*f_corr<< "covariance" <<  GetName() <<std::endl;
  for (int x=0;x<totalbins;x++)
    {
    
    double binwidcorri;
    
    binwidcorri = GetXaxis()->GetBinWidth(x);
        
    
    if (!binwidth)binwidcorri = 1.0;
    
    if (x==0 ||  x==nbins_x+1 ) continue; // Do not print overflow and underflow
    
    for (int this_x=0;this_x<totalbins;this_x++)
      {
      
      
      if  (this_x==0 || this_x==nbins_x+1 ) continue; // Do not print overflow and underflow
      double binwidcorrj;
      
      binwidcorrj = GetXaxis()->GetBinWidth(this_x);
          
      if (!binwidth)binwidcorrj = 1.0;
      if (this_x > 1) *f_corr << ",";
      if (this_x > 1) *f_cov << ",";
      if(!fullprecision){
        *f_cov<<Form("%.2e",covariance_matrix[x][this_x]/binwidcorri/binwidcorrj);
        *f_corr<<Form("%.2e",correlation_matrix[x][this_x]);
      }
      else{
        *f_cov<<Form("%.17e",covariance_matrix[x][this_x]/binwidcorri/binwidcorrj);
        *f_corr<<Form("%.17e",correlation_matrix[x][this_x]);
      }
      // need to include bin widths
      }
    *f_corr<<std::endl;
    *f_cov<<std::endl;
    
  
    }
  f_corr->close();
  f_cov->close();
  delete f_corr;
  delete f_cov;
  //    if (!syserrors){
  //      std::cout << " no systematic errors to consider " << std::endl;
  //      return;
  //    }
  if (!syserrors) return;

  std::ofstream * f_lat = new std::ofstream();
  std::ofstream * f_vert = new std::ofstream();
  f_cov = new std::ofstream(); // reuse the name, sorry
  f_lat->open((directory+name+"_latdump.csv").c_str());
  f_vert->open((directory+name+"_vertdump.csv").c_str());
  f_cov->open((directory+name+"_covdump.csv").c_str());
  std::vector<std::string> vert_errBandNames = GetVertErrorBandNames();
  std::vector<std::string> lat_errBandNames  = GetLatErrorBandNames();
  std::vector<std::string> uncorr_errBandNames  = GetUncorrErrorNames();
  std::vector<std::string> cov_errNames = GetCovMatricesNames();
  *f_vert << " vert " <<  vert_errBandNames.size() <<  std::endl;
  for( std::vector<std::string>::iterator name=vert_errBandNames.begin(); name!=vert_errBandNames.end(); ++name ){
    MnvVertErrorBand* v = GetVertErrorBand(*name);
    unsigned int nunis = v->GetNHists();
    
    for (unsigned int i = 0; i< nunis; i++){
        *f_vert <<GetName() << Form(" %s_%d,",name->c_str(),i);// << std::endl;
      TH1* h = v->GetHist(i);
      for (int j=1;j <=h->GetXaxis()->GetNbins();j++)
        {
        if (j>1) {
          *f_vert << ",";
        }
        double bincor = 1.0;
        if(binwidth) bincor = h->GetXaxis()->GetBinWidth(j);
            if(percentage) bincor = total.GetBinContent(j);
        double frac= (h->GetBinContent(j)/bincor);
        if (!fullprecision){
        *f_vert<<Form("%.2f",frac);
        }
        else{
          *f_vert<<Form("%.17e",frac);
        }
            
        }
      *f_vert << std::endl;
    }
  }
  f_vert->close();
  delete f_vert;
  *f_lat << " Lat " << lat_errBandNames.size() <<  std::endl;
  for( std::vector<std::string>::iterator name=lat_errBandNames.begin(); name!=lat_errBandNames.end(); ++name ){
    MnvLatErrorBand* v = GetLatErrorBand(*name);
    unsigned int nunis = v->GetNHists();
    
    for (unsigned int i = 0; i< nunis; i++){
        *f_lat << GetName() << Form(" %s_%d,",name->c_str(),i) ; //<< std::endl;
      TH1* h = v->GetHist(i);
      for (int j=1;j <=h->GetXaxis()->GetNbins();j++)
        {
        if (j>1) {
          *f_lat << ",";
        }

        double bincor = 1.0;
        if(binwidth) bincor = h->GetXaxis()->GetBinWidth(j);
        if(percentage) bincor = total.GetBinContent(j);
        double frac = (h->GetBinContent(j)/bincor);
        if (!fullprecision){
          *f_lat<<Form("%.2f",frac);
        }
        else{
          *f_lat<<Form("%.17e",frac);
        }
        
        }
      *f_lat << std::endl;
    }
  }
  f_lat->close();
  delete f_lat;
  *f_cov << " covariance " << cov_errNames.size() <<  std::endl;
  for( std::vector<std::string>::iterator name=cov_errNames.begin(); name!=cov_errNames.end(); ++name ){
    TMatrixD  v = GetSysErrorMatrix(*name);
    *f_cov << GetName() << Form(" %s",name->c_str()) << std::endl;
    int nbins_x = GetXaxis()->GetNbins();
    int totalbins=(nbins_x+2);
    
    for (int x=0;x<totalbins;x++)
      {
      double binwidcorri;
        binwidcorri = 1.0;
      if(binwidth) binwidcorri = GetXaxis()->GetBinWidth(x);
      if(percentage) binwidcorri = GetBinContent(x)*scale;
      if (x==0 ||  x==nbins_x+1 ) continue; // Do not print overflow and underflow
      
      for (int this_x=0;this_x<totalbins;this_x++)
        {
        
        if  (this_x==0 || this_x==nbins_x+1 ) continue; // Do not print overflow and underflow
        double binwidcorrj;
            binwidcorrj = 1.0;
        if(binwidth) binwidcorrj = GetXaxis()->GetBinWidth(this_x);
            if(percentage) binwidcorrj = GetBinContent(this_x)*scale;
        if (this_x > 1) *f_cov << ",";
        
        if(!fullprecision){
          *f_cov<<Form("%.2f",v[x][this_x]/binwidcorri/binwidcorrj*scale*scale);
        }
        else{
        *f_cov<<Form("%.17e",v[x][this_x]/binwidcorri/binwidcorrj*scale*scale);
        }
        // need to include bin widths
        }
      *f_cov<<std::endl;
      }
  }
  f_cov->close();
  delete f_cov;
}

#endif
