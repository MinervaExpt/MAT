#ifndef MNV_MnvH2D_cxx
#define MNV_MnvH2D_cxx 1

//#define MNVDBG 1  // Turn this on to get debug covariance matrices
#include "MnvH2D.h"
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;

using namespace PlotUtils;

//==================================================================================
// CONSTRUCTORS
//==================================================================================

//--------------------------------------------------------
// Copy constructors with default normalized bin width
//--------------------------------------------------------
MnvH2D::MnvH2D() : 
TH2D(),
fNormBinWidthX(1.),
fNormBinWidthY(1.)
{ 
  SilentSumw2();
}

MnvH2D::MnvH2D(Double_t normBinWidthX, Double_t normBinWidthY) :
TH2D(),
fNormBinWidthX(normBinWidthX),
fNormBinWidthY(normBinWidthY)
{
  SilentSumw2();
}

MnvH2D::MnvH2D( const TH2D& h2d ) :
TH2D( h2d )
{ 
  fNormBinWidthX = h2d.GetXaxis()->GetBinWidth(1);
  fNormBinWidthY = h2d.GetYaxis()->GetBinWidth(1);
  SilentSumw2();
}

//----------------------------------------------------------------------------------
// Copy constructors with specified normalized bin width for the X and Y projections
//----------------------------------------------------------------------------------
MnvH2D::MnvH2D(const TH2D& h2d, Double_t normBinWidthX, Double_t normBinWidthY) :
TH2D( h2d ),
fNormBinWidthX(normBinWidthX),
fNormBinWidthY(normBinWidthY)
{
  SilentSumw2();
}

//--------------------------------------------------------
// Constructors with default normalization bin width
//--------------------------------------------------------
MnvH2D::MnvH2D( const char* name, const char* title, Int_t nbinsx, const Float_t* xbins, Int_t nbinsy, const Float_t* ybins) :
TH2D( name, title, nbinsx, xbins, nbinsy, ybins)
{ 
  //! Default normlization bin width is the first bin width
  fNormBinWidthX = xbins[1] - xbins[0];
  fNormBinWidthY = ybins[1] - ybins[0];
  SilentSumw2();
}

MnvH2D::MnvH2D( const char* name, const char* title, Int_t nbinsx, const Double_t* xbins, Int_t nbinsy, const Double_t* ybins) :
TH2D( name, title, nbinsx, xbins, nbinsy, ybins)
{ 
  //! Default normlization bin width is the first bin width
  fNormBinWidthX = xbins[1] - xbins[0];
  fNormBinWidthY = ybins[1] - ybins[0];
  SilentSumw2();
}

MnvH2D::MnvH2D( const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup):
TH2D( name, title, nbinsx, xlow, xup, nbinsy, ylow, yup )
{
  //! Default normalization bin width is the constant width of the bins
  fNormBinWidthX = (xup - xlow) / float(nbinsx);
  fNormBinWidthY = (yup - ylow) / float(nbinsy);
  SilentSumw2();
}

//! Deep copy constructor (the real copy)
MnvH2D::MnvH2D( const MnvH2D& h ) :
TH2D( h )
{
  //! Deep copy the variables
  //std::cout << " copy constructor " << std::endl;
  DeepCopy( h );
}

MnvH2D::~MnvH2D()
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
  
}


MnvH2D& MnvH2D::operator=( const MnvH2D& h )
{
  //! If this is me, then no copy is necessary
  if( this == &h )
    return *this;
  
  //call the base class assignment operator
  this->TH2D::operator=(h);
  
  //! Delete and clear all vert and lat error bands
  for( std::map<std::string, MnvVertErrorBand2D*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
    delete it->second;
  fVertErrorBandMap.clear();
  
  for( std::map<std::string, MnvLatErrorBand2D*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
    delete it->second;
  fLatErrorBandMap.clear();
  
  for( std::map<std::string, TMatrixD*>::iterator it = fSysErrorMatrix.begin(); it != fSysErrorMatrix.end(); ++it )
    delete it->second;
  fSysErrorMatrix.clear();
  
  //! Then deep copy the variables
  DeepCopy(h);
  
  return *this;
}

void MnvH2D::DeepCopy( const MnvH2D& h )
{
  //! Set bin norm width
  fNormBinWidthX = h.GetNormBinWidthX();
  fNormBinWidthY = h.GetNormBinWidthY();
  
  //! Copy the vert and lat error bands
  std::vector<std::string> vertNames = h.GetVertErrorBandNames();
  for( std::vector<std::string>::iterator name = vertNames.begin(); name != vertNames.end(); ++name )
    fVertErrorBandMap[*name] = new MnvVertErrorBand2D( *h.GetVertErrorBand(*name) );
  
  std::vector<std::string> latNames = h.GetLatErrorBandNames();
  for( std::vector<std::string>::iterator name = latNames.begin(); name != latNames.end(); ++name )
    fLatErrorBandMap[*name] = new MnvLatErrorBand2D( *h.GetLatErrorBand(*name) );
  
  std::vector<std::string> sysNames = h.GetCovMatricesNames();  // note that I am using the short Covariance list not the big one
  
  for( std::vector<std::string>::iterator name = sysNames.begin(); name != sysNames.end(); ++name ){
#ifdef MNVDBG
    std::cout << " MnvH2D::DeepCopy " << *name << std::endl;
#endif
    fSysErrorMatrix[*name] = new TMatrixD( h.GetSysErrorMatrix(*name) );
  }
  
  
}

MnvH2D* MnvH2D::Clone(const char* name) const
{
  MnvH2D* a = new MnvH2D(*this);
  if (TString(name) != TString("")){
    a->SetName(name);
  }
  
  
  return a;
}

//--------------------------------------------------------
// Constructors with specified normalization bin width
//--------------------------------------------------------
MnvH2D::MnvH2D( const char* name, const char* title, Int_t nbinsx, const Float_t* xbins, Int_t nbinsy, const Float_t* ybins, Double_t normBinWidthX, Double_t normBinWidthY) :
TH2D( name, title, nbinsx, xbins, nbinsy, ybins),
fNormBinWidthX(normBinWidthX),
fNormBinWidthY(normBinWidthY)
{ 
  SilentSumw2();
}

MnvH2D::MnvH2D( const char* name, const char* title, Int_t nbinsx, const Double_t* xbins, Int_t nbinsy, const Double_t* ybins, Double_t normBinWidthX, Double_t normBinWidthY) :
TH2D( name, title, nbinsx, xbins, nbinsy, ybins),
fNormBinWidthX(normBinWidthX),
fNormBinWidthY(normBinWidthY)
{ 
  SilentSumw2();
}

MnvH2D::MnvH2D( const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, Double_t normBinWidthX, Double_t normBinWidthY):
TH2D( name, title, nbinsx, xlow, xup, nbinsy, ylow, yup ),
fNormBinWidthX(normBinWidthX),
fNormBinWidthY(normBinWidthY)
{
  SilentSumw2();
}

//------------------------------------------------------------------------
// Get an MnvH2D which has its bin content and errors normalized to bin width so it looks smooth
//------------------------------------------------------------------------
MnvH2D MnvH2D::GetBinNormalizedCopy( Double_t normBinWidthX /* = fNormBinWidthX */, Double_t normBinWidthY /* = fNormBinWidthY */ ) const
{
  //! If normBinWidthX or normBinWidthXY is negative, use the default bin width for this MnvH2D for the respective axis
  if( normBinWidthX <= 0 )
    normBinWidthX = fNormBinWidthX;
  if( normBinWidthY <= 0 )
    normBinWidthY = fNormBinWidthY;
  
  MnvH2D rval(*this);
  
  if( normBinWidthX > 0 &&  normBinWidthY > 0 )
    rval.Scale( normBinWidthX * normBinWidthY, "width" );
  
  return rval;
}

MnvH2D  MnvH2D::GetAreaNormalizedCopy() const
{
  MnvH2D new_mnvh2d(*this);
  TH2D hmean = TH2D(*this);
  double cv_area=hmean.Integral();
  
  // Scale all the vertical error band histograms to the CV
  std::vector<std::string> vertNames = new_mnvh2d.GetVertErrorBandNames();
  for( unsigned int i = 0; i != vertNames.size(); ++i )
    {
    MnvVertErrorBand2D *errBand = new_mnvh2d.GetVertErrorBand(vertNames[i]);
    errBand->Scale(cv_area/errBand->Integral()); // Scale the CV histogram for the error band
    std::vector<TH2D*> hists = errBand->GetHists();
    for (unsigned int j = 0; j < hists.size(); j++)
      { // Scale the individual universe histograms
        double old_area=hists[j]->Integral();
        hists[j]->Scale(cv_area/old_area);
      }
    }
  // Scale all the lateral error band histograms to the CV
  std::vector<std::string> latNames = new_mnvh2d.GetLatErrorBandNames();
  for( unsigned int i = 0; i != latNames.size(); ++i )
    {
    MnvLatErrorBand2D *errBand = new_mnvh2d.GetLatErrorBand(latNames[i]);
    errBand->Scale(cv_area/errBand->Integral());
    std::vector<TH2D*> hists = errBand->GetHists();
    for (unsigned int j = 0; j < hists.size(); j++)
      { // Scale the individual universe histograms
        double old_area=hists[j]->Integral();
        hists[j]->Scale(cv_area/old_area);
      }
    }
  
  
  return new_mnvh2d;
}

MnvH1D * MnvH2D::Projection(const char* name, bool projectionX, Int_t firstbin, Int_t lastbin, Option_t* option) const
{
  TH1D *cv_px = 0;
  if(projectionX) cv_px=this->TH2D::ProjectionX(name, firstbin, lastbin, option);
  else            cv_px=this->TH2D::ProjectionY(name, firstbin, lastbin, option);
  
  MnvH1D *h_px = new MnvH1D( *cv_px );
  
  // Getting Vertical Names
  std::vector<std::string> vertNames = GetVertErrorBandNames();
  for( unsigned int i = 0; i != vertNames.size(); ++i )
    {
    std::vector<TH1D*> vert_hists;
    const MnvVertErrorBand2D *errBand = GetVertErrorBand(vertNames[i]);
    int nUniverses = errBand->GetNHists();
    for (int j = 0; j != nUniverses; ++j)
      {
      TH2D *h_universe = (TH2D*) dynamic_cast<const TH2D*>( errBand->GetHist(j) );
      //Preserve the cv UserRange for all universes
      if(this->GetXaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetXaxis()->SetRange( this->TH2D::GetXaxis()->GetFirst(), this->TH2D::GetXaxis()->GetLast() );
      if(this->GetYaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetYaxis()->SetRange( this->TH2D::GetYaxis()->GetFirst(), this->TH2D::GetYaxis()->GetLast() );
      TH1D *h_universe_px = 0;
      if(projectionX) h_universe_px=h_universe->ProjectionX( Form("%s_%s_universe%i", name, vertNames[i].c_str(), j),
                                                            firstbin, lastbin, option );
      else            h_universe_px=h_universe->ProjectionY( Form("%s_%s_universe%i", name, vertNames[i].c_str(), j),
                                                            firstbin, lastbin, option );
      vert_hists.push_back(h_universe_px);
      } // end loop over universes
    
    h_px->AddVertErrorBand(vertNames[i], vert_hists);
    
    // Copy the universe weights to the 1D error band
    if ( errBand->GetUnivWgts() ) h_px->GetVertErrorBand( vertNames[i] )->SetUnivWgts( *errBand->GetUnivWgts() );
    
    //cleaning
    for (std::vector<TH1D*>::iterator itHist = vert_hists.begin() ; itHist != vert_hists.end() ; ++ itHist){
      delete *itHist;
    }
    
    } // end loop over vert errors
  
  
  // Getting Lateral Names
  std::vector<std::string> latNames = GetLatErrorBandNames();
  for( unsigned int i = 0; i != latNames.size(); ++i )
    {
    std::vector<TH1D*> lat_hists;
    const MnvLatErrorBand2D *errBand = GetLatErrorBand(latNames[i]);
    int nUniverses = errBand->GetNHists();
    for (int j = 0; j != nUniverses; ++j)
      {
      TH2D *h_universe = (TH2D*) dynamic_cast<const TH2D*>( errBand->GetHist(j) );
      //Preserve the cv UserRange for all universes
      if(this->GetXaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetXaxis()->SetRange( this->TH2D::GetXaxis()->GetFirst(), this->TH2D::GetXaxis()->GetLast() );
      if(this->GetYaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetYaxis()->SetRange( this->TH2D::GetYaxis()->GetFirst(), this->TH2D::GetYaxis()->GetLast() );
      TH1D *h_universe_px = 0;
      if(projectionX) h_universe_px=h_universe->ProjectionX( Form("%s_%s_universe%i", name, latNames[i].c_str(), j),
                                                            firstbin, lastbin, option );
      else            h_universe_px=h_universe->ProjectionY( Form("%s_%s_universe%i", name, latNames[i].c_str(), j),
                                                            firstbin, lastbin, option );
      lat_hists.push_back(h_universe_px);
      } // end loop over universes
    
    h_px->AddLatErrorBand(latNames[i], lat_hists);
    
    // Copy the universe weights to the 1D error band
    if ( errBand->GetUnivWgts() ) h_px->GetLatErrorBand( latNames[i] )->SetUnivWgts( *errBand->GetUnivWgts() );
    
    // cleaning
    for (std::vector<TH1D*>::iterator itHist = lat_hists.begin() ; itHist != lat_hists.end() ; ++ itHist){
      delete *itHist;
    }
    } // end loop over lat error bands
  
  // HMS - implement for covariance matrix
  
  
  std::vector<std::string> sysNames = GetCovMatricesNames();
  
  
  for( unsigned int n = 0; n != sysNames.size(); ++n )
    {
    // the list includes things that are not matrices, only do the matrices
    
    if(! HasErrorMatrix(sysNames[n])) continue;
    
    TMatrixD sys_mat = GetSysErrorMatrix(sysNames[n]);
    
    int nxbins = GetNbinsX()+2;
    int nybins = GetNbinsY()+2;
    TMatrixD newmat;
    if (projectionX){
      newmat.ResizeTo(nxbins,nxbins);
      if (lastbin == -1) lastbin = nybins - 1;  // fix
    }
    else{
      newmat.ResizeTo(nybins,nybins);
      if (lastbin == -1) lastbin = nxbins - 1;
    }
    
    for (int i=0; i < nxbins*nybins; i++){
      int xibin = i%nxbins;
      int yibin = i/nxbins;
      for (int j=0; j < nxbins*nybins; j++){
        int xjbin = j%nxbins;
        int yjbin = j/nxbins;
        if (!projectionX){
          if (xibin < firstbin || xibin > lastbin) continue;
          if (xjbin < firstbin || xjbin > lastbin) continue;
          newmat[yibin][yjbin] += sys_mat[i][j];
        }
        else{
          if (yibin < firstbin || yibin > lastbin) continue;
          if (yjbin < firstbin || yjbin > lastbin) continue;
          
          newmat[xibin][xjbin] += sys_mat[i][j];
        }
      }
    }
    
    h_px->FillSysErrorMatrix(sysNames[n],newmat);
    }
  
  
  return h_px;
}


MnvH1D * MnvH2D::ProjectionX(const char* name /*= "_px"*/, Int_t firstybin /*= 0*/, Int_t lastybin /*= -1*/, Option_t* option /*= ""*/) const
{
  return this->Projection(name, true, firstybin, lastybin, option);
}

MnvH1D * MnvH2D::ProjectionY(const char* name /*= "_px"*/, Int_t firstxbin /*= 0*/, Int_t lastxbin /*= -1*/, Option_t* option /*= ""*/) const
{
  return this->Projection(name, false, firstxbin, lastxbin, option);
}

bool MnvH2D::AddVertErrorBand( const std::string& name, const int nhists /* = -1 */ )
{
  //! Make sure there are no ErrorBands with this name already
  if( HasErrorBand( name ) )
    {
    std::cout << "Warning [MnvH2D::AddVertErrorBand] : There is already an error band with name \"" << name << "\".  Doing nothing." << std::endl;
    return false;
    }
  
  //! Error bands we own have this MnvH2D's name as a prefix
  const std::string errName( std::string(GetName()) + "_" + name );
  
  //! non-positive nhists means you want to use the VertErrorBand's default
  if( nhists > 0 )
    fVertErrorBandMap[name] = new MnvVertErrorBand2D( errName, (TH2D*)this, nhists );
  else
    fVertErrorBandMap[name] = new MnvVertErrorBand2D( errName, (TH2D*)this );
  
  return true;
}

bool MnvH2D::AddVertErrorBand( const std::string& name, const std::vector<TH2D*>& base )
{
  //! Make sure there are no ErrorBands with this name already
  if( HasErrorBand( name ) )
    {
    std::cout << "Warning [MnvH2D::AddVertErrorBand] : There is already an error band with name \"" << name << "\".  Doing nothing." << std::endl;
    return false;
    }
  
  //! Error bands we own have this MnvH1D's name as a prefix
  const std::string errName( std::string(GetName()) + "_" + name );
  
  //!Set the ErrorBand
  fVertErrorBandMap[name] = new MnvVertErrorBand2D( errName, (TH2D*)this, base );
  
  return true;
}

bool MnvH2D::AddVertErrorBandAndFillWithCV( const std::string& name, const unsigned int nhists )
{
  // Make sure there are no ErrorBands with this name already
  if( HasErrorBand( name ) )
    {
    Warning("MnvH2D::AddVertErrorBandAndFillWithCV", Form("There is already an error band with name \"%s\".  Doing nothing.", name.c_str()) );
    return false;
    }
  
  // Make a vector of histos with the CV
  std::vector<TH2D*> histos(nhists, 0);
  for( unsigned int universe=0; universe<nhists; ++universe )
    {
    TH2D* histo = new TH2D( *this );
    histo->SetName( Form( "tmp_universe_%i", universe ) );
    histos[universe] = histo;
    }
  
  // Add the error band and fill it with the vector of histos
  bool ok = this->AddVertErrorBand( name, histos );
  
  // Clean vector of histos
  for( std::vector<TH2D*>::iterator it=histos.begin(); it!=histos.end(); ++it )
    delete *it;
  histos.clear();
  
  return ok;
}

bool MnvH2D::AddLatErrorBand( const std::string& name, const int nhists /* = -1 */ )
{
  //! Make sure there are no ErrorBands with this name already
  if( HasErrorBand( name ) )
    {
    std::cout << "Warning [MnvH2D::AddLatErrorBand] : There is already an error band with name \"" << name << "\".  Doing nothing." << std::endl;
    return false;
    }
  
  //! Error bands we own have this MnvH2D's name as a prefix
  const std::string errName( std::string(GetName()) + "_" + name );
  
  //! non-positive nhists means you want to use the LatErrorBand's default
  if( nhists > 0 )
    fLatErrorBandMap[name] = new MnvLatErrorBand2D( errName, (TH2D*)this, nhists );
  else
    fLatErrorBandMap[name] = new MnvLatErrorBand2D( errName, (TH2D*)this );
  
  return true;
}

bool MnvH2D::AddLatErrorBand( const std::string& name, const std::vector<TH2D*>& base )
{
  //! Make sure there are no ErrorBands with this name already
  if( HasErrorBand( name ) )
    {
    std::cout << "Warning [MnvH2D::AddLatErrorBand] : There is already an error band with name \"" << name << "\".  Doing nothing." << std::endl;
    return false;
    }
  
  //! Error bands we own have this MnvH1D's name as a prefix
  const std::string errName( std::string(GetName()) + "_" + name );
  
  //!Set the ErrorBand
  fLatErrorBandMap[name] = new MnvLatErrorBand2D( errName, (TH2D*)this, base );
  
  return true;
}


bool MnvH2D::AddLatErrorBandAndFillWithCV( const std::string& name, const unsigned int nhists )
{
  // Make sure there are no ErrorBands with this name already
  if( HasErrorBand( name ) )
    {
    Warning("MnvH2D::AddLatErrorBandAndFillWithCV", Form("There is already an error band with name \"%s\".  Doing nothing.", name.c_str()) );
    return false;
    }
  
  // Make a vector of histos with the CV
  std::vector<TH2D*> histos(nhists, 0);
  for( unsigned int universe=0; universe<nhists; ++universe )
    {
    TH2D* histo = new TH2D( *this );
    histo->SetName( Form( "tmp_universe_%i", universe ) );
    histos[universe] = histo;
    }
  
  // Add the error band and fill it with the vector of histos
  bool ok = this->AddLatErrorBand( name, histos );
  
  // Clean vector of histos
  for( std::vector<TH2D*>::iterator it=histos.begin(); it!=histos.end(); ++it )
    delete *it;
  histos.clear();
  
  return ok;
}

//==================================
// Transfer error bands
//==================================
bool MnvH2D::TransferErrorBands( MnvH2D *hist, bool removeFromMe )
{
  bool allOK = true;
  
  std::vector<std::string> names = GetVertErrorBandNames();
  for( std::vector<std::string>::iterator i = names.begin(); i != names.end(); ++i )
    allOK = TransferVertErrorBand( hist, *i, removeFromMe ) && allOK;
  
  names = GetLatErrorBandNames();
  for( std::vector<std::string>::iterator i = names.begin(); i != names.end(); ++i )
    allOK = TransferLatErrorBand( hist, *i, removeFromMe ) && allOK;
  
  names = GetCovMatricesNames();
  for( std::vector<std::string>::iterator i = names.begin(); i != names.end(); ++i )
    allOK = TransferSysErrorMatrix( hist, *i, removeFromMe ) && allOK;
  
  return allOK;
}

bool MnvH2D::TransferVertErrorBand( MnvH2D *hist, const std::string& errName, bool removeFromMe )
{
  MnvVertErrorBand2D *errBand(0);
  if( removeFromMe )
    errBand = PopVertErrorBand( errName );
  else
    {
    MnvVertErrorBand2D *myErrBand = GetVertErrorBand( errName );
    if( 0 != myErrBand )
      errBand = new MnvVertErrorBand2D( *myErrBand );
    }
  
  if( 0 == errBand )
    {
    Error( "MnvH2D::TransferVertErrorBand", Form("Could not find error band %s to transfer.", errName.c_str() ) );
    return false;
    }
  
  //deep scale the error band by ratio CV histos
  TH2D ratio( *( dynamic_cast<TH2D*>(this) ) );
  ratio.Divide( hist, this );
  errBand->MultiplySingle( errBand, &ratio);
  
  bool pushOK = hist->PushErrorBand( errName, errBand );
  
  return pushOK;
}

bool MnvH2D::TransferLatErrorBand( MnvH2D *hist, const std::string& errName, bool removeFromMe )
{
  MnvLatErrorBand2D *errBand(0);
  if( removeFromMe )
    errBand = PopLatErrorBand( errName );
  else
    {
    MnvLatErrorBand2D *myErrBand = GetLatErrorBand( errName );
    if( 0 != myErrBand )
      errBand = new MnvLatErrorBand2D( *myErrBand );
    }
  
  if( 0 == errBand )
    {
    Error( "MnvH2D::TransferLatErrorBand", Form("Could not find error band %s to transfer.", errName.c_str() ) );
    return false;
    }
  
  //deep scale the error band by ratio CV histos
  TH2D ratio( *( dynamic_cast<TH2D*>(this) ) );
  ratio.Divide( hist, this );
  errBand->MultiplySingle( errBand, &ratio);
  
  
  bool pushOK = hist->PushErrorBand( errName, errBand );
  
  return pushOK;
}

bool MnvH2D::TransferSysErrorMatrix( MnvH2D *hist, const std::string& errName, bool removeFromMe )
{
  TMatrixD* errBand;
  if( removeFromMe )
    errBand = new TMatrixD(*PopSysErrorMatrix( errName ));
  else
    {
    TMatrixD myErrBand = TMatrixD(GetSysErrorMatrix( errName ));
    errBand = new TMatrixD( myErrBand );
    }
  
  if( 0 == errBand )
    {
    Error( "MnvH2D::TransferSysErrorMatrix", Form("Could not find error band %s to transfer.", errName.c_str() ) );
    return false;
    }
  
  //deep scale the error band by ratio CV histos
  TH2D ratio( *( dynamic_cast<TH2D*>(this) ) );
  ratio.Divide( hist, this );
  //errBand->MultiplySingle( errBand, &ratio);
  
  for (int i = 0; i < errBand->GetNrows(); i++){
    double hi = ratio.GetBinContent(i);
    for (int j = 0; j < errBand->GetNcols(); j++){
      double hj = ratio.GetBinContent(j);
      errBand[i][j] *= hi*hj;
    }
  }
  std::cout << "MnvH2D::TransferSysErrorMatrix from " << hist->GetName() << " to " << GetName() << " just transferred " << errName << std::endl;
  bool pushOK = hist->PushSysErrorMatrix( errName, errBand );
  
  return pushOK;
}




MnvLatErrorBand2D* MnvH2D::PopLatErrorBand( const std::string& name )
{
  std::map<std::string, MnvLatErrorBand2D*>::iterator i = fLatErrorBandMap.find( name );
  if( i == fLatErrorBandMap.end() )
    {
    Warning( "MnvH2D::PopLatErrorBand", Form("There is no lateral error band with name \"%s\".  Returning NULL.", name.c_str() ));
    return NULL;
    }
  
  //get a pointer to the error band and remove it from the MnvH2D's vector
  MnvLatErrorBand2D* rval = i->second;
  fLatErrorBandMap.erase(i);
  
  return rval;
}



MnvVertErrorBand2D* MnvH2D::PopVertErrorBand( const std::string& name )
{
  std::map<std::string, MnvVertErrorBand2D*>::iterator i = fVertErrorBandMap.find( name );
  if( i == fVertErrorBandMap.end() )
    {
    Warning( "MnvH2D::PopVertErrorBand", Form("There is no vertical error band with name \"%s\".  Returning NULL.", name.c_str() ));
    return NULL;
    }
  
  //get a pointer to the error band and remove it from the MnvH2D's vector
  MnvVertErrorBand2D* rval = i->second;
  fVertErrorBandMap.erase(i);
  
  return rval;
}


TMatrixD* MnvH2D::PopSysErrorMatrix( const std::string& name )
{
  std::map<std::string, TMatrixD*>::iterator i = fSysErrorMatrix.find( name );
  if( i == fSysErrorMatrix.end() )
    {
    Warning( "MnvH2D::PopSysErrorMatrix", Form("There is no systematic error matrix with name \"%s\". Returning NULL.", name.c_str()));
    return NULL;
    }
  
  //get a pointer to the error band and remove it from the MnvH1D's vector
  TMatrixD* rval = i->second;
  fSysErrorMatrix.erase(i);
  
  return rval;
}

bool MnvH2D::PushErrorBand( const std::string& name, MnvVertErrorBand2D* err )
{
  if( HasVertErrorBand(name) )
    {
    Warning( "MnvH2D::PushErrorBand", Form("I already had vert error band %s.  I'm deleting it and adding the new one.", name.c_str() ) );
    delete PopVertErrorBand(name);
    }
  fVertErrorBandMap[name] = err;
  return true;
}

bool MnvH2D::PushErrorBand( const std::string& name, MnvLatErrorBand2D* err )
{ 
  if( HasLatErrorBand(name) )
    {
    Warning( "MnvH2D::PushErrorBand", Form("I already had lat error band %s.  I'm deleting it and adding the new one.", name.c_str() ) );
    delete PopLatErrorBand(name);
    }
  fLatErrorBandMap[name] = err;
  return true;
}

bool MnvH2D::PushSysErrorMatrix( const std::string& name, TMatrixD* err )
{
  if( HasErrorMatrix(name) )
    {
    Warning( "MnvH2D::PushSysErrorMatrix", Form("I already had sys error matrix %s.  I'm deleting it and adding the new one.", name.c_str() ) );
    delete PopSysErrorMatrix(name);
    }
  fSysErrorMatrix[name] = err;
  return true;
}

//! Rename all histograms inside MnvH1D + Error Bands
void MnvH2D::RenameHistosAndErrorBands( const std::string& name )
{
  
  this->SetName( name.c_str() );
  
  std::vector<std::string> vert_errBandNames = this->GetVertErrorBandNames();
  std::vector<std::string> lat_errBandNames  = this->GetLatErrorBandNames();
  for (std::vector<std::string>::iterator itName = vert_errBandNames.begin(); itName != vert_errBandNames.end(); ++itName)
    {
    MnvVertErrorBand2D* tmp_band = this->GetVertErrorBand(*itName);
    std::string band_name = std::string(name + "_" + *itName);
    tmp_band->SetName( band_name.c_str() );
    for (unsigned int i=0; i < tmp_band->GetNHists(); ++i)
      tmp_band->GetHist(i)->SetName( Form("%s_universe%d",band_name.c_str(),i) );
    }
  
  for (std::vector<std::string>::iterator itName = lat_errBandNames.begin(); itName != lat_errBandNames.end(); ++itName)
    {
    MnvLatErrorBand2D* tmp_band = this->GetLatErrorBand(*itName);
    std::string band_name = std::string(name + "_" + *itName);
    tmp_band->SetName( band_name.c_str() );
    for (unsigned int i=0; i < tmp_band->GetNHists(); ++i)
      tmp_band->GetHist(i)->SetName( Form("%s_universe%d",band_name.c_str(),i) );
    }
  
  // no uncorrelated errors for MnvH2D
  
}

//This removes SumW2 from the error band hists.  This is saves about half the space, but
// if you care about the stat error of your error bands, calling this will change that error
void MnvH2D::UnSumw2Universes()
{
  std::vector<std::string> vert_errBandNames = this->GetVertErrorBandNames();

  for (std::vector<std::string>::iterator itName = vert_errBandNames.begin(); itName != vert_errBandNames.end(); ++itName) 
  {
    MnvVertErrorBand2D* band = this->GetVertErrorBand(*itName);
    for (unsigned int i=0; i < band->GetNHists(); ++i){
      band->GetHist(i)->Sumw2(false);
    }
  }

  std::vector<std::string> lat_errBandNames  = this->GetLatErrorBandNames();

  for (std::vector<std::string>::iterator itName = lat_errBandNames.begin(); itName != lat_errBandNames.end(); ++itName) 
  {
    MnvLatErrorBand2D* band = this->GetLatErrorBand(*itName);
    for (unsigned int i=0; i < band->GetNHists(); ++i){
      band->GetHist(i)->Sumw2(false);
    }
  }

}


void MnvH2D::ClearAllErrorBands()
{
  
  ClearSysErrorMatrices( );
  
  // Delete and clear all vert and lat error bands
  for( std::map<std::string, MnvVertErrorBand2D*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
    delete it->second;
  fVertErrorBandMap.clear();
  
  for( std::map<std::string, MnvLatErrorBand2D*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
    delete it->second;
  fLatErrorBandMap.clear();
  
}


bool MnvH2D::FillVertErrorBand( const std::string& name, const double xval, const double yval, const std::vector<double>& weights, const double cvweight /* = 1.0 */, double cvWeightFromMe /*= 1.*/ )
{
  return FillVertErrorBand( name, xval, yval, &(weights[0]), cvweight, cvWeightFromMe );
}

bool MnvH2D::FillVertErrorBand( const std::string& name, const double xval, const double yval, const double * weights, const double cvweight /* = 1.0 */, double cvWeightFromMe /*= 1.*/ )
{
  //! Try to fill a vertical error band
  MnvVertErrorBand2D* vert = GetVertErrorBand( name );
  if( vert )
    return vert->Fill( xval, yval, weights, cvweight, cvWeightFromMe );
  
  std::cout << "Warning [MnvH2D::FillVertErrorBand] : Could not find a vertical error band to fill with name = " << name << std::endl;
  return false;
}

bool MnvH2D::FillVertErrorBand( const std::string& name, const double xval, const double yval, const double weightDown, const double weightUp, const double cvweight  /*= 1.0*/, double cvWeightFromMe /*= 1.*/ )
{
  //! Try to fill a vertical error band
  MnvVertErrorBand2D *vert = GetVertErrorBand( name );
  if( vert )
    return vert->Fill( xval, yval, weightDown, weightUp, cvweight, cvWeightFromMe );
  
  std::cout << "Warning [MnvH2D::FillVertErrorBand] : Could not find a vertical error band to fill with name = " << name << std::endl;
  return false;
}

bool MnvH2D::FillLatErrorBand( const std::string& name, const double xval, const double yval, const std::vector<double>& xshifts, const std::vector<double>& yshifts, const double cvweight  /*= 1.0*/, const bool fillcv /*= true*/, const double* weights /*= 0*/  )
{
  return FillLatErrorBand( name, xval, yval, &(xshifts[0]), &(yshifts[0]), cvweight, fillcv, weights );
}

bool MnvH2D::FillLatErrorBand( const std::string& name, const double xval, const double yval, const double *xshifts, const double *yshifts, const double cvweight  /*= 1.0*/, const bool fillcv /*= true*/, const double* weights /*= 0*/ )
{
  //! Try to fill a lateral error band
  MnvLatErrorBand2D* lat = GetLatErrorBand( name );
  if( lat )
    return lat->Fill( xval, yval, xshifts, yshifts, cvweight, fillcv, weights );
  
  std::cout << "Warning [MnvH2D::FillLatErrorBand] : Could not find a lateral error band to fill with name = " << name << std::endl;
  
  return false;
}

bool MnvH2D::FillLatErrorBand( const std::string& name, const double xval, const double yval, const double xshiftDown, const double xshiftUp, const double yshiftDown, const double yshiftUp, const double cvweight /*= 1.0*/, const bool fillcv /*= true*/ )
{
  //! Try to fill a vertical error band
  MnvLatErrorBand2D *lat = GetLatErrorBand( name );
  if( lat )
    return lat->Fill( xval, yval, xshiftDown, xshiftUp, yshiftDown, yshiftUp, cvweight, fillcv );
  
  std::cout << "Warning [MnvH2D::FillLatErrorBand] : Could not find a lateral error band to fill with name = " << name << std::endl;
  
  return false;
}

bool MnvH2D::FillSysErrorMatrix(const std::string& name, const TMatrixD& matrix){
  
  if (!HasErrorMatrix(name)){
    std::cout << "MnvH2D::FillSysErrorMatrix: creating systematic error matrix in " << GetName() << " " << name << endl;
    fSysErrorMatrix[name] = new TMatrixD(matrix);
    return true;
  }
  std::cout << "MnvH2D::FillSysErrorMatrix: modifying systematic error matrix in " << GetName() << " " << name << endl;
  delete fSysErrorMatrix[name];
  fSysErrorMatrix[name] = new TMatrixD(matrix);
  return true;
}




MnvVertErrorBand2D* MnvH2D::GetVertErrorBand( const std::string& name )
{
  std::map<std::string, MnvVertErrorBand2D*>::iterator i = fVertErrorBandMap.find( name );
  if( i == fVertErrorBandMap.end() )
    {
    std::cout << "Warning [MnvH2D::GetVertErrorBand] : There is no vertical error band with name \"" << name << "\".  Returning NULL." << std::endl;
    return NULL;
    }
  
  return i->second;
}

const MnvVertErrorBand2D* MnvH2D::GetVertErrorBand( const std::string& name ) const
{
  std::map<std::string, MnvVertErrorBand2D*>::const_iterator i = fVertErrorBandMap.find( name );
  if( i == fVertErrorBandMap.end() )
    {
    std::cout << "Warning [MnvH2D::GetVertErrorBand] : There is no vertical error band with name \"" << name << "\".  Returning NULL." << std::endl;
    return NULL;
    }
  
  return i->second;
}

MnvLatErrorBand2D* MnvH2D::GetLatErrorBand( const std::string& name )
{
  std::map<std::string, MnvLatErrorBand2D*>::iterator i = fLatErrorBandMap.find( name );
  if( i == fLatErrorBandMap.end() )
    {
    std::cout << "Warning [MnvH2D::GetLatErrorBand] : There is no lateral error band with name \"" << name << "\".  Returning NULL." << std::endl;
    return NULL;
    }
  
  return i->second;
}

const MnvLatErrorBand2D* MnvH2D::GetLatErrorBand( const std::string& name ) const
{
  std::map<std::string, MnvLatErrorBand2D*>::const_iterator i = fLatErrorBandMap.find( name );
  if( i == fLatErrorBandMap.end() )
    {
    std::cout << "Warning [MnvH2D::GetLatErrorBand] : There is no lateral error band with name \"" << name << "\".  Returning NULL." << std::endl;
    return NULL;
    }
  
  return i->second;
}

bool MnvH2D::HasLatErrorBand( const std::string& name ) const
{
  //! Check the MnvLatErrorBands
  if( fLatErrorBandMap.find( name ) != fLatErrorBandMap.end() )
    return true;
  
  return false;
}

bool MnvH2D::HasVertErrorBand( const std::string& name ) const
{
  //! Check the MnvVertErrorBands
  if( fVertErrorBandMap.find( name ) != fVertErrorBandMap.end() )
    return true;
  
  return false;
}

bool MnvH2D::HasErrorBand( const std::string& name ) const
{
  //! Check the MnvLatErrorBands
  if( HasLatErrorBand( name ) )
    return true;
  
  //! Check the MnvVertErrorBands
  if( HasVertErrorBand( name ) )
    return true;
  
  return false;
}

bool MnvH2D::HasErrorMatrix( const std::string& name ) const
{
  //! Check the fSysErrorMatrix
  if( fSysErrorMatrix.find( name ) != fSysErrorMatrix.end() )
    return true;
  
  return false;
}

std::vector<std::string> MnvH2D::GetErrorBandNames() const
{
  std::vector<std::string> rval;
  for( std::map<std::string, MnvVertErrorBand2D*>::const_iterator i = fVertErrorBandMap.begin(); i != fVertErrorBandMap.end(); ++i )
    rval.push_back( i->first );
  for( std::map<std::string, MnvLatErrorBand2D*>::const_iterator i = fLatErrorBandMap.begin(); i != fLatErrorBandMap.end(); ++i )
    rval.push_back( i->first );
  return rval;
}

std::vector<std::string> MnvH2D::GetVertErrorBandNames() const
{
  std::vector<std::string> rval;
  for( std::map<std::string, MnvVertErrorBand2D*>::const_iterator i = fVertErrorBandMap.begin(); i != fVertErrorBandMap.end(); ++i )
    rval.push_back( i->first );
  return rval;
}

std::vector<std::string> MnvH2D::GetLatErrorBandNames() const
{
  std::vector<std::string> rval;
  for( std::map<std::string, MnvLatErrorBand2D*>::const_iterator i = fLatErrorBandMap.begin(); i != fLatErrorBandMap.end(); ++i )
    rval.push_back( i->first );
  return rval;
}

std::vector<std::string> MnvH2D::GetSysErrorMatricesNames() const
{
  
  std::vector<std::string> rval;
  //Vertical Errors
  for( std::map<std::string, MnvVertErrorBand2D*>::const_iterator i = fVertErrorBandMap.begin(); i != fVertErrorBandMap.end(); ++i )
    rval.push_back( i->first );
  
  //Lateral Errors
  for( std::map<std::string, MnvLatErrorBand2D*>::const_iterator i = fLatErrorBandMap.begin(); i != fLatErrorBandMap.end(); ++i )
    rval.push_back( i->first );
  
  // Special Errors
  for( std::map<std::string, TMatrixD*>::const_iterator i = fSysErrorMatrix.begin(); i != fSysErrorMatrix.end(); ++i )
    if ( !HasEnding(i->first, "_asShape") )
      rval.push_back( i->first );
  
  return rval;
}

std::vector<std::string> MnvH2D::GetCovMatricesNames() const
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
//------------------------------------------------------------------------
// Get a Specific Covariance Matrix from the map
//------------------------------------------------------------------------
TMatrixD MnvH2D::GetSysErrorMatrix(const std::string& name, bool asFrac /*= false*/, bool cov_area_normalize /*= false*/) const
{
  
  if (HasEnding(name,"_asShape") )
    std::cout << "Warning [MnvH2D::GetSysErrorMatrix]: You are calling the error Matrix: " << name <<".\nAssuming the error Band wanted is: " << name.substr(0,name.length()-8) << " with cov_area_normalize = true" << std::endl;
  
  const std::string name_condition = ( (cov_area_normalize) && !(HasEnding(name,"_asShape")) )?  "_asShape" : "";
  const std::string fname = name + name_condition ;
  const std::string errName = HasEnding(fname,"_asShape") ? fname.substr(0,fname.length()-8) : fname;
  
  //! @todo what to do with underflow( bin=0 ) and overflow ( bin=nbins+1 )?
  const int highBinX = GetNbinsX() + 1;
  const int highBinY = GetNbinsY() + 1;
  const int highBin  = GetBin( highBinX, highBinY );
  const int lowBin = 0;
  TMatrixD covmx(highBin+1,highBin+1);

  if ( HasErrorMatrix( fname ) ){
    #ifdef MNVDBG
      std::cout << " MnvH2D::GetSysErrorMatrix getting covariance matrix" << fname << std::endl;
    #endif
    //std::cout << fSysErrorMatrix.find(fname)->second->GetNrows() << std::endl;
    covmx = *(fSysErrorMatrix.find(fname)->second);
  }
  else if( fLatErrorBandMap.find( errName ) != fLatErrorBandMap.end() )
    {
    std::map<std::string, MnvLatErrorBand2D*>::const_iterator it = fLatErrorBandMap.find( errName );
    covmx = it->second->CalcCovMx( ( HasEnding(fname,"_asShape") ) );
    }
  else if( fVertErrorBandMap.find( errName ) != fVertErrorBandMap.end() )
    {
    std::map<std::string, MnvVertErrorBand2D*>::const_iterator it = fVertErrorBandMap.find( errName );
    covmx = it->second->CalcCovMx( ( HasEnding(fname,"_asShape") ) );
    }
  else
    std::cout << "Warning [MnvH2D::GetSysErrorMatrix]: There is no Covariance Matrix with name " << fname << " in " << GetName() << " . Returning an empty Matrix." << std::endl;
  
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

void MnvH2D::ClearSysErrorMatrices()
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



TMatrixD MnvH2D::GetStatErrorMatrix( bool asFrac /* =false */ ) const
{
  const int highBinX = GetNbinsX() + 1;
  const int highBinY = GetNbinsY() + 1;
  const int highBin  = GetBin( highBinX, highBinY );
  const int lowBin = 0;
  TMatrixD covmx(highBin+1,highBin+1);
  //! stat error
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

//--------------------------------------------------------------------------------
// Push a Covariance Matrix to the Systematic/Statistical fSys/StatErrorMatrix map
//--------------------------------------------------------------------------------
bool MnvH2D::PushCovMatrix(
                           const std::string& name,
                           TMatrixD covmx,
                           bool cov_area_normalize /*=false*/)
{
  
  if ( HasEnding(name,"_asShape") )
    {
    std::cout << "Warning [MnvH2D::pushCovMatrix] : You cannot push a Covariance matrix with prefix: \"_asShape\", at the end. Doing Nothing."<< std::endl;
    std::cout << "\nTo Add an \"only shape\" error Matrix, set the last boolean in this constructor to true:\n        PushCovMatrix(\"errorMatrix_name\", TMatrixD errorMatrix, true).\n  "<< std::endl;
    
    return false;
    }
  
  const std::string name_condition = (cov_area_normalize)? "_asShape" : "";
  const std::string fname = name + name_condition ;
  
  const int highBin =  (GetNbinsX() + 2) * (GetNbinsY() + 2);  // account for over- and underflow
  
  if (covmx.GetNrows()!=highBin || covmx.GetNcols()!=highBin )
    {
    std::cout << "Warning [MnvH2D::pushCovMatrix] : The pushed covariance matrix dimensions are incorrect (it should be a " << highBin << "x" << highBin << " matrix, but is actually " << covmx.GetNrows() << "x" << covmx.GetNcols() << "). Doing nothing." <<std::endl;
    if(fStrict){
      assert(0);
    }
    
    return false;
    }
  
  // Make sure there are no Covariance Matrices with this name already
  if( HasErrorMatrix( fname ) )
    {
    std::cout << "Warning [MnvH2D::PushCovMatrix] : There is a Matrix with name \"" << fname << "\" already.  Doing nothing." << std::endl;
    return false;
    }
  else if ( HasErrorBand( name ) )
    {
    std::cout << "Warning [MnvH2D::PushCovMatrix] : " << name << " is already in either a Vertical or Lateral Error" << std::endl;
    return false;
    }
  
  fSysErrorMatrix[fname] = new TMatrixD(covmx);
  
  return true;
}


//------------------------------------------------------------------------
// Modify statistical uncertainty based on studies by DGR related to unfolding and finite MC sizes. See docDB 28992,28899
//------------------------------------------------------------------------

void MnvH2D::ModifyStatisticalUnc(double factor, std::string covmatrixname){

  
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
  for(int i=0;i<this->GetNbinsX()+2;i++){
    for(int j=0;j<this->GetNbinsY()+2;j++){
      this->SetBinError(i,j,this->GetBinError(i,j)*sqrtmod);
    }
  }
}

//------------------------------------------------------------------------
// Get the Total Covariance Matrix
//------------------------------------------------------------------------
TMatrixD MnvH2D::GetTotalErrorMatrix(
                                     bool includeStat /*= true*/,
                                     bool asFrac /*= false*/,
                                     bool cov_area_normalize /*= false*/ ) const
{
  
  //! @todo what to do with underflow( bin=0 ) and overflow ( bin=nbins+1 )?
  const int highBinX = GetNbinsX() + 1;
  const int highBinY = GetNbinsY() + 1;
  const int highBin  = GetBin( highBinX, highBinY );
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

TMatrixD MnvH2D::GetTotalCorrelationMatrix( bool cov_area_normalize /*= false*/ ) const
{
  TMatrixD covmx = GetTotalErrorMatrix(false, false, cov_area_normalize);
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

TH2D MnvH2D::GetTotalError(
                           bool includeStat /* = true */,
                           bool asFrac /* = false */ ,
                           bool cov_area_normalize /*= false */) const
{
  //! Make a copy of this histogram as a TH1D and rename it
  TH2D err( *this );
  err.Reset();
  std::string tmpName( std::string(GetName()) + "_TotalError");
  err.SetName( tmpName.c_str() );
  
  //! @todo what to do with underflow( bin=0 ) and overlow ( bin=nbins+1 )?
  const int highBinX = GetNbinsX() + 1;
  const int highBinY = GetNbinsY() + 1;
  const int highBin  = GetBin( highBinX, highBinY );
  const int lowBin = 0;
  
  //!Get the Total Error Matrix
  TMatrixD errMatrix = GetTotalErrorMatrix(includeStat, asFrac, cov_area_normalize);
  
  for( int iBin = lowBin; iBin <= highBin; ++iBin )
    {
    double derr = errMatrix[iBin][iBin];
    err.SetBinContent( iBin, ( derr > 0 ) ? sqrt(derr): 0. );
    }
  
  return err;
}

//------------------------------------------------------------------------
// Get a histogram of the statistical errors only.
//------------------------------------------------------------------------
TH2D MnvH2D::GetStatError( bool asFrac /* = false */ ) const
{
  //! Make a copy of this histogram as a TH1D and rename it
  TH2D err( *this );
  err.Reset();
  std::string tmpName( std::string(GetName()) + "_StatError");
  err.SetName( tmpName.c_str() );
  
  //! @todo what to do with underflow( bin=0 ) and overlow ( bin=nbins+1 )?
  const int highBinX = GetNbinsX() + 1;
  const int highBinY = GetNbinsY() + 1;
  const int highBin  = GetBin( highBinX, highBinY );
  const int lowBin = 0;
  
  //! stat error
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
  
  return err;
}

TH2D MnvH2D::GetCVHistoWithError( bool includeStat /* = true */ , bool cov_area_normalize /* = false */) const
{
  //! @todo Check the highBin value for here (in MnvH2D, overflow is not set for rval.SetBinError)
  const int highBinX = GetNbinsX() + 1;
  const int highBinY = GetNbinsY() + 1;
  const int highBin  = GetBin( highBinX, highBinY );
  const int lowBin = 0;
  
  //! Get the error band
  TH2D err = GetTotalError( includeStat , false, cov_area_normalize);
  
  //! Create a copy of this histogram and rename it
  TH2D rval( *this );
  std::string tmpName( std::string( GetName() ) + "_CV_WithErr" );
  rval.SetName( tmpName.c_str() );
  
  for( int iBin = lowBin; iBin <= highBin; ++iBin )
    rval.SetBinError( iBin, err.GetBinContent(iBin) );
  
  return rval;
}

TH2D MnvH2D::GetCVHistoWithStatError() const
{
  //! @todo Check the highBin value for here (in MnvH2D, overflow is not set for rval.SetBinError)
  const int highBinX = GetNbinsX() + 1;
  const int highBinY = GetNbinsY() + 1;
  const int highBin  = GetBin( highBinX, highBinY );
  const int lowBin = 0;
  
  //! Get the stat. error band
  TH2D err = GetStatError( false );
  
  //! Create a copy of this histogram and rename it
  TH2D rval( *this );
  std::string tmpName( std::string( GetName() ) + "_CV_WithStatErr" );
  rval.SetName( tmpName.c_str() );
  
  for( int iBin = lowBin; iBin <= highBin; ++iBin )
    rval.SetBinError( iBin, err.GetBinContent(iBin) );
  
  return rval;
}

//======================================================================
// Replacements of ROOT versions of functions
//======================================================================

void MnvH2D::Scale( Double_t c1 /*= 1.*/, Option_t* option /*=""*/, Bool_t allUniv /*=true*/ )
{
  // Scale yourself using TH1D::Scale
  this->TH2D::Scale( c1, option );
  
  if (!allUniv)
    return;
  
  // Scale the vertical error bands
  for( std::map<std::string, MnvVertErrorBand2D*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
    it->second->Scale( c1, option );
  
  // Scale the lateral error bands
  for( std::map<std::string, MnvLatErrorBand2D*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
    it->second->Scale( c1, option );
  
  for( std::map<std::string, TMatrixD*>::const_iterator i = fSysErrorMatrix.begin(); i != fSysErrorMatrix.end(); ++i )
    {
#ifdef MNVDBG
    std::cout << " MnvH2D::Scale Scaling systematic error matrix " << i->first << std::endl;
#endif
    if ( !HasEnding(i->first, "_asShape") ){
      *(i->second) *=c1*c1;
    }
    TMatrixD * m = i->second;
    // implement bin width correction!
    
    if (TString(option).Contains("width")){
      
      int nx = GetNbinsX()+2;
      int ny = GetNbinsY()+2;
      for (int ix = 0; ix < nx; ix++){
        double widix = GetXaxis()->GetBinWidth(ix);
        for (int iy = 0; iy < ny; iy++){
          double widiy = GetYaxis()->GetBinWidth(iy);
          int i =  GetBin(ix,iy);
          for (int jx = 0; jx < nx; jx++){
            double widjx = GetXaxis()->GetBinWidth(jx);
            for (int jy = 0; jy < ny; jy++){
              double widjy = GetYaxis()->GetBinWidth(jy);
              int j = GetBin(jx,jy);
              if (widix*widjx*widiy*widjy != 0.0){
                (*m)[i][j] /=(widix*widiy*widjy*widjx);
                
              }
              else{
                (*m)[i][j] = 0.0;
              } // end else
            } // end jy
          } // end jx
        } // end iy
      } // end ix
    } // end width
#ifdef MNVDBG
    std::cout << " MnvH2D::Scale error matrix scale " <<  (*m)[41][40] << std::endl;
#endif
    } // end loop over syst
}

void MnvH2D::Add( const TH2* h1, const Double_t c1 /*= 1.*/ )
{
  //! Try to cast the input TH2 to a MnvH2D
  const MnvH2D *mnv1 = dynamic_cast<const MnvH2D*>(h1);
  
  if( mnv1 )
    {
    
    //! Add as a TH2D
    this->TH2D::Add( h1, c1 );
    
    //! Call Add for all vertical error bands
    for( std::map<std::string, MnvVertErrorBand2D*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
      {
      const MnvVertErrorBand2D* err1 = mnv1->GetVertErrorBand( it->first );
      if( !err1  )
        {
        Error("Add", Form("Could not add MnvH2Ds because they all don't have the %s MnvVertErrorBand2D", it->first.c_str()) );
        return;
        }
      
      Bool_t ok = it->second->Add( err1, c1 );
      
      if( ! ok )
        {
        Error("Add", Form("Could not add MnvH2Ds because histogram add failed for MnvVertErrorBand2D %s ", it->first.c_str() ) );
        return;
        }
      }//done adding Vert errors
    
    //! Call Add for all lateral error bands
    for( std::map<std::string, MnvLatErrorBand2D*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
      {
      const MnvLatErrorBand2D* err1 = mnv1->GetLatErrorBand( it->first );
      if( !err1  )
        {
        Error("Add", Form("Could not add MnvH2Ds because they all don't have the %s MnvLatErrorBand2D", it->first.c_str()) );
        return;
        }
      
      Bool_t ok = it->second->Add( err1, c1 );
      
      if( ! ok )
        {
        Error("Add", Form("Could not add MnvH2Ds because histogram add failed for MnvLatErrorBand2D %s ", it->first.c_str() ) );
        return;
        }
      }//done adding Lat errors
    
      // HMS add for error matrices.
     
        std::vector<std::string> names = GetCovMatricesNames();
        for( std::vector<std::string>::iterator name=names.begin(); name!=names.end(); ++name )
        {
        if(!mnv1->HasErrorMatrix(*name))
        {
        Error("MnvH2d::Add", Form("Could not add MnvH2Ds because they all don't have the %s Cov Matrix ", (*name).c_str()) );
        return;
        }
          TMatrixD newmatrix = GetSysErrorMatrix(*name);
          TMatrixD addmatrix = mnv1->GetSysErrorMatrix(*name);
          addmatrix*=c1;
          newmatrix += addmatrix;
#ifdef MNVDBG
           std::cout << "MnvH2D::Adding " << mnv1->GetName() << " to " << GetName() << "error matrix from h1 is scaled by c1" << newmatrix[41][40] <<  std::endl;
#endif
          FillSysErrorMatrix(*name,newmatrix);
        
      }
    
    }// end if cast to MnvH2D worked
  else
    {
    Error( "MnvH2D::Add", "Unable to add histogram because it could not be cast to an MnvH2D.  Did nothing." );
    }
  
}

void MnvH2D::Multiply( const MnvH2D* h1, const MnvH2D* h2, Double_t c1 /*= 1*/, Double_t c2 /*= 1*/ )
{
  //! @todo Would love to return a bool here, but we want this Multiply to override TH1's and that is void
  
  //! Call the TH1D Multiply
  this->TH2D::Multiply( (TH2D*)h1, (TH2D*)h2, c1, c2 );
  
  //! Scale the vertical error bands
  for( std::map<std::string, MnvVertErrorBand2D*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
    {
    const MnvVertErrorBand2D* err1 = h1->GetVertErrorBand( it->first );
    const MnvVertErrorBand2D* err2 = h2->GetVertErrorBand( it->first );
    if( !err1 || !err2 )
      {
      Error("Multiply", Form("Could not Multiply  MnvH2Ds because they all don't have the %s MnvVertErrorBand2D", it->first.c_str()) );
      return;
      }
    it->second->Multiply( err1, err2, c1, c2 );
    }
  
  //! Scale the lateral error bands
  for( std::map<std::string, MnvLatErrorBand2D*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
    {
    const MnvLatErrorBand2D* err1 = h1->GetLatErrorBand( it->first );
    const MnvLatErrorBand2D* err2 = h2->GetLatErrorBand( it->first );
    if( !err1 || !err2 )
      {
      Error("Multiply", Form("Could not Multiply MnvH2Ds because they all don't have the %s MnvLatErrorBand2D", it->first.c_str()) );
      return;
      }
    it->second->Multiply( err1, err2, c1, c2 );
    }
  
  //HMS error matrices
  
  for( std::map<std::string, TMatrixD*>::const_iterator i = fSysErrorMatrix.begin(); i != fSysErrorMatrix.end(); ++i ){
    
     
    TMatrixD newmatrix = h1->GetSysErrorMatrix(i->first);
    
    for (int i = 0; i < newmatrix.GetNrows(); i++){
      double hi = h2->GetBinContent(i);
      for (int j = 0; j < newmatrix.GetNcols(); j++){
        double hj = h2->GetBinContent(j);
        newmatrix[i][j] *=  hi*hj*c1*c2*c1*c2;
      }
    }
#ifdef MNVDBG
    std::cout << " MnvH2D::Multiplying " << h1->GetName() << " by " << h2->GetName() << "error matrix from h1 is scaled by h2^2" << newmatrix[41][40] << std::endl;
#endif
    
    FillSysErrorMatrix(i->first,newmatrix);
    
  }
  
  return;
}

void MnvH2D::MultiplySingle( const MnvH2D* h1, const TH2* h2, const Double_t c1 /*= 1*/, const Double_t c2 /*= 1*/ )
{
  // Call the TH1D Multiply
  this->TH2D::Multiply( (TH2D*)h1, h2, c1, c2 );
  
  // Call Multiply for all lateral error bands
  for( std::map<std::string, MnvLatErrorBand2D*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
    {
    const MnvLatErrorBand2D* err1 = h1->GetLatErrorBand( it->first );
    if( !err1  )
      {
      Error("MultiplySingle", Form("Could not multiply MnvH2Ds because they all don't have the %s MnvLatErrorBand", it->first.c_str()) );
      return;
      }
    it->second->MultiplySingle( err1, h2, c1, c2 );
    }
  
  // Call Multiply for all vertical error bands
  for( std::map<std::string, MnvVertErrorBand2D*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
    {
    const MnvVertErrorBand2D* err1 = h1->GetVertErrorBand( it->first );
    if( !err1 )
      {
      Error("MultiplySingle", Form("Could not multiply MnvH2Ds because they all don't have the %s MnvVertErrorBand", it->first.c_str()) );
      return;
      }
    it->second->MultiplySingle( err1, h2, c1, c2 );
    }
  
  for( std::map<std::string, TMatrixD*>::const_iterator i = fSysErrorMatrix.begin(); i != fSysErrorMatrix.end(); ++i ){
   
    TMatrixD newmatrix = h1->GetSysErrorMatrix(i->first);
    
    for (int i = 0; i < newmatrix.GetNrows(); i++){
      double hi = h2->GetBinContent(i);
      for (int j = 0; j < newmatrix.GetNcols(); j++){
        double hj = h2->GetBinContent(j);
        newmatrix[i][j] *=  hi*hj*c1*c2*c1*c2;
      }
    }
#ifdef MNVDBG
     std::cout << " MnvH2D::MultiplySingle " << h1->GetName() << " by " << h2->GetName() << "error matrix from h1 is scaled by h2^2" << newmatrix[41][40] << std::endl;
#endif
    FillSysErrorMatrix(i->first,newmatrix);
    
  }
  
  // Do we have special errors in the systemics - for now leave them alone and use Get and Fill to deal with externally.
  
  
  
  
  
  return;
}

//! this multiplies a TMatrix by h2*h2*c1*c1

/* 
 void MnvH2D::MultiplyTMatrix( TMatrixD* h1, const TH2* h2, const Double_t c1 , const Double_t c2  ){
 int nbinsX = h2->GetNbinsX();
 int nbinsY = h2->GetNbinsY();
 int nRows = h2->GetBin(nbinsX+1,nbinsY+1);
 for (int ix = 0; ix <= nbinsX+1; ix++){
 for int jx = 0; jx <= nbinsY+1; jx++){
 int rowx = h2->GetBin(ix,jx);
 double valx = h2->GetBinContent(ix,jx);
 for (int iy = 0; iy <= nbinsX+1; iy++){
 for (int jy = 0; jy <= nbinsY+1; jy++){
 rowy = h2->Getbin(iy,jy);
 double valy = h2->GetBinContent(iy,jy);
 h1[rowx][rowy]*=c1*c2*valx*valy;
 }
 }
 }
 }
 }
 */



void MnvH2D::Divide( const MnvH2D* h1, const MnvH2D* h2, Double_t c1 /*= 1*/, Double_t c2 /*= 1*/, Option_t* option /*=""*/ )
{
  //! @todo Would love to return a bool here, but we want this Divide to override TH1's and that is void
  
  //! Call the TH2D Divide
  this->TH2D::Divide( (TH2D*)h1, (TH2D*)h2, c1, c2, option );
  
  //! Scale the vertical error bands
  for( std::map<std::string, MnvVertErrorBand2D*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
    {
    const MnvVertErrorBand2D* err1 = h1->GetVertErrorBand( it->first );
    const MnvVertErrorBand2D* err2 = h2->GetVertErrorBand( it->first );
    if( !err1 || !err2 )
      {
      Error("MnvH2D::Divide", Form("Could not divide MnvH2Ds because they all don't have the %s MnvVertErrorBand2D", it->first.c_str()) );
      return;
      }
    //cout << "Dividing error band " << err1->GetName() << " by " << err2->GetName() << endl;
    it->second->Divide( err1, err2, c1, c2, option );
    }
  
  //! Scale the lateral error bands
  for( std::map<std::string, MnvLatErrorBand2D*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
    {
    const MnvLatErrorBand2D* err1 = h1->GetLatErrorBand( it->first );
    const MnvLatErrorBand2D* err2 = h2->GetLatErrorBand( it->first );
    if( !err1 || !err2 )
      {
      Error("MnvH2D::Divide", Form("Could not divide MnvH2Ds because they all don't have the %s MnvLatErrorBand2D", it->first.c_str()) );
      return;
      }
    it->second->Divide( err1, err2, c1, c2, option );
    
    }
  
  for( std::map<std::string, TMatrixD*>::const_iterator i = fSysErrorMatrix.begin(); i != fSysErrorMatrix.end(); ++i ){
    
   
    //std::cout << "i->first " << i->first << std::endl;
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
#ifdef MNVDBG
     std::cout << " MnvH2D::Dividing " << h1->GetName() << " by " << h2->GetName() << " error matrix from h1 is scaled by h2^2" <<  newmatrix[41][40] << std::endl;
#endif
   
    //std::cout << " try to replace the error matrix " << this->GetName() << " " <<  std::endl;
    //this->PopSysErrorMatrix(i->first);
    //std::cout << " what is this " << this->IsA() << std::endl;
    FillSysErrorMatrix(i->first,newmatrix);
    //std::cout << " error matrix is replaced " << std::endl;
    
  }
  //std::cout << " leaving Divide " << h1->GetName() << " " << h2->GetName() << std::endl;
  return;
}

void MnvH2D::DivideSingle( const MnvH2D* h1, const TH2* h2, Double_t c1 /*= 1*/, Double_t c2 /*= 1*/, Option_t* option /*=""*/ )
{
  //! Call the TH2D Divide
  this->TH2D::Divide( (TH2D*)h1, h2, c1, c2, option );
  
  //! Scale the vertical error bands
  for( std::map<std::string, MnvVertErrorBand2D*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
    {
    const MnvVertErrorBand2D* err1 = h1->GetVertErrorBand( it->first );
    if( !err1 )
      {
      Error("MnvH2D::Divide", Form("Could not divide MnvH2Ds because they all don't have the %s MnvVertErrorBand2D", it->first.c_str()) );
      return;
      }
    it->second->DivideSingle( err1, h2, c1, c2, option );
    }
  
  //! Scale the lateral error bands
  for( std::map<std::string, MnvLatErrorBand2D*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
    {
    const MnvLatErrorBand2D* err1 = h1->GetLatErrorBand( it->first );
    if( !err1 )
      {
      Error("MnvH2D::Divide", Form("Could not divide MnvH2Ds because they all don't have the %s MnvLatErrorBand2D", it->first.c_str()) );
      return;
      }
    it->second->DivideSingle( err1, h2, c1, c2, option );
    }
  
  for( std::map<std::string, TMatrixD*>::const_iterator i = fSysErrorMatrix.begin(); i != fSysErrorMatrix.end(); ++i ){
    
   
    
    TMatrixD newmatrix = h1->GetSysErrorMatrix(i->first);
    
    for (int i = 0; i < newmatrix.GetNrows(); i++){
      double hi = h2->GetBinContent(i);
      for (int j = 0; j < newmatrix.GetNcols(); j++){
        double hj = h2->GetBinContent(j);
        if (hj*hj!=0){
          newmatrix[i][j] /=  hi*hj*c2*c2;
          newmatrix[i][j] *= c1*c1;
        }
        else{
          newmatrix[i][j] = 0;
        }
      }
    }
#ifdef MNVDBG
     std::cout << " MnvH2D::Dividing " << h1->GetName() << " by " << h2->GetName() << "error matrix from h1 is scaled by h2^2" << newmatrix[41][40] << std::endl;
#endif
    FillSysErrorMatrix(i->first,newmatrix);
    
  }
  return;
}

//=========================================
// Reset
//=========================================
void MnvH2D::Reset(Option_t *option)
{
  //reset the base class
  this->TH2D::Reset(option);
  
  //reset all vert and lat error bands, but do not remove them
  for( std::map<std::string, MnvLatErrorBand2D*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
    it->second->Reset(option);
  for( std::map<std::string, MnvVertErrorBand2D*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
    it->second->Reset(option);
  
  //if there are any error matrices, clear them
  ClearSysErrorMatrices( );
}

//--------------------------------------------------------
// trivial helper functions
//--------------------------------------------------------
bool MnvH2D::HasEnding (std::string const &fullString, std::string const &ending) const
{
  TString a(fullString);
  return a.EndsWith( ending.c_str() );
}

void MnvH2D::SilentSumw2()
{
  if( 0 == GetSumw2N() )
    Sumw2();
}

void MnvH2D::MnvH2DToCSV(std::string name, std::string directory, double scale, bool fullprecision, bool syserrors, bool percentage,bool binwidth){
    //std::cout << "entering H2DToCsV" << name << " " << (directory+name+".csv").c_str() << std::endl;
    std::ofstream *f_values =new std::ofstream();
    std::ofstream *f_err =new std::ofstream();
    std::ofstream *f_staterr =new std::ofstream();
    std::ofstream *f_syserr =new std::ofstream();
    std::ofstream *f_bins =new std::ofstream();
    std::ofstream *f_corr =new std::ofstream();
    std::ofstream *f_cov =new std::ofstream();
    
    
    
    f_values->open((directory+name+".csv").c_str());
    f_err->open((directory+name+"_errors.csv").c_str());
    f_staterr->open((directory+name+"_staterrors.csv").c_str());
    f_syserr->open((directory+name+"_syserrors.csv").c_str());
    f_bins->open((directory+name+"_bins.csv").c_str());
    f_corr->open((directory+name+"_correlation.csv").c_str());
    f_cov->open((directory+name+"_covariance.csv").c_str());
    
    
    TH2D stat=GetStatError(); //stat error
    TH2D total=GetCVHistoWithError(); // CV with total error
    TH2D err = GetTotalError(true);
    TH2D sys=GetTotalError(false); //sys error only
    //sys.Print("ALL");
  //  *f_bins<<"Bins,"<< GetName()<<std::endl;
    
  /*  *f_values << "Values\t" << GetName();
    *f_err << "err\t" << GetName();
    *f_staterr << "staterr\t"<< GetName();
    *f_syserr << "syserr\t"<< GetName();
    
    *f_values<<std::endl;
    *f_err<<std::endl;
    *f_staterr<<std::endl;
    *f_syserr<<std::endl;
    */
    
    *f_bins<<GetXaxis()->GetBinLowEdge(1)<< ","; //<<std::endl;
    
    for (int x=1;x<=GetXaxis()->GetNbins();x++){
        *f_bins<<GetXaxis()->GetBinUpEdge(x)<< ",";//<<std::endl;
        for (int y=1;y<=GetYaxis()->GetNbins();y++){
            if (y> 1) {
                *f_values<<",";
                *f_err<<",";
                *f_staterr<<",";
                *f_syserr<<",";
            }
            double widcor = 1;
            if (binwidth) widcor=GetXaxis()->GetBinWidth(x)*GetYaxis()->GetBinWidth(y);
            if (!fullprecision){
                
                
                
                *f_values<<Form("%.2f",total.GetBinContent(x,y)/widcor*scale);
                *f_err<<Form("%.2f",err.GetBinContent(x,y)/widcor*scale);
                *f_staterr<<Form("%.2f",stat.GetBinContent(x,y)/widcor*scale);
                *f_syserr<<Form("%.2f",sys.GetBinContent(x,y)/widcor*scale);
                //std::cout << "syscheck" <<  sys.GetBinContent(x,y) << std::endl;
            }
            else{
                
                *f_values<<Form("%.17e",total.GetBinContent(x,y)/widcor*scale);
                *f_err<<Form("%.17e",total.GetBinError(x,y)/widcor*scale);
                *f_staterr<<Form("%.17e",stat.GetBinContent(x,y)/widcor*scale);
                *f_syserr<<Form("%.17e",sys.GetBinContent(x,y)/widcor*scale);
                //std::cout << "syscheck" <<  sys.GetBinContent(x,y) << std::endl;
            }
        }
        *f_values<<std::endl;
        *f_err<<std::endl;
        *f_staterr<<std::endl;
        *f_syserr<<std::endl;
    }
    
    *f_bins<< std::endl;
    for (int y=0;y<=GetYaxis()->GetNbins();y++){
        if (y!=0) *f_bins<< ",";
        *f_bins<<GetYaxis()->GetBinUpEdge(y);
    }
    *f_bins << std::endl;
    
    f_values->close();
    f_err->close();
    f_staterr->close();
    f_syserr->close();
    f_bins->close();
    
    //    TMatrixD correlation_matrix= GetTotalCorrelationMatrix();
    TMatrixD correlation_matrix= GetTotalCorrelationMatrix();
    TMatrixD covariance_matrix= GetTotalErrorMatrix();
    correlation_matrix *= (scale*scale); // scale by factor of 10^41
   // *f_corr << "Covariance" << std::endl;
    
    int nbins_x=GetNbinsX();
    int nbins_y=GetNbinsY();
    int totalbins=(nbins_x+2)*(nbins_y+2);
    
    for (int i=0;i<totalbins;i++)
    {
        int x=i%(nbins_x+2);
        int y=i/(nbins_x+2);
        double binwidcorri;
        
        binwidcorri = GetXaxis()->GetBinWidth(x)*GetYaxis()->GetBinWidth(y);
        if (!binwidth) binwidcorri = 1.0;
        
        
        
        if (x==0 || y==0 || x==nbins_x+1 || y== nbins_y+1) continue; // Do not print overflow and underflow
        
        //*f_corr<< "bin_"<<x<<"_"<<y;
        int first = 0;
        for (int j=0;j<totalbins;j++)
        {
            int this_x=j%(nbins_x+2);
            int this_y=j/(nbins_x+2);
            if  (this_x==0 || this_y==0 || this_x==nbins_x+1 || this_y== nbins_y+1) continue; // Do not print overflow and underflow
            double binwidcorrj;
            if ( first != 0) {
                *f_corr<<",";
                *f_cov<<",";
                
            }
            first ++;
            
            binwidcorrj = GetXaxis()->GetBinWidth(this_x)*GetYaxis()->GetBinWidth(this_y);
            if (binwidcorri*binwidcorrj == 0 ){
            //cout << " i,j " << i << ", " << j << ", " << correlation_matrix[i][j] << ", " << binwidcorri << ", " <<  binwidcorrj << ", " << correlation_matrix[i][j]/binwidcorri/binwidcorrj << endl;
            }
            if (!binwidth) binwidcorrj = 1.0;
            if (!fullprecision){
                *f_cov<<Form("%.2e",covariance_matrix[i][j]/binwidcorri/binwidcorrj);   // need to include bin widths
                *f_corr<<Form("%.2e",correlation_matrix[i][j]);   // do not include bin widths
            }
            else{
                *f_cov<<Form("%.17e",covariance_matrix[i][j]/binwidcorri/binwidcorrj);   // need to include bin widths
                *f_corr<<Form("%.2e",correlation_matrix[i][j]);   // do not include bin widths
            }
        }
        *f_corr<<std::endl;
    *f_cov << std::endl;
    }
    f_corr->close();
    f_cov->close();
    delete f_cov;
    delete f_corr;
    if (!syserrors) return;
  
    
  
  // code below tries to dump a whole lot of individual error bands for diagnosis
  
    std::ofstream * f_vert = new std::ofstream();
    std::ofstream * f_lat = new std::ofstream();
    f_cov = new std::ofstream();
    f_lat->open((directory+name+"_latdump.csv").c_str());
    f_vert->open((directory+name+"_vertdump.csv").c_str());
    f_cov->open((directory+name+"_covdump.csv").c_str());
    
    std::vector<std::string> vert_errBandNames = GetVertErrorBandNames();
    std::vector<std::string> lat_errBandNames  = GetLatErrorBandNames();
    std::vector<std::string> uncorr_errBandNames  = GetUncorrErrorNames();
    
    *f_vert << " vert " <<  vert_errBandNames.size() <<  std::endl;
    for( std::vector<std::string>::iterator name=vert_errBandNames.begin(); name!=vert_errBandNames.end(); ++name ){
        MnvVertErrorBand2D* v = GetVertErrorBand(*name);
        unsigned int nunis = v->GetNHists();
        
        for (unsigned int i = 0; i< nunis; i++){
            *f_vert <<GetName() << Form(" %s_%d,",name->c_str(),i);// << std::endl;
            TH2* h = v->GetHist(i);
            for (int x=1;x <=h->GetXaxis()->GetNbins();x++){
                for (int y = 1; y <= h->GetYaxis()->GetNbins(); y++)
                {
                    if (y>1) {
                        *f_vert << ",";
                    }
                    double bincor = 1.0;
                    if(binwidth) bincor = h->GetXaxis()->GetBinWidth(x)*h->GetYaxis()->GetBinWidth(y);
                    if (percentage) bincor = total.GetBinContent(x,y);
                    double frac= (h->GetBinContent(x,y)/bincor);
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
    }
    f_vert->close();
    delete f_vert;
    *f_lat << " Lat " << lat_errBandNames.size() <<  std::endl;
    for( std::vector<std::string>::iterator name=lat_errBandNames.begin(); name!=lat_errBandNames.end(); ++name ){
        MnvLatErrorBand2D* v = GetLatErrorBand(*name);
        unsigned int nunis = v->GetNHists();
        
        for (unsigned int i = 0; i< nunis; i++){
            *f_lat << GetName() << Form(" %s_%d,",name->c_str(),i); // << std::endl;
            TH2* h = v->GetHist(i);
            for (int x=1;x <=h->GetXaxis()->GetNbins();x++){
                for (int y = 1; y <= h->GetYaxis()->GetNbins(); y++)
                {
                    if (y>1) {
                        *f_lat << ",";
                    }
                    
                    double bincor = 1.0;
                    if(binwidth) bincor = h->GetXaxis()->GetBinWidth(x)*h->GetYaxis()->GetBinWidth(y);
                    if(percentage) bincor = total.GetBinContent(x,y);
                    double frac = (h->GetBinContent(x,y)/bincor);
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
    }
    f_lat->close();
    delete f_lat;
  
    std::vector<std::string> cov_errNames = GetCovMatricesNames();
    *f_cov << " covariance " << cov_errNames.size() <<  std::endl;
    for( std::vector<std::string>::iterator name=cov_errNames.begin(); name!=cov_errNames.end(); ++name ){
        TMatrixD  v = GetSysErrorMatrix(*name);
        *f_cov  << GetName() << Form("%s",name->c_str()) << std::endl;
        //std::cout << "Try a matrix" << v[100][100] << std::endl;
      //v.Print();
        int nbins_x=GetNbinsX();
        int nbins_y=GetNbinsY();
        int totalbins=(nbins_x+2)*(nbins_y+2);
        
        for (int i=0;i<totalbins;i++)
        {
            int x=i%(nbins_x+2);
            int y=i/(nbins_x+2);
            double binwidcorri;
            binwidcorri = 1.0;
            if(binwidth) binwidcorri = GetXaxis()->GetBinWidth(x)*GetYaxis()->GetBinWidth(y);
            //if(percentage) binwidcorri = total.GetBinContent(x,y);
            
            
            
            if (x==0 || y==0 || x==nbins_x+1 || y== nbins_y+1) continue; // Do not print overflow and underflow
            
            //*f_corr<< "bin_"<<x<<"_"<<y;
            int first = 0;
            for (int j=0;j<totalbins;j++)
            {
                int this_x=j%(nbins_x+2);
                int this_y=j/(nbins_x+2);
                if  (this_x==0 || this_y==0 || this_x==nbins_x+1 || this_y== nbins_y+1) continue; // Do not print overflow and underflow
                double binwidcorrj;
                if ( first != 0) {
                    *f_cov<<",";
                    //std::cout << ",\t";
                    
                }
                first ++;
                binwidcorrj = 1.0;
                if(binwidth) binwidcorrj = GetXaxis()->GetBinWidth(this_x)*GetYaxis()->GetBinWidth(this_y);
                //if(percentage) binwidcorrj = total.GetBinContent(this_x,this_y);
            
                if(!fullprecision){
                  if (binwidcorri*binwidcorrj != 0.0){
                    *f_cov<<Form("%.4f",v[i][j]/binwidcorri/binwidcorrj*scale*scale);   // need to include bin widths
                    //std::cout <<Form("%.4f",v[i][j]/binwidcorri/binwidcorrj*scale*scale); // /binwidcorri/binwidcorrj*scale*scale);
                  }
                  else{
                    *f_cov<<Form("%.4f",0.0);
                    //std::cout <<Form("%.4f",0.0);
                  }
                  //std::cout << " test " << i << " " << j << " " << v[i][j] << " " << binwidcorri << " " << binwidcorrj << " " << scale << std::endl;
                }
                else
                {
                 if (binwidcorri*binwidcorrj != 0.0){
                    *f_cov<<Form("%.17e",v[i][j]/binwidcorri/binwidcorrj*scale*scale);
                 }
                 else{
                   *f_cov<<Form("%.2f",0.0);
                 }
            // need to include bin widths
                }
            }
            *f_cov<<std::endl;
            //std::cout << endl;
        }
    }
    f_cov->close();
    delete f_cov;
}

#endif
