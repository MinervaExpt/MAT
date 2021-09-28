#ifndef MNV_MnvH3D_cxx
#define MNV_MnvH3D_cxx 1

#include "PlotUtils/MnvH3D.h"

using namespace PlotUtils;

//==================================================================================
// CONSTRUCTORS
//==================================================================================

//--------------------------------------------------------
// Copy constructors with default normalized bin width
//--------------------------------------------------------
MnvH3D::MnvH3D() : 
  TH3D(),
  fNormBinWidthX(1.),
  fNormBinWidthY(1.),
  fNormBinWidthZ(1.)
{ 
  Sumw2(); 
}

MnvH3D::MnvH3D(Double_t normBinWidthX, Double_t normBinWidthY, Double_t normBinWidthZ) :
    TH3D(),
  fNormBinWidthX(normBinWidthX),
  fNormBinWidthY(normBinWidthY),
  fNormBinWidthZ(normBinWidthZ)
{
  Sumw2();
}

MnvH3D::MnvH3D( const TH3D& h2d ) :
    TH3D( h2d ) 
{ 
  fNormBinWidthX = h2d.GetXaxis()->GetBinWidth(1);
  fNormBinWidthY = h2d.GetYaxis()->GetBinWidth(1);
  fNormBinWidthZ = h2d.GetZaxis()->GetBinWidth(1);
  Sumw2();
}

//----------------------------------------------------------------------------------
// Copy constructors with specified normalized bin width for the X and Y projections
//----------------------------------------------------------------------------------
MnvH3D::MnvH3D(const TH3D& h2d, Double_t normBinWidthX, Double_t normBinWidthY, Double_t normBinWidthZ) :
    TH3D( h2d ),
    fNormBinWidthX(normBinWidthX),
    fNormBinWidthY(normBinWidthY),
    fNormBinWidthZ(normBinWidthZ)
{
  Sumw2();
}

//--------------------------------------------------------
// Constructors with default normalization bin width
//--------------------------------------------------------
MnvH3D::MnvH3D( const char* name, const char* title, Int_t nbinsx, const Float_t* xbins, Int_t nbinsy, const Float_t* ybins, Int_t nbinsz, const Float_t* zbins) :
    TH3D( name, title, nbinsx, xbins, nbinsy, ybins, nbinsz, zbins) 
{
  //! Default normlization bin width is the first bin width
  fNormBinWidthX = xbins[1] - xbins[0];
  fNormBinWidthY = ybins[1] - ybins[0];
  fNormBinWidthZ = zbins[1] - zbins[0];
  Sumw2();
}

MnvH3D::MnvH3D( const char* name, const char* title, Int_t nbinsx, const Double_t* xbins, Int_t nbinsy, const Double_t* ybins, Int_t nbinsz, const Double_t* zbins) :
    TH3D( name, title, nbinsx, xbins, nbinsy, ybins, nbinsz, zbins)
{
  //! Default normlization bin width is the first bin width
  fNormBinWidthX = xbins[1] - xbins[0];
  fNormBinWidthY = ybins[1] - ybins[0];
  fNormBinWidthZ = zbins[1] - zbins[0];
  Sumw2();
}

MnvH3D::MnvH3D( const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, Int_t nbinsz, Double_t zlow, Double_t zup):
    TH3D( name, title, nbinsx, xlow, xup, nbinsy, ylow, yup, nbinsz, zlow, zup )
{ 
  //! Default normalization bin width is the constant width of the bins
  fNormBinWidthX = (xup - xlow) / float(nbinsx);
  fNormBinWidthY = (yup - ylow) / float(nbinsy);
  fNormBinWidthY = (zup - zlow) / float(nbinsz);
  Sumw2();
}

//! Deep copy constructor (the real copy)
MnvH3D::MnvH3D( const MnvH3D& h ) :
    TH3D( h )
{
  //! Deep copy the variables
  DeepCopy( h );
}

MnvH3D& MnvH3D::operator=( const MnvH3D& h )
{
  //! If this is me, then no copy is necessary
  if( this == &h )
    return *this;
    
  //call the base class assignment operator
  this->TH3D::operator=(h);

  //! Delete and clear all vert and lat error bands
  for( std::map<std::string, MnvVertErrorBand3D*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
    delete it->second;
  fVertErrorBandMap.clear();

  for( std::map<std::string, MnvLatErrorBand3D*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
    delete it->second;
  fLatErrorBandMap.clear();

  //! Then deep copy the variables
  DeepCopy(h);

  return *this;
}

void MnvH3D::DeepCopy( const MnvH3D& h )
{
  //! Set bin norm width
  fNormBinWidthX = h.GetNormBinWidthX();
  fNormBinWidthY = h.GetNormBinWidthY();
  fNormBinWidthZ = h.GetNormBinWidthZ();

  //! Copy the vert and lat error bands
  std::vector<std::string> vertNames = h.GetVertErrorBandNames();
  for( std::vector<std::string>::iterator name = vertNames.begin(); name != vertNames.end(); ++name )
    fVertErrorBandMap[*name] = new MnvVertErrorBand3D( *h.GetVertErrorBand(*name) );

  std::vector<std::string> latNames = h.GetLatErrorBandNames();
  for( std::vector<std::string>::iterator name = latNames.begin(); name != latNames.end(); ++name )
    fLatErrorBandMap[*name] = new MnvLatErrorBand3D( *h.GetLatErrorBand(*name) );
}

//--------------------------------------------------------
// Constructors with specified normalization bin width
//--------------------------------------------------------
MnvH3D::MnvH3D( const char* name, const char* title, Int_t nbinsx, const Float_t* xbins, Int_t nbinsy, const Float_t* ybins, Int_t nbinsz, const Float_t* zbins, Double_t normBinWidthX, Double_t normBinWidthY, Double_t normBinWidthZ) :
    TH3D( name, title, nbinsx, xbins, nbinsy, ybins, nbinsz, zbins),
  fNormBinWidthX(normBinWidthX),
  fNormBinWidthY(normBinWidthY),
  fNormBinWidthZ(normBinWidthZ)
{ 
  Sumw2();
}

MnvH3D::MnvH3D( const char* name, const char* title, Int_t nbinsx, const Double_t* xbins, Int_t nbinsy, const Double_t* ybins, Int_t nbinsz, const Double_t* zbins, Double_t normBinWidthX, Double_t normBinWidthY, Double_t normBinWidthZ) :
    TH3D( name, title, nbinsx, xbins, nbinsy, ybins, nbinsz, zbins),
  fNormBinWidthX(normBinWidthX),
  fNormBinWidthY(normBinWidthY),
  fNormBinWidthZ(normBinWidthZ)
{ 
  Sumw2();
}

MnvH3D::MnvH3D( const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, Int_t nbinsz, Double_t zlow, Double_t zup, Double_t normBinWidthX, Double_t normBinWidthY, Double_t normBinWidthZ):
    TH3D( name, title, nbinsx, xlow, xup, nbinsy, ylow, yup, nbinsz, zlow, zup ),
  fNormBinWidthX(normBinWidthX),
  fNormBinWidthY(normBinWidthY),
  fNormBinWidthZ(normBinWidthZ)
{
  Sumw2();
}

//------------------------------------------------------------------------
// Get an MnvH3D which has its bin content and errors normalized to bin width so it looks smooth
//------------------------------------------------------------------------
MnvH3D MnvH3D::GetBinNormalizedCopy( Double_t normBinWidthX /* = fNormBinWidthX */, Double_t normBinWidthY /* = fNormBinWidthY */, Double_t normBinWidthZ /* = fNormBinWidthZ */ ) const
{
  //! If normBinWidthX or normBinWidthXY is negative, use the default bin width for this MnvH3D for the respective axis
  if( normBinWidthX <= 0 )
    normBinWidthX = fNormBinWidthX;
  if( normBinWidthY <= 0 )
    normBinWidthY = fNormBinWidthY;
  if( normBinWidthZ <= 0 )
    normBinWidthZ = fNormBinWidthZ;

  MnvH3D rval(*this);

  if( normBinWidthX > 0 &&  normBinWidthY > 0 && normBinWidthZ > 0 )
    rval.Scale( normBinWidthX * normBinWidthY * normBinWidthZ, "width" );

  return rval;
}

MnvH1D * MnvH3D::ProjectionX(const char* name /*= "_px"*/, Int_t firstybin /*= 0*/, Int_t lastybin /*= -1*/, Int_t firstzbin /*= 0*/, Int_t lastzbin /*= -1*/, Option_t* option /*= ""*/) const
{
  TH1D *cv_px = this->TH3D::ProjectionX(name, firstybin, lastybin, firstzbin, lastzbin, option); 
  MnvH1D *h_px = new MnvH1D( *cv_px );

  //! Getting Vertical Names
  std::vector<std::string> vertNames = GetVertErrorBandNames();
  for( unsigned int i = 0; i != vertNames.size(); ++i )
  {
    std::vector<TH1D*> vert_hists;
    const MnvVertErrorBand3D *errBand = GetVertErrorBand(vertNames[i]);
    int nUniverses = errBand->GetNHists();
    for (int j = 0; j != nUniverses; ++j)
    {
      TH3D *h_universe = (TH3D*) dynamic_cast<const TH3D*>( errBand->GetHist(j) );
      //Preserve the cv UserRange for all universes
      if(this->GetXaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetXaxis()->SetRange( this->TH3D::GetXaxis()->GetFirst(), this->TH3D::GetXaxis()->GetLast() );
      if(this->GetYaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetYaxis()->SetRange( this->TH3D::GetYaxis()->GetFirst(), this->TH3D::GetYaxis()->GetLast() );
      if(this->GetZaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetZaxis()->SetRange( this->TH3D::GetZaxis()->GetFirst(), this->TH3D::GetZaxis()->GetLast() );
      TH1D *h_universe_px = h_universe->ProjectionX( Form("%s_%s_universe%i", name, vertNames[i].c_str(), j), firstybin, lastybin, firstzbin, lastzbin, option );
      vert_hists.push_back(h_universe_px);
    }
    h_px->AddVertErrorBand(vertNames[i], vert_hists);

      //cleaning
    for (std::vector<TH1D*>::iterator itHist = vert_hists.begin() ; itHist != vert_hists.end() ; ++ itHist)
      delete *itHist;
  }

  //! Getting Lateral Names
  std::vector<std::string> latNames = GetLatErrorBandNames();
  for( unsigned int i = 0; i != latNames.size(); ++i )
  {
    std::vector<TH1D*> lat_hists;
    const MnvLatErrorBand3D *errBand = GetLatErrorBand(latNames[i]);
    int nUniverses = errBand->GetNHists();
    for (int j = 0; j != nUniverses; ++j)
    {
      TH3D *h_universe = (TH3D*)dynamic_cast<const TH3D*>( errBand->GetHist(j) );
      //Preserve the cv UserRange for all universes
      if(this->GetXaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetXaxis()->SetRange( this->TH3D::GetXaxis()->GetFirst(), this->TH3D::GetXaxis()->GetLast() );
      if(this->GetYaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetYaxis()->SetRange( this->TH3D::GetYaxis()->GetFirst(), this->TH3D::GetYaxis()->GetLast() );
      if(this->GetZaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetZaxis()->SetRange( this->TH3D::GetZaxis()->GetFirst(), this->TH3D::GetZaxis()->GetLast() );
      TH1D *h_universe_px = h_universe->ProjectionX( Form("%s_%s_universe%i", name, latNames[i].c_str(), j), firstybin, lastybin, firstzbin, lastzbin, option );
      lat_hists.push_back(h_universe_px);
    }
    h_px->AddLatErrorBand(latNames[i], lat_hists);

      //cleaning
    for (std::vector<TH1D*>::iterator itHist = lat_hists.begin() ; itHist != lat_hists.end() ; ++ itHist)
      delete *itHist;
  }

  return h_px;
}

MnvH1D * MnvH3D::ProjectionY(const char* name /*= "_py"*/, Int_t firstxbin /*= 0*/, Int_t lastxbin /*= -1*/, Int_t firstzbin /*= 0*/, Int_t lastzbin /*= -1*/, Option_t* option /*= ""*/) const
{
  TH1D *cv_py = TH3D::ProjectionY(name, firstxbin, lastxbin, firstzbin, lastzbin, option); 
  MnvH1D *h_py = new MnvH1D( *cv_py );
  
  //! Getting Vertical Names
  std::vector<std::string> vertNames = GetVertErrorBandNames();
  for( unsigned int i = 0; i != vertNames.size(); ++i )
  {
    std::vector<TH1D*> vert_hists;
    const MnvVertErrorBand3D *errBand = GetVertErrorBand(vertNames[i]);
    int nUniverses = errBand->GetNHists();
    for (int j = 0; j != nUniverses; ++j)
    {
      TH3D *h_universe = (TH3D*)dynamic_cast<const TH3D*>( errBand->GetHist(j) );
      //Preserve the cv UserRange for all universes
      if(this->GetXaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetXaxis()->SetRange( this->TH3D::GetXaxis()->GetFirst(), this->TH3D::GetXaxis()->GetLast() );
      if(this->GetYaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetYaxis()->SetRange( this->TH3D::GetYaxis()->GetFirst(), this->TH3D::GetYaxis()->GetLast() );
      if(this->GetZaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetZaxis()->SetRange( this->TH3D::GetZaxis()->GetFirst(), this->TH3D::GetZaxis()->GetLast() );
      TH1D *h_universe_py = h_universe->ProjectionY( Form("%s_%s_universe%i", name, vertNames[i].c_str(), j), firstxbin, lastxbin, firstzbin, lastzbin, option );
      vert_hists.push_back(h_universe_py);
    }
    h_py->AddVertErrorBand(vertNames[i], vert_hists);

      //cleaning
    for (std::vector<TH1D*>::iterator itHist = vert_hists.begin() ; itHist != vert_hists.end() ; ++ itHist)
      delete *itHist;
  }
  
  //! Getting Lateral Names
  std::vector<std::string> latNames = GetLatErrorBandNames();
  for( unsigned int i = 0; i != latNames.size(); ++i )
  {
    std::vector<TH1D*> lat_hists;
    const MnvLatErrorBand3D *errBand = GetLatErrorBand(latNames[i]);
    int nUniverses = errBand->GetNHists();
    for (int j = 0; j != nUniverses; ++j)
    {
      TH3D *h_universe = (TH3D*)dynamic_cast<const TH3D*>( errBand->GetHist(j) );
      //Preserve the cv UserRange for all universes
      if(this->GetXaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetXaxis()->SetRange( this->TH3D::GetXaxis()->GetFirst(), this->TH3D::GetXaxis()->GetLast() );
      if(this->GetYaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetYaxis()->SetRange( this->TH3D::GetYaxis()->GetFirst(), this->TH3D::GetYaxis()->GetLast() );
      if(this->GetZaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetZaxis()->SetRange( this->TH3D::GetZaxis()->GetFirst(), this->TH3D::GetZaxis()->GetLast() );
      TH1D *h_universe_py = h_universe->ProjectionY( Form("%s_%s_universe%i", name, latNames[i].c_str(), j), firstxbin, lastxbin, firstzbin, lastzbin, option );
      lat_hists.push_back(h_universe_py);
    }
    h_py->AddLatErrorBand(latNames[i], lat_hists);

      //cleaning
    for (std::vector<TH1D*>::iterator itHist = lat_hists.begin() ; itHist != lat_hists.end() ; ++ itHist)
      delete *itHist;
  }

  return h_py;
}

MnvH1D * MnvH3D::ProjectionZ(const char* name /*= "_pz"*/, Int_t firstxbin /*= 0*/, Int_t lastxbin /*= -1*/, Int_t firstybin /*= 0*/, Int_t lastybin /*= -1*/, Option_t* option /*= ""*/) const
{
  TH1D *cv_pz = TH3D::ProjectionZ(name, firstxbin, lastxbin, firstybin, lastybin, option); 
  MnvH1D *h_pz = new MnvH1D( *cv_pz );
  
  //! Getting Vertical Names
  std::vector<std::string> vertNames = GetVertErrorBandNames();
  for( unsigned int i = 0; i != vertNames.size(); ++i )
  {
    std::vector<TH1D*> vert_hists;
    const MnvVertErrorBand3D *errBand = GetVertErrorBand(vertNames[i]);
    int nUniverses = errBand->GetNHists();
    for (int j = 0; j != nUniverses; ++j)
    {
      TH3D *h_universe = (TH3D*)dynamic_cast<const TH3D*>( errBand->GetHist(j) );
      //Preserve the cv UserRange for all universes
      if(this->GetXaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetXaxis()->SetRange( this->TH3D::GetXaxis()->GetFirst(), this->TH3D::GetXaxis()->GetLast() );
      if(this->GetYaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetYaxis()->SetRange( this->TH3D::GetYaxis()->GetFirst(), this->TH3D::GetYaxis()->GetLast() );
      if(this->GetZaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetZaxis()->SetRange( this->TH3D::GetZaxis()->GetFirst(), this->TH3D::GetZaxis()->GetLast() );
      TH1D *h_universe_pz = h_universe->ProjectionZ( Form("%s_%s_universe%i", name, vertNames[i].c_str(), j), firstxbin, lastxbin, firstybin, lastybin, option );
      vert_hists.push_back(h_universe_pz);
    }
    h_pz->AddVertErrorBand(vertNames[i], vert_hists);

      //cleaning
    for (std::vector<TH1D*>::iterator itHist = vert_hists.begin() ; itHist != vert_hists.end() ; ++ itHist)
      delete *itHist;
  }
  
  //! Getting Lateral Names
  std::vector<std::string> latNames = GetLatErrorBandNames();
  for( unsigned int i = 0; i != latNames.size(); ++i )
  {
    std::vector<TH1D*> lat_hists;
    const MnvLatErrorBand3D *errBand = GetLatErrorBand(latNames[i]);
    int nUniverses = errBand->GetNHists();
    for (int j = 0; j != nUniverses; ++j)
    {
      TH3D *h_universe = (TH3D*)dynamic_cast<const TH3D*>( errBand->GetHist(j) );
      //Preserve the cv UserRange for all universes
      if(this->GetXaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetXaxis()->SetRange( this->TH3D::GetXaxis()->GetFirst(), this->TH3D::GetXaxis()->GetLast() );
      if(this->GetYaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetYaxis()->SetRange( this->TH3D::GetYaxis()->GetFirst(), this->TH3D::GetYaxis()->GetLast() );
      if(this->GetZaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetZaxis()->SetRange( this->TH3D::GetZaxis()->GetFirst(), this->TH3D::GetZaxis()->GetLast() );
      TH1D *h_universe_pz = h_universe->ProjectionZ( Form("%s_%s_universe%i", name, latNames[i].c_str(), j), firstxbin, lastxbin, firstybin, lastybin, option );
      lat_hists.push_back(h_universe_pz);
    }
    h_pz->AddLatErrorBand(latNames[i], lat_hists);

      //cleaning
    for (std::vector<TH1D*>::iterator itHist = lat_hists.begin() ; itHist != lat_hists.end() ; ++ itHist)
      delete *itHist;
  }

  return h_pz;
}

TH1 * MnvH3D::Project3D(Option_t* option /*= "x"*/) const
{

  TH1 *cv_p = TH3D::Project3D(option);
  int dim = cv_p->GetDimension();

  //!Crappy way to fill the universes for MnvH1D/2D
  MnvH1D *h_p1D = NULL;
  MnvH2D *h_p2D = NULL;
  
  if (dim == 1)
  {
    h_p1D = new MnvH1D( *(dynamic_cast<TH1D*>( cv_p )) );
    //! Getting Vertical Names
    std::vector<std::string> vertNames = GetVertErrorBandNames();
    for( unsigned int i = 0; i != vertNames.size(); ++i )
    {
      std::vector<TH1D*> vert_hists;
      const MnvVertErrorBand3D *errBand = GetVertErrorBand(vertNames[i]);
      int nUniverses = errBand->GetNHists();
      for (int j = 0; j != nUniverses; ++j)
      {
        TH3D *h_universe = (TH3D*)dynamic_cast<const TH3D*>( errBand->GetHist(j) );
        //Preserve the cv UserRange for all universes
        if(this->GetXaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetXaxis()->SetRange( this->TH3D::GetXaxis()->GetFirst(), this->TH3D::GetXaxis()->GetLast() );
        if(this->GetYaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetYaxis()->SetRange( this->TH3D::GetYaxis()->GetFirst(), this->TH3D::GetYaxis()->GetLast() );
        if(this->GetZaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetZaxis()->SetRange( this->TH3D::GetZaxis()->GetFirst(), this->TH3D::GetZaxis()->GetLast() );
        TH1D *h_universe_p = dynamic_cast<TH1D*>( h_universe->Project3D( option ) );
        vert_hists.push_back( h_universe_p );
      }
      h_p1D->AddVertErrorBand(vertNames[i], vert_hists);

      //cleaning
      for (std::vector<TH1D*>::iterator itHist = vert_hists.begin() ; itHist != vert_hists.end() ; ++ itHist)
        delete *itHist;
    }
  
  //! Getting Lateral Names
    std::vector<std::string> latNames = GetLatErrorBandNames();
    for( unsigned int i = 0; i != latNames.size(); ++i )
    {
      std::vector<TH1D*> lat_hists;
      const MnvLatErrorBand3D *errBand = GetLatErrorBand(latNames[i]);
      int nUniverses = errBand->GetNHists();
      for (int j = 0; j != nUniverses; ++j)
      {
        TH3D *h_universe = (TH3D*)dynamic_cast<const TH3D*>( errBand->GetHist(j) );
        //Preserve the cv UserRange for all universes
        if(this->GetXaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetXaxis()->SetRange( this->TH3D::GetXaxis()->GetFirst(), this->TH3D::GetXaxis()->GetLast() );
        if(this->GetYaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetYaxis()->SetRange( this->TH3D::GetYaxis()->GetFirst(), this->TH3D::GetYaxis()->GetLast() );
        if(this->GetZaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetZaxis()->SetRange( this->TH3D::GetZaxis()->GetFirst(), this->TH3D::GetZaxis()->GetLast() );
        TH1D *h_universe_p = dynamic_cast<TH1D*>( h_universe->Project3D( option ) );
        lat_hists.push_back( h_universe_p );
      }
      h_p1D->AddLatErrorBand(latNames[i], lat_hists);

      //cleaning
      for (std::vector<TH1D*>::iterator itHist = lat_hists.begin() ; itHist != lat_hists.end() ; ++ itHist)
        delete *itHist;
    }
    return dynamic_cast<TH1*>( h_p1D );
  }//end of dim==1
  else if (dim == 2)
  {
    h_p2D = new MnvH2D( *(dynamic_cast<TH2D*>( cv_p )) );
    //! Getting Vertical Names
    std::vector<std::string> vertNames = GetVertErrorBandNames();
    for( unsigned int i = 0; i != vertNames.size(); ++i )
    {
      std::vector<TH2D*> vert_hists;
      const MnvVertErrorBand3D *errBand = GetVertErrorBand(vertNames[i]);
      int nUniverses = errBand->GetNHists();
      for (int j = 0; j != nUniverses; ++j)
      {
        TH3D *h_universe = (TH3D*)dynamic_cast<const TH3D*>( errBand->GetHist(j) );
      //Preserve the cv UserRange for all universes
        if(this->GetXaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetXaxis()->SetRange( this->TH3D::GetXaxis()->GetFirst(), this->TH3D::GetXaxis()->GetLast() );
        if(this->GetYaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetYaxis()->SetRange( this->TH3D::GetYaxis()->GetFirst(), this->TH3D::GetYaxis()->GetLast() );
        if(this->GetZaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetZaxis()->SetRange( this->TH3D::GetZaxis()->GetFirst(), this->TH3D::GetZaxis()->GetLast() );
        TH2D *h_universe_p = dynamic_cast<TH2D*>( h_universe->Project3D( option ) );
        vert_hists.push_back( h_universe_p );
      }
      h_p2D->AddVertErrorBand(vertNames[i], vert_hists);

      //cleaning
      for (std::vector<TH2D*>::iterator itHist = vert_hists.begin() ; itHist != vert_hists.end() ; ++ itHist)
        delete *itHist;
    }
  
  //! Getting Lateral Names
    std::vector<std::string> latNames = GetLatErrorBandNames();
    for( unsigned int i = 0; i != latNames.size(); ++i )
    {
      std::vector<TH2D*> lat_hists;
      const MnvLatErrorBand3D *errBand = GetLatErrorBand(latNames[i]);
      int nUniverses = errBand->GetNHists();
      for (int j = 0; j != nUniverses; ++j)
      {
        TH3D *h_universe = (TH3D*) dynamic_cast<const TH3D*>( errBand->GetHist(j) );
      //Preserve the cv UserRange for all universes
        if(this->GetXaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetXaxis()->SetRange( this->TH3D::GetXaxis()->GetFirst(), this->TH3D::GetXaxis()->GetLast() );
        if(this->GetYaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetYaxis()->SetRange( this->TH3D::GetYaxis()->GetFirst(), this->TH3D::GetYaxis()->GetLast() );
        if(this->GetZaxis()->TestBit(TAxis::kAxisRange)) h_universe->GetZaxis()->SetRange( this->TH3D::GetZaxis()->GetFirst(), this->TH3D::GetZaxis()->GetLast() );
        TH2D *h_universe_p = dynamic_cast<TH2D*>( h_universe->Project3D( option ) );
        lat_hists.push_back( h_universe_p );
      }
      h_p2D->AddLatErrorBand(latNames[i], lat_hists);

      //cleaning
      for (std::vector<TH2D*>::iterator itHist = lat_hists.begin() ; itHist != lat_hists.end() ; ++ itHist)
        delete *itHist;
    }
    return dynamic_cast<TH1*>( h_p2D );
  }//end of dim==2
  else
  { // just in case
    std::cout<<"[MnvH3D::Project3D]: There was an error getting the projection dimension. Returning NULL pointer"<<std::endl;
    return (TH1*)NULL;
  }

}

bool MnvH3D::AddVertErrorBand( const std::string& name, const int nhists /* = -1 */ )
{
  //! Make sure there are no ErrorBands with this name already
  if( HasErrorBand( name ) )
  {
    std::cout << "Warning [MnvH3D::AddVertErrorBand] : There is already an error band with name \"" << name << "\".  Doing nothing." << std::endl;
    return false;
  }

  //! Error bands we own have this MnvH3D's name as a prefix
  const std::string errName( std::string(GetName()) + "_" + name );

  //! non-positive nhists means you want to use the VertErrorBand's default
  if( nhists > 0 )
    fVertErrorBandMap[name] = new MnvVertErrorBand3D( errName, (TH3D*)this, nhists );
  else
    fVertErrorBandMap[name] = new MnvVertErrorBand3D( errName, (TH3D*)this );

  return true;
}

bool MnvH3D::AddVertErrorBand( const std::string& name, const std::vector<TH3D*>& base )
{
  //! Make sure there are no ErrorBands with this name already
  if( HasErrorBand( name ) )
  {
    std::cout << "Warning [MnvH3D::AddVertErrorBand] : There is already an error band with name \"" << name << "\".  Doing nothing." << std::endl;
    return false;
  }
  
  //! Error bands we own have this MnvH1D's name as a prefix
  const std::string errName( std::string(GetName()) + "_" + name );
  
  //!Set the ErrorBand
  fVertErrorBandMap[name] = new MnvVertErrorBand3D( errName, (TH3D*)this, base );

  return true;
}

bool MnvH3D::AddVertErrorBandAndFillWithCV( const std::string& name, const unsigned int nhists )
{
  // Make sure there are no ErrorBands with this name already
  if( HasErrorBand( name ) )
  {
    Warning("MnvH3D::AddVertErrorBandAndFillWithCV", Form("There is already an error band with name \"%s\".  Doing nothing.", name.c_str()) );
    return false;
  }

  // Make a vector of histos with the CV 
  std::vector<TH3D*> histos(nhists, 0);
  for( unsigned int universe=0; universe<nhists; ++universe )
  {
    TH3D* histo = new TH3D( *this );
    histo->SetName( Form( "tmp_universe_%i", universe ) );
    histos[universe] = histo;
  }

  // Add the error band and fill it with the vector of histos
  bool ok = this->AddVertErrorBand( name, histos );

  // Clean vector of histos
  for( std::vector<TH3D*>::iterator it=histos.begin(); it!=histos.end(); ++it )
    delete *it;
  histos.clear();

  return ok;
}

bool MnvH3D::AddLatErrorBand( const std::string& name, const int nhists /* = -1 */ )
{
  //! Make sure there are no ErrorBands with this name already
  if( HasErrorBand( name ) )
  {
    std::cout << "Warning [MnvH3D::AddLatErrorBand] : There is already an error band with name \"" << name << "\".  Doing nothing." << std::endl;
    return false;
  }

  //! Error bands we own have this MnvH3D's name as a prefix
  const std::string errName( std::string(GetName()) + "_" + name );

  //! non-positive nhists means you want to use the LatErrorBand's default
  if( nhists > 0 )
    fLatErrorBandMap[name] = new MnvLatErrorBand3D( errName, (TH3D*)this, nhists );
  else
    fLatErrorBandMap[name] = new MnvLatErrorBand3D( errName, (TH3D*)this );

  return true;
}

bool MnvH3D::AddLatErrorBand( const std::string& name, const std::vector<TH3D*>& base )
{
  //! Make sure there are no ErrorBands with this name already
  if( HasErrorBand( name ) )
  {
    std::cout << "Warning [MnvH3D::AddLatErrorBand] : There is already an error band with name \"" << name << "\".  Doing nothing." << std::endl;
    return false;
  }
  
  //! Error bands we own have this MnvH1D's name as a prefix
  const std::string errName( std::string(GetName()) + "_" + name );
  
  //!Set the ErrorBand
  fLatErrorBandMap[name] = new MnvLatErrorBand3D( errName, (TH3D*)this, base );

  return true;
}

bool MnvH3D::AddLatErrorBandAndFillWithCV( const std::string& name, const unsigned int nhists )
{
  // Make sure there are no ErrorBands with this name already
  if( HasErrorBand( name ) )
  {
    Warning("MnvH3D::AddLatErrorBandAndFillWithCV", Form("There is already an error band with name \"%s\".  Doing nothing.", name.c_str()) );
    return false;
  }

  // Make a vector of histos with the CV 
  std::vector<TH3D*> histos(nhists, 0);
  for( unsigned int universe=0; universe<nhists; ++universe )
  {
    TH3D* histo = new TH3D( *this );
    histo->SetName( Form( "tmp_universe_%i", universe ) );
    histos[universe] = histo;
  }

  // Add the error band and fill it with the vector of histos
  bool ok = this->AddLatErrorBand( name, histos );

  // Clean vector of histos
  for( std::vector<TH3D*>::iterator it=histos.begin(); it!=histos.end(); ++it ) 
    delete *it;
  histos.clear();

  return ok; 
}

bool MnvH3D::AddMissingErrorBandsAndFillWithCV( const MnvH3D& ref )
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
  }

  return true;

}

void MnvH3D::ClearAllErrorBands()
{

  ClearSysErrorMatrices( );

  // Delete and clear all vert and lat error bands
  for( std::map<std::string, MnvVertErrorBand3D*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
    delete it->second;
  fVertErrorBandMap.clear();

  for( std::map<std::string, MnvLatErrorBand3D*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
    delete it->second;
  fLatErrorBandMap.clear();

}

bool MnvH3D::FillVertErrorBand( const std::string& name, const double xval, const double yval, const double zval, const std::vector<double>& weights, const double cvweight /* = 1.0 */, double cvWeightFromMe /*= 1.*/ )
{
  return FillVertErrorBand( name, xval, yval, zval, &(weights[0]), cvweight, cvWeightFromMe );
}

bool MnvH3D::FillVertErrorBand( const std::string& name, const double xval, const double yval, const double zval, const double * weights, const double cvweight /* = 1.0 */, double cvWeightFromMe /*= 1.*/ )
{
  //! Try to fill a vertical error band
  MnvVertErrorBand3D* vert = GetVertErrorBand( name );
  if( vert )
    return vert->Fill( xval, yval, zval, weights, cvweight, cvWeightFromMe );

  std::cout << "Warning [MnvH3D::FillVertErrorBand] : Could not find a vertical error band to fill with name = " << name << std::endl;
  return false;
}

bool MnvH3D::FillVertErrorBand( const std::string& name, const double xval, const double yval, const double zval, const double weightDown, const double weightUp, const double cvweight  /*= 1.0*/, double cvWeightFromMe /*= 1.*/ )
{
  //! Try to fill a vertical error band
  MnvVertErrorBand3D *vert = GetVertErrorBand( name );
  if( vert )
    return vert->Fill( xval, yval, zval, weightDown, weightUp, cvweight, cvWeightFromMe );

  std::cout << "Warning [MnvH3D::FillVertErrorBand] : Could not find a vertical error band to fill with name = " << name << std::endl;
  return false;
}

bool MnvH3D::FillLatErrorBand( const std::string& name, const double xval, const double yval, const double zval, const std::vector<double>& xshifts, const std::vector<double>& yshifts, const std::vector<double>& zshifts, const double cvweight  /*= 1.0*/, const bool fillcv /*= true*/, const double* weights /*= NULL*/  )
{
  return FillLatErrorBand( name, xval, yval, zval, &(xshifts[0]), &(yshifts[0]), &(zshifts[0]), cvweight, fillcv, weights );
}

bool MnvH3D::FillLatErrorBand( const std::string& name, const double xval, const double yval, const double zval, const double *xshifts, const double *yshifts, const double *zshifts, const double cvweight  /*= 1.0*/, const bool fillcv /*= true*/, const double* weights /*= NULL*/ )
{
  //! Try to fill a lateral error band
  MnvLatErrorBand3D* lat = GetLatErrorBand( name );
  if( lat )
    return lat->Fill( xval, yval, zval, xshifts, yshifts, zshifts, cvweight, fillcv, weights );

  std::cout << "Warning [MnvH3D::FillLatErrorBand] : Could not find a lateral error band to fill with name = " << name << std::endl;

  return false;
}

bool MnvH3D::FillLatErrorBand( const std::string& name, const double xval, const double yval, const double zval, const double xshiftDown, const double xshiftUp, const double yshiftDown, const double yshiftUp, const double zshiftDown, const double zshiftUp, const double cvweight /*= 1.0*/, const bool fillcv /*= true*/ )
{
  //! Try to fill a vertical error band
  MnvLatErrorBand3D *lat = GetLatErrorBand( name );
  if( lat )
    return lat->Fill( xval, yval, zval, xshiftDown, xshiftUp, yshiftDown, yshiftUp, zshiftDown, zshiftUp, cvweight, fillcv );

  std::cout << "Warning [MnvH3D::FillLatErrorBand] : Could not find a lateral error band to fill with name = " << name << std::endl;

  return false;
}


MnvVertErrorBand3D* MnvH3D::GetVertErrorBand( const std::string& name )
{
  std::map<std::string, MnvVertErrorBand3D*>::iterator i = fVertErrorBandMap.find( name );
  if( i == fVertErrorBandMap.end() )
  {
    std::cout << "Warning [MnvH3D::GetVertErrorBand] : There is no vertical error band with name \"" << name << "\".  Returning NULL." << std::endl;
    return NULL;
  }

  return i->second;
}

const MnvVertErrorBand3D* MnvH3D::GetVertErrorBand( const std::string& name ) const
{
  std::map<std::string, MnvVertErrorBand3D*>::const_iterator i = fVertErrorBandMap.find( name );
  if( i == fVertErrorBandMap.end() )
  {
    std::cout << "Warning [MnvH3D::GetVertErrorBand] : There is no vertical error band with name \"" << name << "\".  Returning NULL." << std::endl;
    return NULL;
  }

  return i->second;
}

MnvLatErrorBand3D* MnvH3D::GetLatErrorBand( const std::string& name )
{
  std::map<std::string, MnvLatErrorBand3D*>::iterator i = fLatErrorBandMap.find( name );
  if( i == fLatErrorBandMap.end() )
  {
    std::cout << "Warning [MnvH3D::GetLatErrorBand] : There is no lateral error band with name \"" << name << "\".  Returning NULL." << std::endl;
    return NULL;
  }

  return i->second;
}

const MnvLatErrorBand3D* MnvH3D::GetLatErrorBand( const std::string& name ) const
{
  std::map<std::string, MnvLatErrorBand3D*>::const_iterator i = fLatErrorBandMap.find( name );
  if( i == fLatErrorBandMap.end() )
  {
    std::cout << "Warning [MnvH3D::GetLatErrorBand] : There is no lateral error band with name \"" << name << "\".  Returning NULL." << std::endl;
    return NULL;
  }

  return i->second;
}

bool MnvH3D::HasLatErrorBand( const std::string& name ) const
{
  //! Check the MnvLatErrorBands
  if( fLatErrorBandMap.find( name ) != fLatErrorBandMap.end() )
    return true;

  return false;
}

bool MnvH3D::HasVertErrorBand( const std::string& name ) const
{
  //! Check the MnvVertErrorBands
  if( fVertErrorBandMap.find( name ) != fVertErrorBandMap.end() )
    return true;

  return false;
}

bool MnvH3D::HasErrorBand( const std::string& name ) const
{
  //! Check the MnvLatErrorBands
  if( HasLatErrorBand( name ) )
    return true;

  //! Check the MnvVertErrorBands
  if( HasVertErrorBand( name ) )
    return true;

  return false;
}

bool MnvH3D::HasErrorMatrix( const std::string& name ) const
{
  //! Check the fSysErrorMatrix
  if( fSysErrorMatrix.find( name ) != fSysErrorMatrix.end() )
    return true;

  return false;
}

std::vector<std::string> MnvH3D::GetVertErrorBandNames() const
{
  std::vector<std::string> rval;
  for( std::map<std::string, MnvVertErrorBand3D*>::const_iterator i = fVertErrorBandMap.begin(); i != fVertErrorBandMap.end(); ++i )
    rval.push_back( i->first );
  return rval;
}

std::vector<std::string> MnvH3D::GetLatErrorBandNames() const
{
  std::vector<std::string> rval;
  for( std::map<std::string, MnvLatErrorBand3D*>::const_iterator i = fLatErrorBandMap.begin(); i != fLatErrorBandMap.end(); ++i )
    rval.push_back( i->first );
  return rval;
}

std::vector<std::string> MnvH3D::GetSysErrorMatricesNames() const
{

  std::vector<std::string> rval;
  //Vertical Errors
  for( std::map<std::string, MnvVertErrorBand3D*>::const_iterator i = fVertErrorBandMap.begin(); i != fVertErrorBandMap.end(); ++i )
    rval.push_back( i->first );

  //Lateral Errors
  for( std::map<std::string, MnvLatErrorBand3D*>::const_iterator i = fLatErrorBandMap.begin(); i != fLatErrorBandMap.end(); ++i )
    rval.push_back( i->first );
  /*
  //Special Errors
  for( std::map<std::string, TMatrixD*>::const_iterator i = fSysErrorMatrix.begin(); i != fSysErrorMatrix.end(); ++i )
    if ( !HasEnding(i->first, "_asShape") )
      rval.push_back( i->first );
  */
  return rval;
}

//------------------------------------------------------------------------
// Get a Specific Covariance Matrix from the map
//------------------------------------------------------------------------
TMatrixD MnvH3D::GetSysErrorMatrix(const std::string& name, bool asFrac /*= false*/, bool cov_area_normalize /*= false*/) const
{
  if (HasEnding(name,"_asShape") )
    std::cout << "Warning [MnvH3D::GetSysErrorMatrix]: You are calling the error Matrix: " << name <<".\nAssuming the error Band wanted is: " << name.substr(0,name.length()-8) << " with cov_area_normalize = true" << std::endl;

  const std::string name_condition = ( (cov_area_normalize) && !(HasEnding(name,"_asShape")) )?  "_asShape" : "";
  const std::string fname = name + name_condition ; 
  const std::string errName = HasEnding(fname,"_asShape") ? fname.substr(0,fname.length()-8) : fname;

  //! @todo what to do with underflow( bin=0 ) and overflow ( bin=nbins+1 )?
  const int highBinX = GetNbinsX() + 1;
  const int highBinY = GetNbinsY() + 1;
  const int highBinZ = GetNbinsZ() + 1;
  const int highBin  = GetBin( highBinX, highBinY, highBinZ );
  const int lowBin = 0;
  TMatrixD covmx(highBin+1,highBin+1);

  if ( HasErrorMatrix( fname ) )
    covmx = *(fSysErrorMatrix.find(fname)->second);

  else if( fLatErrorBandMap.find( errName ) != fLatErrorBandMap.end() )
  {
    std::map<std::string, MnvLatErrorBand3D*>::const_iterator it = fLatErrorBandMap.find( errName );
    covmx = it->second->CalcCovMx( ( HasEnding(fname,"_asShape") ) );
  }
  else if( fVertErrorBandMap.find( errName ) != fVertErrorBandMap.end() )
  {
    std::map<std::string, MnvVertErrorBand3D*>::const_iterator it = fVertErrorBandMap.find( errName );
    covmx = it->second->CalcCovMx( ( HasEnding(fname,"_asShape") ) );
  }
  else 
    std::cout << "Warning [MnvH3D::GetSysErrorMatrix]: There is no Covariance Matrix with name " << fname << ".Returning and empty Matrix." << std::endl;

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

void MnvH3D::ClearSysErrorMatrices()
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

TMatrixD MnvH3D::GetStatErrorMatrix( bool asFrac /* =false */ ) const
{
  const int highBinX = GetNbinsX() + 1;
  const int highBinY = GetNbinsY() + 1;
  const int highBinZ = GetNbinsZ() + 1;
  const int highBin  = GetBin( highBinX, highBinY, highBinZ );
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

//------------------------------------------------------------------------
// Get the Total Covariance Matrix
//------------------------------------------------------------------------
TMatrixD MnvH3D::GetTotalErrorMatrix(
    bool includeStat /*= true*/, 
    bool asFrac /*= false*/, 
    bool cov_area_normalize /*= false*/ ) const
{

  //! @todo what to do with underflow( bin=0 ) and overflow ( bin=nbins+1 )?
  const int highBinX = GetNbinsX() + 1;
  const int highBinY = GetNbinsY() + 1;
  const int highBinZ = GetNbinsZ() + 1;
  const int highBin  = GetBin( highBinX, highBinY, highBinZ );
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

TMatrixD MnvH3D::GetTotalCorrelationMatrix( bool cov_area_normalize /*= false*/ ) const
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

TH3D MnvH3D::GetTotalError(
    bool includeStat /* = true */, 
    bool asFrac /* = false */ ,
    bool cov_area_normalize /*= false */) const
{
  //! Make a copy of this histogram as a TH1D and rename it
  TH3D err( *this );
  err.Reset();
  std::string tmpName( std::string(GetName()) + "_TotalError");
  err.SetName( tmpName.c_str() );

  //! @todo what to do with underflow( bin=0 ) and overlow ( bin=nbins+1 )?
  const int highBinX = GetNbinsX() + 1;
  const int highBinY = GetNbinsY() + 1;
  const int highBinZ = GetNbinsZ() + 1;
  const int highBin  = GetBin( highBinX, highBinY, highBinZ );
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
TH3D MnvH3D::GetStatError( bool asFrac /* = false */ ) const
{
  //! Make a copy of this histogram as a TH1D and rename it
  TH3D err( *this );
  err.Reset();
  std::string tmpName( std::string(GetName()) + "_StatError");
  err.SetName( tmpName.c_str() );

  //! @todo what to do with underflow( bin=0 ) and overlow ( bin=nbins+1 )?
  const int highBinX = GetNbinsX() + 1;
  const int highBinY = GetNbinsY() + 1;
  const int highBinZ = GetNbinsZ() + 1;
  const int highBin  = GetBin( highBinX, highBinY, highBinZ );
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

TH3D MnvH3D::GetCVHistoWithError( bool includeStat /* = true */ , bool cov_area_normalize /* = false */) const
{
  //! @todo Check the highBin value for here (in MnvH1D, overflow is not set for rval.SetBinError)
  const int highBinX = GetNbinsX() + 1;
  const int highBinY = GetNbinsY() + 1;
  const int highBinZ = GetNbinsZ() + 1;
  const int highBin  = GetBin( highBinX, highBinY, highBinZ );
  const int lowBin = 0;
  
  //! Get the error band
  TH3D err = GetTotalError( includeStat , false, cov_area_normalize);

  //! Create a copy of this histogram and rename it
  TH3D rval( *this );
  std::string tmpName( std::string( GetName() ) + "_CV_WithErr" );
  rval.SetName( tmpName.c_str() );

  for( int iBin = lowBin; iBin <= highBin; ++iBin )
    rval.SetBinError( iBin, err.GetBinContent(iBin) );

  return rval;
}

TH3D MnvH3D::GetCVHistoWithStatError() const
{
  //! @todo Check the highBin value for here (in MnvH1D, overflow is not set for rval.SetBinError)
  const int highBinX = GetNbinsX() + 1;
  const int highBinY = GetNbinsY() + 1;
  const int highBinZ = GetNbinsZ() + 1;
  const int highBin  = GetBin( highBinX, highBinY, highBinZ );
  const int lowBin = 0;
  
  //! Get the stat. error band
  TH3D err = GetStatError( false );

  //! Create a copy of this histogram and rename it
  TH3D rval( *this );
  std::string tmpName( std::string( GetName() ) + "_CV_WithStatErr" );
  rval.SetName( tmpName.c_str() );

  for( int iBin = lowBin; iBin <= highBin; ++iBin )
    rval.SetBinError( iBin, err.GetBinContent(iBin) );

  return rval;
}

//======================================================================
// Replacements of ROOT versions of functions
//======================================================================

void MnvH3D::Scale( Double_t c1 /*= 1.*/, Option_t* option /*=""*/, Bool_t allUniv /*=true*/)
{
  // Scale yourself using TH1D::Scale
  this->TH3D::Scale( c1, option );
  
  if (!allUniv)
    return;

  // Scale the vertical error bands
  for( std::map<std::string, MnvVertErrorBand3D*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
    it->second->Scale( c1, option );

  // Scale the lateral error bands
  for( std::map<std::string, MnvLatErrorBand3D*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
    it->second->Scale( c1, option );
}

void MnvH3D::Add( const TH3* h1, const Double_t c1 /*= 1.*/ )
{
  //! Try to cast the input TH3 to a MnvH3D
  const MnvH3D *mnv1 = dynamic_cast<const MnvH3D*>(h1);

  if( mnv1 )
  {

    //! Add as a TH3D
    this->TH3D::Add( h1, c1 );

    //! Call Add for all vertical error bands
    for( std::map<std::string, MnvVertErrorBand3D*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
    {
      const MnvVertErrorBand3D* err1 = mnv1->GetVertErrorBand( it->first );
      if( !err1  )
      {
        Error("Add", Form("Could not add MnvH3Ds because they all don't have the %s MnvVertErrorBand3D", it->first.c_str()) );
        return;
      }

      Bool_t ok = it->second->Add( err1, c1 );

      if( ! ok )
      {
        Error("Add", Form("Could not add MnvH3Ds because histogram add failed for MnvVertErrorBand3D %s ", it->first.c_str() ) );
        return;
      }
    }//done adding Vert errors

    //! Call Add for all lateral error bands
    for( std::map<std::string, MnvLatErrorBand3D*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
    {
      const MnvLatErrorBand3D* err1 = mnv1->GetLatErrorBand( it->first );
      if( !err1  )
      {
        Error("Add", Form("Could not add MnvH3Ds because they all don't have the %s MnvLatErrorBand3D", it->first.c_str()) );
        return;
      }

      Bool_t ok = it->second->Add( err1, c1 );

      if( ! ok )
      {
        Error("Add", Form("Could not add MnvH3Ds because histogram add failed for MnvLatErrorBand3D %s ", it->first.c_str() ) );
        return;
      }
    }//done adding Lat errors

  }// end if cast to MnvH3D worked
  else
  {
    Error( "MnvH3D::Add", "Unable to add histogram because it could not be cast to an MnvH3D.  Did nothing." );
  }

}

void MnvH3D::Multiply( const MnvH3D* h1, const MnvH3D* h2, Double_t c1 /*= 1*/, Double_t c2 /*= 1*/ )
{
  //! @todo Would love to return a bool here, but we want this Multiply to override TH1's and that is void

  //! Call the TH1D Multiply 
  this->TH3D::Multiply( (TH3D*)h1, (TH3D*)h2, c1, c2 );

  //! Scale the vertical error bands
  for( std::map<std::string, MnvVertErrorBand3D*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
  {
    const MnvVertErrorBand3D* err1 = h1->GetVertErrorBand( it->first );
    const MnvVertErrorBand3D* err2 = h2->GetVertErrorBand( it->first );
    if( !err1 || !err2 )
    {
      Error("Multiply", Form("Could not divide MnvH3Ds because they all don't have the %s MnvVertErrorBand3D", it->first.c_str()) );
      return;
    }
    it->second->Multiply( err1, err2, c1, c2 );
  }
  
  //! Scale the lateral error bands
  for( std::map<std::string, MnvLatErrorBand3D*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
  {
    const MnvLatErrorBand3D* err1 = h1->GetLatErrorBand( it->first );
    const MnvLatErrorBand3D* err2 = h2->GetLatErrorBand( it->first );
    if( !err1 || !err2 )
    {
      Error("Multiply", Form("Could not divide MnvH3Ds because they all don't have the %s MnvLatErrorBand3D", it->first.c_str()) );
      return;
    }
    it->second->Multiply( err1, err2, c1, c2 );
  }

  return;
}

void MnvH3D::Divide( const MnvH3D* h1, const MnvH3D* h2, Double_t c1 /*= 1*/, Double_t c2 /*= 1*/, Option_t* option /*=""*/ )
{
  //! @todo Would love to return a bool here, but we want this Divide to override TH1's and that is void

  //! Call the TH1D Divide
  this->TH3D::Divide( (TH3D*)h1, (TH3D*)h2, c1, c2, option );

  //! Scale the vertical error bands
  for( std::map<std::string, MnvVertErrorBand3D*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
  {
    const MnvVertErrorBand3D* err1 = h1->GetVertErrorBand( it->first );
    const MnvVertErrorBand3D* err2 = h2->GetVertErrorBand( it->first );
    if( !err1 || !err2 )
    {
      Error("Divide", Form("Could not divide MnvH3Ds because they all don't have the %s MnvVertErrorBand3D", it->first.c_str()) );
      return;
    }
    it->second->Divide( err1, err2, c1, c2, option );
  }
  
  //! Scale the lateral error bands
  for( std::map<std::string, MnvLatErrorBand3D*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
  {
    const MnvLatErrorBand3D* err1 = h1->GetLatErrorBand( it->first );
    const MnvLatErrorBand3D* err2 = h2->GetLatErrorBand( it->first );
    if( !err1 || !err2 )
    {
      Error("Divide", Form("Could not divide MnvH3Ds because they all don't have the %s MnvLatErrorBand3D", it->first.c_str()) );
      return;
    }
    it->second->Divide( err1, err2, c1, c2, option );
  }

  return;
}

void MnvH3D::DivideSingle( const MnvH3D* h1, const TH3* h2, Double_t c1 /*= 1*/, Double_t c2 /*= 1*/, Option_t* option /*=""*/ )
{
  //! Call the TH1D Divide
  this->TH3D::Divide( (TH3D*)h1, h2, c1, c2, option );

  //! Scale the vertical error bands
  for( std::map<std::string, MnvVertErrorBand3D*>::iterator it = fVertErrorBandMap.begin(); it != fVertErrorBandMap.end(); ++it )
  {
    const MnvVertErrorBand3D* err1 = h1->GetVertErrorBand( it->first );
    if( !err1 )
    {
      Error("Divide", Form("Could not divide MnvH3Ds because they all don't have the %s MnvVertErrorBand3D", it->first.c_str()) );
      return;
    }
    it->second->DivideSingle( err1, h2, c1, c2, option );
  }
  
  //! Scale the lateral error bands
  for( std::map<std::string, MnvLatErrorBand3D*>::iterator it = fLatErrorBandMap.begin(); it != fLatErrorBandMap.end(); ++it )
  {
    const MnvLatErrorBand3D* err1 = h1->GetLatErrorBand( it->first );
    if( !err1 )
    {
      Error("Divide", Form("Could not divide MnvH3Ds because they all don't have the %s MnvLatErrorBand3D", it->first.c_str()) );
      return;
    }
    it->second->DivideSingle( err1, h2, c1, c2, option );
  }

  return;
}



//--------------------------------------------------------
// trivial helper functions
//--------------------------------------------------------
bool MnvH3D::HasEnding (std::string const &fullString, std::string const &ending) const
{
  TString a(fullString);
  return a.EndsWith( ending.c_str() );
}

#endif
