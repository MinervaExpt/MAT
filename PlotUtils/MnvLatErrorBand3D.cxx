#ifndef MNV_MnvLatErrorBand3D_cxx
#define MNV_MnvLatErrorBand3D_cxx 1

#include "PlotUtils/MnvLatErrorBand3D.h"
#include "PlotUtils/HistogramUtils.h"
#include <algorithm>

using namespace PlotUtils;


MnvLatErrorBand3D::MnvLatErrorBand3D( const std::string& name, const TH3D* base, const unsigned int nHists /* = 1000 */ ) :
  TH3D( *base )
{
  SetName( name.c_str() );
  SetTitle( name.c_str() );

  fNHists = nHists; 
  char tmpName[256];

  //set the good colors
  if( fGoodColors.size() == 0 )
  {
    fGoodColors.push_back( 2 );
    fGoodColors.push_back( 4 );
    fGoodColors.push_back( 6 );
    fGoodColors.push_back( 8 );
    fGoodColors.push_back( 9 );
    for( int i = 20; i < 50; i++ )
      fGoodColors.push_back( i );
  }

  for( unsigned int i = 0; i < fNHists; i++ )
  {
    sprintf(tmpName, "%s_universe%d", name.c_str(), i );
    TH3D *tmp = new TH3D( *base );
    tmp->Reset();
    tmp->SetName(tmpName);

    //give the universe histos a style and color
    tmp->SetLineColor( fGoodColors[ i % fGoodColors.size() ] );
    tmp->SetLineStyle( i % 10 + 1 );

    fHists.push_back( tmp );
  }

  if( nHists == 1 )
    fUseSpreadError = true;
  else
    fUseSpreadError = false;

}

MnvLatErrorBand3D::MnvLatErrorBand3D( const std::string& name, const TH3D* base, const std::vector<TH3D*>& hists ) :
TH3D (*base)
{

  SetName( name.c_str() );
  SetTitle( name.c_str() );

  fNHists = hists.size();
  char tmpName[256];

  //set the good colors
  if( fGoodColors.size() == 0 )
  {
    fGoodColors.push_back( 2 );
    fGoodColors.push_back( 4 );
    fGoodColors.push_back( 6 );
    fGoodColors.push_back( 8 );
    fGoodColors.push_back( 9 );
    for( int i = 20; i < 50; i++ )
      fGoodColors.push_back( i );
  }

  std::vector<TH3D*>::const_iterator it = hists.begin();
  int it_pos = 0;
  for( ; it != hists.end(); ++it, ++it_pos )
  {
    sprintf(tmpName, "%s_universe%d", name.c_str(), it_pos );
    TH3D *tmp = new TH3D( **it );
    tmp->SetName(tmpName);

    //give the universe histos a style and color
    tmp->SetLineColor( fGoodColors[ it_pos % fGoodColors.size() ] );
    tmp->SetLineStyle( it_pos % 10 + 1 );

    fHists.push_back( tmp );
  }

  if( fNHists == 1 )
    fUseSpreadError = true;
  else
    fUseSpreadError = false;

}

MnvLatErrorBand3D::MnvLatErrorBand3D( const MnvLatErrorBand3D& h ) :
    TH3D( h )
{
  //!Deep copy the variables
  DeepCopy( h );
}

MnvLatErrorBand3D& MnvLatErrorBand3D::operator=( const MnvLatErrorBand3D& h ) 
{
  //! If this is me, no copy is needed
  if( this == &h )
    return *this;

  //! Call the base class's assignment
  TH3D::operator=(h);

  //! Delete and clear the hists vector
  for( unsigned int i = 0; i < h.fNHists; ++i )
    delete fHists[i];
  fHists.clear();
  fGoodColors.clear();

  DeepCopy( h );

  return *this;
}

void MnvLatErrorBand3D::DeepCopy( const MnvLatErrorBand3D& h )
{
  fUseSpreadError = h.GetUseSpreadError();
  fNHists = h.fNHists;
  for( unsigned int i = 0; i < fNHists; ++i )
    fHists.push_back( new TH3D(*h.GetHist(i)) );

  //set the good colors
  if( fGoodColors.size() == 0 )
  {
    fGoodColors.push_back( 2 );
    fGoodColors.push_back( 4 );
    fGoodColors.push_back( 6 );
    fGoodColors.push_back( 8 );
    fGoodColors.push_back( 9 );
    for( int i = 20; i < 50; i++ )
      fGoodColors.push_back( i );
  }
}

bool MnvLatErrorBand3D::Fill( const double xval, const double yval, const double zval, const double *xshifts, const double *yshifts, const double *zshifts, const double cvweight /*= 1.0*/, const bool fillcv /*= true*/, const double* weights /*= NULL*/ )
{
  //! Fill the CV hist with the CV weight and value
  if( fillcv ) 
    this->TH3D::Fill( xval, yval, zval, cvweight );

  //! Use FindBin to get the CV bin 
  //const int cvbin = FindBin( xval, yval, zval );

  //! Add to the bin content of all the universes
  //! Assume that the bin we will fill is close to the bin at the central value
  //!@todo find a clever way to find the bin based on the cvbin information

  for( unsigned int i = 0; i != fNHists; ++i )
  {
    //!@todo Check consistency of this condition
    if( MnvHist::IsNotPhysicalShift( xshifts[i] ) || MnvHist::IsNotPhysicalShift( yshifts[i] ) || MnvHist::IsNotPhysicalShift( zshifts[i] ) )
      continue;

    const double x_shiftVal = xval + xshifts[i];
    const double y_shiftVal = yval + yshifts[i];
    const double z_shiftVal = zval + zshifts[i];
    //!@todo find a clever way to find the bin based on the cvbin information
    int bin = FindBin( x_shiftVal, y_shiftVal, z_shiftVal );
    
    if( 0==weights )
    {
      fHists[i]->AddBinContent( bin, cvweight );
    }
    else
    {
      fHists[i]->AddBinContent( bin, cvweight*weights[i] );
    }

  }

  return true;
}

bool MnvLatErrorBand3D::Fill( const double xval, const double yval, const double zval, const double xshiftDown, const double xshiftUp, const double yshiftDown, const double yshiftUp, const double zshiftDown, const double zshiftUp, const double cvweight /*= 1.0*/, const bool fillcv /*= true*/ )
{
  //! Throw an exception if nUniverses is not 2
  if( fNHists != 2 )
  {
    std::cout << "ERROR [MnvLatErrorBand3D::Fill] - You used the specialized 2 shifts version, but you do not have 2 universes." << std::endl;
    throw 1;
  }

  //! put the weights in an array and use the standard Fill
  double xshifts[2] = {xshiftDown, xshiftUp};
  double yshifts[2] = {yshiftDown, yshiftUp};
  double zshifts[2] = {zshiftDown, zshiftUp};

  return Fill( xval, yval, zval, xshifts, yshifts, zshifts, cvweight, fillcv);
}

TH3D MnvLatErrorBand3D::GetErrorBand(bool asFrac /*= false*/ , bool cov_area_normalize /*= false*/) const
{
  TH3D errBand( *this );
  errBand.Reset();
  const int lowBin = 0;
  const int highBinX = GetNbinsX() + 1;
  const int highBinY = GetNbinsY() + 1;
  const int highBinZ = GetNbinsZ() + 1;
  const int highBin  = GetBin( highBinX, highBinY, highBinZ );
  
  TMatrixD covmx(highBin+1,highBin+1);
  covmx = CalcCovMx(cov_area_normalize,asFrac);

  for( int i = lowBin; i <= highBin; ++i )
  {
    double err = (covmx[i][i]>0.)? sqrt( covmx[i][i] ): 0.; //Protect against odd sqrt(0) = NaN errors.

    errBand.SetBinContent( i, err );
    errBand.SetBinError(i, 0.);
  }

  return errBand;
}

const TH3D *MnvLatErrorBand3D::GetHist( unsigned int i ) const
{
  if( i >= fNHists )
  {
    Error("GetHist", Form("Cannot return histogram of universe %d because this object has %d universes.", i, fNHists) );
    return NULL;
  }
  return fHists[i];
}

TH3D *MnvLatErrorBand3D::GetHist( unsigned int i )
{
  if( i >= fNHists )
  {
    Error("GetHist", Form("Cannot return histogram of universe %d because this object has %d universes.", i, fNHists) );
    return NULL;
  }
  return fHists[i];
}

TMatrixD MnvLatErrorBand3D::CalcCovMx(bool area_normalize, bool asFrac) const
{
  //Calculating the Mean
  TH3D hmean = TH3D(*this);

  // if there's more than one varied universe, calculate the mean; if not, use the CV
  // histo as the 'mean'
  if(fNHists > 1)
    hmean.Reset();

  //!Storing Area Normalization Factors for the many universes
  //! @todo Need to Check this! 
  std::vector<double> normFactors;

  for( unsigned int j = 0; j < fNHists; ++j ) {
    if (area_normalize)
    {
      double area_scale = fHists[j]->Integral();
      normFactors.push_back(area_scale/Integral());

      if (fHists[j]->Integral()!=0) //just in case
      {
        fHists[j]->Scale(Integral()/area_scale);
      }
    }
    if(fNHists>1)
      hmean.Add(fHists[j]);
  }

  if(fNHists>1)
    hmean.Scale(1.0/(double) fNHists);

  // Calculating Covariance Matrix
  //! @todo what to do with underflow( bin=0 ) and overlow ( bin=nbins+1 )?
  const int lowBin = 0; // considering underflow bin
  const int highBin = hmean.GetBin( hmean.GetNbinsX()+1, hmean.GetNbinsY()+1, hmean.GetNbinsZ()+1 );// considering under/overflow

  TMatrixD covmx(highBin+1, highBin+1);

  if (fUseSpreadError)
  {
    std::vector< std::vector<double> > binValsVec(highBin+1);
    //calculate maximum and minum values for every bin
    for( int i = lowBin; i <= highBin; ++i )
    {
      std::vector<double> binVals;
      for( unsigned int j = 0; j < fNHists; ++j )
      {
        const double val = fHists[j]->GetBinContent(i);
        binVals.push_back( val );
      }
      //get the CV value for this bin
      const double cv = GetBinContent(i);
      binVals.push_back( cv );
      sort( binVals.begin(), binVals.end() );

      binValsVec.at(i) =  binVals ;
    }

    //! For spread errors in only 1 universe take the full max spread
    //! For spread errors in less than 10 universes take 1/2 the max spread
    //! For spread errors with more than 10 universes, use the interquartile spread
    for( int i = lowBin; i <= highBin; ++i )
    {
      for( int k = i; k <= highBin; ++k )
      {
        if ( fNHists == 1 )
          covmx[i][k] = ( (binValsVec.at(i)).back()- (binValsVec.at(i)).front() ) * ( (binValsVec.at(k)).back()- (binValsVec.at(k)).front() );
        else if ( fNHists < 10 )
          covmx[i][k] = ( (binValsVec.at(i)).back() - (binValsVec.at(i)).front() ) * ( (binValsVec.at(k)).back() - (binValsVec.at(k)).front() ) / 4. ;
        else
          covmx[i][k] = MnvHist::GetInterquartileRange( binValsVec.at(i) ) * MnvHist::GetInterquartileRange( binValsVec.at(k) ) * pow(MnvHist::InterquartileRangeToSigma,2);

        covmx[k][i] = covmx[i][k];
      }
    }
  }
  else
  {
    //Calculating Covariance
    for( unsigned int j = 0; j < fNHists; ++j ) 
    {
      for( int i = lowBin; i <= highBin; ++i )
      {
        double xi=fHists[j]->GetBinContent(i);
        double ximean=hmean.GetBinContent(i);
        for( int k = i; k <= highBin; ++k )
        {
          double xk=fHists[j]->GetBinContent(k);
          double xkmean=hmean.GetBinContent(k);
          covmx[i][k] +=(xi-ximean)*(xk-xkmean);
        }
      }
    }

    for( int i = lowBin; i <= highBin; ++i )
    {
      for( int k = i; k <= highBin; ++k )
      {
        covmx[i][k]/=((double)fNHists);
        covmx[k][i]=covmx[i][k];
      }
    }
  }

  if (asFrac)
  {
    for( int i = lowBin; i <= highBin; ++i )
    {
      for( int k = i; k <= highBin; ++k )
      {
        //double sign = covmx[i][k]>0 ? 1. : -1. ;
        //Gettting the the CV value for bin i
        const double cv_i = GetBinContent(i);
        const double cv_k = GetBinContent(k);

        //covmx[i][k]=sign * sqrt(TMath::Abs(covmx[i][k])/(cv_i * cv_k));
        covmx[i][k]= ( cv_i != 0. && cv_k != 0. ) ? covmx[i][k]/(cv_i * cv_k) : 0.;
        covmx[k][i]=covmx[i][k];
      }
    }
  }

  //! Getting fHists back to normal
  //! @todo Need to Check this
  if (area_normalize)
  {
    for( unsigned int j = 0; j < fNHists; ++j ) {
      if (fHists[j]->Integral()!=0)
      {
        fHists[j]->Scale(normFactors[j]);
      }
    }
  }
  return covmx;
}

Bool_t MnvLatErrorBand3D::Add( const MnvLatErrorBand3D* h1, const Double_t c1 /*= 1.*/ )
{
  //! Check that we all have the same number of universes.
  if( h1->GetNHists() != this->GetNHists() )
  {
    Error("Add", "Attempt to Add with different numbers of universes" );
    return kFALSE;
  }

  //! Call Add on the CVHists
  this->TH3D::Add( h1, c1 );

  //! Call Add for all universes
  for( unsigned int iHist = 0; iHist != fNHists; ++iHist )
    fHists[iHist]->Add( h1->GetHist(iHist), c1 );

  return true;
}

Bool_t MnvLatErrorBand3D::Multiply( const MnvLatErrorBand3D* h1, const MnvLatErrorBand3D* h2, Double_t c1 /*= 1*/, Double_t c2 /*= 1*/ )
{
  //! Check that we all have the same number of universes.
  if( h1->GetNHists() != h2->GetNHists() || h1->GetNHists() != this->GetNHists() )
  {
    Error("Multiply", "Attempt to divide by different numbers of universes" );
    return kFALSE;
  }

  //! Call Multiply on the CVHists
  this->TH3D::Multiply( h1, h2, c1, c2 );

  //! Call Multiply for all universes
  for( unsigned int iHist = 0; iHist != fNHists; ++iHist )
    fHists[iHist]->Multiply( h1->GetHist(iHist), h2->GetHist(iHist), c1, c2 );

  return true;
}

Bool_t MnvLatErrorBand3D::DivideSingle( const MnvLatErrorBand3D* h1, const TH3* h2, Double_t c1 /*= 1*/, Double_t c2 /*= 1*/, Option_t* option /*=""*/ )
{
  //! Check that we all have the same number of universes.
  if( h1->GetNHists() != this->GetNHists() )
  {
    Error("DivideSingle", "Attempt to divide by different numbers of universes" );
    return kFALSE;
  }

  //! Call Divide on the CVHists
  this->TH3D::Divide( (TH3D*)h1, h2, c1, c2, option);

  //! Call Divide for all universes
  for( unsigned int iHist = 0; iHist != fNHists; ++iHist )
    fHists[iHist]->Divide( h1->GetHist(iHist), h2, c1, c2, option );

  return true;
}

Bool_t MnvLatErrorBand3D::Divide( const MnvLatErrorBand3D* h1, const MnvLatErrorBand3D* h2, Double_t c1 /*= 1*/, Double_t c2 /*= 1*/, Option_t* option /*=""*/ )
{
  //! Check that we all have the same number of universes.
  if( h1->GetNHists() != h2->GetNHists() || h1->GetNHists() != this->GetNHists() )
  {
    Error("Divide", "Attempt to divide by different numbers of universes" );
    return kFALSE;
  }

  //! Call Divide on the CVHists
  this->TH3D::Divide( h1, h2, c1, c2, option);

  //! Call Divide for all universes
  for( unsigned int iHist = 0; iHist != fNHists; ++iHist )
    fHists[iHist]->Divide( h1->GetHist(iHist), h2->GetHist(iHist), c1, c2, option );

  return true;
}

void MnvLatErrorBand3D::Scale( Double_t c1 /*= 1.*/, Option_t *option /*= ""*/, Bool_t allUniv /*= true*/ )
{ 
  // Scale the CVHist
  this->TH3D::Scale( c1, option );

  if (!allUniv)
    return;

  // Scale all universes
  for( unsigned int iHist = 0; iHist != fNHists; ++iHist )
    fHists[iHist]->Scale( c1, option );
}

#endif
