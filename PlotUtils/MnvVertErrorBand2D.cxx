#ifndef MNV_MnvVertErrorBand2D_cxx
#define MNV_MnvVertErrorBand2D_cxx 1

#include "MnvVertErrorBand2D.h"
#include "HistogramUtils.h"
#include "Exceptions.h"
#include <algorithm>

using namespace MAT;


MnvVertErrorBand2D::MnvVertErrorBand2D( const std::string& name, const TH2D* base, const unsigned int nHists /* = 1000 */ ) :
  TH2D( *base )
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
    TH2D *tmp = new TH2D( *base );
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

MnvVertErrorBand2D::MnvVertErrorBand2D( const std::string& name, const TH2D* base, const std::vector<TH2D*>& hists ) :
TH2D (*base)
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

  std::vector<TH2D*>::const_iterator it = hists.begin();
  int it_pos = 0;
  for( ; it != hists.end(); ++it, ++it_pos )
  {
    sprintf(tmpName, "%s_universe%d", name.c_str(), it_pos );
    TH2D *tmp = new TH2D( **it );
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

MnvVertErrorBand2D::MnvVertErrorBand2D( const MnvVertErrorBand2D& h ) :
    TH2D( h )
{
  //!Deep copy the variables
  DeepCopy( h );
}

MnvVertErrorBand2D::~MnvVertErrorBand2D()
{
  while (fHists.size() > 0)
  {
    delete *(fHists.rbegin());
    fHists.pop_back();
  }
}

MnvVertErrorBand2D& MnvVertErrorBand2D::operator=( const MnvVertErrorBand2D& h ) 
{
  //! If this is me, no copy is needed
  if( this == &h )
    return *this;

  //! Call the base class's assignment
  TH2D::operator=(h);

  //! Delete and clear the hists vector
  for( unsigned int i = 0; i < h.fNHists; ++i )
    delete fHists[i];
  fHists.clear();
  fGoodColors.clear();

  DeepCopy( h );

  return *this;
}

void MnvVertErrorBand2D::DeepCopy( const MnvVertErrorBand2D& h )
{
  fUseSpreadError = h.GetUseSpreadError();
  fNHists = h.fNHists;
  for( unsigned int i = 0; i < fNHists; ++i )
    fHists.push_back( new TH2D(*h.GetHist(i)) );

  if (h.GetUnivWgts())
    fUnivWgts = *(h.GetUnivWgts());

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

bool MnvVertErrorBand2D::Fill( const double xval, const double yval, const double *weights, const double cvweight, double cvWeightFromMe )
{
  //! Fill the CV hist with the CV weight and value
  int cvbin = this->TH2D::Fill( xval, yval, cvweight );

  //! If cvbin is -1, it means we Filled the under or overflow, so use FindFind to get the right one
  //! @note FindBin is fast for under/overlow because there is no search
  if( cvbin == -1 )
    cvbin = FindBin( xval, yval );

  //! Add bin content to the bin for all the universes using their weights.
  //! Note that all universes will be filled in the same bin as the CV hist.
  const double applyWeight = cvweight / cvWeightFromMe;
  for( unsigned int i = 0; i != fNHists; ++i ) 
  {
      fHists[i]->AddBinContent( cvbin, weights[i]*applyWeight );
  }

  return true;
}

bool MnvVertErrorBand2D::Fill( const double xval, const double yval, const double weightDown, const double weightUp, const double cvweight, double cvWeightFromMe )
{
  //! Throw an exception if there are nUniverses is not 2
  if( fNHists != 2 )
  {
    std::cout << "ERROR [MnvVertErrorBand2D::Fill] - You used the specialized 2 shifts version, but you do not have 2 universes." << std::endl;
    throw 1;
  }

  //! put the weights in an array and use the standard Fill
  double weights[2] = {weightDown, weightUp};

  return Fill( xval, yval, weights, cvweight, cvWeightFromMe );
}

TH2D MnvVertErrorBand2D::GetErrorBand(bool asFrac /*= false*/ , bool cov_area_normalize /*= false*/) const
{
  TH2D errBand( *this );
  errBand.Reset();
  const int lowBin = 0;
  const int highBinX = GetNbinsX() + 1;
  const int highBinY = GetNbinsY() + 1;
  const int highBin  = GetBin( highBinX, highBinY );
  
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

const TH2D *MnvVertErrorBand2D::GetHist( unsigned int i ) const
{
  if( i >= fNHists )
  {
    Error("GetHist", Form("Cannot return histogram of universe %d because this object has %d universes.", i, fNHists) );
    return NULL;
  }
  return fHists[i];
}

TH2D *MnvVertErrorBand2D::GetHist( unsigned int i )
{
  if( i >= fNHists )
  {
    Error("GetHist", Form("Cannot return histogram of universe %d because this object has %d universes.", i, fNHists) );
    return NULL;
  }
  return fHists[i];
}

double MnvVertErrorBand2D::GetUnivWgt(unsigned int uid) const
{
  if (fUnivWgts.size() == 0)
    return -1;

  if (uid >= fUnivWgts.size())
    throw BadUnivWgtError("Requested universe index out of range");

  return fUnivWgts[uid];
}

const std::vector<double> * MnvVertErrorBand2D::GetUnivWgts() const
{
  if (fUnivWgts.size() == 0)
    return NULL;

  return &fUnivWgts;
}

void MnvVertErrorBand2D::SetUnivWgt(unsigned int uid, double wgt)
{
  // initialize all the elements with weight 1 if any aren't already initialized
  if (fUnivWgts.size() == 0)
    fUnivWgts = std::vector<double>(fNHists, 1.0);

  if (uid >= fUnivWgts.size())
    throw BadUnivWgtError("Requested universe index out of range");

  fUnivWgts.at(uid) = wgt;
}

void MnvVertErrorBand2D::SetUnivWgts(const std::vector<double>& wgts)
{
  if (wgts.size() != fNHists)
    throw BadUnivWgtError("Attempt to initialize weight list with different size vector than number of universes");

  fUnivWgts = wgts;
}

TMatrixD MnvVertErrorBand2D::CalcCovMx(bool area_normalize, bool asFrac) const
{
  //Calculating the Mean
  TH2D hmean = TH2D(*this);

  // if there's more than one varied universe, calculate the mean; if not, use the CV
  // histo as the 'mean'
  if(fNHists > 1)
  hmean.Reset();

  //!Storing Area Normalization Factors for the many universes
  //! @todo Need to Check this! 
  std::vector<double> normFactors;

  double totalWgtSum = 0;
  for( unsigned int hist_idx = 0; hist_idx < fNHists; ++hist_idx )
  {
    if (area_normalize)
    {
      double area_scale = fHists[hist_idx]->Integral();
      normFactors.push_back(area_scale/Integral());

      if( 0 < area_scale ) //just in case
      {
        fHists[hist_idx]->Scale(Integral()/area_scale);
      }
    }

    double univWgt = 1.0;
    if (fUnivWgts.size() > 0)
      univWgt = fUnivWgts[hist_idx];

    totalWgtSum += univWgt;

    if(fNHists>1)
      hmean.Add(fHists[hist_idx], univWgt);
  }

  if(fNHists>1)
    hmean.Scale(1.0/totalWgtSum);

  // Calculating Covariance Matrix
  //! @todo what to do with underflow( bin=0 ) and overlow ( bin=nbins+1 )?
  const int lowBin = 0; // considering underflow bin
  const int highBin = hmean.GetBin( hmean.GetNbinsX()+1, hmean.GetNbinsY()+1 ); // considering under/overflow

  TMatrixD covmx(highBin+1, highBin+1);
  
  if (fUseSpreadError)
  {
    if (fUnivWgts.size() > 0)
      throw NoWgtdSpreadError("Can't calculate spread errors for weighted universes!");

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
    for( unsigned int hist_idx = 0; hist_idx < fNHists; ++hist_idx )
    {
      double univWgt = 1.0;
      if (fUnivWgts.size() > 0)
        univWgt = fUnivWgts[hist_idx];

      // first calculate the weighted sum of squares of the deviations
      // (note that if no weights are specified, we use 1 for each univ)
      for( int i = lowBin; i <= highBin; ++i )
      {
        double xi=fHists[hist_idx]->GetBinContent(i);
        double ximean=hmean.GetBinContent(i);
        for( int k = i; k <= highBin; ++k )
        {
          double xk=fHists[hist_idx]->GetBinContent(k);
          double xkmean=hmean.GetBinContent(k);
          covmx[i][k] += univWgt * (xi-ximean)*(xk-xkmean);
        }
      }
    }

    // now normalize
    for( int i = lowBin; i <= highBin; ++i )
    {
      for( int k = i; k <= highBin; ++k )
      {
        covmx[i][k] /= totalWgtSum;
        covmx[k][i] = covmx[i][k];  // covariance matrices are symmetric
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

Bool_t MnvVertErrorBand2D::Add( const MnvVertErrorBand2D* h1, const Double_t c1 /*= 1.*/ )
{
  //! Check that we all have the same number of universes.
  if( h1->GetNHists() != this->GetNHists() )
  {
    Error("Add", "Attempt to Add with different numbers of universes" );
    return kFALSE;
  }

  //! Call Add on the CVHists
  this->TH2D::Add( h1, c1 );

  //! Call Add for all universes
  for( unsigned int iHist = 0; iHist != fNHists; ++iHist )
    fHists[iHist]->Add( h1->GetHist(iHist), c1 );

  return true;
}

Bool_t MnvVertErrorBand2D::Multiply( const MnvVertErrorBand2D* h1, const MnvVertErrorBand2D* h2, Double_t c1 /*= 1*/, Double_t c2 /*= 1*/ )
{
  //! Check that we all have the same number of universes.
  if( h1->GetNHists() != h2->GetNHists() || h1->GetNHists() != this->GetNHists() )
  {
    Error("Multiply", "Attempt to divide by different numbers of universes" );
    return kFALSE;
  }

  //! Call Multiply on the CVHists
  //! @note root documentation says this function returns a bool but its void in our version
  this->TH2D::Multiply( h1, h2, c1, c2 );

  //! Call Multiply for all universes
  for( unsigned int iHist = 0; iHist != fNHists; ++iHist )
    fHists[iHist]->Multiply( h1->GetHist(iHist), h2->GetHist(iHist), c1, c2 );

  return true;
}

Bool_t MnvVertErrorBand2D::MultiplySingle( const MnvVertErrorBand2D* h1, const TH2* h2, Double_t c1 /*= 1*/, Double_t c2 /*= 1*/ )
{
  // Check that we all have the same number of universes.
  if( h1->GetNHists() != this->GetNHists() )
  {
    Error("MultiplySingle", "Attempt to multiply by different numbers of universes" );
    return kFALSE; 
  }

  // Call Divide on the CVHists
  this->TH2D::Multiply( (TH2D*)h1, h2, c1, c2 );

  // Call Multiply for all universes
  for( unsigned int iHist = 0; iHist != fNHists; ++iHist )
    fHists[iHist]->Multiply( h1->GetHist(iHist), h2, c1, c2 );

  return true;
}

Bool_t MnvVertErrorBand2D::DivideSingle( const MnvVertErrorBand2D* h1, const TH2* h2, Double_t c1 /*= 1*/, Double_t c2 /*= 1*/, Option_t* option /*=""*/ )
{
  //! Check that we all have the same number of universes.
  if( h1->GetNHists() != this->GetNHists() )
  {
    Error("DivideSingle", "Attempt to divide by different numbers of universes" );
    return kFALSE;
  }

  //! Call Divide on the CVHists
  //! @note root documentation says this function returns a bool but its void in our version
  this->TH2D::Divide( (TH2D*)h1, h2, c1, c2, option);

  //! Call Divide for all universes
  for( unsigned int iHist = 0; iHist != fNHists; ++iHist )
    fHists[iHist]->Divide( h1->GetHist(iHist), h2, c1, c2, option );

  return true;
}

Bool_t MnvVertErrorBand2D::Divide( const MnvVertErrorBand2D* h1, const MnvVertErrorBand2D* h2, Double_t c1 /*= 1*/, Double_t c2 /*= 1*/, Option_t* option /*=""*/ )
{
  //! Check that we all have the same number of universes.
  if( h1->GetNHists() != h2->GetNHists() || h1->GetNHists() != this->GetNHists() )
  {
    Error("Divide", "Attempt to divide by different numbers of universes" );
    return kFALSE;
  }

  //! Call Divide on the CVHists
  //! @note root documentation for recent versions says TH1D::Divide returns a bool but its void in our current version of ROOT (5.30)
  this->TH2D::Divide( h1, h2, c1, c2, option);

  //! Call Divide for all universes
  for( unsigned int iHist = 0; iHist != fNHists; ++iHist )
    fHists[iHist]->Divide( h1->GetHist(iHist), h2->GetHist(iHist), c1, c2, option );

  return true;
}

void MnvVertErrorBand2D::Reset( Option_t* option /* = "" */ )
{
  //Reset the base
  this->TH2D::Reset(option);

  //reset all universes
  for( std::vector<TH2D*>::iterator i = fHists.begin(); i != fHists.end(); ++i )
    (*i)->Reset(option);
}


void MnvVertErrorBand2D::Scale( Double_t c1 /*= 1.*/, Option_t *option /*= ""*/, Bool_t allUniv /*= true */ )
{ 
  // Scale the CVHist
  this->TH2D::Scale( c1, option );

  if (!allUniv)
    return;

  // Scale all universes
  for( unsigned int iHist = 0; iHist != fNHists; ++iHist )
    fHists[iHist]->Scale( c1, option );
}





#endif
