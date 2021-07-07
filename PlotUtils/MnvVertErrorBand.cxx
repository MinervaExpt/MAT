#ifndef MNV_MnvVertErrorBand_cxx
#define MNV_MnvVertErrorBand_cxx 1

#include "MnvVertErrorBand.h"
#include "HistogramUtils.h"
#include "Exceptions.h"
#include <algorithm>
#include <cassert>

using namespace PlotUtils;


MnvVertErrorBand::MnvVertErrorBand( const std::string& name, const TH1D* base, const unsigned int nHists /* = 1000 */ ) :
  TH1D( *base )
{
  SetName( name.c_str() );
  SetTitle( name.c_str() );

  SetBit(kMustCleanup, false);

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
    TH1D *tmp = new TH1D( *base );
    tmp->SetBit(kMustCleanup, false);
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

MnvVertErrorBand::MnvVertErrorBand( const std::string& name, const TH1D* base, const std::vector<TH1D*>& hists ) :
TH1D (*base)
{

  SetName( name.c_str() );
  SetTitle( name.c_str() );

  SetBit(kMustCleanup, false);

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

  std::vector<TH1D*>::const_iterator it = hists.begin();
  int it_pos = 0;
  for( ; it != hists.end(); ++it, ++it_pos )
  {
    sprintf(tmpName, "%s_universe%d", name.c_str(), it_pos );
    TH1D *tmp = new TH1D( **it );
    tmp->SetBit(kMustCleanup, false);
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


MnvVertErrorBand::MnvVertErrorBand( const MnvVertErrorBand& h ) :
  TH1D( h )
{
  //!Deep copy the variables
  DeepCopy( h );
}

MnvVertErrorBand::~MnvVertErrorBand()
{
  while (fHists.size() > 0)
  {
    delete *(fHists.rbegin());
    fHists.pop_back();
  }
}

MnvVertErrorBand& MnvVertErrorBand::operator=( const MnvVertErrorBand& h ) 
{
  //! If this is me, no copy is needed
  if( this == &h )
    return *this;

  //! Call the base class's assignment
  TH1D::operator=(h);

  //! Delete and clear the hists vector
  for( unsigned int i = 0; i < h.fNHists; ++i )
    delete fHists[i];
  fHists.clear();
  fGoodColors.clear();

  DeepCopy( h );

  return *this;
}

void MnvVertErrorBand::DeepCopy( const MnvVertErrorBand& h )
{
  fUseSpreadError = h.GetUseSpreadError();
  fNHists = h.fNHists;
  for( unsigned int i = 0; i < fNHists; ++i )
    fHists.push_back( new TH1D(*h.GetHist(i)) );

  if (h.GetUnivWgts())
    fUnivWgts = *(h.GetUnivWgts());

  //! meh.
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

const TH1D *MnvVertErrorBand::GetHist( unsigned int i ) const
{
  if( i >= fNHists )
  {
    Error("GetHist", Form("Cannot return histogram of universe %d because this object has %d universes.", i, fNHists) );
    return NULL;
  }
  return fHists[i];
}

TH1D *MnvVertErrorBand::GetHist( unsigned int i )
{
  if( i >= fNHists )
  {
    Error("GetHist", Form("Cannot return histogram of universe %d because this object has %d universes.", i, fNHists) );
    return NULL;
  }
  return fHists[i];
}


bool MnvVertErrorBand::Fill( const double val, const double *weights, const double cvweight /* = 1.0 */, const double cvweightFromMe /* = 1. */)
{
  //! Fill the CV hist with the CV weight and value
  int cvbin = this->TH1D::Fill( val, cvweight );

  //! If cvbin is -1, it means we Filled the under or overflow, so use FindFind to get the right one
  //! @note FindBin is fast for under/overlow because there is no search
  if( cvbin == -1 )
    cvbin = FindBin( val );

  //! Add bin content to the bin for all the universes using their weights.
  //! Note that all universes will be filled in the same bin as the CV hist.
  const double applyWeight = cvweight / cvweightFromMe;
  for( unsigned int i = 0; i != fNHists; ++i ) 
  {
    const double wgtU = weights[i]*applyWeight;
    fHists[i]->AddBinContent( cvbin, wgtU );

    const double err = fHists[i]->GetBinError(cvbin);
    const double newerr2 = err*err + wgtU*wgtU; 
    const double newerr = (0.<newerr2) ? sqrt(newerr2) : 0.;
    fHists[i]->SetBinError( cvbin, newerr );
  }

  return true;
}

bool MnvVertErrorBand::Fill( const double val, const double weightDown, const double weightUp, const double cvweight /*= 1.0*/, const double cvWeightFromMe /* = 1. */ )
{
  //! Throw an exception if there are nUniverses is not 2
  if( fNHists != 2 )
  {
    std::cout << "ERROR [MnvVertErrorBand::Fill] - You used the specialized 2 shifts version, but you do not have 2 universes." << std::endl;
    throw 1;
  }

  //! put the weights in an array and use the standard Fill
  double weights[2] = {weightDown, weightUp};

  return Fill( val, weights, cvweight, cvWeightFromMe );
}


TH1D MnvVertErrorBand::GetErrorBand( bool asFrac /* = false */ , bool cov_area_normalize /* = false */) const
{
  //! @todo what to do with underflow( bin=0 ) and overlow ( bin=nbins+1 )?
  TH1D errBand( *this );
  errBand.Reset();
  const int lowBin = 0;
  const int highBin = GetNbinsX() + 1;
  TMatrixD* tmp_covmx;
  try{
    tmp_covmx = new TMatrixD(CalcCovMx(cov_area_normalize,asFrac));
  } catch(ErrBandEmptyCVError e){
    std::cerr << e.what() << std::endl;
    std::cerr << "ERROR MnvVertErrorBand::GetErrorBand: Unable to calculate error band " << this->GetName() << std::endl;
    exit(1);
  }

  TMatrixD covmx = *tmp_covmx;

  for( int i = lowBin; i <= highBin; ++i )
  {
    double err = (covmx[i][i]>0.)? sqrt( covmx[i][i] ): 0.; //Protect against odd sqrt(0) = NaN errors.

    errBand.SetBinContent( i, err );
    errBand.SetBinError(i, 0.);
  }

  return errBand;
}

TMatrixD MnvVertErrorBand::CalcCovMx(bool area_normalize /* = false */ , bool asFrac /* = false */ ) const
{
  //If this band's CV histo is empty, and we want to use the spread error then we should not be here.
  if (fUseSpreadError && GetEntries() == 0){
    std::cerr << "ERROR From MnvVertErrorBand::CalcCovMx. Problem with " << GetName() << std::endl;
    throw ErrBandEmptyCVError("Vertical error band's CV histo is empty! Cannot Calculate Spread Error.");
  }

  //Calculating the Mean
  TH1D hmean = TH1D(*this);

  // if there's more than one varied universe, calculate the mean; if not, use the CV
  // histo as the 'mean'
  if(fNHists > 1)
    hmean.Reset();

  //!Storing Area Normalization Factors for the many universes
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

    //    std::cout << " universe " << hist_idx << " in error band " << this->GetName() << " has weight: " << univWgt << std::endl;
    totalWgtSum += univWgt;

    if(fNHists>1)
      hmean.Add(fHists[hist_idx], univWgt);
  }

  if(fNHists>1)
    hmean.Scale(1.0/totalWgtSum);

  // Calculating Covariance Matrix
  //! @todo what to do with underflow( bin=0 ) and overlow ( bin=nbins+1 )?
  const int lowBin = 0; // considering underflow bin
  const int highBin = hmean.GetNbinsX()+1; //considering overflow bin

  TMatrixD covmx(highBin+1, highBin+1); // this is #bins + 2(bins for under&overflow)

  if (fUseSpreadError)
  {
    if (fUnivWgts.size() > 0)
      throw NoWgtdSpreadError("Can't calculate spread errors for weighted universes!");

    std::vector< std::vector<double> > binValsVec(highBin+1);
    //calculate maximum and minum values for every bin
    for( int i = lowBin; i <= highBin; ++i )
    {
      std::vector<double> binVals;
      for( unsigned int hist_idx = 0; hist_idx < fNHists; ++hist_idx )
      {
        const double val = fHists[hist_idx]->GetBinContent(i);
        if( isnan(val) )
          Warning( "MnvVertErrorBand::CalcCovMx", Form( "%s is trying to add nan val in bin %d,%d", GetName(), i, hist_idx ) );
        else
          binVals.push_back( val );
      }
      //get the CV value for this bin
      const double cv = GetBinContent(i);
      binVals.push_back( cv );
      sort( binVals.begin(), binVals.end() );

      binValsVec.at(i) =  binVals ;
    }

    // For spread errors in only 1 universe take the full max spread
    // For spread errors in less than 10 universes take 1/2 the max spread
    // For spread errors with more than 10 universes, use the interquartile spread
    for( int i = lowBin; i <= highBin; ++i )
    {
      for( int k = i; k <= highBin; ++k )
      {
        if ( fNHists == 1 )
          covmx[i][k] = ( binValsVec.at(i).back()- binValsVec.at(i).front() ) * ( binValsVec.at(k).back()- binValsVec.at(k).front() );
        else if ( fNHists < 10 )
          covmx[i][k] = ( binValsVec.at(i).back() - binValsVec.at(i).front() ) * ( binValsVec.at(k).back() - binValsVec.at(k).front() ) / 4. ;
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
    //    std::cout << " normalizing by sum of weights: " << totalWgtSum << std::endl;
    // now normalize
    for( int i = lowBin; i <= highBin; ++i )
    {
      for( int k = i; k <= highBin; ++k )
      {
        covmx[i][k] /= totalWgtSum;
        covmx[k][i] = covmx[i][k];  // covariance matrices are symmetric
      }
    }
    //    covmx.Print();
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
  if (area_normalize)
  {
    for( unsigned int hist_idx = 0; hist_idx < fNHists; ++hist_idx ) {
      if (fHists[hist_idx]->Integral()!=0)
      {
        fHists[hist_idx]->Scale(normFactors[hist_idx]);
      }
    }
  }
  return covmx;
}

TMatrixD MnvVertErrorBand::CalcCorrMx(bool area_normalize /* = false */) const
{
  //! Getting covariance matrix
  TMatrixD* tmp_covmx;
  try{
    tmp_covmx = new TMatrixD(CalcCovMx(area_normalize));
  } catch(ErrBandEmptyCVError e){
    std::cerr << e.what() << std::endl;
    std::cerr << "ERROR MnvVertErrorBand::CalcCorrMx: Unable to calc corr matrix for " << this->GetName() << std::endl;
    exit(1);
  }

  TMatrixD covmx = *tmp_covmx;

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

void MnvVertErrorBand::DrawAll( const char *option /* = "" */, bool drawCV /* = false */, bool area_normalize /* = false */, double normBinWidth /* = 0.0 */ ) const
{

  //! make a copy of each universe
  std::vector<TH1D*> histsCopy;
  TH1D *tallest(NULL);
  double maxVal = 0.;
  for( unsigned int i = 0; i < fNHists; ++i )
  {
    TH1D* histCopy = (TH1D*)fHists[i]->Clone( Form( "%s_tmp", fHists[i]->GetName() ) );

    //! area normalize universe hist if desired
    if( area_normalize && histCopy->Integral()!=0 )
    {
      histCopy->Scale( Integral( 1, GetNbinsX()) / histCopy->Integral( 1, histCopy->GetNbinsX() ) );
    }

    //! bin width normalize universe hist if desired
    if( normBinWidth > 0 )
    {
      histCopy->Scale( normBinWidth, "width" );
    }

    //! copy the CVHist axes attributes to each copy so they propagetes to the draw and tallest finding
    GetXaxis()->Copy( *histCopy->GetXaxis() );
    GetYaxis()->Copy( *histCopy->GetYaxis() );

    //! find the tallest universe hist
    if( histCopy->GetBinContent( histCopy->GetMaximumBin() ) > maxVal )
    {
      tallest = histCopy;
      maxVal  = histCopy->GetBinContent( histCopy->GetMaximumBin() );
    }

    histCopy->SetMinimum( 0 );
    histCopy->SetFillColor( 0 );
    histsCopy.push_back( histCopy );
  }

  //! make a copy of cv
  TH1D* cvcopy = (TH1D*)Clone( Form( "%s_tmp", GetName() ) );
  cvcopy->SetMinimum( 0 );
  cvcopy->SetLineWidth( 2.0 );
  cvcopy->SetFillColor( 0 );
  cvcopy->SetLineColor( kBlack );

  //! bin width normalize cv if desired
  if( normBinWidth > 0 )
  {
    cvcopy->Scale( normBinWidth, "width" );
  }

  if( ! tallest )
  {
    std::cout << "ERROR [MnvVertErrorBand::DrawAll] : Could not find any histograms with entries.  Nothing to draw." << std::endl;
    std::cout << "      Perhaps you didn't fill or fill entirely in the under/overlow?" << std::endl;
    return;
  }

  //! all but the first get the SAME option
  char optionSAME[156];
  sprintf( optionSAME, "SAME %s", option );

  if( drawCV )
  {
    //! draw the CV hist if desired
    if( cvcopy->GetBinContent( cvcopy->GetMaximumBin() ) > maxVal )
    {
      //! if CVHist is the tallest, draw it and mark it as tallest
      tallest = cvcopy;
      tallest->Draw( option );//note not DrawCopy because it will be redrawn below
    }
    else
    {
      //! otherwise, draw the tallest, then the CV hist
      tallest->DrawCopy( option );
      cvcopy->Draw( optionSAME ); //note not DrawCopy because it will be redrawn below
    }
  }
  else
  {
    //! if not drawing CV hist, just draw the tallest now
    tallest->DrawCopy();
  }

  //! now draw the rest
  for( unsigned int i = 0; i < fNHists; ++i )
  {
    if( tallest == fHists[i] )
      continue;
    histsCopy[i]->DrawCopy( optionSAME );
  }

  //! re-draw CV so it is visible
  if( drawCV ){
    cvcopy->DrawCopy( optionSAME );
  }

  //! clean up
//  for( unsigned int i = 0; i < fNHists; ++i ){
//    delete histsCopy[i];
//  }
//  delete cvcopy;

}


double MnvVertErrorBand::GetUnivWgt(unsigned int uid) const
{
  if (fUnivWgts.size() == 0)
    return -1;

  if (uid >= fUnivWgts.size())
    throw BadUnivWgtError("Requested universe index out of range");

  return fUnivWgts[uid];
}

const std::vector<double> * MnvVertErrorBand::GetUnivWgts() const
{
  if (fUnivWgts.size() == 0)
    return NULL;

  return &fUnivWgts;
}

void MnvVertErrorBand::SetUnivWgt(unsigned int uid, double wgt)
{
  // initialize all the elements with weight 1 if any aren't already initialized
  if (fUnivWgts.size() == 0)
    fUnivWgts = std::vector<double>(fNHists, 1.0);

  if (uid >= fUnivWgts.size())
    throw BadUnivWgtError("Requested universe index out of range");

  fUnivWgts.at(uid) = wgt;
}

void MnvVertErrorBand::SetUnivWgts(const std::vector<double>& wgts)
{
  if (wgts.size() != fNHists)
    throw BadUnivWgtError("Attempt to initialize weight list with different size vector than number of universes");

  fUnivWgts = wgts;
}


void MnvVertErrorBand::Scale( Double_t c1 /*= 1.*/, Option_t *option /*= ""*/, Bool_t allUniv /* = true */ )
{ 
  //! Scale the CVHist
  this->TH1D::Scale( c1, option );

  if (!allUniv)
    return;
    
  //! Scale all universes
  for( unsigned int iHist = 0; iHist != fNHists; ++iHist )
    fHists[iHist]->Scale( c1, option );
}

Bool_t MnvVertErrorBand::Divide( const MnvVertErrorBand* h1, const MnvVertErrorBand* h2, Double_t c1 /*= 1*/, Double_t c2 /*= 1*/, Option_t* option /*=""*/ )
{
  //! Check that we all have the same number of universes.
  if( h1->GetNHists() != h2->GetNHists() || h1->GetNHists() != this->GetNHists() )
  {
    const char* err = Form(
        "Attempt to divide by different number of universes.\n"
        "Clue about which error band is the problem: %s.", GetName());
      std::cout << " # of error bands " << h1->GetNHists() << " " << h2->GetNHists() << " " << this->GetNHists() <<  std::endl;
    Error("Divide", err);
    return kFALSE;
  }

  //! Call Divide on the CVHists
  //! @note root documentation for recent versions says TH1D::Divide returns a bool but its void in our current version of ROOT (5.30)
  this->TH1D::Divide( h1, h2, c1, c2, option);

  //! Call Divide for all universes
  for( unsigned int iHist = 0; iHist != fNHists; ++iHist )
    fHists[iHist]->Divide( h1->GetHist(iHist), h2->GetHist(iHist), c1, c2, option );

  return true;
}

Bool_t MnvVertErrorBand::Divide( TF1 *f1, Double_t c1  )
{
  //! Call Divide on the CVHists
  //! @note root documentation for recent versions says TH1D::Divide returns a bool but its void in our current version of ROOT (5.30)
  this->TH1D::Divide(f1, c1);

  //! Call Divide for all universes
  for( unsigned int iHist = 0; iHist != fNHists; ++iHist )
    fHists[iHist]->Divide(f1, c1);

  return true;
}

Bool_t MnvVertErrorBand::DivideSingle( const MnvVertErrorBand* h1, const TH1* h2, Double_t c1 /*= 1*/, Double_t c2 /*= 1*/, Option_t* option /*=""*/ )
{
  //! Check that we all have the same number of universes.
  if( h1->GetNHists() != this->GetNHists() )
  {
    Error("DivideSingle", "Attempt to divide by different numbers of universes" );
    return kFALSE;
  }

  //! Call Divide on the CVHists
  //! @note root documentation says this function returns a bool but its void in our version
  this->TH1D::Divide( (TH1D*)h1, h2, c1, c2, option);

  //! Call Divide for all universes
  for( unsigned int iHist = 0; iHist != fNHists; ++iHist )
    fHists[iHist]->Divide( h1->GetHist(iHist), h2, c1, c2, option );

  return true;
}


Bool_t MnvVertErrorBand::Multiply( const MnvVertErrorBand* h1, const MnvVertErrorBand* h2, Double_t c1 /*= 1*/, Double_t c2 /*= 1*/ )
{
  //! Check that we all have the same number of universes.
  if( h1->GetNHists() != h2->GetNHists() || h1->GetNHists() != this->GetNHists() )
  {
    Error("Multiply", "Attempt to divide by different numbers of universes" );
    return kFALSE;
  }

  //! Call Multiply on the CVHists
  //! @note root documentation says this function returns a bool but its void in our version
  this->TH1D::Multiply( h1, h2, c1, c2 );

  //! Call Multiply for all universes
  for( unsigned int iHist = 0; iHist != fNHists; ++iHist )
    fHists[iHist]->Multiply( h1->GetHist(iHist), h2->GetHist(iHist), c1, c2 );

  return true;
}


Bool_t MnvVertErrorBand::MultiplySingle( const MnvVertErrorBand* h1, const TH1* h2, Double_t c1 /*= 1*/, Double_t c2 /*= 1*/ )
{
  // Check that we all have the same number of universes.
  if( h1->GetNHists() != this->GetNHists() )
  {
    Error("MultiplySingle", "Attempt to multiply by different numbers of universes" );
    return kFALSE; 
  }

  // Call Divide on the CVHists
  this->TH1D::Multiply( (TH1D*)h1, h2, c1, c2 );

  // Call Multiply for all universes
  for( unsigned int iHist = 0; iHist != fNHists; ++iHist )
    fHists[iHist]->Multiply( h1->GetHist(iHist), h2, c1, c2 );

  return true;
}


Bool_t MnvVertErrorBand::AddSingle( const TH1* h1, const Double_t c1 /*= 1.*/ )
{
  //add to CV
  this->TH1D::Add( h1, c1 );

  //add to all universes
  for( unsigned int iHist = 0; iHist != fNHists; ++iHist )
    fHists[iHist]->Add( h1, c1 );

  return true;
}


Bool_t MnvVertErrorBand::Add( const MnvVertErrorBand* h1, const Double_t c1 /*= 1.*/ )
{
  //! Check that we all have the same number of universes.
  if( h1->GetNHists() != this->GetNHists() )
  {
    Error("Add", "Attempt to Add with different numbers of universes" );
    return kFALSE;
  }

  //! Call Add on the CVHists
  this->TH1D::Add( h1, c1 );

  //! Call Add for all universes
  for( unsigned int iHist = 0; iHist != fNHists; ++iHist )
    fHists[iHist]->Add( h1->GetHist(iHist), c1 );

  return true;
}

MnvVertErrorBand * MnvVertErrorBand::Rebin( const Int_t ngroup /*= 2*/, const Double_t * xbins /*= 0*/  )
{
  MnvVertErrorBand * rval;
  MnvVertErrorBand * clone = NULL;
  if ( xbins )
  {
    clone = new MnvVertErrorBand(*this);
    rval = clone;
  }
  else
    rval = this;

  // Call Rebin on the CVHist
  TH1 * newH;
  if (clone)
    newH = this->TH1D::Rebin( ngroup, "tmp_vert_eb", xbins );
  else
    newH = this->TH1D::Rebin(ngroup);

  if (newH == NULL)
    return NULL;
  if (clone)
  {
    dynamic_cast<TH1D*>(newH)->Copy(*clone);  // for some reason TH1::Copy is protected, but TH1D::Copy isn't
    clone->SetName(this->GetName());
    delete newH;
  }

  //we only need the rebin warning once
  const int oldVerbosity = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kError;

  // Call Rebin for all universes
  for( unsigned int iHist = 0; iHist != fNHists; ++iHist )
  {
    if (clone)
      newH = fHists[iHist]->Rebin( ngroup, "tmp", xbins );
    else
      newH = fHists[iHist]->Rebin(ngroup);
    if (newH == NULL)
      return NULL;
    if (clone)
    {
      dynamic_cast<TH1D*>(newH)->Copy(*(clone->GetHist(iHist)));
      clone->GetHist(iHist)->SetName(fHists[iHist]->GetName());
      delete newH;
    }
  }

  gErrorIgnoreLevel = oldVerbosity;

  return rval;
}

void MnvVertErrorBand::Reset( Option_t* option /* = "" */ )
{
  //Reset the base
  this->TH1D::Reset(option);

  //reset all universes
  for( std::vector<TH1D*>::iterator i = fHists.begin(); i != fHists.end(); ++i )
    (*i)->Reset(option);
}

void MnvVertErrorBand::SetBit( UInt_t f, Bool_t set)
{
  //Set the base class bit
  this->TH1D::SetBit(f,set);

  //set all universes' bits
  for( std::vector<TH1D*>::iterator i = fHists.begin(); i != fHists.end(); ++i )
    (*i)->SetBit(f,set);
}

#endif
