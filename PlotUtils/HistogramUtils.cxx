//#############################################################################
//
// HistogramFunctions.cpp - these functions are useful for working with ROOT
//                          histograms.  For example, dividing two histograms 
//                          and assigning binomial errors, or averaging over
//                          a dimension to reduce the dimensionality by 1. 
//
//
// D. Schmitz  October 28, 2005
//
//#############################################################################

#ifndef HISTOGRAMUTILS_cxx
#define HISTOGRAMUTILS_cxx

#include "HistogramUtils.h"
#include <TMath.h>
#include <algorithm>

using namespace MAT;

// initialize the static members here
const double MnvHist::NotPhysicalShiftNumber = -12345678.87654321;
const double MnvHist::AutoAxisLimit = -1111;  
const double MnvHist::InterquartileRangeToSigma = 1.0 / 1.34896;


bool MnvHist::IsNotPhysicalShift( double shift ){ 
    return abs(NotPhysicalShiftNumber - shift) == 0;
}

bool MnvHist::IsAutoAxisLimit( double x ){ 
  return fabs(AutoAxisLimit - x) < 1E-6; 
}



//===========================
double MnvHist::GetInterquartileRange( const std::vector<double>& vec ){
  if( vec.empty() ) {
    Warning( "MnvHist::GetInterquartileRange", "Attempt to get inter-quartile range for an empty vector.  Returning 0." );
    return 0.;
  }

  //make a copy so we can sort
  std::vector<double> vecCopy = vec;
  sort( vecCopy.begin(), vecCopy.end() );
  size_t firstQ  = vecCopy.size() / 4;
  size_t thrirdQ = 3 * ( vecCopy.size() / 4 );
  return vecCopy[thrirdQ] - vecCopy[firstQ];
}


//=============================================================================
// divide3D( )
//
// note that -99 is the flag for a 0 denominator. we need to distinguish this 
// case from the 0 / 4 case, for example.
// errors are binomial errors!!
//=============================================================================
TH3D* MnvHist::divide3D( const TH3D *num, const TH3D *den ){

  TH3D *result = (TH3D*)num->Clone();

  for( int i = 1; i <= num->GetNbinsX(); i++ )
  {
    for( int j = 1; j <= num->GetNbinsY(); j++ )
    {
      for( int k = 1; k <= num->GetNbinsZ(); k++ )
      {
        if(  den->GetBinContent(i,j,k) != 0 )
        {
          result->SetBinContent(i,j,k, num->GetBinContent(i,j,k)/den->GetBinContent(i,j,k) );
          double error = sqrt( fabs((num->GetBinContent(i,j,k) * ( 1.0 - num->GetBinContent(i,j,k) / den->GetBinContent(i,j,k) )))) / den->GetBinContent(i,j,k);
          result->SetBinError(i,j,k,error);
          //cout << i << ", " << j << ", " << k << " = " << num->GetBinContent(i,j,k) << " / " << den->GetBinContent(i,j,k) << " = " 
          //<< result->GetBinContent(i,j,k) << " +- " << error << endl; 
        }
        else
        {
          result->SetBinContent(i,j,k,-99);
          result->SetBinError(i,j,k,1);
        }
      }
    }
  }

  return result;

}


//=============================================================================
// divide2D( )
//
// note that -99 is the flag for a 0 denominator. we need to distinguish this 
// case from the 0 / 4 case, for example.
// errors are binomial errors!!
//=============================================================================
TH2D* MnvHist::divide2D( const TH2D *num, const TH2D *den ){

  TH2D *result = (TH2D*)num->Clone();

  for( int i = 1; i <= num->GetNbinsX(); i++ )
  {
    for( int j = 1; j <= num->GetNbinsY(); j++ )
    {
      if(  den->GetBinContent(i,j) != 0 )
      {
        result->SetBinContent(i,j, num->GetBinContent(i,j)/den->GetBinContent(i,j) );
        double error = sqrt( fabs((num->GetBinContent(i,j) * ( 1.0 - num->GetBinContent(i,j) / den->GetBinContent(i,j) )))) / den->GetBinContent(i,j);
        result->SetBinError(i,j,error);
        //cout << i << ", " << j << " = " << num->GetBinContent(i,j) << " / " << den->GetBinContent(i,j) << " = " 
        //<< result->GetBinContent(i,j) << " +- " << error << endl; 
      }
      else
      {
        result->SetBinContent(i,j,-99);
        result->SetBinError(i,j,1);
      }
    }
  }

  return result;

}


//=============================================================================
// divide1D( )
//
// note that -99 is the flag for a 0 denominator. we need to distinguish this 
// case from the 0 / 4 case, for example.
// errors are binomial errors!!
//=============================================================================
MnvH1D* MnvHist::divide1D( const MnvH1D *num, const MnvH1D *den ){

  MnvH1D *result = (MnvH1D*)num->Clone();

  for( int i = 1; i <= num->GetNbinsX(); i++ )
  {
    if(  den->GetBinContent(i) != 0 )
    {
      result->SetBinContent(i, num->GetBinContent(i)/den->GetBinContent(i) );
      double error = sqrt( fabs((num->GetBinContent(i) * ( 1.0 - num->GetBinContent(i) / den->GetBinContent(i) )))) / den->GetBinContent(i);
      result->SetBinError(i,error);
      //cout << i << " = " << num->GetBinContent(i) << " / " << den->GetBinContent(i) << " = " 
      //<< result->GetBinContent(i) << " +- " << error << endl; 
    }
    else
    {
      result->SetBinContent(i,-99);
      result->SetBinError(i,1);
    }
  }

  return result;

}



//=============================================================================
// average3D_to_2D()
//=============================================================================
TH2D* MnvHist::average3D_to_2D( const TH3D *object, const char *axis1, const char *axis2 ){

  int xbins=0, ybins=0, zbins=0;
  double* x_bins(0);
  double* y_bins(0);
  double* z_bins(0);
  char axis3[1];
  char axes[3];

  if( (0 == strcmp(axis1,"1") && 0 == strcmp(axis2,"2")) || (0 == strcmp(axis1,"2") && 0 == strcmp(axis2,"1")) )
    sprintf( axis3, "3" );
  if( (0 == strcmp(axis1,"1") && 0 == strcmp(axis2,"3")) || (0 == strcmp(axis1,"3") && 0 == strcmp(axis2,"1")) )
    sprintf( axis3, "2" );
  if( (0 == strcmp(axis1,"2") && 0 == strcmp(axis2,"3")) || (0 == strcmp(axis1,"3") && 0 == strcmp(axis2,"2")) )
    sprintf( axis3, "1" );

  sprintf( axes, "%s%s%s", axis1, axis2, axis3 );

  if( 0 == strcmp(axes,"123") )
  { 
    xbins = object->GetNbinsX();
    x_bins = new double[xbins+1];
    for( int i = 0; i <= xbins; i++ )
      x_bins[i] = object->GetXaxis()->GetBinLowEdge(i+1);
    ybins = object->GetNbinsY();
    y_bins = new double[ybins+1];
    for( int i = 0; i <= ybins; i++ )
      y_bins[i] = object->GetYaxis()->GetBinLowEdge(i+1);
    zbins = object->GetNbinsZ();
    z_bins = new double[zbins+1];
    for( int i = 0; i <= zbins; i++ )
      z_bins[i] = object->GetZaxis()->GetBinLowEdge(i+1);
  }
  else if( 0 == strcmp(axes,"132") )
  { 
    xbins = object->GetNbinsX();
    x_bins = new double[xbins+1];
    for( int i = 0; i <= xbins; i++ )
      x_bins[i] = object->GetXaxis()->GetBinLowEdge(i+1);
    ybins = object->GetNbinsZ();
    y_bins = new double[ybins+1];
    for( int i = 0; i <= ybins; i++ )
      y_bins[i] = object->GetZaxis()->GetBinLowEdge(i+1);
    zbins = object->GetNbinsY();
    z_bins = new double[zbins+1];
    for( int i = 0; i <= zbins; i++ )
      z_bins[i] = object->GetYaxis()->GetBinLowEdge(i+1);
  }
  else if( 0 == strcmp(axes,"213") )
  {
    xbins = object->GetNbinsY();
    x_bins = new double[xbins+1];
    for( int i = 0; i <= xbins; i++ )
      x_bins[i] = object->GetYaxis()->GetBinLowEdge(i+1);
    ybins = object->GetNbinsX();
    y_bins = new double[ybins+1];
    for( int i = 0; i <= ybins; i++ )
      y_bins[i] = object->GetXaxis()->GetBinLowEdge(i+1);
    zbins = object->GetNbinsZ();
    z_bins = new double[zbins+1];
    for( int i = 0; i <= zbins; i++ )
      z_bins[i] = object->GetZaxis()->GetBinLowEdge(i+1);
  }
  else if( 0 == strcmp(axes,"231") )
  { 
    xbins = object->GetNbinsY();
    x_bins = new double[xbins+1];
    for( int i = 0; i <= xbins; i++ )
      x_bins[i] = object->GetYaxis()->GetBinLowEdge(i+1);
    ybins = object->GetNbinsZ();
    y_bins = new double[ybins+1];
    for( int i = 0; i <= ybins; i++ )
      y_bins[i] = object->GetZaxis()->GetBinLowEdge(i+1);
    zbins = object->GetNbinsX();
    z_bins = new double[zbins+1];
    for( int i = 0; i <= zbins; i++ )
      z_bins[i] = object->GetXaxis()->GetBinLowEdge(i+1);
  }
  else if( 0 == strcmp(axes,"312") )
  {
    xbins = object->GetNbinsZ();
    x_bins = new double[xbins+1];
    for( int i = 0; i <= xbins; i++ )
      x_bins[i] = object->GetZaxis()->GetBinLowEdge(i+1);
    ybins = object->GetNbinsX();
    y_bins = new double[ybins+1];
    for( int i = 0; i <= ybins; i++ )
      y_bins[i] = object->GetXaxis()->GetBinLowEdge(i+1);
    zbins = object->GetNbinsY();
    z_bins = new double[zbins+1];
    for( int i = 0; i <= zbins; i++ )
      z_bins[i] = object->GetYaxis()->GetBinLowEdge(i+1);
  }
  else if( 0 == strcmp(axes,"321") )
  { 
    xbins = object->GetNbinsZ();
    x_bins = new double[xbins+1];
    for( int i = 0; i <= xbins; i++ )
      x_bins[i] = object->GetZaxis()->GetBinLowEdge(i+1);
    ybins = object->GetNbinsY();
    y_bins = new double[ybins+1];
    for( int i = 0; i <= ybins; i++ )
      y_bins[i] = object->GetYaxis()->GetBinLowEdge(i+1);
    zbins = object->GetNbinsX();
    z_bins = new double[zbins+1];
    for( int i = 0; i <= zbins; i++ )
      z_bins[i] = object->GetXaxis()->GetBinLowEdge(i+1);
  }
  else
  {
    cout << " Unhandled axis settings.  here' a NULL pointer." << endl;
    return NULL;
  }


  TH2D *result = new TH2D("result", "result", xbins, x_bins, ybins, y_bins );

  double sum = 0;
  int count = 0;
  double errorSum2=0;

  for( int i = 1; i <= xbins; i++ )
  {
    for( int j = 1; j <= ybins; j++ )
    {
      count = 0;
      sum = 0;
      errorSum2 = 0;
      for( int k = 1; k <= zbins; k++ )
      {
        if( 0 == strcmp(axes,"123") ){
          if( object->GetBinContent(i,j,k) != -99 ){
            sum += object->GetBinContent(i,j,k);
            count++;
            errorSum2 += pow(object->GetBinError(i,j,k),2);
          }
        }else if( 0 == strcmp(axes,"132") ){
          if( object->GetBinContent(i,k,j) != -99 ){
            sum += object->GetBinContent(i,k,j);
            count++;
            errorSum2 += pow(object->GetBinError(i,k,j),2);
          }
        }else if( 0 == strcmp(axes,"213") ){
          if( object->GetBinContent(j,i,k) != -99 ){
            sum += object->GetBinContent(j,i,k);
            count++;
            errorSum2 += pow(object->GetBinError(j,i,k),2);
          }
        }else if( 0 == strcmp(axes,"231") ){
          if( object->GetBinContent(k,i,j) != -99 ){
            sum += object->GetBinContent(k,i,j);
            count++;
            errorSum2 += pow(object->GetBinError(k,i,j),2);
          }
        }else if( 0 == strcmp(axes,"312") ){
          if( object->GetBinContent(j,k,i) != -99 ){
            sum += object->GetBinContent(j,k,i);
            count++;
            errorSum2 += pow(object->GetBinError(j,k,i),2);
          }
        }else if( 0 == strcmp(axes,"321") ){
          if( object->GetBinContent(k,j,i) != -99 ){
            sum += object->GetBinContent(k,j,i);
            count++;
            errorSum2 += pow(object->GetBinError(k,j,i),2);
          }
        }
      }
      if( count != 0 )
      {
        result->SetBinContent(i,j,sum/(double)count);
        double error = sqrt( errorSum2 ) / count;
        result->SetBinError(i,j, error );
        //cout << i << ", " << j << " = " << sum << " / " << count << " = " 
        //<< result->GetBinContent(i,j) << " +- " << error << endl;
      }
      else
      {
        result->SetBinContent(i,j,0);
        result->SetBinError(i,j,0);
      }
    }
  }

  return result;
}


//=============================================================================
// average2D_to_1D()
//=============================================================================
MnvH1D* MnvHist::average2D_to_1D( const TH2D* object, const char *axis1 ){

  int xbins=0, ybins=0;
  double* x_bins(0);
  double* y_bins(0);
  char axis2[1];

  if( 0 == strcmp(axis1,"1") )
  {
    sprintf( axis2, "2" );
    xbins = object->GetNbinsX();
    x_bins = new double[xbins+1];
    for( int i = 0; i <= xbins; i++ )
      x_bins[i] = object->GetXaxis()->GetBinLowEdge(i+1);
    ybins = object->GetNbinsY();
    y_bins = new double[ybins+1];
    for( int i = 0; i <= ybins; i++ )
      y_bins[i] = object->GetYaxis()->GetBinLowEdge(i+1);
  }
  else
  {
    sprintf( axis2, "1" );
    xbins = object->GetNbinsY();
    x_bins = new double[xbins+1];
    for( int i = 0; i <= xbins; i++ )
      x_bins[i] = object->GetYaxis()->GetBinLowEdge(i+1);
    ybins = object->GetNbinsX();
    y_bins = new double[ybins+1];
    for( int i = 0; i <= ybins; i++ )
      y_bins[i] = object->GetXaxis()->GetBinLowEdge(i+1);
  }

  MnvH1D *result = new MnvH1D( "result", "result", xbins, x_bins );

  double sum = 0;
  int count = 0;
  double errorSum2=0;

  for( int i = 1; i <= xbins; i++ )
  {
    sum = 0;
    count = 0;
    errorSum2 = 0;
    for( int j = 1; j <= ybins; j++ )
    {
      if( 0 == strcmp(axis1,"1") ){
        if( object->GetBinContent(i,j) != -99 )
        {
          sum += object->GetBinContent(i,j);
          count++;
          errorSum2 += pow(object->GetBinError(i,j),2);
        }
      }else{
        if( object->GetBinContent(j,i) != -99 )
        {
          sum += object->GetBinContent(j,i);
          count++;
          errorSum2 += pow(object->GetBinError(j,i),2);
        }
      }
    }
    if( count != 0 )
    {
      result->SetBinContent(i,sum/(double)count);
      double error = sqrt( errorSum2 ) / count;
      result->SetBinError(i,error );
      //cout << i << " = " << sum << " / " << count << " = " 
      //<< result->GetBinContent(i) << " +- " << error << endl;
    }
    else
    {
      result->SetBinContent(i,-99);
      result->SetBinError(i,1);
    }
  }

  return result;

}


//=============================================================================
// integrate2D_to_1D()
//=============================================================================
MnvH1D* MnvHist::integrate2D_to_1D( const TH2D* object, const char *axis1, const double mult, const char *option ){

  int xbins=0, ybins=0;
  double* x_bins(0);
  double* y_bins(0);
  char axis2[1];

  if( 0 == strcmp(axis1,"1") )
  {
    sprintf( axis2, "2" );
    xbins = object->GetNbinsX();
    x_bins = new double[xbins+1];
    for( int i = 0; i <= xbins; i++ )
      x_bins[i] = object->GetXaxis()->GetBinLowEdge(i+1);
    ybins = object->GetNbinsY();
    y_bins = new double[ybins+1];
    for( int i = 0; i <= ybins; i++ )
      y_bins[i] = object->GetYaxis()->GetBinLowEdge(i+1);
  }
  else
  {
    sprintf( axis2, "1" );
    xbins = object->GetNbinsY();
    x_bins = new double[xbins+1];
    for( int i = 0; i <= xbins; i++ )
      x_bins[i] = object->GetYaxis()->GetBinLowEdge(i+1);
    ybins = object->GetNbinsX();
    y_bins = new double[ybins+1];
    for( int i = 0; i <= ybins; i++ )
      y_bins[i] = object->GetXaxis()->GetBinLowEdge(i+1);
  }

  MnvH1D *result = new MnvH1D( "result", "result", xbins, x_bins );

  double sum = 0;
  int count = 0;
  double errorSum2=0;

  for( int i = 1; i <= xbins; i++ )
  {
    sum = 0;
    count = 0;
    errorSum2 = 0;
    for( int j = 1; j <= ybins; j++ )
    {
      if( 0 == strcmp(axis1,"1") ){
        if( object->GetBinContent(i,j) != -99 )
        {
          if( strcmp( option, "omega" ) == 0 )
          {
            sum += object->GetBinContent(i,j) * 
              (cos(object->GetYaxis()->GetBinLowEdge(j)) - cos(object->GetYaxis()->GetBinUpEdge(j))) * mult;
            count++;
            errorSum2 += pow(object->GetBinError(i,j),2) *
              (cos(object->GetYaxis()->GetBinLowEdge(j)) - cos(object->GetYaxis()->GetBinUpEdge(j))) * mult;
          }
          else
          {
            sum += object->GetBinContent(i,j) * object->GetYaxis()->GetBinWidth(j) * mult;
            count++;
            errorSum2 += pow(object->GetBinError(i,j),2);
          }
        }
      }else{
        if( object->GetBinContent(j,i) != -99 )
        {
          sum += object->GetBinContent(j,i) * object->GetXaxis()->GetBinWidth(j) * mult;
          count++;
          errorSum2 += pow(object->GetBinError(j,i),2) * object->GetXaxis()->GetBinWidth(j) * mult;;
        }
      }
    }
    if( count != 0 )
    {
      result->SetBinContent(i,sum);
      double error = sqrt( errorSum2 );
      result->SetBinError(i,error );
      //cout << i << " = " << sum << " / " << count << " = " 
      //<< result->GetBinContent(i) << " +- " << error << endl;
    }
    else
    {
      result->SetBinContent(i,-99);
      result->SetBinError(i,1);
    }
  }

  return result;

}


//=============================================================================
// integrate2D()
//=============================================================================
double MnvHist::integrate2D( const TH2D* object, const double mult, const char *option ){

  double sum = 0;

  for( int i = 1; i <= object->GetNbinsX(); i++ )
  {
    for( int j = 1; j <= object->GetNbinsY(); j++ )
    {
      if( strcmp( option, "omega" ) == 0 )
      {
        sum += object->GetBinContent(i,j) * 
          object->GetXaxis()->GetBinWidth(i) *
          (cos(object->GetYaxis()->GetBinLowEdge(j)) - cos(object->GetYaxis()->GetBinUpEdge(j))) * mult;
      }
      else
      {
        sum += object->GetBinContent(i,j) * 
          object->GetXaxis()->GetBinWidth(i) * 
          object->GetYaxis()->GetBinWidth(j) * mult;
      }
    }
  }

  return sum;

}

//==============================================================
// print formatted matrix
//==============================================================
void MnvHist::printHisto( TH2D *hist, string name ){

  double sumX;

  cout << endl << name << ":" << endl;

  cout << setw(10) << "Y  X->  " << flush;
  for( int i = 0; i <= hist->GetNbinsX()+1; i++ ) {
    cout << setprecision(3) << setw(8) << hist->GetXaxis()->GetBinLowEdge(i) << " " << flush;
  }
  
  cout << endl;
  cout << "----------" << flush;
  for( int i = 0; i <= hist->GetNbinsX()+1; i++ ) {
    cout << "---------" << flush;
  }
  cout << endl;
  for( int j = 0; j <= hist->GetNbinsY()+1; j++ ) {
    cout << setprecision(3) << setw(8) << hist->GetYaxis()->GetBinLowEdge(j) << " | " << flush;
    sumX = 0;
    for(int i = 0; i <= hist->GetNbinsX()+1; i++) {
      cout << setprecision(3) << setw(8) << hist->GetBinContent(i,j) << " " << flush;
      sumX += hist->GetBinContent(i,j);
    }
    cout << " = " << setprecision(3) << sumX << flush;
    cout << endl;
  }
  cout << endl;
  
}

//==============================================================
// print formatted matrix
//==============================================================
void MnvHist::printMatrix( TMatrix matrix, string name ){

  double sumX;

  cout << endl << name << ":" << endl;

  for( int i = 0; i < matrix.GetNrows(); i++ ) {
    sumX = 0;
    for(int j = 0; j < matrix.GetNcols(); j++) {
      if( fabs(matrix[i][j]) > 1e-6 ){
	cout << setprecision(3) << setw(9) << matrix[i][j] << " " << flush;
	sumX += matrix[i][j];
      }
      else
	cout << setprecision(3) << setw(9) << 0 << " " << flush;
    }
    cout << " = " << setprecision(3) << sumX << flush;
    cout << endl;
  }
  cout << endl;
  
}


int MnvHist::AddInQuadrature( TH1* a, const TH1* b){

  //check axes
  if( a->GetNbinsX() != b->GetNbinsX() ) {
    Error("MnvHist::AddInQuadrature", "Histogram axes not compatible.  Doing nothing." );
    return 1;
  }

  //add all bins in quadrature (including under/overflow)
  //also add errors accordingly
  int firstBin = 0;
  int lastBin  = a->GetNbinsX()+1;
  for( int iBin = firstBin; iBin <= lastBin; ++iBin )
  {
    const double aVal = a->GetBinContent(iBin);
    const double bVal = b->GetBinContent(iBin);
    const double aErr = a->GetBinError(iBin);
    const double bErr = b->GetBinError(iBin);

    const double val = sqrt( aVal*aVal + bVal*bVal );
    const double err = ( 1E-8 < fabs(val) ) ? sqrt( pow(aVal*aErr,2) + pow(bVal*bErr,2) ) / val : 0.;

    a->SetBinContent(iBin, val);
    a->SetBinError(iBin, err);
  }
  return 0;
}

TMatrixD MnvHist::GetErrorsAsMatrix( const TH1D *h )
{
  const int highBin = h->GetNbinsX() + 1;
  const int lowBin = 0;
  TMatrixD covmx(highBin+1,highBin+1);

  // stat error
  for (int iBin=lowBin; iBin<=highBin; ++iBin )
    covmx[iBin][iBin] = TMath::Power( h->GetBinError(iBin), 2 );

  return covmx;
}

#endif

//#############################################################################
//
// END 
//
//#############################################################################
