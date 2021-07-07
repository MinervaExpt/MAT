#include <iostream>
#include <string>
#include <vector>
#include "TRandom3.h"
#include "TCanvas.h"
#include "TFile.h"

#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvApplication.h"
#include "PlotUtils/MnvPlotter.h"

using namespace std;
using namespace PlotUtils;

int tryToWrite()
{
  cout << "Create a new MnvH1D.  Also create a partner." << endl;
  MnvH1D* a = new MnvH1D("a","This baby is going to disk", 10, 0., 10. );
  MnvH1D* b = new MnvH1D("a","This baby is also going to disk", 10, 0., 10. );

  cout << "Add a vertical error A wth many universes" << endl;
  a->AddVertErrorBand("ErrA" );
  b->AddVertErrorBand("ErrA" );

  cout << "Add a lateral error B with 2 universes (+/-1 sigma) to use max spread" << endl;
  a->AddLatErrorBand("ErrB" );
  b->AddLatErrorBand("ErrB" );

  cout << "Add a lateral sys error C to use max spread" << endl;
  a->AddLatErrorBand("ErrC", 5); //5 universes
  b->AddLatErrorBand("ErrC", 5); //5 universes

  //set ErrC to get its errors from the spread method manually, although spread is the default for less than 10 universes
  a->GetLatErrorBand("ErrC")->SetUseSpreadError(true);
  cout << "      ErrC has " << a->GetLatErrorBand("ErrC")->GetNHists() << " hists and UseSpreadError = " << a->GetLatErrorBand("ErrC")->GetUseSpreadError() << endl;

  // get the vector of error band names like this
  vector<string> errNames = a->GetErrorBandNames();
  cout << "MnvErrorBand names: " << endl;
  for( vector<string>::iterator it = errNames.begin(); it != errNames.end(); ++it )
    cout << "         name = " << *it << endl;

  //fill 1000 universes of nonsense weights
  //these should be branches in the AnaTuple.  e.g. mc_wgt_Flux_Tertiary
  TRandom3 r;
  vector<double> weightsA;
  for( int i = 0; i < 1000; ++i )
  {
    weightsA.push_back( r.Gaus(1.,.05) ); //vertical shifts use weights centered around a CV value of 1 (same weight)
  }

  //B has asymmetrical errors.  15% up and 10% down
  const double bWeightSigmaUp   =  .15;
  const double bWeightSigmaDown = -.10;

  //only do 5 weights for C since we are using max spread.  pick something obvious here...
  vector<double> weightsC;
  for( int i = 0; i < 5; ++i )
    weightsC.push_back( -.1 + i*0.05 );

  //the lat vectors store the +/- shift from the fill value and may depend on the fill value
  vector<double> weightsLatB, weightsLatC;
  weightsLatB.reserve( 2 );
  weightsLatC.reserve(weightsC.size() );

  vector<double> weightsLatB2, weightsLatC2;
  weightsLatB2.reserve( 2 );
  weightsLatC2.reserve(weightsC.size() );

  double cvweight = 1.;
  for( int i = 0; i < 100000; ++i )
  {
    double val = r.Gaus(5,2);
    while( val < 0 )
      val = r.Gaus(5,2);

    double val2 = r.Gaus(5,2);
    while( val2 < 0 )
      val2 = r.Gaus(5,2);

    //fill the +/-1 sigma shifts
    weightsLatB[0] = val * bWeightSigmaDown;
    weightsLatB[1] = val * bWeightSigmaUp;

    weightsLatB2[0] = val2 * bWeightSigmaDown;
    weightsLatB2[1] = val2 * bWeightSigmaUp;


    for( unsigned int j = 0; j != weightsC.size(); ++j )
    {
      weightsLatC[j] = val * weightsC[j];
      weightsLatC2[j] = val2 * weightsC[j];
    }


    //fill MnvH1D like you would a normal histogram
    a->Fill(val);
    b->Fill(val2);

    //fill the error bands of the MnvH1D by name
    a->FillVertErrorBand( "ErrA", val, weightsA, cvweight );
    a->FillLatErrorBand( "ErrB", val, weightsLatB, cvweight );
    a->FillLatErrorBand( "ErrC", val, weightsLatC, cvweight );


    b->FillVertErrorBand( "ErrA", val2, weightsA, cvweight );
    b->FillLatErrorBand( "ErrB", val2, weightsLatB2, cvweight );
    b->FillLatErrorBand( "ErrC", val2, weightsLatC2, cvweight );
  }


  TH1D totErr = a->GetTotalError(true); //total error.  sys + stat
  TH1D totSysErr = a->GetTotalError(false); //total sys error.  no stat.
  for( int i = 0; i < a->GetNbinsX(); ++i )
  {
    cout << "bin " << i << endl;
    cout << "     val = " << a->GetBinContent(i) << endl;
    cout << "     totErr = " << totErr.GetBinContent(i) << endl;
    cout << "     totSysErr = " << totSysErr.GetBinContent(i) << endl;
  }

  TCanvas c;
  c.Divide(2,2);
  c.cd(1);
  a->Draw();
  c.cd(2);
  a->GetVertErrorBand("ErrA")->GetErrorBand(true).DrawCopy(); //get the error band from a spread in A's universes as a fraction of CV
  c.cd(3);
  a->GetLatErrorBand("ErrB")->GetErrorBand(true).DrawCopy();
  c.cd(4);
  a->GetLatErrorBand("ErrC")->GetErrorBand(true).DrawCopy();
  c.Print("canDuringWrite.png","png");

  //you can use the namespace MnvPlot to get MinervaPlottingUtils
  //MultiPrint prints a canvas using this name and a comma separated list of types.
  MnvPlotter mnvPlotter;
  mnvPlotter.MultiPrint( &c, "canDuringWrite", "eps,gif" );


  //draw all the universes of ErrC
  TCanvas cAll( "canDrawAllUniverses", "canDrawAllUniverses");
  a->GetLatErrorBand("ErrC")->DrawAll("HIST", true);
  mnvPlotter.MultiPrint( &cAll );

  //get the final result (CV + errror) and draw
  TCanvas cCVErr("canDrawCVWithError", "canDrawCVWithError");
  TH1D cvWithErr = a->GetCVHistoWithError();
  cvWithErr.Draw();
  mnvPlotter.MultiPrint( &cCVErr );

  //open a TFile and write the MnvH1D and its MnvSysErrorBands to disk
  TFile f("wroteSomething.root","recreate");
  a->Write("MyHist"); //the name in the file is MyHist
  f.Close();

  TFile f2("wroteSomethingElse.root","recreate");
  b->Write("MyHist");
  f2.Close();

  return 0;
}

int main() {
  PlotUtils::Initialize();

  return tryToWrite();

}
