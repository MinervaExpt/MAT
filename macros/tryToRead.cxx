#include <iostream>
#include <string>
#include <vector>
#include "TRandom3.h"
#include "TCanvas.h"
#include "TFile.h"

#include "PlotUtils/MnvApplication.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvPlotter.h"

using namespace std;
using namespace PlotUtils;

int tryToRead()
{
  MnvPlotter mnvPlotter;

  //open the file and get the MnvH1D
  TFile *f = new TFile("wroteSomething.root");
  MnvH1D* a = (MnvH1D*)f->Get("MyHist");

  TFile *f2 = new TFile("wroteSomethingElse.root");
  MnvH1D* b = (MnvH1D*)f2->Get("MyHist");


  //Lets make a clone and divide
  MnvH1D* divHist = (MnvH1D*)b->Clone("MyHists_Divided");
  vector<string> errNames = divHist->GetErrorBandNames();
  for( vector<string>::iterator i = errNames.begin(); i != errNames.end(); ++i )
    cout << "   ErrorBand : " << *i << endl;

  divHist->Divide(a, b); //can also scale Divide(a,b,c1,c2)

  TCanvas cD("cD","Look a division",1280,900);
  cD.Divide(3);
  cD.cd(1);
  a->GetCVHistoWithError().Draw();
  cD.cd(2);
  b->GetCVHistoWithError().Draw();
  cD.cd(3);
  divHist->GetCVHistoWithError().Draw();
  mnvPlotter.MultiPrint( &cD, "MnvH1D_Divide" ); //specify a print name optionally

  //All of the MnvSysErrorBands are still in the file as well
  TH1D totErr = a->GetTotalError(true);
  TH1D totSysErr = a->GetTotalError(false);

  TH1D totErr2 = b->GetTotalError(true);
  TH1D totSysErr2 = b->GetTotalError(false);

  for( int i = 0; i < a->GetNbinsX(); ++i )
  {
    cout << "bin " << i << endl;
    double aval =  a->GetBinContent(i);
    double bval =  b->GetBinContent(i);
    double divval =  divHist->GetBinContent(i);
    double abval = bval > 0 ? aval / bval : 0.;
    cout << "      a = " << aval << ", b = " << bval << ", divided = " << divval << ", a/b = " << abval << endl;
  }


  TCanvas c( "canDuringRead", "canDuringRead" );
  c.Divide(2,2);
  c.cd(1);
  a->Draw();
  c.cd(2);
  a->GetVertErrorBand("ErrA")->GetErrorBand().DrawCopy();
  c.cd(3);
  a->GetLatErrorBand("ErrB")->GetErrorBand().DrawCopy();
  c.cd(4);
  a->GetLatErrorBand("ErrC")->GetErrorBand().DrawCopy();
  mnvPlotter.MultiPrint( &c );

  //draw all the universes of ErrC
  TCanvas cAll("canDrawAllUniverses", "canDrawAllUniverses");
  a->GetLatErrorBand("ErrC")->DrawAll("HIST", true);
  mnvPlotter.MultiPrint( &cAll );


  return 0;
}

int main() {
  PlotUtils::Initialize();
  return tryToRead();

}
