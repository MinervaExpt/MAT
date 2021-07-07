#include <PlotUtils/MnvH1D.h>
#include <TString.h>
#include <iostream>
#include <TH1.h>
#include<TCanvas.h>

using PlotUtils::MnvH1D;
using namespace std;

void ValidationforTrung() {

  TFile *f = new TFile("output.root");
  f->ls();

  h1 = (PlotUtils::MnvH1D*)f->Get("flux_E_cvweighted");
  vector<string> vertNames = h1->GetVertErrorBandNames();
  cout << "vertNames" << vertNames[] << endl;
 
  TH1F *h1_compare = new TH1F("h1_compare","h1_compare",100,0,100);
  
  for(vector<string>::iterator name = vertNames.begin();
      name != vertNames.end(); ++name){
    string errName = *name;
    cout <<"the name is "<<errName;
    int n_universe = h1->GetVertErrorBand(errName)->GetNHists();
    cout << n_universe << endl;
    // for (int i =0; i < h1->GetNbinsX(); ++i){
      vector<TH1D*> h_univ = h1->GetVertErrorBand(errName)->GetHists();  
      int size = h_univ.size();
      TH1D* mean_univ = (TH1D*)h_univ[0]->Clone("mean_universe");
      mean_univ->Reset();
      for(int j = 0;j<size;j++)mean_univ->Add(h_univ[j]);
      mean_univ->Scale(1.0/size);
      for(int i = 0;i<100;i++){
	//h_univ[]->GetMean();
	// h1->GetBinContent(i);
	//cout << "the mean of the universe in 1 bin is " <<  i<< mean_univ->GetBinContent(i) << endl;
	 mean_cont = mean_univ->GetBinContent(i);
	if(mean_cont==0)mean_cont=1.0;
	cout << "the central value is" <<  i << " "<< h1->GetCVHistoWithStatError().GetBinContent(i)/mean_cont << endl;
	double value = h1->GetCVHistoWithStatError().GetBinContent(i)/mean_cont;
	h1_compare->Fill(i,value);
    }

  }
  
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->cd();
  h1_compare->Draw();
  h1_compare->GetYaxis()->SetRangeUser(0.98,1.03);
  c1->Update();
}
  
  	








