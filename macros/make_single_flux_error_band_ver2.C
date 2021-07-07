#include <iostream>
#include <string>

#include <TStyle.h>
#include <TFile.h>
#include <TColor.h>
#include <TRandom1.h>
#include <TCanvas.h>
#include <TString.h>
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TDecompChol.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TSystem.h"
#include <PlotUtils/MnvH1D.h>


using PlotUtils::MnvH1D;
using PlotUtils::MnvVertErrorBand;

void make_single_flux_error_band_ver2(const std::string& filename,
                                 const std::string& name="sample_Lownu_MC_Enu")
{
    TMatrixD GetTransposeChMatrix(TMatrixD hcov);
    const int n_universe = 1000;
    
    TFile* file = new TFile(filename.c_str(), "read");
    MnvH1D* hist = static_cast<MnvH1D*>(file->Get(name.c_str()));
    assert(hist);

    if (hist->HasVertErrorBand("AKGY Model")) hist->PopVertErrorBand("AKGY Model");
    TMatrixD chmatrix = hist->GetTotalErrorMatrix();    
    //TMatrixD norm_chmatrix=chmatrix/cov_sum;
    TMatrixD transposed_chmatrix = GetTransposeChMatrix(chmatrix);
    
    transposed_chmatrix.Print();    

    
    TH1D* result = (TH1D*) hist->GetCVHistoWithError().Clone("result");
    assert(result);
/*
    gStyle->SetErrorX(0.5);
    gStyle->SetOptStat(0000);
    result->SetFillColor(TColor::GetColor("#ff82ab"));
    result->SetMarkerStyle(kDot);
    result->DrawCopy("E2");
    result->Draw("SAME AXIS");
    result->SetFillColor(0);
    result->Draw("SAME HIST");
*/
    for (int bin = 1; bin < result->GetNbinsX(); ++bin) {
        printf("\t %5d    %4.1f - %4.1f %10.4f %10.4f\n",
               bin,
               result->GetBinLowEdge(bin),
               result->GetBinLowEdge(bin) + result->GetBinWidth(bin),
               result->GetBinContent(bin),
               result->GetBinError(bin)/result->GetBinContent(bin));
    }
    std::vector<TRandom1*> generators;
    TRandom1 master_generator(1234567890);
    for (int bin = 1; bin <= result->GetNbinsX(); ++bin) {
        int seed = master_generator.Uniform(1e7,1e9);
        printf("\t bin %2d seed %10d\n", bin, seed);
        generators.push_back(new TRandom1(seed));
    }


    std::vector<TH1D*> universe_histograms;
    for (int i = 0; i < n_universe; ++i) {
        TH1D* h_univ = (TH1D*) result->Clone(Form("flux_universe_%3d", i));
        h_univ->Reset();
        universe_histograms.push_back(h_univ);
    }
    std::vector<TVectorD>rand_vect;
    TVectorD cv_nom(result->GetNbinsX()-1);
    for(int i=1;i<result->GetNbinsX();i++)cv_nom[i-1]=result->GetBinContent(i);
    
    for(int i=0;i<n_universe;++i){
       TVectorD z(transposed_chmatrix.GetNrows());
       for(int j=0;j<transposed_chmatrix.GetNrows();++j){
       TRandom1 *generator = generators[j];
       z(j)=generator->Gaus();
       }
       TVectorD x(transposed_chmatrix*z);
      // std::cout<<x.GetNrows()<<" "<<cv_nom.GetNrows()<<" "<<z.GetNrows()<<std::endl;
       TVectorD newpar = cv_nom+x;
       rand_vect.push_back(newpar);    
    }
   // std::cout<<"number of bins "<<result->GetNbinsX()<<" vector size "<< rand_vect[0].GetNrows()<<std::endl;
    for(int i =0;i<rand_vect.size();i++){
        TVectorD temp_vect= rand_vect[i];
        TH1D *h_univ = universe_histograms.at(i);
	for(int j=0;j<result->GetNbinsX();j++){

	  double temp_cont = result->GetBinContent(j);
	  if(j>0)temp_cont = temp_vect[j-1];
	  h_univ->SetBinContent(j,temp_cont);
	
	}
    
    }
    /*for (int bin = 1; bin <= result->GetNbinsX(); ++bin) {
        double cv = result->GetBinContent(bin);
        double sigma = result->GetBinError(bin);

        TRandom1* generator = generators[bin-1];
        assert(generator);
        for (int i = 0; i < n_universe; ++i) {
            double random_normal = generator->Gaus();
            double random_shift = sigma * random_normal;

            TH1D* h_univ = universe_histograms.at(i);
            h_univ->SetBinContent(bin, cv + random_shift);
            
        }
    }*/

    
    
#ifdef old_universe_code
    
    for (int i = 0; i < n_universe; ++i) { 
        TH1D* h_univ = universe_histograms[i];
        double random_normal = normal_vec[i];
        for (int bin = 1; bin <= result->GetNbinsX(); ++bin) {
            double error = result->GetBinError(bin);
            double random_shift = error * random_normal;
            double cv = result->GetBinContent(bin);
            h_univ->SetBinContent(bin, cv + random_shift);
        }
    }
    
#endif
    
    MnvH1D* output = new MnvH1D(*hist);
    output->SetName("beamfit_reweight_function");
    output->ClearAllErrorBands();
    output->AddVertErrorBand("Flux", universe_histograms);

    TH1D* reproduced = (TH1D*) output->GetCVHistoWithError().Clone("reproduced");
    assert(reproduced);

    std::cout << "Trying to reproduce the total flux error" << std::endl;
    for (int bin = 1; bin < result->GetNbinsX(); ++bin) {
        printf("\t %5d    %4.1f - %4.1f %10.4f %10.4f\n",
               bin,
               reproduced->GetBinLowEdge(bin),
               reproduced->GetBinLowEdge(bin) + reproduced->GetBinWidth(bin),
               reproduced->GetBinContent(bin),
               reproduced->GetBinError(bin)/reproduced->GetBinContent(bin));
    }
/*
    TCanvas* c2 = new TCanvas("c2", "reproduced");

    reproduced->SetFillColor(TColor::GetColor("#ff82ab"));
    reproduced->SetMarkerStyle(kDot);
    reproduced->DrawCopy("E2");
    reproduced->Draw("SAME AXIS");
    reproduced->SetFillColor(0);
    reproduced->Draw("SAME HIST");

*/
    
    TFile* new_file = new TFile("new_file.root", "recreate");
    output->Write();

    new_file->Close();
    
}

TMatrixD GetTransposeChMatrix(TMatrixD old_covmatrix){
//do all the BS here....dont destroy the main function.....
    TH2D *hcov = new TH2D(old_covmatrix);
    //I really dont know why this works and not the straightforward method....
    const int nbinsx = hcov->GetNbinsX();
    const int nbinsy = hcov->GetNbinsY();
    int newsize= nbinsx;
    for(int i=nbinsx-1;i>=0;--i){
    	double element = hcov->GetBinContent(i,i);
	if(element!=0.0){
	   newsize=i-1;
	   break;
	 }
    }
    std::cout<<"New size is "<<newsize<<std::endl;
    TMatrixD *cov_matrix = new TMatrixD(nbinsx+2,nbinsy+2,hcov->GetArray(),"F");
    //now we need to determine the diagonal element where they are 0
    TMatrixD covmatrix(newsize,newsize);
    for(int i=0;i!=newsize;++i){
      for(int j=0;j!=newsize;++j){
        covmatrix(i,j) = hcov->GetBinContent(i+2,j+2);      
        }
      }
    TMatrixDEigen coveigen(covmatrix);
    
    const double tol=0.0;
    TDecompChol decomp(covmatrix,tol);//0.0 is the tolerance...
    std::cout<<decomp.Decompose()<<std::endl;
    TMatrixD U = decomp.GetU();
    coveigen.GetEigenValuesRe().Print();
   // std::cout<<decomp.Decompose()<<" Decomposition status "<<std::endl;
   if(!decomp.Decompose()){
     std::cout<<"Decomposition failed....*&^&*$@#$%#$#@#@##$$%%$%$$......."<<std::endl;
     coveigen.GetEigenValuesRe().Print();
   }
  // std::cout<<"U matrix is"<<std::endl;
  // U.Print();
   TMatrixD UT(TMatrixD::kTransposed,U);
   return UT;
}
