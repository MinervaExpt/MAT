#include <iostream>
#include <string>

#include <TStyle.h>
#include <TFile.h>
#include <TColor.h>
#include <TRandom1.h>
#include <TCanvas.h>
#include <TString.h>

#include <PlotUtils/MnvH1D.h>


using PlotUtils::MnvH1D;
using PlotUtils::MnvVertErrorBand;

void make_single_flux_error_band(const std::string& filename,
                                 const std::string& name="sample_Lownu_MC_Enu")
{

    const int n_universe = 1000;
    
    TFile* file = new TFile(filename.c_str(), "read");
    MnvH1D* hist = static_cast<MnvH1D*>(file->Get(name.c_str()));
    assert(hist);

    if (hist->HasVertErrorBand("AKGY Model")) hist->PopVertErrorBand("AKGY Model");
    
    TH1D* result = (TH1D*) hist->GetCVHistoWithError().Clone("result");
    assert(result);

    gStyle->SetErrorX(0.5);
    gStyle->SetOptStat(0000);
    result->SetFillColor(TColor::GetColor("#ff82ab"));
    result->SetMarkerStyle(kDot);
    result->DrawCopy("E2");
    result->Draw("SAME AXIS");
    result->SetFillColor(0);
    result->Draw("SAME HIST");

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

    for (int bin = 1; bin <= result->GetNbinsX(); ++bin) {
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
    }

    
    
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

    TCanvas* c2 = new TCanvas("c2", "reproduced");

    reproduced->SetFillColor(TColor::GetColor("#ff82ab"));
    reproduced->SetMarkerStyle(kDot);
    reproduced->DrawCopy("E2");
    reproduced->Draw("SAME AXIS");
    reproduced->SetFillColor(0);
    reproduced->Draw("SAME HIST");


    
    TFile* new_file = new TFile("new_file.root", "recreate");
    output->Write();

    new_file->Close();
    
}
