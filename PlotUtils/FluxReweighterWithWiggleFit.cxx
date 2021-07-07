#include <iostream>
#include <cstdlib>
#undef NDEBUG
#include <cassert>

#include <TSystem.h>
#include <TString.h>

#include "MnvH1D.h"
#include "MnvVertErrorBand.h"
#include "FluxReweighter.h"
#include "FluxReweighterWithWiggleFit.h"

using PlotUtils::FluxReweighter;
using PlotUtils::MnvH1D;
using PlotUtils::MnvVertErrorBand;

PlotUtils::FluxReweighterWithWiggleFit::FluxReweighterWithWiggleFit(int nuPDG,
                                                                    bool applyNuEConstraint,
                                                                    enum EPlaylist playlist,
                                                                    enum EFluxVersion fluxVersion,
                                                                    enum EG4NumiVersion g4NumiVersion)
    : FluxReweighter(nuPDG, applyNuEConstraint, playlist, fluxVersion, g4NumiVersion),
      m_wiggleFitReweightHist(0)
{

    assert((nuPDG == 12 || nuPDG == 14 || nuPDG == -12 || nuPDG == -14) && "Neutrino types not supported");
    
    m_wiggleFitReweightHist = GetFluxWiggleMnvH1D();

    MnvVertErrorBand* original_flux_errorband = m_fluxReweightNu->PopVertErrorBand("Flux");
    
    m_fluxReweightNu->ClearAllErrorBands();
    m_fluxReweightNu->AddMissingErrorBandsAndFillWithCV(*m_wiggleFitReweightHist);

    m_fluxReweightNu = MultiplyHists(m_fluxReweightNu, m_wiggleFitReweightHist);

    MnvVertErrorBand* modified_flux_errorband = m_fluxReweightNu->GetVertErrorBand("Flux");
    unsigned int n_universe = modified_flux_errorband->GetNHists();
    assert(original_flux_errorband->GetNHists() == n_universe);

    const double max_enu = m_wiggleFitReweightHist->GetXaxis()->GetXmax();    
    for (unsigned int univ = 0; univ < n_universe; ++univ) {
        TH1D* h_univ_original = original_flux_errorband->GetHists()[univ];
        TH1D* h_univ_modified = modified_flux_errorband->GetHists()[univ];

        for (int bin = 1; bin <= h_univ_original->GetNbinsX(); ++bin) {
            double content_original = h_univ_original->GetBinContent(bin);
            double Enu = h_univ_modified->GetBinCenter(bin);

            if (Enu > max_enu) h_univ_modified->SetBinContent(bin, content_original);

        }

    }
    
}

PlotUtils::FluxReweighterWithWiggleFit::FluxReweighterWithWiggleFit(int nuPDG,
                                                                    bool applyNuEConstraint,
                                                                    std::string playlist_str,
                                                                    enum EFluxVersion fluxVersion,
                                                                    enum EG4NumiVersion g4NumiVersion)
: FluxReweighter(nuPDG, applyNuEConstraint, playlist_str, fluxVersion, g4NumiVersion),
m_wiggleFitReweightHist(0)
{
    
    assert((nuPDG == 12 || nuPDG == 14 || nuPDG == -12 || nuPDG == -14) && "Neutrino types not supported");
    
    m_wiggleFitReweightHist = GetFluxWiggleMnvH1D();
    
    MnvVertErrorBand* original_flux_errorband = m_fluxReweightNu->PopVertErrorBand("Flux");
    
    m_fluxReweightNu->ClearAllErrorBands();
    m_fluxReweightNu->AddMissingErrorBandsAndFillWithCV(*m_wiggleFitReweightHist);
    
    m_fluxReweightNu = MultiplyHists(m_fluxReweightNu, m_wiggleFitReweightHist);
    
    MnvVertErrorBand* modified_flux_errorband = m_fluxReweightNu->GetVertErrorBand("Flux");
    unsigned int n_universe = modified_flux_errorband->GetNHists();
    assert(original_flux_errorband->GetNHists() == n_universe);
    
    const double max_enu = m_wiggleFitReweightHist->GetXaxis()->GetXmax();
    for (unsigned int univ = 0; univ < n_universe; ++univ) {
        TH1D* h_univ_original = original_flux_errorband->GetHists()[univ];
        TH1D* h_univ_modified = modified_flux_errorband->GetHists()[univ];
        
        for (int bin = 1; bin <= h_univ_original->GetNbinsX(); ++bin) {
            double content_original = h_univ_original->GetBinContent(bin);
            double Enu = h_univ_modified->GetBinCenter(bin);
            
            if (Enu > max_enu) h_univ_modified->SetBinContent(bin, content_original);
            
        }
        
    }
    
}





PlotUtils::FluxReweighterWithWiggleFit::~FluxReweighterWithWiggleFit()
{
    delete m_wiggleFitReweightHist;
}


double PlotUtils::FluxReweighterWithWiggleFit::GetFluxCVWeight( double Enu, int nuPDG )
{
    return FluxReweighter::GetFluxCVWeight(Enu, nuPDG);
}

std::vector<double> PlotUtils::FluxReweighterWithWiggleFit::GetFluxErrorWeights( double Enu, int nuPDG )
{
    return FluxReweighter::GetFluxErrorWeights(Enu, nuPDG);
}

    
double PlotUtils::FluxReweighterWithWiggleFit::GetFluxErrorWeight( double Enu, int nuPDG, unsigned int universe )
{
    const double max_enu = m_wiggleFitReweightHist->GetXaxis()->GetXmax();

    double weight = 1.0;

    if (Enu < max_enu) weight = FluxReweighter::GetFluxErrorWeight(Enu, nuPDG, universe);

    return weight;
    
}


MnvH1D* PlotUtils::FluxReweighterWithWiggleFit::GetFluxWiggleMnvH1D()
{
    const char* plotutils=gSystem->Getenv("PLOTUTILSROOT");
    if(!plotutils || plotutils[0] == '\0'){
        std::cout << "$PLOTUTILSROOT is not set. Can't find flux histograms" << std::endl;
        std::exit(1);
    }
    
    const char* filename = Form("%s/data/dataMCWiggle/reweight_function_from_wiggle_fit.root", plotutils);

    return GetMnvH1D(filename, "beamfit_reweight_function"); // put amits format
}

MnvH1D* PlotUtils::FluxReweighterWithWiggleFit::MultiplyHists(MnvH1D* h1, MnvH1D* mnvh1d_weight)
{

    MnvH1D* result = new MnvH1D(*h1);
    
    const double max_enu = mnvh1d_weight->GetXaxis()->GetXmax();
        // central-value histogram
    for (int bin = 1; bin <= h1->GetNbinsX(); ++bin) {
        double content = h1->GetBinContent(bin);
        double Enu = h1->GetBinCenter(bin);
        double weight = 1.0;
        int weight_bin = mnvh1d_weight->FindBin(Enu);
        if (Enu < max_enu) weight = mnvh1d_weight->GetBinContent(weight_bin);
        result->SetBinContent(bin, content * weight);
    }

        // synchronize error band cv with MnvH1D cv
    MnvVertErrorBand* flux_error_band = result->GetVertErrorBand("Flux");
    for (int bin = 1; bin <= result->GetNbinsX(); ++bin) {
        double content = result->GetBinContent(bin);
        flux_error_band->SetBinContent(bin, content);
    }
    
        // universe histograms
    std::vector<std::string> vertNames = result->GetVertErrorBandNames();
    for (std::vector<std::string>::iterator name = vertNames.begin();
         name != vertNames.end(); ++name) {
        std::string errName = *name;
        const unsigned int n_universe = result->GetVertErrorBand(errName)->GetNHists();
        for (unsigned int i = 0; i < n_universe; ++i) {
            TH1D* h_univ        = result->GetVertErrorBand(errName)->GetHists().at(i);
            TH1D* h_univ_weight = mnvh1d_weight->GetVertErrorBand(errName)->GetHists().at(i);
               
            for (int bin = 1; bin <= h_univ->GetNbinsX(); ++bin) {
                double content = h_univ->GetBinContent(bin);
                double Enu = h_univ->GetBinCenter(bin);
                double weight = 1.0;
                int weight_bin = h_univ_weight->FindBin(Enu);
                if (Enu < max_enu) weight = h_univ_weight->GetBinContent(weight_bin);
                h_univ->SetBinContent(bin, content * weight);
            }
  
  
        }

    }

    return result;
    
}
