#include <iostream>
#include <numeric>

#include <TCanvas.h>
#include <TRandom1.h>
#include <TFile.h>

#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/FluxReweighter.h>
#include <PlotUtils/FluxReweighterWithWiggleFit.h>


using namespace PlotUtils;

void print(MnvH1D* h)
{
    return; // not doing anything now
    
    printf("\n");
    TH1D herr = h->GetTotalError(false,true,false);
    int N = h->GetNbinsX();
    for (int bin = 1; bin <= std::min(100,N); ++bin) {
        printf("\t %5d %10.4e %10.4f\n", bin,
               h->GetBinContent(bin),
               herr.GetBinContent(bin));
    }
}

void test_flux_reweighter_wiggle()
{

    
    FluxReweighter* frw = new FluxReweighterWithWiggleFit(14, false,
                                             FluxReweighter::minervame1A,
                                             FluxReweighter::gen2thin,
                                             FluxReweighter::g4numiv6);


    MnvH1D* gen_flux = frw->GetFluxGenerated(14);
    MnvH1D* reweighted_flux = frw->GetFluxReweighted(14);
    MnvH1D* reweight_function = frw->GetFluxWiggleMnvH1D();
    

    print(gen_flux);

    print(reweighted_flux);

    print(reweight_function);

    bool not_stat = false;
    bool as_frac = true;
    bool not_cov_area_norm = false;


    TFile f("output.root", "recreate");

    gen_flux->GetCVHistoWithError().Write("flux_gen");
    reweighted_flux->GetCVHistoWithError().Write("flux_reweighted");
    reweight_function->GetCVHistoWithError().Write("reweight_function");
    reweighted_flux->GetTotalError(not_stat, as_frac, not_cov_area_norm).Write("flux_reweighted_frac_err");
    reweight_function->GetTotalError(not_stat, as_frac, not_cov_area_norm).Write("reweight_function_frac_err");
    
    f.Close();
}
