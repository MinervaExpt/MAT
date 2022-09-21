#include "PlotUtils/FluxReweighter.h"

#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvFluxConstraint.h"
#include "PlotUtils/MnvVertErrorBand.h"

#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TDecompChol.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1D.h"
#include "TObjString.h"
#include "TDirectory.h"

#include <iostream>
#include <cstdlib> // For exit()
#include <vector>
#include <string>
#include <algorithm> // tolower()

namespace {
  inline std::vector<TH1D*> GetVector(PlotUtils::MnvH1D* h)  {
    return std::vector<TH1D*>();
  }

    inline std::vector<TH2D*> GetVector(PlotUtils::MnvH2D* h)  {
      return std::vector<TH2D*>();
    }

  inline TH1D* NewHist(const TH1D& h) {
    return new TH1D (h);
  }

  inline TH2D* NewHist(const TH2D& h) {
    return new TH2D (h);
  }

  PlotUtils::FluxReweighter* _flux_reweighter(std::string plist,int nu_pdg,
                                              bool use_nuE_constraint,
                                              int n_flux_universes) {
    if (plist == "minerva1") {
      return new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv5,
                                         n_flux_universes  );
    }
    else if (plist == "minerva5") {
      return new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv5,
                                         n_flux_universes  );
    }
    else if (plist == "minerva7") {
      return new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv5,
                                         n_flux_universes  );
    }
    else if (plist == "minerva9") {
      return new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv5,
                                         n_flux_universes  );
    }
    else if (plist == "minerva13" || plist == "minerva13C" || plist == "minerva13c"|| plist == "minerva13E"|| plist == "minerva13e" ) {
      return new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         "minerva13",
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv5,
                                         n_flux_universes  );
    }
    else if (plist == "2p2h") {
      return new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         "minerva2p2h",
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv5,
                                         n_flux_universes  );
    }
    else if(plist == "minervame1a" || plist == "minervame1A"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );
    }
    else if(plist == "minervame1b" || plist == "minervame1B"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );

    }
    else if(plist == "minervame1c" || plist == "minervame1C"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );

    }
    else if(plist == "minervame1d" || plist == "minervame1D"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );

    }
    else if(plist == "minervame1e" || plist == "minervame1E"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );

    }
    else if(plist == "minervame1f" || plist == "minervame1F"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );

    }
    else if(plist == "minervame1g" || plist == "minervame1G"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );

    }
    else if(plist == "minervame1h" || plist == "minervame1H"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );

    }
    else if(plist == "minervame1l" || plist == "minervame1L"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );
    }
    else if(plist == "minervame1m" || plist == "minervame1M"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );
    }
    else if(plist == "minervame1n" || plist == "minervame1N"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );

    }
    else if(plist == "minervame1o" || plist == "minervame1O"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );
    }
    else if(plist == "minervame1p" || plist == "minervame1P"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );

    }
    else if (plist == "minervame1d1m1nweightedave") {
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );

    }
    else if(plist == "minervame5a" || plist == "minervame5A"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );

    }
    else if(plist == "minervame6a" || plist == "minervame6A"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );
    }
    else if(plist == "minervame6b" || plist == "minervame6B"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );
    }
    else if(plist == "minervame6c" || plist == "minervame6C"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );
    }
    else if(plist == "minervame6d" || plist == "minervame6D"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );
    }
    else if(plist == "minervame6e" || plist == "minervame6E"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );
    }
    else if(plist == "minervame6f" || plist == "minervame6F"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );
    }
    else if(plist == "minervame6g" || plist == "minervame6G"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );
    }
    else if(plist == "minervame6h" || plist == "minervame6H"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );
    }
    else if(plist == "minervame6i" || plist == "minervame6I"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );
    }
    else if(plist == "minervame6j" || plist == "minervame6J"){
      return
          new PlotUtils::FluxReweighter( nu_pdg, use_nuE_constraint,
                                         plist,
                                         PlotUtils::FluxReweighter::gen2thin,
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         n_flux_universes );
    }


    else {
      std::cerr << "ERROR flux_reweighter playlist string unrecognized "
                << plist << std::endl;
      std::exit(1);
    }
  }
}

namespace PlotUtils
{
  //============================================================================
  FluxReweighter::FluxReweighter(int nuPDG,
      bool applyNuEConstraint,
      enum EPlaylist playlist,
      enum EFluxVersion fluxVersion,
      enum EG4NumiVersion g4NumiVersion,
      int nUniverses/*= 200*/)
  {
    // cache the current directory
    TDirectory* original_dir = gDirectory->CurrentDirectory();

    // for the LE nUniverses can't be greater than 100
    if (LeOrMe(playlist)==1&&nUniverses>100) {
      nUniverses = 100;
    }

    // These member variables need to be set here because they're used downstream within the constructor
    m_applyNuEConstraint = applyNuEConstraint;
    m_nFluxUniverses = nUniverses;
    m_fluxErrorName = "Flux";

    m_fluxReweightNu    = GetFluxMnvH1D(abs(nuPDG),    playlist, fluxVersion, g4NumiVersion);
    m_fluxReweightNubar = GetFluxMnvH1D(-1*abs(nuPDG), playlist, fluxVersion, g4NumiVersion);

    d_Playlist = playlist;
    if (applyNuEConstraint) {
      std::cout << "Applying flux constraint to nu histogram" << std::endl;
      m_fluxReweightNu = Constrainer().ConstrainHisto<MnvH1D, MnvVertErrorBand>( m_fluxReweightNu );
      std::cout << "Applying flux constraint to nubar histogram" << std::endl;
      m_fluxReweightNubar = Constrainer().ConstrainHisto<MnvH1D, MnvVertErrorBand>( m_fluxReweightNubar );
    }

    // Downstream methods require access to the original set of universe weights, so we will propagate an un-truncated
    // copy of the flux histograms for local use, but not to be passed on to the user
  //  m_fluxReweightNu_Ref = m_fluxReweightNu->Clone();
  //  m_fluxReweightNubar_Ref = m_fluxReweightNubar->Clone();

    if(nuPDG>0)m_fluxReweightNu_Ref = m_fluxReweightNu->Clone();
    else{

    m_fluxReweightNu_Ref = m_fluxReweightNubar->Clone();

    }
    if (nUniverses < 1000) {
      TruncateNumberOfFluxUniverses( m_fluxReweightNu , nUniverses );
      TruncateNumberOfFluxUniverses( m_fluxReweightNubar , nUniverses );
    }

    m_fluxGenNu    = GetFluxMnvH1D(abs(nuPDG),    playlist, fluxVersion, g4NumiVersion, true);
    m_fluxGenNubar = GetFluxMnvH1D(-1*abs(nuPDG), playlist, fluxVersion, g4NumiVersion, true);
    m_reweightMELowNuDataMC = GetMELowNuMnvH1D();

    SetFluxSysMnvH1D(nuPDG,fluxVersion);

    m_MELowNuDataMCRatioSpline = GetSpline(m_reweightMELowNuDataMC);

    // return to the original directory
    original_dir->cd();
  }


  //============================================================================
  FluxReweighter::FluxReweighter(int nuPDG,
      bool applyNuEConstraint,
      std::string playlist_str,
      enum EFluxVersion fluxVersion,
      enum EG4NumiVersion g4NumiVersion,
      int nUniverses/*= 200*/)
  {
    // cache the current directory
    TDirectory* original_dir = gDirectory->CurrentDirectory();

    // These member variables need to be set here because they're used
    // downstream within the constructor
    m_applyNuEConstraint = applyNuEConstraint;
    m_nFluxUniverses = nUniverses;
    m_fluxErrorName = "Flux";

    FluxReweighter::EPlaylist playlist = GetPlaylistEnum(playlist_str);
    d_Playlist = playlist;

    // for the LE nUniverses can't be greater than 100
    if (LeOrMe(playlist)==1&&nUniverses>100) {
      nUniverses = 100;
    }

    m_fluxReweightNu    = GetFluxMnvH1D(abs(nuPDG),    playlist, fluxVersion, g4NumiVersion);
    m_fluxReweightNubar = GetFluxMnvH1D(-1*abs(nuPDG), playlist, fluxVersion, g4NumiVersion);

    if (applyNuEConstraint) {
      //std::cout << "Applying flux constraint to nu histogram" << std::endl;
      m_fluxReweightNu =
          Constrainer().ConstrainHisto<MnvH1D, MnvVertErrorBand>( m_fluxReweightNu );
      //std::cout << "Applying flux constraint to nubar histogram" << std::endl;
      m_fluxReweightNubar =
          Constrainer().ConstrainHisto<MnvH1D, MnvVertErrorBand>( m_fluxReweightNubar );
    }

    // Downstream methods require access to the original set of universe
    // weights, so we will propagate an un-truncated copy of the flux
    // histograms for local use, but not to be passed on to the user
    //m_fluxReweightNu_Ref = m_fluxReweightNu->Clone();
    //m_fluxReweightNubar_Ref = m_fluxReweightNubar->Clone();
     if(nuPDG>0)m_fluxReweightNu_Ref = m_fluxReweightNu->Clone();
    else{

    m_fluxReweightNu_Ref = m_fluxReweightNubar->Clone();

    }

    if (nUniverses < 1000) {
      TruncateNumberOfFluxUniverses( m_fluxReweightNu,    nUniverses );
      TruncateNumberOfFluxUniverses( m_fluxReweightNubar, nUniverses );
    }

    m_fluxGenNu    = GetFluxMnvH1D(abs(nuPDG),    playlist, fluxVersion, g4NumiVersion, true);
    m_fluxGenNubar = GetFluxMnvH1D(-1*abs(nuPDG), playlist, fluxVersion, g4NumiVersion, true);
    m_reweightMELowNuDataMC = GetMELowNuMnvH1D();

    SetFluxSysMnvH1D(nuPDG,fluxVersion);

    m_MELowNuDataMCRatioSpline = GetSpline(m_reweightMELowNuDataMC);

    // return to the original directory
    original_dir->cd();
  }


  //============================================================================
  FluxReweighter::FluxReweighter(MnvH1D* fluxGenNu,
                                 MnvH1D* fluxGenNubar,
                                 MnvH1D* fluxReweightNu,
                                 MnvH1D* fluxReweightNubar,
                                 TSpline3* MELowNuDataMCRatioSpline,
                                 bool applyNuEConstraint)
    : m_fluxGenNu(fluxGenNu),
      m_fluxGenNubar(fluxGenNubar),
      m_fluxReweightNu(fluxReweightNu),
      m_fluxReweightNubar(fluxReweightNubar),
      m_MELowNuDataMCRatioSpline(MELowNuDataMCRatioSpline)
  {
    // cache the current directory
    TDirectory* original_dir = gDirectory->CurrentDirectory();

    if (applyNuEConstraint) {
      std::cout << "Applying flux constraint to nu histogram" << std::endl;
      m_fluxReweightNu =
          Constrainer().ConstrainHisto<MnvH1D, MnvVertErrorBand>( m_fluxReweightNu );
      std::cout << "Applying flux constraint to nubar histogram" << std::endl;
      m_fluxReweightNubar =
          Constrainer().ConstrainHisto<MnvH1D, MnvVertErrorBand>( m_fluxReweightNubar );
    }
    m_applyNuEConstraint = applyNuEConstraint;
    m_nFluxUniverses = 1000;
    m_fluxErrorName = "Flux";

    // return to the original directory
    original_dir->cd();
  }


  //============================================================================
  FluxReweighter::~FluxReweighter()
  {
    delete m_fluxGenNu;
    delete m_fluxGenNubar;
    delete m_fluxReweightNu;
    delete m_fluxReweightNubar;
    delete m_fluxReweightNu_Ref;
    //delete m_fluxReweightNubar_Ref;
    delete m_MELowNuDataMCRatioSpline;
    delete m_reweightMELowNuDataMC;
  }


  //============================================================================
  MnvHistoConstrainer& FluxReweighter::Constrainer()
  {
    static MnvHistoConstrainer* _constrainer=0;
    if (!_constrainer) {
      _constrainer=new MnvHistoConstrainer;

      const char* mparam=gSystem->Getenv("MPARAMFILES");
      if (!mparam || !strlen(mparam)) {
        std::cout << "$MPARAMFILES is not set. Can't find universe weights for "
                  << "nu+e constraint" << std::endl;
        std::exit(1);
      }
      std::string filename(mparam);
      int iplaylist = LeOrMe(d_Playlist);
      // Timebeing 1 = LE; 2 = ME (ALL Playlist)
      if (iplaylist == 1) {
        filename+="/FluxConstraints/nu+e_eroica.txt";
      }
      //Using the combined fhc+rhc+imd constraint
      else if (iplaylist == 2 or iplaylist==3) {
	filename+="/FluxConstraints/sorted_NuEConstraint_FHC_RHC_IMD.txt";
      }

      //This is the nu+e method with fhc measurement
      /*      
      else if (iplaylist == 2 || iplaylist==3) {
	filename+="/FluxConstraints/nu+e_ME_spectrum.txt";
      }
      */
      
      //This is the IMD method
      /*
      else if (iplaylist == 2 or iplaylist==3) {
	std::cout << "I am running with the IMD constraint" << std::endl;
      	filename+="/FluxConstraints/IMDweights_unordered.txt";
      }
      */
      //This is the nu+e method with rhc measurement
      /*
      else if (iplaylist ==3){
        filename +="/FluxConstraints/sorted_nu+e_RHC_Final.txt";

      }
      */
      _constrainer->LoadConstraint( "nu+e", filename );

      _constrainer->SetSpectatorCorrectionStrategy("Flux_BeamFocus", MnvHistoConstrainer::PRESERVE_FRACTIONAL_ERR);
      _constrainer->SetSpectatorCorrectionStrategy("ppfx1_Total", MnvHistoConstrainer::PRESERVE_FRACTIONAL_ERR);
    }
    return *_constrainer;
  }


  //======================================================================
  MnvH1D* FluxReweighter::GetMnvH1D(TString filename, TString histname)
  {
    TFile f(filename);
    if (f.IsZombie()) {
      std::cout << "Can't open flux file with name " << filename << std::endl;
      std::exit(1);
    }
    MnvH1D* ret=(MnvH1D*)f.Get(histname);
    if (!ret) {
      std::cout << "Can't find flux histogram with name " << histname << " in file " << filename << std::endl;
      std::exit(1);
    }
    ret->SetDirectory(0);
    return ret;
  }

  /*
  //============================================================================
  MnvH1D* FluxReweighter::GetTargetFluxMnvH1D(int nuPDG,
                                              std::string tar_mat,
                                              std::string project_dir)
  {

    if( nuPDG != 14 ) {
      std::cout << "FRW: Target fluxes have been generated for muon neutrinos only" <<std::endl;
      std::exit(1);
    }

    if( !m_applyNuEConstraint ) {
      std::cout << "FRW: Target fluxes have been generated for nu-e constraints only" <<std::endl;
      std::exit(1);
    }

    const char* plotutils=gSystem->Getenv("PLOTUTILSROOT");
    if (!plotutils || !strlen(plotutils)) {
      std::cout << "$PLOTUTILSROOT is not set. Can't find flux histograms" << std::endl;
      std::exit(1);
    }

    //NOTE: Target fluxes are made with gen2thin and g4numiv6 - 9/7/2021

    TString tracker_filename = TString::Format("%s/data/flux_daisy/%s/flux_with_errors/flux_tracker.root",
                                                plotutils,
                                                project_dir.c_str(),
                                                tar_mat.c_str());

    TString tarfilename = TString::Format("%s/data/flux_daisy/%s/flux_with_errors/flux_%s.root",
                                           plotutils,
                                           project_dir.c_str(),
                                           tar_mat.c_str());

    TString histName = TString::Format("flux");

    MnvH1D* tracker_flux = GetMnvH1D(tracker_filename, histName);
    MnvH1D* tar_flux     = GetMnvH1D(tarfilename, histName);

    //Will use ratio of target/tracker to supply a weight to make the correct flux universes

    MnvH1D* h_flux       = (MnvH1D*)m_fluxReweightNu->Clone(Form("flux_%s",tar_mat.c_str()));
    MnvH1D* ratio_flux   = (MnvH1D*)m_fluxReweightNu->Clone(Form("tmp_flux_ratio_%s",tar_mat.c_str())); 
    ratio_flux->ClearAllErrorBands();
    ratio_flux->Reset();

    //Create a ratio of targets/tracker
    for( int iBin = 0; iBin < ratio_flux->GetNbinsX()+2; ++iBin )
    {
      double binCenter = ratio_flux->GetBinCenter(iBin);

      int iBinTracker  = tracker_flux->FindBin(binCenter);
      int iBinTarget   = tar_flux->FindBin(binCenter);

      double content_tracker = tracker_flux->GetBinContent(iBinTracker);
      double content_target  = tar_flux->GetBinContent(iBinTarget);

      ratio_flux->SetBinContent( iBin, content_tracker > 0 ? content_target/content_tracker : 1.0 );
      ratio_flux->SetBinError(   iBin, 0.0 );
    }

    ratio_flux->AddMissingErrorBandsAndFillWithCV( *h_flux );
    h_flux->Multiply( h_flux, ratio_flux );

    //std::vector< std::string > vert_error_names = h_flux->GetVertErrorBandNames();
    //for( auto &name : vert_error_names ) {
    //  if( strcmp( name.c_str(), "Flux" ) != 0 ) h_flux->PopVertErrorBand(name);
    //}

    //std::vector< std::string > lat_error_names  = h_flux->GetLatErrorBandNames();
    //for( auto &name : lat_error_names ) {
    //  if( strcmp( name.c_str(), "Flux" ) != 0 ) h_flux->PopLatErrorBand(name);
    //}

    return h_flux;

  }

  //============================================================================
  MnvH1D* FluxReweighter::GetDaisyParamMnvH1D(int nuPDG,
                                              std::string tar_mat,
                                              std::string project_dir)
  {

    if( nuPDG != 14 ) {
      std::cout << "FRW: Daisy reweights have been generated for muon neutrinos only" <<std::endl;
      std::exit(1);
    }

    if( !m_applyNuEConstraint ) {
      std::cout << "FRW: Daisy reweights have been generated for nu-e constraints only" <<std::endl;
      std::exit(1);
    }

    const char* plotutils=gSystem->Getenv("PLOTUTILSROOT");
    if (!plotutils || !strlen(plotutils)) {
      std::cout << "$PLOTUTILSROOT is not set. Can't find daisy reweights" << std::endl;
      std::exit(1);
    }

    MnvH1D* h_param;

    if( strcmp( tar_mat.c_str(), "tracker" ) == 0 ) {
      h_param = new MnvH1D("tracker_param_hist", "tracker_param_hist;Flux Bin;Weight", 12, 0., 12.);
      for( int iBin = 1; iBin <= 12; ++iBin ) {
        h_param->SetBinContent(iBin, 1);
        h_param->SetBinError(iBin, 0);
      }
    }
    else{
      TString filename = TString::Format("%s/data/flux_daisy/%s/out_%s_000100.root",
                                           plotutils,
                                           project_dir.c_str(),
                                           tar_mat.c_str());

      TString histName = TString::Format("param_hist");

      //std::cout<<"FRW: Getting daisy reweight parameters"<<std::endl;
      h_param = GetMnvH1D(filename, histName);
    }
    return h_param;    

  }
  */

//============================================================================
// assumes files in
// /data/flux_daisy/neutrinos AND /data/flux_daisy/antineutrinos
// NOTE that it was decided not to use daisy reweight fot water target
MnvH1D* FluxReweighter::GetTargetFluxMnvH1D(int nuPDG,
                                              std::string tar_mat,
                                              std::string project_dir)
{
  
    if( !m_applyNuEConstraint ) {
      std::cout << "FRW: Target fluxes have been generated for FHC + RHC nu-e + IMD constraints only" <<std::endl;
      std::exit(1);
    }
	
    const char* plotutils=gSystem->Getenv("PLOTUTILSROOT");
    if (!plotutils || !strlen(plotutils)) {
      std::cout << "$PLOTUTILSROOT is not set. Can't find flux histograms" << std::endl;
      std::exit(1);
    }

    //NOTE: Target fluxes are made with gen2thin and g4numiv6 - 17/02/2022
	
    // define file variables
    TString tracker_filename;
    TString tarfilename;
    TString histName;
    MnvH1D* tracker_flux;
    MnvH1D* tar_flux;
    MnvH1D* h_flux;
    MnvH1D* ratio_flux;

    if( nuPDG == 14 ) {
      std::cout << "FRW: Using target fluxes for muon NEUTRINOS" <<std::endl;
	  
      tracker_filename = TString::Format("%s/data/flux_daisy/neutrinos/%s/flux_with_errors/flux_tracker.root",
                                                  plotutils,
                                                  project_dir.c_str(),
                                                  tar_mat.c_str());

      tarfilename = TString::Format("%s/data/flux_daisy/neutrinos/%s/flux_with_errors/flux_%s.root",
                                             plotutils,
                                             project_dir.c_str(),
                                             tar_mat.c_str());
      /*
											 
	    // Old fluxes that have been generated for FHC nu-e constraint only										 
      TString tracker_filename = TString::Format("%s/data/flux_daisy/%s/flux_with_errors/flux_tracker.root",
                                                          plotutils,
                                                          project_dir.c_str(),
                                                          tar_mat.c_str());

      TString tarfilename = TString::Format("%s/data/flux_daisy/%s/flux_with_errors/flux_%s.root",
                                                     plotutils,
                                                     project_dir.c_str(),
                                                     tar_mat.c_str());
											 
											 
	    */
											 
      histName = TString::Format("flux");

      tracker_flux = GetMnvH1D(tracker_filename, histName);
      tar_flux     = GetMnvH1D(tarfilename, histName);

      //Will use ratio of target/tracker to supply a weight to make the correct flux universes

      h_flux       = (MnvH1D*)m_fluxReweightNu->Clone(Form("flux_%s",tar_mat.c_str()));
      ratio_flux   = (MnvH1D*)m_fluxReweightNu->Clone(Form("tmp_flux_ratio_%s",tar_mat.c_str())); 
      
    }
	
    else if( nuPDG == -14 ) {
      std::cout << "FRW: Using target fluxes for muon ANTINEUTRINOS" <<std::endl;
	  
      tracker_filename = TString::Format("%s/data/flux_daisy/antineutrinos/%s/flux_with_errors/flux_tracker.root",
                                                  plotutils,
                                                  project_dir.c_str(),
                                                  tar_mat.c_str());

      tarfilename = TString::Format("%s/data/flux_daisy/antineutrinos/%s/flux_with_errors/flux_%s.root",
                                             plotutils,
                                             project_dir.c_str(),
                                             tar_mat.c_str());
											 
      histName = TString::Format("flux");

      tracker_flux = GetMnvH1D(tracker_filename, histName);
      tar_flux     = GetMnvH1D(tarfilename, histName);

      //Will use ratio of target/tracker to supply a weight to make the correct flux universes

      h_flux       = (MnvH1D*)m_fluxReweightNubar->Clone(Form("flux_%s",tar_mat.c_str()));
      ratio_flux   = (MnvH1D*)m_fluxReweightNubar->Clone(Form("tmp_flux_ratio_%s",tar_mat.c_str())); 
      
    }
	
    else {
      std::cout << "FRW: Cannot use target fluxes for given PDG (not generated). " <<std::endl;
      
    }


    ratio_flux->ClearAllErrorBands();
    ratio_flux->Reset();
	

    //Create a ratio of targets/tracker
    for( int iBin = 0; iBin < ratio_flux->GetNbinsX()+2; ++iBin )
    {
      double binCenter = ratio_flux->GetBinCenter(iBin);

      int iBinTracker  = tracker_flux->FindBin(binCenter);
      int iBinTarget   = tar_flux->FindBin(binCenter);

      double content_tracker = tracker_flux->GetBinContent(iBinTracker);
      double content_target  = tar_flux->GetBinContent(iBinTarget);

      ratio_flux->SetBinContent( iBin, content_tracker > 0 ? content_target/content_tracker : 1.0 );
      ratio_flux->SetBinError(   iBin, 0.0 );
    }

    ratio_flux->AddMissingErrorBandsAndFillWithCV( *h_flux );
    h_flux->Multiply( h_flux, ratio_flux );

    //std::vector< std::string > vert_error_names = h_flux->GetVertErrorBandNames();
    //for( auto &name : vert_error_names ) {
    //  if( strcmp( name.c_str(), "Flux" ) != 0 ) h_flux->PopVertErrorBand(name);
    //}

    //std::vector< std::string > lat_error_names  = h_flux->GetLatErrorBandNames();
    //for( auto &name : lat_error_names ) {
    //  if( strcmp( name.c_str(), "Flux" ) != 0 ) h_flux->PopLatErrorBand(name);
    //}

    return h_flux;

  }

//============================================================================
  MnvH1D* FluxReweighter::GetDaisyParamMnvH1D(int nuPDG,
                                              std::string tar_mat,
                                              std::string project_dir)
  {
    
	  if( !m_applyNuEConstraint ) {
      std::cout << "FRW: Target fluxes have been generated for FHC + RHC nu-e + IMD constraints only" <<std::endl;
      std::exit(1);
    }

    const char* plotutils=gSystem->Getenv("PLOTUTILSROOT");
    if (!plotutils || !strlen(plotutils)) {
      std::cout << "$PLOTUTILSROOT is not set. Can't find daisy reweights" << std::endl;
      std::exit(1);
    }

    MnvH1D* h_param;

    if( strcmp( tar_mat.c_str(), "tracker" ) == 0 ) {
      h_param = new MnvH1D("tracker_param_hist", "tracker_param_hist;Flux Bin;Weight", 12, 0., 12.);
      for( int iBin = 1; iBin <= 12; ++iBin ) {
        h_param->SetBinContent(iBin, 1);
        h_param->SetBinError(iBin, 0);
      }
    }
    else{
	  TString filename;
	  
	  if( nuPDG == 14 ){
	      std::cout << "Using tracker flux weights for muon NEUTRINOS" <<std::endl;
	      filename = TString::Format("%s/data/flux_daisy/neutrinos/%s/out_%s_000100.root",
	                                           plotutils,
	                                           project_dir.c_str(),
	                                           tar_mat.c_str());

	  }
	  
	  else if ( nuPDG == -14 ){
		  std::cout << "Using tracker flux weights for muon ANTINEUTRINOS" <<std::endl;
	      filename = TString::Format("%s/data/flux_daisy/antineutrinos/%s/out_%s_000100.root",
	                                           plotutils,
	                                           project_dir.c_str(),
	  										   tar_mat.c_str());
		  
	  }
	  else{
		  std::cout << "Tracker flux weights were not generated for given PDG. " <<std::endl;
		  std::exit(1);
	  	
	  }
		

      TString histName = TString::Format("param_hist");

      //std::cout<<"FRW: Getting daisy reweight parameters"<<std::endl;
      h_param = GetMnvH1D(filename, histName);
    }
    return h_param;    

  }
  
  //============================================================================
  MnvH1D* FluxReweighter::GetFluxMnvH1D(int nuPDG,
      enum EPlaylist playlist,
      enum EFluxVersion fluxVersion,
      enum EG4NumiVersion g4NumiVersion,
      bool useGen  /*=false*/)
  {
    const char* plotutils=gSystem->Getenv("PLOTUTILSROOT");
    if (!plotutils || !strlen(plotutils)) {
      std::cout << "$PLOTUTILSROOT is not set. Can't find flux histograms" << std::endl;
      std::exit(1);
    }

    if(m_applyNuEConstraint&&IsCustomFlux(fluxVersion)){
      std::cout << "Application of the nu-e constraint is incompapitable with using a non-PPFX flux!! I'm exiting." << std::endl;
      std::exit(1);
    }

    const char* constraintString;
    int iplaylist = LeOrMe(playlist);
    // Timebeing 1 = LE; 2 = ME (ALL Playlist)
    if (iplaylist==1) {
      constraintString =  "";
    }
    //This is the nu+e
    if(iplaylist==2 or iplaylist==3)constraintString = m_applyNuEConstraint ? "_rearrangedUniverses" : "";
    //This is the IMD method
    //    if(iplaylist==2 or iplaylist==3)constraintString = m_applyNuEConstraint ? "" : "";

    TString filename;
    const char* histName;

    // If using lowNu or highNu flux, use different prescription to fetch it
    // For the reweighted flux
    if(IsCustomFlux(fluxVersion)&&!useGen){
      histName = fluxForSystematicsHistName(fluxVersion);
      filename = TString::Format("%s/data/flux/sys/%s.root",
            plotutils,
            histName);
    }
    // For the generated flux
    else if(IsCustomFlux(fluxVersion)&&useGen){
      histName = "genFlux_coarseBinning_MnvH1D";
      filename = TString::Format("%s/data/flux/sys/genFlux_coarseBinning_MnvH1D.root",
            plotutils);
    }
    // Default prescription for assembling file and hist names
    else{
      const char* fluxString;
      fluxString = useGen ? g4NumiVersionString(g4NumiVersion) : fluxVersionString(fluxVersion);

      filename = TString::Format("%s/data/flux/flux-%s-pdg%d-%s%s.root",
            plotutils,
            fluxString,
            nuPDG,
            playlistString(playlist),
            constraintString);

      histName = useGen ? "flux_E_unweighted" : "flux_E_cvweighted";
    }

    std::cout << "FRW: this is the flux file I'm using: " << filename << std::endl;
    return GetMnvH1D(filename, histName);
  }
  //============================================================================

  MnvH1D* FluxReweighter::GetMELowNuMnvH1D()
  {
    const char* plotutils=gSystem->Getenv("PLOTUTILSROOT");
    if (!plotutils || !strlen(plotutils)) {
      std::cout << "$PLOTUTILSROOT is not set. Can't find flux histograms" << std::endl;
      std::exit(1);
    }

    TString filename(TString::Format("%s/data/lowNuReweight/lowNuReweight.root",
          plotutils));
    return GetMnvH1D(filename, "ratio");
  }


  //============================================================================
  TSpline3* FluxReweighter::GetSpline(MnvH1D *h)
  {
    const char* plotutils=gSystem->Getenv("PLOTUTILSROOT");
    if (!plotutils || !strlen(plotutils)) {
      std::cout << "$PLOTUTILSROOT is not set. Can't find flux histograms" << std::endl;
      std::exit(1);
    }

    TSpline3 *spl = new TSpline3(h);
    return spl;
  }


  //============================================================================
  double FluxReweighter::GetFluxCVWeight( double Enu, int nuPDG )
  {
    MnvH1D* genFlux = nuPDG > 0 ? m_fluxGenNu      : m_fluxGenNubar;
    MnvH1D* newFlux = nuPDG > 0 ? m_fluxReweightNu : m_fluxReweightNubar;

    double num = newFlux->Interpolate( Enu );
    double dem = genFlux->Interpolate( Enu );
    if(Enu>75.){//Condition to avoid the large slope at large Enu which causes crazy weights. Slope in that region is 10^-22 or so
      num = newFlux->GetBinContent(newFlux->FindBin(Enu));
      dem = genFlux->GetBinContent(genFlux->FindBin(Enu));
    }
    if (num ==0 ||dem ==0) {
      return 1;
    }
    return  num/dem;
  }


  //============================================================================
  // Suggested alternate name: FRW::GetFluxSystematicUniverseWeight
  double FluxReweighter::GetSysUniFluxWeightCorrection( double Enu,
                                                        int nuPDG,
                                                        std::string sys_name,
                                                        int universe)
  // This method returns the multiplicative factor that should modify the flux weight
  // in the evaluation of some systematic universes. Currently the use cases are
  // the evaluation in FHC of the muon-momentum-shifted universes for most users and
  // the high-nu flux analysis
  // ----------
  // POSSIBLY USEFUL TO FUTURE ANALYSTS: This method assumes the user
  // uses the same names for their systematic universes as the available
  // fluxes. To accomodate a different naming scheme, the user would have
  // to add a method that maps user names onto the "standard" naming
  // scheme. The best idea I have for how that would work is that FRW
  // would own a map, and the user would interactively add (or overwrite)
  // entries in that map to reflect their naming scheme. RDF 2020-01-10
  // ----------
  {
    // FluxCorrection not valid for nue
    if (nuPDG == -12 ) return 1;
    if (nuPDG == 12 ) return 1;

    //We also don't want the wrong sign fluxes to pick this correction... -David L. putting in Amit's changes from CVS commit 1.106-8
    int iPlaylist = LeOrMe(d_Playlist);
    if (iPlaylist==2 && nuPDG==-14) return 1;
    if (iPlaylist==3 && nuPDG==14) return 1;

    // FluxCorrection not valid for RHC numu (yet?) It is now!
    //if (nuPDG == -14 ) return 1; Commented out now that RHC is valid. -David L. putting in Amit's changes from CVS commit 1.106-8
    // Correction for Amit's Emu fluxes not useful above 15 GeV
    if (sys_name == "Muon_Energy" && Enu > 15) return 1;
    if (sys_name == "Muon_Energy_MINOS" && Enu > 15) return 1;

    TH1D* cvFlux = m_fluxSystematicsMnvH1D;

    TH1D* sysFlux = (TH1D*)m_fluxSystematicsMnvH1D->GetVertErrorBand(sys_name.c_str())->GetHist(universe);

    /* Commented out now that RHC is valid. Replaced by the above lines. -David L. putting in Amit's changes from CVS commit 1.106-8
    TH1D* cvFlux = nuPDG > 0 ?
      m_fluxSystematicsMnvH1D :
      NULL; // We may have a RHC case in the future

    TH1D* sysFlux = nuPDG > 0 ?
        (TH1D*)m_fluxSystematicsMnvH1D->GetVertErrorBand(sys_name.c_str())->GetHist(universe) :
        NULL; // We may have a RHC case in the future
    */

    double num = sysFlux->Interpolate( Enu );
    double dem = cvFlux->Interpolate( Enu );

    if (num == 0 || dem == 0) return 1;

    return num/dem;
  }


  //======================================================================
  void FluxReweighter::SetFluxSysMnvH1D( int nuPDG, enum EFluxVersion fluxVersion)
  {
    const char* plotutils=gSystem->Getenv("PLOTUTILSROOT");
    bool isNu = nuPDG>0?true:false; //Kept around in Amit CVS commit 1.108 for correlated muon fluxes. -David L.
    const char* fluxName = fluxForSystematicsHistName(fluxVersion);

    //if(nuPDG>0){ Commented out so that RHC correlations in fact grab the necessary file/histo. -David L.
    //std::cout << "Creating alternate sysmap using the " << fluxName << " flux." << std::endl;
    TFile myfile(TString::Format("%s/data/flux/sys/%s.root",plotutils,fluxName));
    m_fluxSystematicsMnvH1D = (MnvH1D*)myfile.Get(TString::Format("%s",fluxName));
    // assert will provide useful debug info only when the code is compiled in debug mode
    assert(m_fluxSystematicsMnvH1D != NULL && "Failed to find a histogram of that name in FluxReweighter!  Maybe its name changed?");
    m_fluxSystematicsMnvH1D->SetDirectory(0);
    //} Commented out so that RHC correlations in fact grab the necessary file/histo. -David L.

  }


  //============================================================================
  double FluxReweighter::GetMELowNuFluxWeight( double Enu )
  {
    return m_MELowNuDataMCRatioSpline->Eval(Enu);
  }


  //============================================================================
  MnvH1D* FluxReweighter::GetRebinnedFluxGenerated(int nuPDG,
                                                   MnvH1D* template_hist)
  {
    MnvH1D* h_flux = this->GetFluxGenerated(nuPDG);
        //h_flux->SaveAs("fluxchecker_wowiggle.root");

    MnvH1D* h_flux_rebinned =
        (MnvH1D*)template_hist->Clone("generatedflux_rebinned");
    h_flux_rebinned->ClearAllErrorBands();
    h_flux_rebinned->Reset();

    TH1D* tmp_flux_cv     =
        new TH1D(h_flux->GetCVHistoWithStatError());
    TH1D* tmp_template_cv =
        new TH1D(h_flux_rebinned->GetCVHistoWithStatError());

    FluxReweighter::RebinFluxHist(tmp_flux_cv,tmp_template_cv);

    //CV first
    for(int i=0;i<h_flux_rebinned->GetNbinsX()+2;i++)
      h_flux_rebinned->SetBinContent(i,tmp_template_cv->GetBinContent(i));

    //clean my mess
    delete tmp_flux_cv;
    delete tmp_template_cv;

    //Now Flux Universes // DON'T Assume Flux is the only error
    std::vector<std::string> vertNames = h_flux->GetVertErrorBandNames();
    for(unsigned int k=0; k<vertNames.size(); ++k ) {
      MnvVertErrorBand *errBand = h_flux->GetVertErrorBand( vertNames[k] );
      const int universes = errBand->GetNHists();
      std::vector<TH1D*> vert_hists;
      for(int u=0;u<universes;++u) {
        TH1D* tmp_flux = new TH1D(*errBand->GetHist( u ));
        TH1D* tmp_template = new TH1D(h_flux_rebinned->GetCVHistoWithStatError());
        tmp_template->SetName(Form("Flux_rebinned_universe_%d",u));
        RebinFluxHist(tmp_flux,tmp_template);
        vert_hists.push_back(tmp_template);
      }
      h_flux_rebinned->AddVertErrorBand( vertNames[k],vert_hists);
      //clean my mess
      for(std::vector<TH1D*>::iterator itHist = vert_hists.begin(); itHist != vert_hists.end(); ++itHist)
        delete *itHist;
    }
    //AddMissingError with CV
    h_flux_rebinned->AddMissingErrorBandsAndFillWithCV(*template_hist);
    //h_flux_rebinned->SaveAs("Enu_Flux_nukecc.root");

    return h_flux_rebinned;
  }


  //============================================================================
  MnvH1D* FluxReweighter::GetRebinnedFluxReweighted(int nuPDG,
                                                    MnvH1D* template_hist)
  {
    MnvH1D* h_flux = this->GetFluxReweighted(nuPDG);
    h_flux->SaveAs("fluxchecker_wowiggle.root");

    MnvH1D* h_flux_rebinned =
        (MnvH1D*)template_hist->Clone("reweightedflux_rebinned");
    h_flux_rebinned->ClearAllErrorBands();
    h_flux_rebinned->Reset();

    TH1D* tmp_flux_cv     =
        new TH1D(h_flux->GetCVHistoWithStatError());
    TH1D* tmp_template_cv =
        new TH1D(h_flux_rebinned->GetCVHistoWithStatError());

    FluxReweighter::RebinFluxHist(tmp_flux_cv,tmp_template_cv);

    //CV first
    for(int i=0;i<h_flux_rebinned->GetNbinsX()+2;i++)
      h_flux_rebinned->SetBinContent(i,tmp_template_cv->GetBinContent(i));

    //clean my mess
    delete tmp_flux_cv;
    delete tmp_template_cv;

    //Now Flux Universes // DON'T Assume Flux is the only error
    std::vector<std::string> vertNames = h_flux->GetVertErrorBandNames();
    for(unsigned int k=0; k<vertNames.size(); ++k ) {
      MnvVertErrorBand *errBand = h_flux->GetVertErrorBand( vertNames[k] );
      const int universes = errBand->GetNHists();
      std::vector<TH1D*> vert_hists;
      for(int u=0;u<universes;++u) {
        TH1D* tmp_flux = new TH1D(*errBand->GetHist( u ));
        TH1D* tmp_template = new TH1D(h_flux_rebinned->GetCVHistoWithStatError());
        tmp_template->SetName(Form("Flux_rebinned_universe_%d",u));
        RebinFluxHist(tmp_flux,tmp_template);
        vert_hists.push_back(tmp_template);
      }
      h_flux_rebinned->AddVertErrorBand( vertNames[k],vert_hists);
      //clean my mess
      for(std::vector<TH1D*>::iterator itHist = vert_hists.begin();
          itHist != vert_hists.end(); ++itHist)
        delete *itHist;
    }
    //AddMissingError with CV
    h_flux_rebinned->AddMissingErrorBandsAndFillWithCV(*template_hist);
    //h_flux_rebinned->SaveAs("Enu_Flux_nukecc.root");

    // This method will fill the constraint weights in the Flux error band if appropriate
    CheckAndFixFluxErrorBand( h_flux_rebinned );

    return h_flux_rebinned;
  }


  //============================================================================
  MnvH1D* FluxReweighter::GetRebinnedFluxReweighted_FromInputFlux(MnvH1D* input_flux,
                                                                  MnvH1D* template_hist)
  {
    MnvH1D* h_flux = (MnvH1D*)input_flux->Clone("input_flux"); //this->GetFluxReweighted(nuPDG);

    MnvH1D* h_flux_rebinned =
        (MnvH1D*)template_hist->Clone("reweightedflux_rebinned");
    h_flux_rebinned->ClearAllErrorBands();
    h_flux_rebinned->Reset();

    TH1D* tmp_flux_cv     =
        new TH1D(h_flux->GetCVHistoWithStatError());
    TH1D* tmp_template_cv =
        new TH1D(h_flux_rebinned->GetCVHistoWithStatError());

    FluxReweighter::RebinFluxHist(tmp_flux_cv,tmp_template_cv);

    //CV first
    for(int i=0;i<h_flux_rebinned->GetNbinsX()+2;i++)
      h_flux_rebinned->SetBinContent(i,tmp_template_cv->GetBinContent(i));

    //clean my mess
    delete tmp_flux_cv;
    delete tmp_template_cv;

    //Now Flux Universes // DON'T Assume Flux is the only error
    std::vector<std::string> vertNames = h_flux->GetVertErrorBandNames();
    for(unsigned int k=0; k<vertNames.size(); ++k ) {
      MnvVertErrorBand *errBand = h_flux->GetVertErrorBand( vertNames[k] );
      const int universes = errBand->GetNHists();
      std::vector<TH1D*> vert_hists;
      for(int u=0;u<universes;++u) {
        TH1D* tmp_flux = new TH1D(*errBand->GetHist( u ));
        TH1D* tmp_template = new TH1D(h_flux_rebinned->GetCVHistoWithStatError());
        tmp_template->SetName(Form("Flux_rebinned_universe_%d",u));
        RebinFluxHist(tmp_flux,tmp_template);
        vert_hists.push_back(tmp_template);
      }
      h_flux_rebinned->AddVertErrorBand( vertNames[k],vert_hists);
      //clean my mess
      for(std::vector<TH1D*>::iterator itHist = vert_hists.begin();
          itHist != vert_hists.end(); ++itHist)
        delete *itHist;
    }

    h_flux_rebinned->AddMissingErrorBandsAndFillWithCV(*template_hist);

    //if (m_applyNuEConstraint) {
    //  std::cout << "Applying flux constraint to rebinned histogram" << std::endl;
    //  h_flux_rebinned = Constrainer().ConstrainHisto<MnvH1D, MnvVertErrorBand>( h_flux_rebinned );
    //}

    return h_flux_rebinned;
  }

  //============================================================================
  TH1D* FluxReweighter::GetRebinnedFluxReweighted_FromInputFlux(TH1D* input_flux,
                                                                TH1D* template_hist)

  {
    TH1D* h_flux = (TH1D*)input_flux->Clone("input_flux"); //this->GetFluxReweighted(nuPDG);
    TH1D* h_flux_rebinned = (TH1D*)template_hist->Clone("reweightedflux_rebinned");
    //h_flux_rebinned->ClearAllErrorBands();
    h_flux_rebinned->Reset();

    TH1D* tmp_flux_cv     =(TH1D*)h_flux->Clone();
    //new TH1D(h_flux_rebinned->GetCVHistoWithStatError());
    TH1D* tmp_template_cv = (TH1D*)h_flux_rebinned->Clone();
    //  new TH1D(h_flux_rebinned->GetCVHistoWithStatError());

    FluxReweighter::RebinFluxHist(tmp_flux_cv,tmp_template_cv);
    //CV first
    for(int i=0;i<h_flux_rebinned->GetNbinsX()+2;i++)
      h_flux_rebinned->SetBinContent(i,tmp_template_cv->GetBinContent(i));

    //clean my mess
    delete tmp_flux_cv;
    delete tmp_template_cv;


    //if (m_applyNuEConstraint) {
    //  std::cout << "Applying flux constraint to rebinned histogram" << std::endl; 
    //  h_flux_rebinned = Constrainer().ConstrainHisto<MnvH1D, MnvVertErrorBand>( h_flux_rebinned );
    //}

    return h_flux_rebinned;
  }

  //============================================================================
  void FluxReweighter::RebinFluxHist(TH1D* h_flux, TH1D*&h_rebinned_flux)
  {
    //strategy is to recale orig by bin width (undo bin width normalization)
    //then combine bins and then rescale by binwidth again
    TH1D* scaler = (TH1D*)h_flux->Clone("fluxscaler");
    TH1D* flux_cv = (TH1D*)h_flux->Clone("fluxcvtomod");
    for(int i=1;i<scaler->GetNbinsX();i++)
      scaler->SetBinContent(i,scaler->GetBinWidth(i));
    flux_cv->Multiply(scaler);//undid bin width normalization
    std::vector<double>rebinned_flux_bin_edges;
    for(int i=1;i<h_rebinned_flux->GetNbinsX()+2;i++) { //need low edge of overflow (high edge of last bin)
      rebinned_flux_bin_edges.push_back(h_rebinned_flux->GetBinLowEdge(i));
    }
    for(unsigned int i=0;i<rebinned_flux_bin_edges.size();i++) { //need low edge of overflow (high edge of last bin)
      //std::cout<<i<<"  "<<h_rebinned_flux->GetBinLowEdge(i)<<std::endl;
      //rebinned_flux_bin_edges.push_back(h_rebinned_flux->GetBinLowEdge(i));
      //std::cout<<i<<"  "<<rebinned_flux_bin_edges[i]<<std::endl;
    }
    h_rebinned_flux = (TH1D*)flux_cv->Rebin(rebinned_flux_bin_edges.size()-1,
                                            "Fluxrebinned",
                                            &rebinned_flux_bin_edges[0]);
    h_rebinned_flux->Scale(1.0,"width");//And redo bin width norm
  }

  //============================================================================
  template<class MnvHistoType>
  MnvHistoType* FluxReweighter::GetIntegratedFluxReweighted(int nuPDG,
                                                      MnvHistoType* template_hist,
                                                      double min_energy,
                                                      double max_energy,
                                                      bool use_muon_correlations)
  {

    // This is the hist that will be retured to the user. Start from a clean slate
    MnvHistoType* h_flux_integrated = (MnvHistoType*)template_hist->Clone("reweightedflux_integrated");
    h_flux_integrated->ClearAllErrorBands();
    h_flux_integrated->Reset();

    // The default flux prescription is PPFX plus correlated fluxes for
    // "Muon_Energy" or "Muon_Energy_MINOS" systematic universes
    if(m_useStandardFlux){

      MnvH1D* h_flux_ppfx = this->GetFluxReweighted(nuPDG);

      //Get bins to integrate over
      int ppfx_b_min = h_flux_ppfx->FindBin( min_energy );
      int ppfx_b_max = h_flux_ppfx->FindBin( max_energy );
      //get integral
      double int_flux_cv = h_flux_ppfx->Integral( ppfx_b_min, ppfx_b_max , "width" );
      // The correlated muon fluxes only provide predictions below 15 GeV, which makes correcting
      // the CV PPFX integrate flux a little more complicated
      int ppfx_b_15GeV = h_flux_ppfx->FindBin( 15. );
      double int_flux_cv_below15GeV;
      double int_flux_cv_above15GeV;
      if(min_energy>15.0){int_flux_cv_below15GeV = 0.;}
      else{int_flux_cv_below15GeV = h_flux_ppfx->Integral( ppfx_b_min, ppfx_b_15GeV, "width" );}
      if(max_energy<=15.0){int_flux_cv_above15GeV = 0.;}
      else{int_flux_cv_above15GeV = h_flux_ppfx->Integral( ppfx_b_15GeV+1, ppfx_b_max , "width" );}
      //CV first
      for(int i=0;i<h_flux_integrated->GetSize();i++)
        h_flux_integrated->SetBinContent(i,int_flux_cv);

      // Integrate flux in every "Flux" systematic universe
      // To use a different PPFX error band, modify FRW initialization
      MnvVertErrorBand *errBand = h_flux_ppfx->GetVertErrorBand(m_fluxErrorName);
      const int universes = errBand->GetNHists();
      auto flux_sys_hists = GetVector(h_flux_integrated);
      for(int u=0;u<universes;++u) {
        TH1D* tmp_flux = new TH1D(*errBand->GetHist( u ));
        auto tmp_template = h_flux_integrated->GetCVHistoWithStatError();
        tmp_template.SetName(Form("Flux_integrated_universe_%d",u));
        double flux_uni = tmp_flux->Integral(ppfx_b_min,ppfx_b_max,"width");
        for(int i=0;i<h_flux_integrated->GetSize();i++)
          tmp_template.SetBinContent(i,flux_uni);
        flux_sys_hists.push_back(NewHist(tmp_template));
      }
      // Push the constructed error band into the return hist
      h_flux_integrated->AddVertErrorBand(m_fluxErrorName,flux_sys_hists);

      // Integrate flux in special correlated flux universes. This is a correction on
      // the CV integrated flux calculated using separately produced systematic and CV
      // fluxes. This matches the implementation of FRW::GetSysUniFluxWeightCorrection
      if(use_muon_correlations){

        // The correlated flux is constructed relative to a different CV
        TH1D* tmp_flux_special_cv = m_fluxSystematicsMnvH1D;
        //Get bins to integrate over
        int special_b_min = tmp_flux_special_cv->FindBin( min_energy );
        int special_b_max = tmp_flux_special_cv->FindBin( max_energy );
        //get integral
        double int_flux_special_cv = tmp_flux_special_cv->Integral(special_b_min,special_b_max,"width");

        // Don't assume which convention the analyzer is using
        std::vector<std::string> errorBandsToLoop;
        if(template_hist->HasErrorBand("Muon_Energy")){errorBandsToLoop.push_back("Muon_Energy");}
        if(template_hist->HasErrorBand("Muon_Energy_MINOS")){errorBandsToLoop.push_back("Muon_Energy_MINOS");}

        for(unsigned int k=0; k<errorBandsToLoop.size(); ++k ) {
          auto correlated_flux_hists = GetVector(h_flux_integrated);
          MnvVertErrorBand *errBand = m_fluxSystematicsMnvH1D->GetVertErrorBand( errorBandsToLoop[k] );
          for(int u=0;u<2;++u) { // There are two systematic variations corresponding to +/- 1 sigma
            // Fetch correlated flux from separate container
            TH1D* tmp_flux_special_sys = (TH1D*)errBand->GetHist(u);
            auto tmp_template = h_flux_integrated->GetCVHistoWithStatError();
            std::string tempString = Form("Flux_integrated_%s_universe_%d",errorBandsToLoop[k].c_str(),u);
            tmp_template.SetName(tempString.c_str());
            double int_flux_special_sys = tmp_flux_special_sys->Integral(special_b_min,special_b_max,"width");
            for(int i=0;i<h_flux_integrated->GetSize();i++)
              tmp_template.SetBinContent(i,(int_flux_special_sys/int_flux_special_cv)*int_flux_cv_below15GeV+int_flux_cv_above15GeV);
            correlated_flux_hists.push_back(NewHist(tmp_template));
          }
          // Push the constructed error band into the return hist.
          // Figure out if the analyzer is using a LatErrorBand or VertErrorBand :/
          if(template_hist->HasLatErrorBand(errorBandsToLoop[k])){
            h_flux_integrated->AddLatErrorBand(errorBandsToLoop[k],correlated_flux_hists);
          }
          else{
            h_flux_integrated->AddVertErrorBand(errorBandsToLoop[k],correlated_flux_hists);
          }
          // When MnvLatErrorBand goes away, replace the above (7) lines with the below line
          //h_flux_integrated->AddVertErrorBand(errorBandsToLoop[k],correlated_flux_hists);
        }

      }

    }

    // If the low-nu or high-nu flux is being used, there is a correlated flux
    // for many universes of the analysis
    else{

      // The correlated flux is constructed relative to a different CV
      MnvH1D* h_flux_special = m_fluxSystematicsMnvH1D;
      //Get bins to integrate over
      int special_b_min = h_flux_special->FindBin( min_energy );
      int special_b_max = h_flux_special->FindBin( max_energy );
      //get integral
      double int_flux_special_cv = h_flux_special->Integral(special_b_min,special_b_max,"width");

      //CV first
      for(int i=0;i<h_flux_integrated->GetSize();i++)
        h_flux_integrated->SetBinContent(i,int_flux_special_cv);

      // Loop over user's vertical error bands
      std::vector<std::string> vertErrorBandNames = template_hist->GetVertErrorBandNames();
      for(unsigned int k=0; k<vertErrorBandNames.size(); ++k ) {
        if(!h_flux_special->HasVertErrorBand(vertErrorBandNames[k])){continue;}
        auto correlated_flux_hists = GetVector(h_flux_integrated);
        //std::vector<THND*> correlated_flux_hists;
        MnvVertErrorBand *errBand = h_flux_special->GetVertErrorBand( vertErrorBandNames[k] );
        const int universes = errBand->GetNHists();
        for(int u=0; u<universes; ++u ) {
          // Fetch correlated flux from separate container
          TH1D* tmp_flux_special_sys = (TH1D*)errBand->GetHist(u);
          auto tmp_template = h_flux_integrated->GetCVHistoWithStatError();
          std::string tempString = Form("Flux_integrated_%s_universe_%d",vertErrorBandNames[k].c_str(),u);
          tmp_template.SetName(tempString.c_str());
          double int_flux_special_sys = tmp_flux_special_sys->Integral(special_b_min,special_b_max,"width");
          for(int i=0;i<h_flux_integrated->GetSize();i++)
            tmp_template.SetBinContent(i,int_flux_special_sys);
          correlated_flux_hists.push_back(NewHist(tmp_template));
        }
        // Push the constructed error band into the return hist
        h_flux_integrated->AddVertErrorBand(vertErrorBandNames[k],correlated_flux_hists);
      }

      // This block can go away when MnvLatErrorBand becomes deprecated
      // Loop over user's lateral error bands
      std::vector<std::string> latErrorBandNames = template_hist->GetLatErrorBandNames();
      for(unsigned int k=0; k<latErrorBandNames.size(); ++k ) {
        if(!h_flux_special->HasVertErrorBand(latErrorBandNames[k])){continue;}
        auto correlated_flux_hists = GetVector(h_flux_integrated);
        //std::vector<THND*> correlated_flux_hists;
        MnvVertErrorBand *errBand = h_flux_special->GetVertErrorBand( latErrorBandNames[k] );
        const int universes = errBand->GetNHists();
        for(int u=0; u<universes; ++u ) {
          // Fetch correlated flux from separate container
          TH1D* tmp_flux_special_sys = (TH1D*)errBand->GetHist(u);
          auto tmp_template = h_flux_integrated->GetCVHistoWithStatError();
          std::string tempString = Form("Flux_integrated_%s_universe_%d",latErrorBandNames[k].c_str(),u);
          tmp_template.SetName(tempString.c_str());
          double int_flux_special_sys = tmp_flux_special_sys->Integral(special_b_min,special_b_max,"width");
          for(int i=0;i<h_flux_integrated->GetSize();i++)
            tmp_template.SetBinContent(i,int_flux_special_sys);
          correlated_flux_hists.push_back(NewHist(tmp_template));
        }
        // Push the constructed error band into the return hist
        h_flux_integrated->AddLatErrorBand(latErrorBandNames[k],correlated_flux_hists);
      }

    }
    // Fill CV for any systematic universes which don't have a non-CV flux prediction
    h_flux_integrated->AddMissingErrorBandsAndFillWithCV(*template_hist);
    this->CheckAndFixFluxErrorBand(h_flux_integrated);
    return h_flux_integrated;

  }
  template MnvH1D* FluxReweighter::GetIntegratedFluxReweighted<MnvH1D>( int nuPDG,
                                                                 MnvH1D* template_hist,
                                                                 double min_energy,
                                                                 double max_energy,
                                                                 bool use_muon_correlations);
  template MnvH2D* FluxReweighter::GetIntegratedFluxReweighted<MnvH2D>( int nuPDG,
                                                                 MnvH2D* template_hist,
                                                                 double min_energy,
                                                                 double max_energy,
                                                                 bool use_muon_correlations);



  //============================================================================
  template<class MnvHistoType>
  MnvHistoType* FluxReweighter::GetIntegratedFluxReweighted_FromInputFlux(MnvH1D* input_flux,
                                                                          MnvHistoType* template_hist,
                                                                          double min_energy,
                                                                          double max_energy)
  {
    MnvH1D* h_flux = (MnvH1D*)input_flux->Clone("input_flux");

    MnvHistoType* h_flux_integrated =
        (MnvHistoType*)template_hist->Clone("reweightedflux_integrated");

    h_flux_integrated->ClearAllErrorBands();
    h_flux_integrated->Reset();
    //Get bins to integrate over
    int b_min = h_flux->FindBin( min_energy );
    int b_max = h_flux->FindBin( max_energy );

    //get integral
    double flux_cv = h_flux->Integral( b_min, b_max , "width" );
    //CV first
    for(int i=0;i<h_flux_integrated->GetSize();i++)
      h_flux_integrated->SetBinContent(i,flux_cv);

    h_flux_integrated->AddMissingErrorBandsAndFillWithCV(*template_hist);
    //h_flux_integrated->SaveAs("x_Flux_XSecLooper.root");
    return h_flux_integrated;
  }
  template MnvH1D* FluxReweighter::GetIntegratedFluxReweighted_FromInputFlux<MnvH1D>(MnvH1D* input_flux,
                                                                                     MnvH1D* template_hist,
                                                                                     double min_energy,
                                                                                     double max_energy);
  template MnvH2D* FluxReweighter::GetIntegratedFluxReweighted_FromInputFlux<MnvH2D>(MnvH1D* input_flux,
                                                                                     MnvH2D* template_hist,
                                                                                     double min_energy,
                                                                                     double max_energy);


  //============================================================================
  MnvH1D* FluxReweighter::GetIntegratedTargetFlux(int nuPDG, 
                                                  std::string tar_mat, 
                                                  MnvH1D* template_hist, 
                                                  double min_energy, 
                                                  double max_energy, 
                                                  std::string project_dir)
  {
    MnvH1D* h_flux = GetTargetFluxMnvH1D(nuPDG, tar_mat, project_dir);
    MnvH1D* h_flux_integrated = GetIntegratedFluxReweighted_FromInputFlux( h_flux, template_hist, min_energy, max_energy );
    
    //Put in the correct flux universes
    if( !h_flux_integrated->HasVertErrorBand("Flux") ) return h_flux_integrated;

    MnvVertErrorBand* flux_int_band = h_flux_integrated->PopVertErrorBand("Flux");
    MnvVertErrorBand* flux_band     = h_flux->GetVertErrorBand("Flux");

    const int nUni = flux_int_band->GetNHists();
    if( flux_int_band->GetNHists() > flux_band->GetNHists() ) 
    {
      std::cout<<"FRW: Target flux band only has "<<flux_band->GetNHists()<<" uni while hist has "<<flux_int_band->GetNHists()<<" uni"<<std::endl;
      std::exit(1);
    }

    const int b_min = h_flux->FindBin( min_energy );
    const int b_max = h_flux->FindBin( max_energy );
    std::vector< TH1D* > flux_int_hists; 
    for( int iUni = 0; iUni < nUni; ++iUni )
    {
      TH1D* tmp_flux_uni = (TH1D*)flux_int_band->GetHist(iUni)->Clone(Form("Custom_flux_integrated_universe_%d",iUni));
      tmp_flux_uni->Clear();

      double flux_int = flux_band->GetHist(iUni)->Integral( b_min, b_max, "width" );
      for( int iBin = 0; iBin < tmp_flux_uni->GetSize(); ++iBin ) tmp_flux_uni->SetBinContent( iBin, flux_int );
      flux_int_hists.push_back(tmp_flux_uni);
    }

    h_flux_integrated->AddVertErrorBand("Flux",flux_int_hists);

    return h_flux_integrated;

  }

  MnvH2D* FluxReweighter::GetIntegratedTargetFlux(int nuPDG, 
                                                        std::string tar_mat, 
                                                        MnvH2D* template_hist, 
                                                        double min_energy, 
                                                        double max_energy, 
                                                        std::string project_dir)
  {
    MnvH1D* h_flux = GetTargetFluxMnvH1D(nuPDG, tar_mat, project_dir);
    MnvH2D* h_flux_integrated = GetIntegratedFluxReweighted_FromInputFlux( h_flux, template_hist, min_energy, max_energy );
    
    //Put in the correct flux universes (if the hist has the Flux vert error band)
    if( !h_flux_integrated->HasVertErrorBand("Flux") ) return h_flux_integrated;

    MnvVertErrorBand2D* flux_int_band = h_flux_integrated->PopVertErrorBand("Flux");
    MnvVertErrorBand* flux_band       = h_flux->GetVertErrorBand("Flux");

    const int nUni = flux_int_band->GetNHists();
    if( flux_int_band->GetNHists() > flux_band->GetNHists() ) 
    {
      std::cout<<"FRW: Target flux band only has "<<flux_band->GetNHists()<<" uni while hist has "<<flux_int_band->GetNHists()<<" uni"<<std::endl;
      std::exit(1);
    }

    const int b_min = h_flux->FindBin( min_energy );
    const int b_max = h_flux->FindBin( max_energy );
    std::vector< TH2D* > flux_int_hists; 
    for( int iUni = 0; iUni < nUni; ++iUni )
    {
      TH2D* tmp_flux_uni = (TH2D*)flux_int_band->GetHist(iUni)->Clone(Form("Custom_flux_integrated_universe_%d",iUni));
      tmp_flux_uni->Clear();

      double flux_int = flux_band->GetHist(iUni)->Integral( b_min, b_max, "width" );
      for( int iBin = 0; iBin < tmp_flux_uni->GetSize(); ++iBin ) tmp_flux_uni->SetBinContent( iBin, flux_int );
      flux_int_hists.push_back(tmp_flux_uni);
    }

    h_flux_integrated->AddVertErrorBand("Flux",flux_int_hists);

    return h_flux_integrated;
  }
  //============================================================================
  MnvH1D* FluxReweighter::GetReweightedDaisySum(int nuPDG, 
                                                std::string tar_mat, 
                                                std::map<int, MnvH1D*> daisy_eff_hists, 
                                                std::string project_dir)
  {
    MnvH1D* h_param = GetDaisyParamMnvH1D( nuPDG, tar_mat, project_dir );

    const int nBins = h_param->GetNbinsX();
    //Make sure there are the correct daisy keys
    for( int iDaisy = 0; iDaisy < nBins; ++iDaisy ) { 
      if( daisy_eff_hists.find(iDaisy) == daisy_eff_hists.end() ) {
        std::cout<<"Daisy eff hists missing daisy "<<iDaisy<<std::endl;
        std::exit(1);
      }
    }
    MnvH1D* daisy_sum = (MnvH1D*)daisy_eff_hists[0]->Clone(Form("%s_daisySum",daisy_eff_hists[0]->GetName()));
    daisy_sum->SetDirectory(0);
    daisy_sum->Reset();
    
    for( int iDaisy = 0; iDaisy < nBins; ++iDaisy ) {
      daisy_sum->Add( daisy_eff_hists[iDaisy], h_param->GetBinContent(iDaisy+1) );
    }

    //Insert error bands from the param hist
    std::vector< std::string > vert_error_names = h_param->GetVertErrorBandNames();
    for( auto& name : vert_error_names ) {
      MnvVertErrorBand* error_band = h_param->GetVertErrorBand(name);
      const int nHists = error_band->GetNHists();
      daisy_sum->AddVertErrorBand(name, nHists);

      for( uint iHist = 0; iHist < error_band->GetNHists(); ++iHist){
        for( int iDaisy = 0; iDaisy < nBins; ++iDaisy ) {
          daisy_sum->GetVertErrorBand(name)->GetHist(iHist)->Add( daisy_eff_hists[iDaisy], 
                          error_band->GetHist(iHist)->GetBinContent(iDaisy+1) );
        }
      }
    }
    return daisy_sum;
  }

  MnvH2D* FluxReweighter::GetReweightedDaisySum(int nuPDG, 
                                                std::string tar_mat, 
                                                std::map<int, MnvH2D*> daisy_eff_hists, 
                                                std::string project_dir)
  {
    MnvH1D* h_param = GetDaisyParamMnvH1D( nuPDG, tar_mat, project_dir );

    const int nBins = h_param->GetNbinsX();
    //Make sure there are the correct daisy keys
    for( int iDaisy = 0; iDaisy < nBins; ++iDaisy ) { 
      if( daisy_eff_hists.find(iDaisy) == daisy_eff_hists.end() ) {
        std::cout<<"Daisy eff hists missing daisy "<<iDaisy<<std::endl;
        std::exit(1);
      }
    }
    MnvH2D* daisy_sum = (MnvH2D*)daisy_eff_hists[0]->Clone(Form("%s_daisySum",daisy_eff_hists[0]->GetName()));
    daisy_sum->SetDirectory(0);
    daisy_sum->Reset();
    
    for( int iDaisy = 0; iDaisy < nBins; ++iDaisy ) {
      daisy_sum->Add( daisy_eff_hists[iDaisy], h_param->GetBinContent(iDaisy+1) );
    }

    //Insert error bands from the param hist
    std::vector< std::string > vert_error_names = h_param->GetVertErrorBandNames();
    for( auto& name : vert_error_names ) {
      MnvVertErrorBand* error_band = h_param->GetVertErrorBand(name);
      const int nHists = error_band->GetNHists();
      daisy_sum->AddVertErrorBand(name, nHists);

      for( uint iHist = 0; iHist < error_band->GetNHists(); ++iHist){
        for( int iDaisy = 0; iDaisy < nBins; ++iDaisy ) {
          daisy_sum->GetVertErrorBand(name)->GetHist(iHist)->Add( daisy_eff_hists[iDaisy], 
                          error_band->GetHist(iHist)->GetBinContent(iDaisy+1) );
        }
      }
    }
    return daisy_sum;
  }
  //============================================================================
  double FluxReweighter::GetFluxErrorWeight( double Enu, int nuPDG,
                                             unsigned int universe )
  {
    if ( universe >= m_nFluxUniverses )
    {
      std::cout << Form( "Flux universe %d requested.  FluxReweighter has universes 0-%d",
                          universe, m_nFluxUniverses - 1 )
                << std::endl;
      std::exit(1);
    }
    MnvH1D* newFlux = nuPDG > 0 ? m_fluxReweightNu : m_fluxReweightNubar;
    double univFluxVal =
        newFlux->GetVertErrorBand( m_fluxErrorName )->GetHist( universe )->Interpolate( Enu );


    double num = univFluxVal;
    double dem = newFlux->Interpolate( Enu );

    if(Enu>75.){//Condition to avoid the large slope at large Enu which causes crazy weights. Slope in that region is 10^-22 or so
      num = newFlux->GetVertErrorBand( m_fluxErrorName )->GetHist( universe )->GetBinContent(newFlux->FindBin(Enu));
      dem = newFlux->GetBinContent(newFlux->FindBin(Enu));
    }
    if (num ==0 ||dem ==0) {
      //      std::cout << " I have the following num and denominator " << num << "\t" << dem << std::endl;
      return 1;
    }
    return  num/dem;
  }


  //============================================================================
  std::vector<double> FluxReweighter::GetFluxErrorWeights(double Enu, int nuPDG)
  {
    std::vector<double> weights;
    for ( unsigned int iuniv = 0; iuniv < m_nFluxUniverses; ++iuniv )
    {
      weights.push_back( GetFluxErrorWeight(Enu, nuPDG, iuniv) );
    }
    return weights;
  }


  //============================================================================
  template<class MnvHistoType>
    void FluxReweighter::AddFluxErrorBand( MnvHistoType* h )
    {
      if ( !h->AddVertErrorBand( m_fluxErrorName, m_nFluxUniverses ) )
      {
        std::cout << Form( "Could not add flux error band to histogram \"%s\"",
                           h->GetName() ) << std::endl;
        std::exit(1);
      }
      PropagateNuEConstraintWeights( h );
    }
  template void FluxReweighter::AddFluxErrorBand<MnvH1D>( MnvH1D* h );
  template void FluxReweighter::AddFluxErrorBand<MnvH2D>( MnvH2D* h );


  //============================================================================
  template<class MnvHistoType>
    bool FluxReweighter::CheckFluxErrorBand( MnvHistoType* h )
    {
      if (h->GetVertErrorBand(m_fluxErrorName)->GetUnivWgts()) return true; //There are weights
      else return false;
    }
  template bool FluxReweighter::CheckFluxErrorBand<MnvH1D>( MnvH1D* h );
  template bool FluxReweighter::CheckFluxErrorBand<MnvH2D>( MnvH2D* h );

  //============================================================================
  template<class MnvHistoType>
  void FluxReweighter::CheckAndFixFluxErrorBand( MnvHistoType* h )
    {
      if (!h->GetVertErrorBand(m_fluxErrorName)->GetUnivWgts())
        PropagateNuEConstraintWeights( h );
    }
  template void FluxReweighter::CheckAndFixFluxErrorBand<MnvH1D>( MnvH1D* h );
  template void FluxReweighter::CheckAndFixFluxErrorBand<MnvH2D>( MnvH2D* h );


  //============================================================================
  void FluxReweighter::TruncateNumberOfFluxUniverses( MnvH1D* h, int nUniverses)
  {
    MnvVertErrorBand *poppedFluxErrorBand = h->PopVertErrorBand("Flux");
    std::vector<TH1D*> fluxUniverses = poppedFluxErrorBand->GetHists();
    fluxUniverses.resize(nUniverses);
    h->AddVertErrorBand("Flux",fluxUniverses);
    CheckAndFixFluxErrorBand(h);
  }


  //============================================================================
  void FluxReweighter::FillFluxErrorBand( MnvH1D* h, double val,
                                          double cvweight, double Enu,
                                          int nuPDG )
  {
    if ( !IsFluxErrorBandOK( h ) ) std::exit(1);
    std::vector<double> fluxErrorWeights = GetFluxErrorWeights(Enu, nuPDG);
    double cvWeightFromMe = 1.0;
    if ( !h->FillVertErrorBand( m_fluxErrorName, val, fluxErrorWeights,
                                cvweight, cvWeightFromMe ) )
    {
      std::cout << Form( "Could not fill flux error band for MnvH1D \"%s\"",
                         h->GetName() )
                << std::endl;
      std::exit(1);
    }
    PropagateNuEConstraintWeights( h );
  }


  //===========================================================================
  void FluxReweighter::FillFluxErrorBand( MnvH2D* h, double xval, double yval,
                                          double cvweight, double Enu,
                                          int nuPDG )
  {
    if ( !IsFluxErrorBandOK( h ) ) std::exit(1);
    std::vector<double> fluxErrorWeights = GetFluxErrorWeights(Enu, nuPDG);
    double cvWeightFromMe = 1.0;
    if ( !h->FillVertErrorBand( m_fluxErrorName, xval, yval, fluxErrorWeights,
                                cvweight, cvWeightFromMe ) )
    {
      std::cout << Form( "Could not fill flux error band for MnvH2D \"%s\"",
                         h->GetName() )
                << std::endl;
      std::exit(1);
    }
    PropagateNuEConstraintWeights( h );
  }


  //============================================================================
  template<class MnvHistoType>
    bool FluxReweighter::IsFluxErrorBandOK( MnvHistoType* h )
    {
      if ( !h->HasVertErrorBand( m_fluxErrorName ) )
      {
        std::cout << Form( "Flux error band not found for histogram \"%s\"",
                           h->GetName() )
                  << std::endl;
        return false;
      }
      else if ( h->GetVertErrorBand( m_fluxErrorName )->GetNHists() !=
                m_nFluxUniverses )
      {
        std::cout << Form( "Histogram \"%s\" has %d flux universes. Exactly %d required",
                           h->GetName(),
                           h->GetVertErrorBand( m_fluxErrorName )->GetNHists(),
                           m_nFluxUniverses )
                  << std::endl;
        return false;
      }
      else return true;
    }
  template bool FluxReweighter::IsFluxErrorBandOK<MnvH1D>( MnvH1D* h );
  template bool FluxReweighter::IsFluxErrorBandOK<MnvH2D>( MnvH2D* h );


  //============================================================================
  template<class MnvHistoType>
    void FluxReweighter::PropagateNuEConstraintWeights( MnvHistoType* h )
    {
      if ( m_applyNuEConstraint &&
           fabs( h->GetVertErrorBand( m_fluxErrorName )->GetUnivWgt(0) -
             m_fluxReweightNu_Ref->GetVertErrorBand( m_fluxErrorName )->GetUnivWgt(0) ) > 1.0e-6 )
      {
        //std::cout << Form( "Propagating v+e constraint weights to flux error band of histogram \"%s\"",
        //                     h->GetName() )
        //            << std::endl;
        for ( unsigned int iuniv = 0; iuniv < m_nFluxUniverses; ++iuniv )
        {
          h->GetVertErrorBand( m_fluxErrorName )->SetUnivWgt( iuniv,
              m_fluxReweightNu_Ref->GetVertErrorBand( m_fluxErrorName )->GetUnivWgt( iuniv ) );
        }
      }
      else {
        //std::cout << "I'm inside PropagateNuEConstraintWeights, "
        //          << "but I didn't satisfy the condition to apply the "
        //          << "weights..."
        //          << std::endl;
      }
    }
  template void FluxReweighter::PropagateNuEConstraintWeights<MnvH1D>( MnvH1D* h );
  template void FluxReweighter::PropagateNuEConstraintWeights<MnvH2D>( MnvH2D* h );


  //============================================================================
  const char* FluxReweighter::g4NumiVersionString(EG4NumiVersion v)
  {
    switch(v) {
      case g4numiv5: // This also handles "eroica"
        return "g4numiv5";
      case g4numiv6:
        return "g4numiv6";
      default:
        std::cerr << "No g4numi version " << v << std::endl;
        std::exit(1);
    }
  }


  //======================================================================
  const char* FluxReweighter::fluxForSystematicsHistName(EFluxVersion fluxVersion)
  {
    int iPlaylist = LeOrMe(d_Playlist); // Adding in Amit's CVS commit 1.106-8 changes -David L.
    switch(fluxVersion){
      case lowNu:
        m_useStandardFlux = false;
        return "lowNuFlux_MnvH1D";
      case highNu:
        m_useStandardFlux = false;
        return "highNuFlux_MnvH1D";
      case ppfxDebug:
        m_useStandardFlux = false;
        return "ppfxDebugFlux_MnvH1D";
      case comboPPFXHighNu:
	m_useStandardFlux = false;
	return "comboPPFXHighNuFlux_MnvH1D";

      default:
        m_useStandardFlux = true;
	if (iPlaylist==3) return "MuonBar_Energy_MnvH1D"; // Adding in Amit's CVS commit 1.106-8 changes -David L.
	else return "Muon_Energy_MnvH1D";
    }
  }


  //======================================================================
  bool FluxReweighter::IsCustomFlux(EFluxVersion fluxVersion)
  {
    switch(fluxVersion){
      case lowNu:
        return true;
      case highNu:
        return true;
      case ppfxDebug:
        return true;
      case comboPPFXHighNu:
	return true;
      default:
        return false;
    }
  }


  //======================================================================
  const char* FluxReweighter::playlistString(EPlaylist p)
  {
    switch(p) {
      case minerva1:
        return "minerva1";
      case minerva5:
        return "minerva5";
      case minerva7:
        return "minerva1";
      case minerva9:
        return "minerva1";
      case minerva13:
        return "minerva13";
      case minerva2p2h:
        return "minerva13";
      case minervaLE_FHC:
        return "minervaLE-FHC";
      case minervame1A:
        return "minervame1D";
      case minervame1B:
        return "minervame1D";
      case minervame1C:
        return "minervame1D";
      case minervame1D:
        return "minervame1D";
      case minervame1E:
        return "minervame1D";
      case minervame1F:
        return "minervame1D";
      case minervame1N:
        return "minervame1N";
      case minervame1O:
        return "minervame1N";
      case minervame1P:
        return "minervame1N";
      case minervame1G:
        return "minervame1M";
      case minervame1L:
        return "minervame1M";
      case minervame1M:
        return "minervame1M";
      case minervame1D1M1NWeightedAve:
        return "minervame1D1M1NWeightedAve";
      case minervame5A:
        return "minervame5A";
      case minervame6A:
        return "minervame6A";
      case minervame6B:
        return "minervame6B";
      case minervame6C:
        return "minervame6C";
      case minervame6D:
        return "minervame6D";
      case minervame6E:
        return "minervame6E";
      case minervame6F:
        return "minervame6F";
      case minervame6G:
      	return "minervame6G";
      case minervame6H:
      	return "minervame6H";
      case minervame6I:
      	return "minervame6I";
      case minervame6J:
      	return "minervame6J";
      default:
        std::cerr << "No playlist enum " << p << std::endl;
        std::exit(1);
    }
  }


  //============================================================================
  int FluxReweighter::LeOrMe(EPlaylist p) {
    /*
     *
     * Given a playlist enum it returns an integer value
     * corresponding to Playlist of LE or ME.
     * All LE should return 1
     * ME will return 2
     * This is to keep a future provision if the
     * weight factors to be implimented playlist by playlist
     * or in a group set different integers other than 1, 2
     */
    switch(p) {
      case minerva1:
        return 1;
      case minerva5:
        return 1;
      case minerva13:
        return 1;
      case minervaLE_FHC:
        return 1;
      case minervame1A:
        return 2;
      case minervame1B:
        return 2;
      case minervame1C:
        return 2;
      case minervame1D:
        return 2;
      case minervame1E:
        return 2;
      case minervame1F:
        return 2;
      case minervame1G:
        return 2;
      case minervame1L:
        return 2;
      case minervame1M:
        return 2;
      case minervame1N:
        return 2;
      case minervame1O:
        return 2;
      case minervame1P:
        return 2;
      case minervame1D1M1NWeightedAve:
        return 2;
      case minervame5A:
        return 3;
      case minervame6A:
        return 3;
      case minervame6B:
      	return 3;
      case minervame6C:
      	return 3;
      case minervame6D:
      	return 3;
      case minervame6E:
        return 3;
      case minervame6F:
      	return 3;
      case minervame6G:
      	return 3;
      case minervame6H:
      	return 3;
      case minervame6I:
      	return 3;
      case minervame6J:
      	return 3;
      default:
        std::cerr << "No playlist enum " << p << std::endl;
        std::exit(1);
    }
  }



  //============================================================================
  FluxReweighter::EPlaylist FluxReweighter::GetPlaylistEnum(std::string& playlist)
  {
    // case-independent
    std::transform(playlist.begin(), playlist.end(), playlist.begin(), ::tolower);

    //std::cout << "Looking for a flux reweighter with the playlist "
    //          << playlist << std::endl;

    if ( playlist.compare("minerva1")== 0 )
      return minerva1;
    if( playlist.compare("minerva5")== 0 )
      return minerva5;
    if ( playlist.compare("minerva7")== 0 )
      return minerva1;
    if ( playlist.compare("minerva9")== 0 )
      return minerva1;
    if ( playlist.compare("minerva13")== 0 )
      return minerva13;
    if ( playlist.compare("minerva13c")== 0 )
      return minerva13;
    if ( playlist.compare("minerva13e")== 0 )
      return minerva13;
    if ( playlist.compare("minerva2p2h")== 0 )
      return minerva13;
    if( playlist.compare("minervame1a")== 0)
      return minervame1D;
    if( playlist.compare("minervame1b")== 0)
      return minervame1D;
    if( playlist.compare("minervame1c")== 0)
      return minervame1D;
    if( playlist.compare("minervame1d")== 0)
      return minervame1D;
    if( playlist.compare("minervame1e")== 0)
      return minervame1D;
    if( playlist.compare("minervame1f")== 0)
      return minervame1D;
    if( playlist.compare("minervame1g")== 0)
      return minervame1M;
    if( playlist.compare("minervame1l")== 0)
      return minervame1M;
    if( playlist.compare("minervame1m")== 0)
      return minervame1M;
    if( playlist.compare("minervame1n")== 0)
      return minervame1N;
    if( playlist.compare("minervame1o")== 0)
      return minervame1N;
    if( playlist.compare("minervame1p")== 0)
      return minervame1N;
    if ( playlist.compare("minervame1d1m1nweightedave")== 0)
      return minervame1D1M1NWeightedAve;
    if( playlist.compare("minervame5a")== 0)
      return minervame5A;
    if( playlist.compare("minervame6a")== 0)
      return minervame6A;
    if( playlist.compare("minervame6b")== 0)
      return minervame6A;
    if( playlist.compare("minervame6c")== 0)
      return minervame6A;
    if( playlist.compare("minervame6d")== 0)
      return minervame6A;
    if( playlist.compare("minervame6e")== 0)
      return minervame6A;
    if( playlist.compare("minervame6f")== 0)
      return minervame6A;
    if( playlist.compare("minervame6g")== 0)
      return minervame6A;
    if( playlist.compare("minervame6h")== 0)
      return minervame6A;
    if( playlist.compare("minervame6i")== 0)
      return minervame6A;
    if( playlist.compare("minervame6j")== 0)
      return minervame6A;
    std::cerr << "No playlist enum " << playlist << std::endl;
    std::exit(1);
  }


  //============================================================================
  const char* FluxReweighter::fluxVersionString(EFluxVersion f)
  {
    switch(f) {
      case gen1:
        return "gen1";
      case gen2thick:
        return "gen2thick";
      case gen2thin:
        return "gen2thin";
      case lowNu:
        return "lowNu";
      case highNu:
        return "highNu";
      case ppfxDebug:
        return "ppfxDebug";
      case comboPPFXHighNu:
        return "comboPPFXHighNu";
      default:
        std::cerr << "No flux version " << f << std::endl;
        std::exit(1);
    }
  }

  FluxReweighter& flux_reweighter(std::string plist, int nu_pdg,
                                  bool use_nuE_constraint, int n_flux_universes /*= 200*/) {
    static std::map<std::string, std::map<int, PlotUtils::FluxReweighter*> > FRW_map;
    std::transform(plist.begin(), plist.end(), plist.begin(), ::tolower);
    if (!FRW_map[plist][abs(nu_pdg)]) {
      FRW_map[plist][abs(nu_pdg)] = _flux_reweighter(plist, nu_pdg,
                                                     use_nuE_constraint,n_flux_universes);
    }
    return *(FRW_map[plist][abs(nu_pdg)]);
  }

} // end namespace PlotUtils

// Local Variables:
// c-basic-offset: 2
// End:
