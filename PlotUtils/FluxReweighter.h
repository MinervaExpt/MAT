#ifndef FLUXREWEIGHTER_H
#define FLUXREWEIGHTER_H

#include "TString.h"
#include "TSpline.h"
#include "TH1D.h"
#include "TMap.h"

#include <vector>
#include <map>

namespace PlotUtils
{
  class MnvH1D;
  class MnvH2D;
  class MnvHistoConstrainer;

    /*! \brief Calculate flux weights for events using histograms of
     *   the reweighted flux and the flux as generated
     *
     *   Example usage:
     *
     *   \code
     *       FluxReweighter* frw = new FluxReweighter(14, true,
     *                                                FluxReweighter::minerva13,
     *                                                FluxReweighter::gen2thin,
     *                                                FluxReweighter::g4numiv5);
     *
     *       // IMPORTANT: Use FluxReweighter to add the flux error band to your MnvH*Ds
     *       // This will propagate the v+e constraint weights to your flux error band
     *       frw->AddFluxErrorBand( hist );
     *
     *       // Loop over MC events
     *       for(int i=0; i<ch->GetEntries()/10; ++i){
     *           ch->GetEntry(i); // Get the event we want
     *           double Enu = ch->mc_incomingE * 1.0e-3; // true neutrino energy (GeV)
     *           int nuPDG = mc->mc_incoming; // neutrino PDG code
     *
     *           // Get the central value weight
     *           double cvWeight = frw->GetFluxCVWeight( Enu, nuPDG );
     *           // Fill your histogram
     *           hist->Fill( recoQ2, cvWeight );
     *
     *           // Use FluxReweighter to fill your histogram's flux error band
     *           frw->FillFluxErrorBand( hist, recoQ2, cvWeight, Enu, nuPDG );
     *           // Alternate method:
     *           // std::vector<double> weights = frw->GetFluxErrorWeights( enu, nuPDG );
     *           // hist->FillVertErrorBand( "Flux", recoQ2, weights, cvWeight ); // cvWeightFromMe = 1.0
     *
     *       }
     *   \endcode
     *
     *   NB: you just need one FluxReweighter per job, not one per
     *   event (which will be slow and might leak memory)
     *
     *   The histograms live in $PLOTUTILSROOT/data/flux with a
     *   hopefully obvious naming scheme. When a new flux or g4numi
     *   version is introduced, someone will have to update the
     *   relevant enum and the function that converts it to a string.
     *
     *   The histograms are produced by
     *   Ana/Flux/python/compute_flux.py, as driven by
     *   Ana/Flux/scripts/make_flux_histograms.sh. That in turn relies
     *   on ntuples being produced for the new flux in all the
     *   relevant playlists.
     */

  class FluxReweighter
  {
    public:

      //! Playlist for which there's a well-defined flux. Includes
      //! "LE" and "ME" as "playlists"
      enum EPlaylist
      {
        minerva1,
        minerva5,
        minerva7,
        minerva9,
        minerva13,
        minerva2p2h,
        minervaLE_FHC,
        minervame1A,
        minervame1B,
        minervame1C,
        minervame1D,
        minervame1E,
        minervame1F,
        minervame1G,
        minervame1L,
        minervame1M,
        minervame1N,
        minervame1O,
        minervame1P,
        minervame1D1M1NWeightedAve,
        minervame5A,
        minervame6A,
        minervame6B,
        minervame6C,
        minervame6D,
        minervame6E,
        minervame6F,
        minervame6G,
        minervame6H,
        minervame6I,
        minervame6J 
      };

      //! A version of the flux reweighting
      enum EFluxVersion
      {
        gen1,
        gen2thick,
        gen2thin,
        lowNu,
        highNu,
        ppfxDebug,
        comboPPFXHighNu,
      };

      //! A version of g4numi used for generating events
      enum EG4NumiVersion
      {
        g4numiv5=0,
        eroica=0, // Eroica is processed with g4numi v5
        g4numiv6=1
      };

      //! Create flux reweighter with specified flux
      //! parameters. We'll get neutrino and antineutrino histograms
      //! for nuPDG and -1*nuPDG
      FluxReweighter(int nuPDG,
                     bool applyNuEConstraint,
                     enum EPlaylist playlist,
                     enum EFluxVersion fluxVersion,
                     enum EG4NumiVersion g4NumiVersion,
                     int nUniverses = 200);

      //! Create flux reweighter with specified flux
      //! parameters. We'll get neutrino and antineutrino histograms
      //! for nuPDG and -1*nuPDG
      FluxReweighter(int nuPDG,
                     bool applyNuEConstraint,
                     std::string playlist_str,
                     enum EFluxVersion fluxVersion,
                     enum EG4NumiVersion g4NumiVersion,
                     int nUniverses = 200);

      //! Create flux reweighter with specified histograms. You
      //! probably don't want this unless you know what you're doing
      FluxReweighter(MnvH1D* fluxGenNu,
                     MnvH1D* fluxGenNubar,
                     MnvH1D* fluxReweightNu,
                     MnvH1D* fluxReweightNubar,
                     TSpline3* MELowNuDataMCRatioSpline,
                     bool applyNuEConstraint);


      virtual ~FluxReweighter();

      //! Get the central value flux weight.  Enu is true neutrino energy (GeV).
      virtual double GetFluxCVWeight( double Enu, int nuPDG );

      //! Get the flux weight correction for a systematic universe.
      virtual double GetSysUniFluxWeightCorrection(double Enu, int nuPDG,
                                                   std::string sys_name, int universe);

      //! Get the flux error weight (GeV).
      double GetMELowNuFluxWeight( double Enu );

      //! Get the flux error weights for all flux universes.  Enu is true neutrino energy (GeV).
      virtual std::vector<double> GetFluxErrorWeights( double Enu, int nuPDG );

      //! Get the flux error weight for a single flux universe.  Enu is true neutrino energy (GeV).
      virtual double GetFluxErrorWeight(double Enu, int nuPDG,
                                        unsigned int universe);

      //! Add the flux error band to a MnvH1D, MnvH2D
      //! If the v+e constraint is applied, the constraint weights are propagated to the flux error band
      template<class MnvHistoType>
	  void AddFluxErrorBand( MnvHistoType* h );
      
      template<class MnvHistoType>
	  bool CheckFluxErrorBand( MnvHistoType* h );
      
      template<class MnvHistoType>
	  void CheckAndFixFluxErrorBand( MnvHistoType* h );

      //! Reduce the 1000 universes of the flux to the first 100
      void TruncateNumberOfFluxUniverses( MnvH1D* h , int nUniverses );

      //! Fill the flux error band of an MnvH1D, MnvH2D
      //! The v+e constraint weights are propagated to the input histogram's flux error band as needed
      void FillFluxErrorBand( MnvH1D* h, double val, double cvweight, double Enu, int nuPDG );

      void FillFluxErrorBand( MnvH2D* h, double xval, double yval, double cvweight, double Enu, int nuPDG );

      MnvH1D* GetFluxGenerated(int nuPDG) {return nuPDG>0 ? m_fluxGenNu : m_fluxGenNubar;}
      MnvH1D* GetFluxReweighted(int nuPDG) {return nuPDG>0 ? m_fluxReweightNu : m_fluxReweightNubar;}

      //! Get a rebinned version of the weighted flux using the input histogram binning. Returns an MnvH1D with all flux universes populated and any other error bands filled with CV
      MnvH1D* GetRebinnedFluxGenerated(int nuPDG, MnvH1D* template_hist);
      MnvH1D* GetRebinnedFluxReweighted(int nuPDG, MnvH1D* template_hist);
      MnvH1D* GetRebinnedFluxReweighted_FromInputFlux(MnvH1D* input_flux, MnvH1D* template_hist);

      //! Get a rebinned version of the unweighted flux using the input histogram binning. Returns an MnvH1D with all flux universes populated and any other error bands filled with CV
      //	MnvH1D* GetRebinnedFluxGenerated(int nuPDg, MnvH1D* template_hist);

      //! Get a histogram filled with integrals of weighted flux universes and cv values for other errors
      //! RDF 05-2020: Added support for special flux systematic universes
      template<class MnvHistoType>
      MnvHistoType* GetIntegratedFluxReweighted(int nuPDG, MnvHistoType* template_hist, double min_energy, double max_energy,bool use_correlations = true);
      //MnvH1D* GetIntegratedFluxReweighted(int nuPDG, MnvH1D* template_hist, double min_energy, double max_energy, bool use_correlations = true);
      
      template<class MnvHistoType>
      MnvHistoType* GetIntegratedFluxReweighted_FromInputFlux(MnvH1D* input_flux, MnvHistoType* template_hist, double min_energy, double max_energy);
      //MnvH1D* GetIntegratedFluxReweighted_FromInputFlux(MnvH1D* input_flux, MnvH1D* template_hist, double min_energy, double max_energy);
      
      MnvH1D* GetIntegratedTargetFlux(int nuPDG, std::string tar_mat, MnvH1D* template_hist, double min_energy, double max_energy, std::string project_dir = "targets_2345_temp");
      MnvH2D* GetIntegratedTargetFlux(int nuPDG, std::string tar_mat, MnvH2D* template_hist, double min_energy, double max_energy, std::string project_dir = "targets_2345_temp");

      //Get flux developed for the targets.  Currently only accepts nu-e constrained, neutrino mode only)
      MnvH1D* GetTargetFluxMnvH1D(int nuPDG, std::string tar_mat, std::string project_dir = "targets_2345_temp");

      //Add eff corrected tracker daisy histograms together so that the flux matches that of the specified target
      //The map key is the daisy number (the method to get this is in TargetUtils::GetDaisyPetal)
      MnvH1D* GetReweightedDaisySum(int nuPDG, std::string tar_mat, std::map<int, MnvH1D*> daisy_eff_hists, std::string project_dir = "targets_2345_temp");

      MnvH2D* GetReweightedDaisySum(int nuPDG, std::string tar_mat, std::map<int, MnvH2D*> daisy_eff_hists, std::string project_dir = "targets_2345_temp");

      //Get the parameter files for daisy reweight 
      MnvH1D* GetDaisyParamMnvH1D(int nuPDG, std::string tar_mat, std::string project_dir = "targets_2345_temp");

      //! Get a histogram filled with unweighted flux universes and cv values for other errors
      //MnvH1D* GetIntegratedFluxGenerated(int nuPDG, MnvH1D* template_hist);

      EPlaylist GetPlaylistEnum(std::string& playlist);


    protected:

      //! Get a singleton MnvHistoConstrainer
      MnvHistoConstrainer& Constrainer();

      //! Get MnvH1D with name histname from file filename
      MnvH1D* GetMnvH1D(TString filename, TString histname);

      //! Get the MnvH1D of the reweighted flux for the given parameters
      MnvH1D* GetFluxMnvH1D(int nuPDG,
                            enum EPlaylist playlist,
                            enum EFluxVersion fluxVersion,
                            enum EG4NumiVersion g4NumiVersion,
                            bool useGen = false /*use reweighted flux by default*/);

      MnvH1D* GetMELowNuMnvH1D();

      //! Fetch the MnVH1D with alternate flux predictions for systematic universes
      void SetFluxSysMnvH1D(int nuPDG, enum EFluxVersion fluxVersion);

      TSpline3* GetSpline(MnvH1D*);

      void RebinFluxHist(TH1D* h_flux, TH1D*&h_rebinned_flux);

      template<class MnvHistoType>
        bool IsFluxErrorBandOK( MnvHistoType* h );

      template<class MnvHistoType>
        void PropagateNuEConstraintWeights( MnvHistoType* h );

      const char* g4NumiVersionString(EG4NumiVersion);
      const char* fluxForSystematicsHistName(EFluxVersion);
      bool IsCustomFlux(EFluxVersion);
      const char* fluxVersionString(EFluxVersion);
      const char* playlistString(EPlaylist);
      EPlaylist PlaylistString;
      EPlaylist d_Playlist;

      //! for flux playlist identification 
      int LeOrMe(EPlaylist);  

      //! The histogram of the flux as generated, for neutrinos
      MnvH1D* m_fluxGenNu;
      //! The histogram of the flux as generated, for antineutrinos
      MnvH1D* m_fluxGenNubar;
      //! The histogram of the flux to be reweighted to, for neutrinos
      MnvH1D* m_fluxReweightNu;
      //! plus a copy that won't be truncated
      MnvH1D* m_fluxReweightNu_Ref;
      //! The histogram of the flux to be reweighted to, for antineutrinos,
      MnvH1D* m_fluxReweightNubar;
      //! plus a copy that won't be truncated
      //MnvH1D* m_fluxReweightNubar_Ref;
      //! The histogram of the low nu reweight from data/MC.  Valid for ME only
      MnvH1D* m_reweightMELowNuDataMC;
      //! The map of systematic universe fluxes.
      std::map<std::string,int> m_fluxSysMapNu;
      //! This will hold the MnvH1D with alternate flux predictions for systematic universes
      MnvH1D* m_fluxSystematicsMnvH1D;
      //! This will indicate whether the standard (PPFX) flux is being used or one of
      //! low/high-nu flux is being used
      bool m_useStandardFlux;


      TSpline3* m_MELowNuDataMCRatioSpline;

      //! Toggle for applying the v+e flux constraint
      bool m_applyNuEConstraint;
      //! Number of flux universes
      unsigned int m_nFluxUniverses;
      //! Total flux error band name
      std::string m_fluxErrorName;
  public:
    void AddFluxErrorBand(PlotUtils::MnvH1D* h) {
      AddFluxErrorBand<PlotUtils::MnvH1D>(h);
    }
    void AddFluxErrorBand(PlotUtils::MnvH2D* h) {
      AddFluxErrorBand<PlotUtils::MnvH2D>(h);
    }
    MnvH1D* GetIntegratedFluxReweighted( int nuPDG,
                                         MnvH1D* template_hist,
                                         double min_energy,
                                         double max_energy,
                                         bool use_muon_correlations) {
      return GetIntegratedFluxReweighted<MnvH1D>(nuPDG,template_hist,min_energy,max_energy,use_muon_correlations);
    }
    MnvH2D* GetIntegratedFluxReweighted( int nuPDG,
                                         MnvH2D* template_hist,
                                         double min_energy,
                                         double max_energy,
                                         bool use_muon_correlations) {
      return GetIntegratedFluxReweighted<MnvH2D>(nuPDG,template_hist,min_energy,max_energy,use_muon_correlations);
    }
  };

  FluxReweighter& flux_reweighter(std::string plist, int nu_pdg, bool use_nuE_constraint, int n_flux_universes = 200);
};


#endif


// Local Variables:
// c-basic-offset: 4
// End:
