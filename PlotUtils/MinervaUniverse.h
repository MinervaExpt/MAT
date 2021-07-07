#ifndef MINERVAUNIVERSE_H
#define MINERVAUNIVERSE_H

#include "BaseUniverse.h"
#include "GeantHadronSystematics.cxx"       // PlotUtils::weight_hadron (MnvHadronReweighter)
#include "FluxReweighter.h"                 // PlotUtils::flux_reweighter
#include "GenieSystematics.cxx"             // PlotUtils::IsNonResPi
#include "MinosMuonEfficiencyCorrection.h"  // PlotUtils::MinosMuonEfficiencyCorrection
#include "MnvNormalization.h"               // PlotUtils::MnvNormalizer
#include "MnvTuneSystematics.cxx"  // PlotUtils::Get<RPA/2p2h/NonResPi/LowQ2Pi>Weight
#include "NSFDefaults.h"
#include "PlotUtilsPhysicalConstants.h"
#include "TVector3.h"          // Needed by SystCalcs/TruthFunctions.h
#include "weightCoherentPi.h"  // PlotUtils::weight_coherent
#include "weight_fsi.h"        // PlotUtils::weight_fsi
#include "weightMK.h"          // PlotUtils::weight_mk

namespace PlotUtils {
class MinervaUniverse : public PlotUtils::BaseUniverse {
 public:
  MinervaUniverse(){PlotUtils::BaseUniverse();};
  // CTOR
  MinervaUniverse(PlotUtils::TreeWrapper* chw, double nsigma = 0)
      : PlotUtils::BaseUniverse(chw, nsigma){};

  //! Playlists
  static std::string GetPlaylist();
  static bool IsPlaylistME(std::string playlist);
  static bool SetPlaylist(std::string playlist);
  //!  what tree is this?
protected:
  static std::string m_treename;
public:
  static void SetTreeName(std::string treename);//{m_TreeName=treename;};
  static std::string GetTreeName();//{return m_TreeName;};
  
  //! MINERvA Software Processing
  static bool SetProcessingEroica();
  static bool IsProcessingNX();

  //! Horn current
  static bool isFHC() { return m_is_FHC; }
  static bool isRHC() { return m_is_RHC; }
  static bool SetHornCurrent(std::string& playlist);

  //! Analysis neutrino identity -- used for flux weights
  static int GetAnalysisNuPDG();
  static bool SetAnalysisNuPDG(int nu_pdg);

  // When this is on, NRP weights are applied to the CV weight
  static bool UseDeuteriumGeniePiTune();
  static bool SetDeuteriumGeniePiTune(bool use_tune);
  static bool UseNonResPiReweight();
  static bool SetNonResPiReweight(bool use_reweight);
  static bool UseZExpansionFaReweight();
  static bool SetZExpansionFaReweight(bool use_reweight);
  //
  static bool UseNuEConstraint();
  static bool SetNuEConstraint(bool use_constraint);

  //
  static int GetNFluxUniverses();
  static bool SetNFluxUniverses(int n_flux_universes);

  //! MnvHadronReweighter
  static bool SetMHRWeightFilename( std::string process_name, std::string directory = Form("%s/data/mhrwKineRenorm",getenv("PLOTUTILSROOT")) );
  static bool SetReadoutVolume( std::string readout_vol );
  static bool SetReadoutVolume( double minZ, double maxZ, double apothem );
  static bool SetMHRWeightNeutronCVReweight( bool useNeutronCVReweight );
  static bool SetMHRWeightElastics( bool useElastics );
  bool SetupMHRWeighter() const;

  // These should probably move, but I'm not sure where, yet.
  virtual double GetBeamAngleOffsetX() const;
  virtual double GetBeamAngleOffsetY() const;

  virtual double GetBatchPOT() const;

 protected:
  static bool _is_muon_momentum_cv_offset_set;  // = false;
  static bool _has_muon_error_been_thrown;      // = false;
  static double m_muon_momentum_cv_offset;      // = 0.0;
  
  //GeantHadronSystematics need these
  static bool m_mhrw_load_renorm;               // = true;
  static bool m_use_mhrw_neutronCV_reweight;    //
  static bool m_use_mhrw_elastics;
  static double m_readout_vol_minZ;
  static double m_readout_vol_maxZ;
  static double m_readout_vol_apothem;
  static std::string m_readout_volname;
  static std::string m_mhrw_process_name;
  static std::string m_mhrw_directory;

  
  
  //proton KE cut variables...
  static double m_true_proton_ke_cut_central; //= NSFDefaults::TrueProtonKECutCentral;
  static double m_reco_proton_ke_cut_central; // = NSFDefaults::RecoProtonKECutCentral;

 private:
  static std::string m_playlist;
  static std::string m_processing;
  
  static bool m_is_RHC;
  static bool m_is_FHC;
  static int m_analysis_nu_pdg;
  static bool m_use_nuE_constraint;
  static int m_n_flux_universes;
  static bool m_use_D2_pion_genie_tune;
  static bool m_use_nonResPi_reweight;
  static bool m_use_zExpansionFa_reweight;

};
}  // namespace PlotUtils
#endif  // MINERVAUNIVERSE
