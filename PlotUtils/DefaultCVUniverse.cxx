
#ifndef DEFAULTCVUNIVERSE_cxx
#define DEFAULTCVUNIVERSE_cxx

#include "DefaultCVUniverse.h"
#include "MinosMuonEfficiencyCorrection.h" 
#include "MnvNormalization.h"
#include "PlotUtilsPhysicalConstants.h"
#include "FluxSystematics.cxx" // PlotUtils::flux_reweighter
#include "MnvTuneSystematics.cxx" // PlotUtils::Get<RPA/2p2h/NonResPi/LowQ2Pi>Weight
#include "GenieSystematics.cxx" //PlotUtils::IsNonResPi
#include "TVector3.h" // Using this to get theta from a 3d vector
#include "weight_fsi.h"
#include "weightCoherentPi.h"
#include "weightMK.h"
#include <algorithm>  //::to_lower

using namespace PlotUtils;

//Flags for outputtting warnings when the user tries to reset a variable
namespace {
  bool _is_nue_constraint_set          = false;
  bool _is_nonResPi_reweight_set       = false;
  bool _is_zExpansionFa_reweight_set   = false;
  bool _is_playlist_set                = false;
  bool _is_horn_current_set            = false;
  bool _is_analysis_nu_pdg_set         = false;
  bool _is_muon_momentum_cv_offset_set = false;
  bool _has_muon_error_been_thrown     = false;
}

 

bool DefaultCVUniverse::m_is_FHC                    = false;
bool DefaultCVUniverse::m_is_RHC                    = false;
bool DefaultCVUniverse::m_use_nuE_constraint        = false;
bool DefaultCVUniverse::m_use_nonResPi_reweight     = false;
bool DefaultCVUniverse::m_use_zExpansionFa_reweight = false;
bool DefaultCVUniverse::m_is_truth                  = false;
std::string DefaultCVUniverse::m_playlist           = "minervame1a";
std::string DefaultCVUniverse::m_processing         = "NX";
int DefaultCVUniverse::m_analysis_nu_pdg            = 14;
int DefaultCVUniverse::m_n_flux_universes           = 200;
double DefaultCVUniverse::m_muon_momentum_cv_offset = 0.0;
double DefaultCVUniverse::m_true_proton_ke_cut_central = NSFDefaults::TrueProtonKECutCentral;
double DefaultCVUniverse::m_reco_proton_ke_cut_central = NSFDefaults::RecoProtonKECutCentral;
//Recoil Energy that needs to be reimplemented by USER in their CVUnivrse//


// CTOR
  DefaultCVUniverse::DefaultCVUniverse(PlotUtils::TreeWrapper* chw,
                                       double nsigma )
    : m_chw(chw), m_nsigma(nsigma), m_entry(-1)
  {}


// Get Weights
  double DefaultCVUniverse::GetGenieWeight() const {
    return UseNonResPiReweight() && PlotUtils::IsNonResPi(*this) ? PlotUtils::kNonResPiWeight : 1.;
  }


  double DefaultCVUniverse::GetFluxAndCVWeight(double Enu /*GeV*/,
                                               int nu_pdg) const {
    if (Enu == -99.)   Enu    = GetDouble("mc_incomingE")*1e-3;
    // For LE, electron-neutrino fluxes aren't available, so we should force nu_pdg to have
    // absolute value 14. The +/- doesn't matter, because FRW sets up both, anyways. This
    // may change in the future if electron-neutrino LE fluxes become available.
    if (nu_pdg == -99 && !IsPlaylistME(DefaultCVUniverse::GetPlaylist())) nu_pdg = 14;
    else if (nu_pdg == -99) nu_pdg = GetInt("mc_incoming");
    return PlotUtils::flux_reweighter(
        GetPlaylist(), nu_pdg,UseNuEConstraint(),GetNFluxUniverses()).GetFluxCVWeight(Enu, nu_pdg);
  }

  
  double DefaultCVUniverse::GetRPAWeight( ) const {
    const int variation = 0; // CV
    return PlotUtils::GetRPAWeight(*this, DefaultCVUniverse::Getq0True()/1000 /* GeV */, DefaultCVUniverse::Getq3True()/1000 /* GeV */, variation,DefaultCVUniverse::IsProcessingNX());
  }


  double DefaultCVUniverse::GetLowRecoil2p2hWeight( ) const {
    const int variation = 0; // CV
    return PlotUtils::GetLowRecoil2p2hWeight(*this, DefaultCVUniverse::Getq0True()/1000 /* GeV */, DefaultCVUniverse::Getq3True()/1000 /* GeV */, variation);
  }


  double DefaultCVUniverse::GetLowQ2PiWeight( std::string channel ) const {
    int variation = 0; // CV
    if (!PlotUtils::IsCCRes(*this))
      return 1.;
    else
      return PlotUtils::weight_lowq2pi().getWeight(DefaultCVUniverse::GetQ2True()*1e-6 /*GeV^2*/, channel, variation);
  }


  double DefaultCVUniverse::GetCoherentPiWeight(double thpi_true /*deg*/,
                                                double tpi_true /*GeV*/) const {
    if (GetInt("mc_intType") != 4) return 1.;
    assert(tpi_true > 0. && "GetCoherentPiWeight failed with tpi < 0.");
    assert(thpi_true > 0. && "GetCoherentPiWeight failed with thpi < 0.");
    return PlotUtils::weight_coherent().get_combined_weight(thpi_true, tpi_true);
  }

  double DefaultCVUniverse::GetMKWeight() const { 
    return PlotUtils::weight_mk().getWeight(m_chw, m_entry);
  }


  double DefaultCVUniverse::GetMinosEfficiencyWeight() const {
    
    if( m_is_truth ) return 1. ;//No efficiency reweigting for truth events 
    if( !isFHC() && !isRHC() ) return 1.;//LE or nonstandard playlist.  Don't have away of dealing with this yet
    if(IsPlaylistME( DefaultCVUniverse::GetPlaylist() )){ //The correction factors are different between ME and LE
       double pmu = GetPmuMinos()/1000; //GetCorrection expects GeV
       return PlotUtils::MinosMuonEfficiencyCorrection::Get( isFHC() ).GetCorrection(pmu,GetBatchPOT(),isFHC()); 
    }
    else{ //Assume if not ME, then it's LE
      double pmu = GetPmuMinos(); //MnVnormalizer GetCorrection expects MeV	
      static PlotUtils::MnvNormalizer mnvNormalizer = PlotUtils::MnvNormalizer("Eroica",GetPlaylist());
      return mnvNormalizer.GetCorrection(pmu);
    }
  }

  double DefaultCVUniverse::GetTargetMassWeight() const {
    return 1.0;  
  }

  double DefaultCVUniverse::GetFSIWeight( int iWeight ) const {
    static PlotUtils::weight_fsi weight_FSI;
    weight_FSI.UseTrackingThreshold();
    weight_FSI.calcWeights(GetInt("mc_incoming"),
                           GetInt("mc_primaryLepton"),
                           GetInt("mc_charm"),
                           GetInt("mc_intType"),
                           GetInt("mc_targetA"),
                           GetInt("mc_targetZ"),
                           GetInt("mc_resID"),
                           GetInt("mc_er_nPart"),
                           GetVecInt("mc_er_ID"),
                           GetVecInt("mc_er_status"),
                           GetVecInt("mc_er_FD"),
                           GetVecInt("mc_er_LD"),
                           GetVecInt("mc_er_mother"),
                           GetVecDouble("mc_er_Px"),
                           GetVecDouble("mc_er_Py"),
                           GetVecDouble("mc_er_Pz"),
                           GetVecDouble("mc_er_E"),
                           m_entry);
    if( iWeight == 0 ) return weight_FSI.GetElasticWeight(1)*weight_FSI.GetAbsorptionWeight(); 
    if( iWeight == 1 ) return weight_FSI.GetElasticWeight(1); 
    if( iWeight == 2 ) return weight_FSI.GetAbsorptionWeight(); 
    return 1;
  }


  double DefaultCVUniverse::GetGeantHadronWeight() const {
    //TODO Need to put the neutron CV reweight here.  If that's the case, we need to instantiate a MnvHadronReweighter here, WITH the correct fiducial volume
    return 1;
  }
 
  double DefaultCVUniverse::GetMichelEfficiencyWeight() const {
    return 1;
  }
 
  bool DefaultCVUniverse::UseNuEConstraint() {
    if (!_is_nue_constraint_set) {
      std::cout << "WARNING: YOU ARE USING DEFAULT CVU NUE CONSTRAINT (false). "
                << "PLEASE SET THE VARIABLE BY "
                << "CVU::SetNuEConstraint(use_nue_constraint) "
                << "BEFORE USING UNIVERSES." << std::endl;
    }
    return m_use_nuE_constraint; 
  }


  bool DefaultCVUniverse::SetNuEConstraint(bool use_constraint) {
    if (_is_nue_constraint_set) {
      std::cout << "WARNING: YOU ATTEMPTED SETTING NUE CONSTRAINT A SECOND TIME. "
                << "THIS IS NOT ALLOWED FOR CONSISTENCY." << std::endl;
      return false;
    }
    else {
      m_use_nuE_constraint = use_constraint;
      _is_nue_constraint_set = true;
    }
    return true;
  }

  bool DefaultCVUniverse::UseNonResPiReweight() {
    if (!_is_nonResPi_reweight_set) {
      std::cout << "WARNING: YOU ARE USING DEFAULT NONRESPI REWEIGHT (false). "
                << "PLEASE SET THE VARIABLE BY "
                << "CVU::SetNonResPiReweight(use_nonResPi_reweight) "
                << "BEFORE USING UNIVERSES." << std::endl;
    }
    return m_use_nonResPi_reweight; 
  }


  bool DefaultCVUniverse::SetNonResPiReweight(bool use_reweight) {
    if (_is_nonResPi_reweight_set) {
      std::cout << "WARNING: YOU ATTEMPTED SETTING NONRESPI REWEIGHT A SECOND TIME. "
                << "THIS IS NOT ALLOWED FOR CONSISTENCY." << std::endl;
      return false;
    }
    else {
      m_use_nonResPi_reweight = use_reweight;
      _is_nonResPi_reweight_set = true;
    }
    return true;
  }

  bool DefaultCVUniverse::UseZExpansionFaReweight() {
    // Debug message that can be activated if/when this feature becomes more widely used
    //if (!_is_zExpansionFa_reweight_set) {
    //  std::cout << "WARNING: YOU ARE USING DEFAULT ZEXPANSION REWEIGHT (false). "
    //            << "PLEASE SET THE VARIABLE BY "
    //            << "CVU::SetZExpansionFaReweight(use_zExpansionFa_reweight) "
    //            << "BEFORE USING UNIVERSES." << std::endl;
    //}
    return m_use_zExpansionFa_reweight; 
  }


  bool DefaultCVUniverse::SetZExpansionFaReweight(bool use_reweight) {
    if (_is_zExpansionFa_reweight_set) {
      std::cout << "WARNING: YOU ATTEMPTED SETTING ZEXPANSION REWEIGHT A SECOND TIME. "
                << "THIS IS NOT ALLOWED FOR CONSISTENCY." << std::endl;
      return false;
    }
    else {
      m_use_zExpansionFa_reweight = use_reweight;
      _is_zExpansionFa_reweight_set = true;
    }
    return true;
  }
  std::string DefaultCVUniverse::GetPlaylist() {
    if (!_is_playlist_set) {
      std::cout << "WARNING: YOU ARE USING DEFAULT CVU PLAYLIST (minervame1a). "
                << "PLEASE SET THE VARIABLE BY CVU::SetPlaylist(playlist) "
                << "BEFORE USING UNIVERSES." <<std::endl;
    }
    return m_playlist; 
  }


  bool DefaultCVUniverse::IsPlaylistME( std::string playlist ) {
    if ( playlist.find("me") == std::string::npos &&
         playlist.find("ME") == std::string::npos ) {
      return false;
    }
    else {
      return true;
    }
  }


  bool DefaultCVUniverse::SetPlaylist( std::string playlist ) {
    if (_is_playlist_set) {
      std::cout << "WARNING: YOU ATTEMPTED SETTING THE PLAYLIST A SECOND TIME. "
                << "THIS IS NOT ALLOWED FOR CONSISTENCY." << std::endl;
      return false;
    }
    else {
      m_playlist       = playlist;
      _is_playlist_set = true;
      if( !SetHornCurrent( playlist ) ) return false;
    }
    return true;
  }


  bool DefaultCVUniverse::SetProcessingEroica() {
    if (m_processing == "Eroica") {
      std::cout << "WARNING: THE PROCESSING IS ALREADY SET AS EROICA. "
                << "EXECUTING THIS AGAIN IS NOT ALLOWED FOR CONSISTENCY." << std::endl;
      return false;
    }
    else {
      m_processing = "Eroica";
    }
    return true;
  }


  bool DefaultCVUniverse::IsProcessingNX() {
    if ( m_processing == "NX" ){
      return true;
    }
    else if( m_processing == "Eroica"){
      return false;
    }
    else{
      std::cout << "WARNING: YOU HAVE AN UNRECOGNIZED PROCESSING SET!" << std::endl;
      return false;
    }
  }

  bool DefaultCVUniverse::SetHornCurrent( std::string& playlist ) {
    
    if (_is_horn_current_set) {
      std::cout << "WARNING: YOU ATTEMPTED SETTING THE HORN CURRENT A SECOND TIME. "
                << "You really shouldn't be trying to do this anyway"
                << "Horn current is set when you input a playlist"
                << "THIS IS NOT ALLOWED FOR CONSISTENCY." << std::endl;
      return false;
    }
    else {
      // case-independent
      std::transform(playlist.begin(), playlist.end(), playlist.begin(), ::tolower);
      if( playlist.compare("minerva1")    == 0 ||
          playlist.compare("minerva7")    == 0 ||
          playlist.compare("minerva9")    == 0 ||
          playlist.compare("minerva13")   == 0 ||
          playlist.compare("minerva13c")  == 0 ||
          playlist.compare("minerva13e")  == 0 ||
          playlist.compare("minervame1a") == 0 ||
          playlist.compare("minervame1b") == 0 ||
          playlist.compare("minervame1c") == 0 ||
          playlist.compare("minervame1d") == 0 ||
          playlist.compare("minervame1e") == 0 ||
          playlist.compare("minervame1f") == 0 ||
          playlist.compare("minervame1g") == 0 ||
          playlist.compare("minervame1l") == 0 ||
          playlist.compare("minervame1m") == 0 ||
          playlist.compare("minervame1n") == 0 ||
          playlist.compare("minervame1o") == 0 ||
          playlist.compare("minervame1p") == 0  ) m_is_FHC = true;

      if( playlist.compare("minerva5")    == 0 ||
          playlist.compare("minervame5a") == 0 ||
          playlist.compare("minervame6a") == 0 ||
          playlist.compare("minervame6b") == 0 ||
          playlist.compare("minervame6c") == 0 ||
          playlist.compare("minervame6d") == 0 ||
          playlist.compare("minervame6e") == 0 ||
          playlist.compare("minervame6g") == 0  ) m_is_RHC = true;

      if( !m_is_FHC && !m_is_RHC ) std::cout<<"You have entered a nonstandard playlist. I can't `SetHornCurrent`, and therefore I can't `SetPlaylist`!"<<std::endl;
      _is_horn_current_set = true;
    }
    return true;
  }


  int DefaultCVUniverse::GetAnalysisNuPDG() {
    if (!_is_analysis_nu_pdg_set) {
      std::cout << "WARNING: YOU ARE USING DEFAULT CVU NU PDG (14). "
                << "PLEASE SET THE VARIABLE BY CVU::SetAnalysisNuPDG(pdg) "
                << "BEFORE USING UNIVERSES." << std::endl;
    }
    return m_analysis_nu_pdg;
  }


  bool DefaultCVUniverse::SetAnalysisNuPDG( int nu_pdg ) {
    if (_is_analysis_nu_pdg_set) {
      std::cout << "WARNING: YOU ATTEMPTED SETTING THE NU PDG A SECOND TIME. THIS IS NOT ALLOWED FOR CONSISTENCY." << std::endl;
      return false;
    } else {
      m_analysis_nu_pdg = nu_pdg;
      _is_analysis_nu_pdg_set = true;
    }
    return true;
  }


  // This functionality should *not* be used by standard users. This will allow data 
  // preservation users in the future to offset the muon energy scale in data when it 
  // is no longer feasible to reproduce ntuples. So, the warning is inverted in that
  // the user will be warned only if they _do_ use the method to specify an offset.
  // n.b. 'offset' refers to CV only, and 'shift' refers to systematic universes.
  double DefaultCVUniverse::GetCVMuonMomentumOffset() {
    if (_is_muon_momentum_cv_offset_set && !_has_muon_error_been_thrown) {
      std::cout << "WARNING: YOU HAVE SET A NONZERO CV MUON MOMENTUM SHIFT. "
                << "IF YOU DON'T KNOW WHAT THIS MEANS, "
                << "YOU SHOULDN'T BE DOING THIS! AND PROBABLY ROB IS THE ONLY "
                << "PERSON WHO SHOULD BE DOING THIS." <<std::endl;
      _has_muon_error_been_thrown = true;
    }
    return m_muon_momentum_cv_offset;
  }


  bool DefaultCVUniverse::SetCVMuonMomentumOffset( double muon_momentum_cv_offset ) {
    if (_is_muon_momentum_cv_offset_set) {
      std::cout << "WARNING: YOU ATTEMPTED SETTING THE MUON MOMENTUM CV SHIFT "
                << "A SECOND TIME. THIS IS NOT ALLOWED FOR CONSISTENCY." 
                << std::endl;
      return false;
    }
    else {
      m_muon_momentum_cv_offset = muon_momentum_cv_offset;
      _is_muon_momentum_cv_offset_set = true;
    }
    return true;
  }

  int DefaultCVUniverse::GetNFluxUniverses() {
    return m_n_flux_universes; 
  }

  bool DefaultCVUniverse::SetNFluxUniverses( int n_flux_universes ) {
    m_n_flux_universes = n_flux_universes;
    return true;
  }

// Specific-use-case Getters 
  double DefaultCVUniverse::GetBatchPOT() const {
    double batch_pot    = 0;

    // Anushree bash pot calculation
    double numi_pot = GetDouble("numi_pot");
    bool is_me = DefaultCVUniverse::IsPlaylistME( DefaultCVUniverse::GetPlaylist() );
    int batch_structure;
    int reco_vertex_batch;
    if (is_me) {
      batch_structure = GetInt("batch_structure");
      reco_vertex_batch = GetInt("reco_vertex_batch");
    }
    else {
      batch_structure = 0;
    }

    if (batch_structure == 0 ) {
      batch_pot = numi_pot / 6. ;
    }
    if (batch_structure == 1 ) {
      if (reco_vertex_batch<3) { batch_pot = numi_pot / 4. ;}
      else { batch_pot = numi_pot / 8. ;}
    }
    if (batch_structure == 2 ) {
      if (reco_vertex_batch<5) { batch_pot = numi_pot / 5. ;}
      else { batch_pot = numi_pot / 10. ;}
    }
    if (batch_structure == 3) {
      batch_pot = numi_pot / 6. ;
    }
    if (batch_structure == -1 ) {
      //std::cout << "There isn't a batch structure for this gate, "
      //          << "so I'm not going to guess the batch POT if the batch structure is 0. Entry = "
      //          << m_entry << std::endl;
      batch_pot = numi_pot / 6. ;
    }
    return batch_pot;    
  }

  double DefaultCVUniverse::GetEnuTrue() const { 
    return GetDouble("mc_incomingE"); 
  }

  double DefaultCVUniverse::GetElepTrue() const { 
    return GetVecElem("mc_primFSLepton",3); 
  }

  double DefaultCVUniverse::GetPlepTrue() const { 
    TVector3 p3lep( GetVecElem("mc_primFSLepton",0), GetVecElem("mc_primFSLepton",1), GetVecElem("mc_primFSLepton",2) );
    return p3lep.Mag();  
  }

  double DefaultCVUniverse::GetThetalepTrue() const {
    TVector3 p3lep( GetVecElem("mc_primFSLepton",0), GetVecElem("mc_primFSLepton",1), GetVecElem("mc_primFSLepton",2) );
    p3lep.RotateX(MinervaUnits::numi_beam_angle_rad);
    return p3lep.Theta();
  }

  double DefaultCVUniverse::GetPhilepTrue() const {
    TVector3 p3lep( GetVecElem("mc_primFSLepton",0), GetVecElem("mc_primFSLepton",1), GetVecElem("mc_primFSLepton",2) );
    p3lep.RotateX(MinervaUnits::numi_beam_angle_rad);
    return p3lep.Phi();
  }

  double DefaultCVUniverse::GetQ2True() const {
    return GetDouble("mc_Q2");
  }

  double DefaultCVUniverse::Getq0True() const { 
    return calcq0(DefaultCVUniverse::GetEnuTrue(), DefaultCVUniverse::GetElepTrue()); 
  }

  double DefaultCVUniverse::Getq3True() const { 
    return calcq3(DefaultCVUniverse::GetQ2True(),  DefaultCVUniverse::GetEnuTrue(), DefaultCVUniverse::GetElepTrue());
  }

  int DefaultCVUniverse::GetTargetZTrue() const { 
    return GetInt("mc_targetZ");
  }

  double DefaultCVUniverse::GetVertexZTrue() const { 
    return GetVecElem("mc_vtx",2);
  }

  //Recoil Energy
  //  Analyzer should specify what their recoil energy is (in MeV) in CVUniverse)
  double DefaultCVUniverse::GetRecoilEnergy() const {
     return GetCalRecoilEnergy() + GetNonCalRecoilEnergy();
  }

  double DefaultCVUniverse::GetCalRecoilEnergy() const {
    std::cout<<"GetCalRecoilEnergy() should be implemented in CVUniverse to"<<std::endl;
    std::cout<<"     return all recoil energy found calorimetrically (MeV)"<<std::endl;
    return -999.99;
  }
  
  double DefaultCVUniverse::GetNonCalRecoilEnergy() const {
    std::cout<<"GetNonCalRecoilEnergy() should be implemented in CVUniverse to"<<std::endl;
    std::cout<<"     return all recoil energy not found calorimetrically like dEdX (MeV)"<<std::endl;
    return -999.99;
  }
  
  // Vertex Z position -- required for muon momentum systematics with recalculation of shifts
  double DefaultCVUniverse::GetVertexZ() const {
    double vertex_Z = GetVecElem("vtx",2);
    return vertex_Z;
  }

// set and get the proton KE cut central value for the DefaultCVUniverse, the variation is set when systematics map is set up
  
  bool DefaultCVUniverse::SetTrueProtonKECutCentral(double cut=NSFDefaults::TrueProtonKECutCentral){
    m_true_proton_ke_cut_central = cut;
    return true;
  }

  double DefaultCVUniverse::GetTrueProtonKECutCentral(){
   return m_true_proton_ke_cut_central;
 }

// the method you use to get the cut value which will be varied by TrueProtonKECutSystematics

  double DefaultCVUniverse::GetTrueProtonKECut() const{
   return m_true_proton_ke_cut_central;
 }

  bool DefaultCVUniverse::SetRecoProtonKECutCentral(double cut=NSFDefaults::RecoProtonKECutCentral){
    m_reco_proton_ke_cut_central = cut;
    return true;
  }

  double DefaultCVUniverse::GetRecoProtonKECutCentral(){
   return m_reco_proton_ke_cut_central;
 }

// the method you use to get the cut value which will be varied by TrueProtonKECutSystematics

  double DefaultCVUniverse::GetRecoProtonKECut() const{
   return m_reco_proton_ke_cut_central;
 }


  double DefaultCVUniverse::GetEmu() const {
    double m_mu = MinervaUnits::M_mu;
    double pmu  = GetPmu();
    return sqrt(m_mu*m_mu + pmu*pmu);
  }


  double DefaultCVUniverse::GetPmu() const {

    double Pmu_minerva = GetPmuMinerva();
    double Pmu_minos = GetPmuMinos();

    return Pmu_minerva + Pmu_minos;

  }

  // Get total muon momentum w/o CV shift
  double DefaultCVUniverse::GetPmu_nominal() const {
    
   const std::string lepton_branch = GetAnaToolName() + "_leptonE"; // AnaTool 
    //const std::string lepton_branch = "Muon_leptonE";
    double px = GetVecElem(lepton_branch.c_str(), 0);
    double py = GetVecElem(lepton_branch.c_str(), 1);
    double pz = GetVecElem(lepton_branch.c_str(), 2);
    double total_p_nominal = sqrt(px*px + py*py + pz*pz);
   
    return total_p_nominal;
 
  }

  double DefaultCVUniverse::GetPmuMinerva() const {
    // this method can be overridden by systematics universes.
    // note that it includes a mechanism for introducing an overall 
    // multiplicative factor
    
    double total_p_nominal = GetPmu_nominal();
    double minos_p_nominal = GetPmuMinos_nominal();
    double minerva_p = total_p_nominal - minos_p_nominal;

    return minerva_p;

  }

  double DefaultCVUniverse::GetPmuMinos() const {
    // this method can be overridden by systematics universes.
    // note that it includes a mechanism for introducing an overall 
    // multiplicative factor

    double minos_p_nominal = GetPmuMinos_nominal();
    
    // this CV offset is set by CVUniverse.SetCVMuonMomentumOffset method
    double momentum_offset = GetCVMuonMomentumOffset();
    if (momentum_offset == 0.0) {
      return minos_p_nominal;
    }
    else {
      double minos_p_shifted = minos_p_nominal * (1.+momentum_offset);
      return minos_p_shifted;
    }
  }

  double DefaultCVUniverse::GetPmuMinos_nominal() const {

//const std::string minos_momentum_branch = "Muon_minos_trk_p";
    const std::string minos_momentum_branch = GetAnaToolName() + "_minos_trk_p";
    double minos_p_nominal = GetDouble(minos_momentum_branch.c_str());

    return minos_p_nominal;

  }

//! override these to implement beam angle universes

  double DefaultCVUniverse::GetBeamAngleOffsetX() const{
    return NSFDefaults::beamThetaX_Central;
  }

  double DefaultCVUniverse::GetBeamAngleOffsetY() const{
    return NSFDefaults::beamThetaY_Central;
  }

  // the default offset is zero - but this allows you to add resolution and possibly change the angle

  double DefaultCVUniverse::GetThetaXmu() const {
    const std::string thetaBranch = "muon_thetaX";
    return GetDouble(thetaBranch.c_str())+GetBeamAngleOffsetX();
  }

  double DefaultCVUniverse::GetThetaYmu() const {
    const std::string thetaBranch = "muon_thetaY";
    return GetDouble(thetaBranch.c_str())+GetBeamAngleOffsetY();
  }

    
  //! Theta is a derived quantity
  double DefaultCVUniverse::GetThetamu() const {
    return theta3D(GetThetaXmu(),GetThetaYmu());
  }

  double DefaultCVUniverse::GetPhimu() const {
    return phi3D(GetThetaXmu(),GetThetaYmu());
  }


//! 4 vector
  ROOT::Math::PxPyPzEVector DefaultCVUniverse::GetMuon4V() const {
    double p      = GetPmu();
    double thetax = GetThetaXmu();
    double thetay = GetThetaYmu();
    double theta  = GetThetamu();
    double px     = p*std::sin(thetax);
    double py     = GetPmu()*std::sin(thetay);
    double pz     = GetPmu()*std::cos(theta);
    double E      = GetEmu();
    return     ROOT::Math::PxPyPzEVector(px,py,pz,E);
  }
    

//! Generic Branch Getters
  int DefaultCVUniverse::GetInt(const char* name) const {
    return m_chw->GetInt(name, m_entry); 
  }


  bool DefaultCVUniverse::GetBool(const char* name) const {
    return m_chw->GetValue(name, m_entry);
  }


  double DefaultCVUniverse::GetDouble(const char* name) const {
    return m_chw->GetValue(name, m_entry); 
  }


  int DefaultCVUniverse::GetVecElemInt(const char* name, const int i) const { 
    return m_chw->GetValue(name, m_entry, i);
  }


  double DefaultCVUniverse::GetVecElem(const char* name, const int i) const {
    return m_chw->GetValue(name, m_entry, i); 
  }

  double DefaultCVUniverse::GetVecElem(const char* name, 
                                       const int i,
                                       const int j) const {
    std::vector<std::vector<double>> vec = GetVecOfVecDouble(name);
    return vec[i][j];
  }

  std::vector<int> DefaultCVUniverse::GetVecInt(const char* name) const {
    return m_chw->GetValueVector<int>(name, m_entry);
  }

  std::vector<double> DefaultCVUniverse::GetVecDouble(const char* name) const {
    return m_chw->GetValueVector<double>(name, m_entry);
  }

  std::vector<std::vector<double> > DefaultCVUniverse::GetVecOfVecDouble(const char* name) const {
    return m_chw->GetValueNDVector<std::vector<std::vector<double>>> (name,m_entry);
  }

// helper functions
  double DefaultCVUniverse::phi3D(double thetaX,double thetaY) const {
    return std::atan2(std::tan(thetaY),std::tan(thetaX));
  }

  double DefaultCVUniverse::theta3D(double thetaX, double thetaY) const {
    double sec_thetaX = 1./std::cos(thetaX);
    double sec_thetaY = 1./std::cos(thetaY);
    double inter = std::sqrt(1./(sec_thetaX*sec_thetaX + sec_thetaY*sec_thetaY -1.));
    double angle = std::acos(inter);

    if(thetaX > 3.14159265/2 || thetaY > 3.14159265/2 ) return 3.14159265 - angle;
    else return angle;
  }

  double DefaultCVUniverse::calcq0(const double Enu, const double Elep) const {
    return Enu-Elep;
  }

  double DefaultCVUniverse::calcq3(const double Q2, const double Enu, const double Elep) const { 
    return sqrt(Q2 + pow(Enu - Elep,2.0)); 
  }

#endif // DEFAULTCVUNIVERSE_cxx
