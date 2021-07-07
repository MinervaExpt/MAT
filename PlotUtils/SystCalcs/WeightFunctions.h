#ifndef WEIGHTFUNCTIONS_H
#define WEIGHTFUNCTIONS_H

#include "TruthFunctions.h"
#include "PlotUtils/NSFDefaults.h" 

// Get Weights
virtual double GetGenieWeight() const {
  double nonResPiWgt = UseNonResPiReweight() && PlotUtils::IsNonResPi(*this)
                     ? PlotUtils::kNonResPiWeight : 1.;
  double deutWgt = UseDeuteriumGeniePiTune() && PlotUtils::IsCCRes(*this) ? 
                   ( PlotUtils::GetGenieParReweight(*this,"truth_genie_wgt_MaRES", 
                                                  NSFDefaults::DEUTERIUM_MaRES,
                                                  NSFDefaults::GENIE_MaRES, 
                                                  NSFDefaults::GENIE_MaRES_1Sig ) * 
                    NSFDefaults::DEUTERIUM_RES_NORM ) : 1.;
  return nonResPiWgt * deutWgt;
}

virtual double GetRPAWeight() const {
  const int variation = 0;  // CV
  return PlotUtils::GetRPAWeight(*this, Getq0True() / 1000 /* GeV */,
                                 Getq3True() / 1000 /* GeV */, variation,
                                 IsProcessingNX());
}

virtual double GetLowRecoil2p2hWeight() const {
  const int variation = 0;  // CV
  return PlotUtils::GetLowRecoil2p2hWeight(*this, Getq0True() / 1000 /* GeV */,
                                           Getq3True() / 1000 /* GeV */,
                                           variation);
}

virtual double GetLowQ2PiWeight(std::string channel) const {
  int variation = 0;  // CV
  if (!PlotUtils::IsCCRes(*this))
    return 1.;
  else
    return PlotUtils::weight_lowq2pi().getWeight(GetQ2True() * 1e-6 /*GeV^2*/,
                                                 channel, variation);
}

virtual double GetCoherentPiWeight(double thpi_true /*deg*/,
                                   double tpi_true /*GeV*/) const {
  if (GetInt("mc_intType") != 4) return 1.;
  assert(tpi_true > 0. && "GetCoherentPiWeight failed with tpi < 0.");
  assert(thpi_true > 0. && "GetCoherentPiWeight failed with thpi < 0.");
  return PlotUtils::weight_coherent().get_combined_weight(thpi_true, tpi_true);
}

virtual double GetFluxAndCVWeight(double Enu = -99. /*GeV*/,
                                  int nu_pdg = -99.) const {
  if (Enu == -99.) Enu = GetDouble("mc_incomingE") * 1e-3;
  // For LE, electron-neutrino fluxes aren't available, so we should force
  // nu_pdg to have absolute value 14. The +/- doesn't matter, because FRW sets
  // up both, anyways. This may change in the future if electron-neutrino LE
  // fluxes become available.
  if (nu_pdg == -99 && !IsPlaylistME(GetPlaylist()))
    nu_pdg = 14;
  else if (nu_pdg == -99)
    nu_pdg = GetInt("mc_incoming");
  return PlotUtils::flux_reweighter(GetPlaylist(), nu_pdg, UseNuEConstraint(),
                                    GetNFluxUniverses())
      .GetFluxCVWeight(Enu, nu_pdg);
}

virtual double GetTargetMassWeight() const {
  return 1.0;  
}

virtual double GetFSIWeight(int iWeight) const {
  static PlotUtils::weight_fsi weight_FSI;
  weight_FSI.UseTrackingThreshold();
  weight_FSI.calcWeights(
      GetInt("mc_incoming"), GetInt("mc_primaryLepton"), GetInt("mc_charm"),
      GetInt("mc_intType"), GetInt("mc_targetA"), GetInt("mc_targetZ"),
      GetInt("mc_resID"), GetInt("mc_er_nPart"), GetVecInt("mc_er_ID"),
      GetVecInt("mc_er_status"), GetVecInt("mc_er_FD"), GetVecInt("mc_er_LD"),
      GetVecInt("mc_er_mother"), GetVecDouble("mc_er_Px"),
      GetVecDouble("mc_er_Py"), GetVecDouble("mc_er_Pz"),
      GetVecDouble("mc_er_E"), m_entry);
  if (iWeight == 0)
    return weight_FSI.GetElasticWeight(1) * weight_FSI.GetAbsorptionWeight();
  if (iWeight == 1) return weight_FSI.GetElasticWeight(1);
  if (iWeight == 2) return weight_FSI.GetAbsorptionWeight();
  return 1;
}

virtual double GetMKWeight() const { 
  return PlotUtils::weight_mk().getWeight(m_chw, m_entry);
}

virtual double GetGeantHadronWeight() const {
  if (m_is_truth) return 1.;  // No efficiency reweighting for truth events
  //Don't want to spend time creating thing weighter if you're not using this
  if( m_use_mhrw_neutronCV_reweight )
  {
    SetupMHRWeighter();
    return PlotUtils::weight_hadron<PlotUtils::TreeWrapper*>().reweightNeutronCV( *this );
  }
  else return 1.;
}

#endif  // WEIGHTFUNCTIONS
