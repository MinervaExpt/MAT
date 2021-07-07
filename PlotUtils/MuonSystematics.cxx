#ifndef MUONSYSTEMATICS_CXX
#define MUONSYSTEMATICS_CXX

#include "MuonSystematics.h"
#include <iostream>

// Helper functions -- get Weighters, containers of systematics universes
namespace PlotUtils {
  //=================================================================================
  // Minerva muon-momentum-shifted universe 
  //=================================================================================
  template <class T>
  std::vector<T*> GetMinervaMuonSystematics(typename T::config_t chain ) {
    std::vector<T*> ret;

    ret.push_back(new PlotUtils::MuonUniverseMinerva<T>(chain, -1.));
    ret.push_back(new PlotUtils::MuonUniverseMinerva<T>(chain, 1.));

    return ret;
  }

  template <class T>
  std::map< std::string, std::vector<T*> > GetMinervaMuonSystematicsMap( typename T::config_t chain ) {
    std::map< std::string, std::vector<T*> > ret;
    
    ret["Muon_Energy_MINERvA"].push_back(new PlotUtils::MuonUniverseMinerva<T>(chain, -1.));
    ret["Muon_Energy_MINERvA"].push_back(new PlotUtils::MuonUniverseMinerva<T>(chain, 1.));

    return ret;
  }


  //=================================================================================
  // Minos muon-momentum-shifted universe 
  //=================================================================================
  template <class T>
  std::vector<T*> GetMinosMuonSystematics(typename T::config_t chain ) {
    std::vector<T*> ret;

    ret.push_back(new PlotUtils::MuonUniverseMinos<T>(chain, -1.));
    ret.push_back(new PlotUtils::MuonUniverseMinos<T>(chain, 1.));

    return ret;
  }

  template <class T>
  std::map< std::string, std::vector<T*> > GetMinosMuonSystematicsMap( typename T::config_t chain ) {
    std::map< std::string, std::vector<T*> > ret;
    
    ret["Muon_Energy_MINOS"].push_back(new PlotUtils::MuonUniverseMinos<T>(chain, -1.));
    ret["Muon_Energy_MINOS"].push_back(new PlotUtils::MuonUniverseMinos<T>(chain, 1.));

    return ret;
  }
}


// Class Definitions
//=================================================================================
// Minerva muon-momentum-shifted universe 
//=================================================================================
  // Constructor
  template<typename T>
  PlotUtils::MuonUniverseMinerva<T>::MuonUniverseMinerva(
      typename T::config_t chw, double nsigma) : T(chw, nsigma) {}
  
  template<typename T>
  double PlotUtils::MuonUniverseMinerva<T>::GetMuonMomentumShiftMinerva() const {
  
    double muon_shift_minerva;
  
    //--------------------------------------------------
    // Uncertainty on Minerva component of muon momentum
    //--------------------------------------------------
  
    double vertex_Z = T::GetVertexZ();
    double dEdxUncertainty = 0.0;
    double materialAssayUncertainty = 0.0;
    // Use a different correction in the nuclear targets
    if(vertex_Z<NSFDefaults::nuclearTargetZEnd){
      dEdxUncertainty = NSFDefaults::dEdxUncertainty;
      materialAssayUncertainty = NSFDefaults::MaterialAssayUncertainty;
    }
    else{
      dEdxUncertainty = NSFDefaults::dEdxUncertaintyNoNuke;
      materialAssayUncertainty = NSFDefaults::MaterialAssayUncertaintyNoNuke;
    } 
    muon_shift_minerva = sqrt(
      dEdxUncertainty*dEdxUncertainty +
      materialAssayUncertainty*materialAssayUncertainty );
   
    double shift_val_minerva = T::m_nsigma * muon_shift_minerva;
    
    return shift_val_minerva;
  }
  
  template<typename T>
  double PlotUtils::MuonUniverseMinerva<T>::GetPmuMinerva() const {
  
    double shift_val = GetMuonMomentumShiftMinerva();
  
    return shift_val+T::GetPmuMinerva();
  }
  
  template<typename T>
  std::string PlotUtils::MuonUniverseMinerva<T>::ShortName() const { return "Muon_Energy_MINERvA"; }
  
  
  template<typename T>
  std::string PlotUtils::MuonUniverseMinerva<T>::LatexName() const { return "Muon Energy MINERvA"; }

//=================================================================================
// Minos muon-momentum-shifted universe
//=================================================================================
  // Constructor
  template<typename T>
  PlotUtils::MuonUniverseMinos<T>::MuonUniverseMinos(
      typename T::config_t chw, double nsigma) : T(chw, nsigma) {}
  
  template<typename T>
  double PlotUtils::MuonUniverseMinos<T>::GetMuonMomentumShiftMinos() const {
  
    double muon_shift_minos;
  
    //--------------------------------------------------
    // Uncertainty on Minos components of muon momentum
    //--------------------------------------------------
  
    // Was Minos component of muon momentum reconstructed by range or curvature
    std::string minosRangeCurveBool = T::GetAnaToolName()+"_minos_used_curvature";
    bool minosRecoByCurve = T::GetBool(minosRangeCurveBool.c_str());
  
    // Sys shifts are made relative to the total (minos component of the) muon momentum
    double minosCVMomentum = T::GetPmuMinos();
    // Uncertainty due to range on Minos component of muon momentum
    double minosRangeUncertainty = minosCVMomentum*NSFDefaults::MinosMuonPRange_Err;
    // Uncertainty due to curvature on Minos component of muon momentum
    double minosCurveUncertainty = 0.0;
    double localCurveUncertainty = 0.0;
    if(minosRecoByCurve){
      localCurveUncertainty = minosCVMomentum>1000.0?NSFDefaults::MinosMuonHighPCurvature_Err:NSFDefaults::MinosMuonLowPCurvature_Err;
      minosCurveUncertainty = minosCVMomentum*localCurveUncertainty;
    }
  
    muon_shift_minos = sqrt(
        minosRangeUncertainty*minosRangeUncertainty +
        minosCurveUncertainty*minosCurveUncertainty );
  
    double shift_val_minos = T::m_nsigma * muon_shift_minos;
    
    return shift_val_minos;
  }
  
  template<typename T>
  double PlotUtils::MuonUniverseMinos<T>::GetPmuMinos() const {
  
    double shift_val = GetMuonMomentumShiftMinos();
  
    return shift_val+T::GetPmuMinos();
  }
  
  template <class T>
  double PlotUtils::MuonUniverseMinos<T>::GetFluxAndCVWeight(double Enu,
                                                             int nu_pdg) const {
    // If this playlist is not ME, `GetFluxAndCVWeight` shouldn't change
    if (!T::IsPlaylistME(T::GetPlaylist())) return T::GetFluxAndCVWeight();
    // If this playlist is ME RHC, `GetFluxAndCVWeight` shouldn't change
    // This will change if/when RHC correlated fluxes become available
    if (!T::isFHC()) return T::GetFluxAndCVWeight();

    if (Enu == -99.)   Enu    = T::GetDouble("mc_incomingE")*1e-3;
    if (nu_pdg == -99) nu_pdg = T::GetInt("mc_incoming");
    int universe = T::m_nsigma < 0 ? 0 : 1;
    // This universes's flux weight is modified because the flux is strongly correlated with the muon energy scale.
    // The fluxes used for these systematic universes were extracted as part of Amit's wiggle studies.
    double wgtMod = PlotUtils::flux_reweighter(
      T::GetPlaylist(),nu_pdg,T::UseNuEConstraint(),T::GetNFluxUniverses()).GetSysUniFluxWeightCorrection(Enu,nu_pdg,"Muon_Energy",universe);
    return wgtMod*T::GetFluxAndCVWeight();
  }

  template <class T>
  double PlotUtils::MuonUniverseMinos<T>::GetWeightRatioToCV() const {
    // If this playlist is not ME, `GetFluxAndCVWeight` shouldn't change
    if (!T::IsPlaylistME(T::GetPlaylist())) return 1;
    // If this playlist is ME RHC, `GetFluxAndCVWeight` shouldn't change
    // This will change if/when RHC correlated fluxes become available
    if (!T::isFHC()) return 1;

    const double Enu    = T::GetDouble("mc_incomingE")*1e-3;
    const int nu_pdg = T::GetInt("mc_incoming");
    int universe = T::m_nsigma < 0 ? 0 : 1;
    // This universes's flux weight is modified because the flux is strongly correlated with the muon energy scale.
    // The fluxes used for these systematic universes were extracted as part of Amit's wiggle studies.
    double wgtMod = PlotUtils::flux_reweighter(
      T::GetPlaylist(),nu_pdg,T::UseNuEConstraint(),T::GetNFluxUniverses()).GetSysUniFluxWeightCorrection(Enu,nu_pdg,"Muon_Energy",universe);
    return wgtMod;
  }

  template<typename T>
  std::string PlotUtils::MuonUniverseMinos<T>::ShortName() const { return "Muon_Energy_MINOS"; }
  
  
  template<typename T>
  std::string PlotUtils::MuonUniverseMinos<T>::LatexName() const { return "Muon Energy MINOS"; }


#endif // MUONSYSTEMATICS_CXX
