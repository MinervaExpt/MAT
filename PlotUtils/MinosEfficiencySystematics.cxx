#ifndef MINOSEFFICIENCYSYSTEMATICS_CXX
#define MINOSEFFICIENCYSYSTEMATICS_CXX

#include "MinosMuonEfficiencyCorrection.h"
#include "MnvNormalization.h"
#include "MinosEfficiencySystematics.h"
#include <iostream>

// Helper functions -- get Weighters, containers of systematics universes
namespace PlotUtils {

  template <class T>
  std::vector<T*> GetMinosEfficiencySystematics(typename T::config_t chain ) {
    std::vector<T*> ret;
    
    ret.push_back(new PlotUtils::MinosEfficiencyUniverse<T>(chain, -1.));
    ret.push_back(new PlotUtils::MinosEfficiencyUniverse<T>(chain, 1.));

    return ret;
  }


  template <class T>
  std::map< std::string, std::vector<T*> > GetMinosEfficiencySystematicsMap(
      typename T::config_t chain ) {
    std::map< std::string, std::vector<T*> > ret;
    
    ret["MinosEfficiency"].push_back(new PlotUtils::MinosEfficiencyUniverse<T>(chain, -1.));
    ret["MinosEfficiency"].push_back(new PlotUtils::MinosEfficiencyUniverse<T>(chain, 1.));

    return ret;
  }

}


// Class Definitions
// Constructor
template<typename T>
PlotUtils::MinosEfficiencyUniverse<T>::MinosEfficiencyUniverse(
    typename T::config_t chw, double nsigma) : T(chw, nsigma) {}


template<typename T>
double PlotUtils::MinosEfficiencyUniverse<T>::GetMinosEfficiencyWeight() const {

  if( T::IsTruth() ) return 1;  
  if( !T::isFHC() && !T::isRHC() ) return 1; //LE or non standard palylist
  double cv_wgt = T::GetMinosEfficiencyWeight();

  // hopefully temporary
  double pmu = T::GetPmuMinos()/1000;          // GetCorrectionErr takes GeV
  double thetamu = T::GetThetamu()*180/3.14159;// GetCorrectionErr takes angle, not radians
  double wgt_shift;

  if(T::IsPlaylistME( T::GetPlaylist() )){ //The correction factors are different between ME and LE
    wgt_shift = PlotUtils::MinosMuonEfficiencyCorrection::Get( T::isFHC() ).GetCorrectionErr(
                           pmu, thetamu, 0, T::isFHC() );
  }
  else{ 
       pmu=pmu*1000;// Normalizer LE wants MeV	
       static PlotUtils::MnvNormalizer mnvNormalizer = PlotUtils::MnvNormalizer("Eroica",T::GetPlaylist());
       wgt_shift = mnvNormalizer.GetCorrectionErr(pmu);
  }

  return  cv_wgt + T::m_nsigma*wgt_shift;
  //The GetCorrectionErr assigns a constant error to batch_pot.  If that
  //changes, you'll have to change that 0 to BatchPOT(), but until then, i'm
  //keeping it a constant. 
}

template <typename T>
double  PlotUtils::MinosEfficiencyUniverse<T>::GetWeightRatioToCV() const {
  //In these 2 situations, the CV weight is also 1.
  if( T::IsTruth() ) return 1;
  if( !T::isFHC() && !T::isRHC() ) return 1; //LE or non standard palylist

  //Calculate the CV weight
  //TODO: Maybe don't duplicate this code?  Seems like it's pretty minimal code duplication to me.
  double cv_wgt;
  if (T::IsPlaylistME(T::GetPlaylist())) {  // The correction factors are different
                                      // between ME and LE
    double pmu = T::GetPmuMinos() / 1000;  // GetCorrection expects GeV
    cv_wgt = PlotUtils::MinosMuonEfficiencyCorrection::Get(T::isFHC()).GetCorrection(
        pmu, T::GetBatchPOT(), T::isFHC());
  } else {                       // Assume if not ME, then it's LE
    double pmu = T::GetPmuMinos();  // MnVnormalizer GetCorrection expects MeV
#ifndef __CINT__
    static PlotUtils::MnvNormalizer mnvNormalizer =
        PlotUtils::MnvNormalizer("Eroica", T::GetPlaylist());
#endif // __CINT__
    cv_wgt = mnvNormalizer.GetCorrection(pmu);
  }

  // hopefully temporary
  double pmu = T::GetPmuMinos()/1000;          // GetCorrectionErr takes GeV
  double thetamu = T::GetThetamu()*180/3.14159;// GetCorrectionErr takes angle, not radians
  double wgt_shift;

  if(T::IsPlaylistME( T::GetPlaylist() )){ //The correction factors are different between ME and LE
    wgt_shift = PlotUtils::MinosMuonEfficiencyCorrection::Get( T::isFHC() ).GetCorrectionErr(
                           pmu, thetamu, 0, T::isFHC() );
  }
  else{
       pmu=pmu*1000;// Normalizer LE wants MeV  
       static PlotUtils::MnvNormalizer mnvNormalizer = PlotUtils::MnvNormalizer("Eroica",T::GetPlaylist());
       wgt_shift = mnvNormalizer.GetCorrectionErr(pmu);
  }

  return 1 + T::m_nsigma*wgt_shift/cv_wgt;
  //The GetCorrectionErr assigns a constant error to batch_pot.  If that
  //changes, you'll have to change that 0 to BatchPOT(), but until then, i'm
  //keeping it a constant. 
}

template<typename T>
std::string PlotUtils::MinosEfficiencyUniverse<T>::ShortName() const { return "MINOS_Reconstruction_Efficiency"; }


template<typename T>
std::string PlotUtils::MinosEfficiencyUniverse<T>::LatexName() const { return "MINOS Reconstruction Efficiency"; }


#endif // MINOSEFFICIENCYSYSTEMATICS_CXX
