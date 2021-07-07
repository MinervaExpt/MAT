#ifndef MNVTUNESYSTEMATICS_CXX
#define MNVTUNESYSTEMATICS_CXX

#include "weight_2p2h.h"
#include "weightRPA.h"
#include "weightLowQ2Pi.h"
#include "MnvTuneSystematics.h"

using namespace PlotUtils;

// Helper functions -- get Weighters, containers of systematics universes
namespace PlotUtils{
  //=================================================================================
  // 2p2h
  //=================================================================================
  
    // Calculate the weight given an event and variation
    template<typename T>
    double GetLowRecoil2p2hWeight(const T& universe, double q0, double q3, int variation) {
      //variations
      //0 : CV value
      //1 : nn+pp pair only
      //2 : np pair only
      //3 : qe 1p1h variation.
      //4 : wgt = 1
      if(variation>=4 /*|| !do2p2h*/ ){
        return 1.0;
      }
      if(universe.GetInt("mc_intType")!=1 && universe.GetInt("mc_intType")!=8){
        return 1.0; //if given invalid variations or the user doesn't want to apply
                    // this modification or the MC event is not CCQE or MEC
      }

      if(universe.GetInt("mc_targetZ")<2) {
	return 1.0; //There is no 2p2h on hydrogen. The QE variation is not appropriate for QE on hydrogen.
      }
      bool applyOn2p2h = false;
      bool applyOn1p1h = false;

      if(variation==0 || variation==1 || variation==2) applyOn2p2h = true;//CV, nn+pp pair only, np pair only fits
      else if (variation ==3) applyOn1p1h = true;
      if(universe.GetInt("mc_intType")==1 && !applyOn1p1h) return 1.0; //if CCQE and don't apply 1p1h, don't apply weights
      if(universe.GetInt("mc_intType")==8 && !applyOn2p2h) return 1.0; //if MEC and don't apply 2p2h, don't apply weights
      if(universe.GetInt("mc_intType")==8 && applyOn1p1h)  return 1.0; //if MEC and don't apply 1p1h, don't apply weights
      bool isnnorpp = false;
      bool isnp = false;
      //now target analysis
      int target = universe.GetInt("mc_targetNucleon");
      if(target-2000000200==0 || target-2000000200==2) isnnorpp = true;
      if(target-2000000200==1) isnp = true;

      if(variation==1 && !isnnorpp) return 1.0;//variation 1 is for nn/pp only interactions
      if(variation==2 && !isnp) return 1.0;//variation 2 is for np only interactions
      double ret = 1.0;
      if(variation==0)      ret=PlotUtils::weight_2p2h_cv().getWeight(q0,q3); // pass as GeV
      else if(variation==1) ret=PlotUtils::weight_2p2h_nn().getWeight(q0,q3); // pass as GeV
      else if(variation==2) ret=PlotUtils::weight_2p2h_np().getWeight(q0,q3); // pass as GeV
      else if(variation==3) ret=PlotUtils::weight_2p2h_qe().getWeight(q0,q3); // pass as GeV
      else std::cout <<"Should not have gotten here. GetLowRecoil2p2hWeight."<< std::endl;
      return ret;
    }


    // Get 2p2h Systematics Containers
    template <typename T>
    std::vector<T*> Get2p2hSystematics(typename T::config_t chain ) { 
      std::vector<T*> ret;
      const double nsigma = 1.;
      const unsigned int n_variations = 3;
      for (unsigned int i = 1; i <= n_variations; ++i)
        ret.push_back(new PlotUtils::Universe2p2h<T>(chain, nsigma, i));
      return ret;
    }


    template <typename T>
    std::map< std::string, std::vector<T*> > Get2p2hSystematicsMap(typename T::config_t chain ) {
      std::map< std::string, std::vector<T*> > ret;
      const double nsigma = 1.;
      const unsigned int n_variations = 3;
      for (unsigned int i = 1; i <= n_variations; ++i)
        ret["2p2h"].push_back(new PlotUtils::Universe2p2h<T>(chain, nsigma, i));
      return ret;
    }


  //=================================================================================
  // RPA
  //=================================================================================
    template<typename T>
    double GetRPAWeight(const T& universe, double q0, double q3, int variation, bool useNX) {
      double ret = 1.0;
      //variations
      //0 : CV value
      //1 : hq2pos
      //2 : hq2neg
      //3 : lq2pos
      //4 : lq2neg
      if(variation > 4/*|| !doRPA*/){
        return 1.0;
      }

      if(universe.GetInt("mc_intType")!=1) return 1.0;
      if(universe.GetInt("mc_targetZ")<6) return 1.0;
      if(variation==0) ret=PlotUtils::weightRPA_cv_and_var(useNX).getWeight(q0,q3);
      else if(variation==1) ret=PlotUtils::weightRPA_cv_and_var(useNX).getWeightHighQ2(q0,q3,1);
      else if(variation==2) ret=PlotUtils::weightRPA_cv_and_var(useNX).getWeightHighQ2(q0,q3,-1);
      else if(variation==3) ret=PlotUtils::weightRPA_cv_and_var(useNX).getWeightLowQ2(q0,q3,1);
      else if(variation==4) ret=PlotUtils::weightRPA_cv_and_var(useNX).getWeightLowQ2(q0,q3,-1);
      else 
        throw std::runtime_error("RPAUniverse::GetRPAWeight: invalid variation");

      return ret;
    }



    template <typename T>
    std::map< std::string, std::vector<T*> > GetRPASystematicsMap(typename T::config_t chain ) {
      std::map< std::string, std::vector<T*> > ret;
      const double nsigma = 1.;
      ret["HighQ2"].push_back(new PlotUtils::RPAUniverse<T>(chain, nsigma, 1, "HighQ2"));
      ret["HighQ2"].push_back(new PlotUtils::RPAUniverse<T>(chain, nsigma, 2, "HighQ2"));
      ret["LowQ2"].push_back(new PlotUtils::RPAUniverse<T>(chain, nsigma, 3, "LowQ2"));
      ret["LowQ2"].push_back(new PlotUtils::RPAUniverse<T>(chain, nsigma, 4, "LowQ2"));
      return ret;
    }


  //=================================================================================
  // LowQ2Pi
  //=================================================================================

    template <typename T>
    bool IsCCRes(const T& universe) {
      bool is_ccres = universe.GetInt("mc_intType") == 2  // Res
                      &&
                      universe.GetInt("mc_current") == 1; // CC
      return is_ccres;
    }


    template <typename T>
    std::vector<T*> GetLowQ2PiSystematics(typename T::config_t chain ) {
      std::vector<T*> ret;
      ret.push_back(new PlotUtils::LowQ2PionUniverse<T>(chain, -1));
      ret.push_back(new PlotUtils::LowQ2PionUniverse<T>(chain, +1));
      return ret;
    }

    template <typename T>
    std::map< std::string, std::vector<T*> > GetLowQ2PiSystematicsMap(typename T::config_t chain ) {
      std::map< std::string, std::vector<T*> > ret;
      ret["LowQ2Pi"].push_back(new PlotUtils::LowQ2PionUniverse<T>(chain, -1));
      ret["LowQ2Pi"].push_back(new PlotUtils::LowQ2PionUniverse<T>(chain, +1));
      return ret;
    }
}


// Class Definitions
//=================================================================================
// 2p2h
//=================================================================================
  // Constructor
  template<typename T>
  Universe2p2h<T>::Universe2p2h(typename T::config_t chw, double nsigma, int variation)
    : T(chw, nsigma),
      m_variation(variation)
  {}


  template<typename T>
  double Universe2p2h<T>::GetLowRecoil2p2hWeight() const {
    return PlotUtils::GetLowRecoil2p2hWeight(*this, T::Getq0True()/1000, T::Getq3True()/1000, m_variation);
  }

  template <typename T>
  double Universe2p2h<T>::GetWeightRatioToCV() const {
    return PlotUtils::GetLowRecoil2p2hWeight(*this, T::Getq0True()/1000, T::Getq3True()/1000, m_variation) / PlotUtils::GetLowRecoil2p2hWeight(*this, T::Getq0True()/1000, T::Getq3True()/1000, 0); //Variation 0 is the CV
  }


  template<typename T>
  void Universe2p2h<T>::SetVariation(int i){ m_variation = i; }


  template<typename T>
  std::string Universe2p2h<T>::ShortName() const { return "Low_Recoil_2p2h_Tune"; }


  template<typename T>
  std::string Universe2p2h<T>::LatexName() const { return "Low Recoil 2p2h Tune"; }


//=================================================================================
// RPA
//=================================================================================
  // Constructor
  template<typename T>
  RPAUniverse<T>::RPAUniverse(typename T::config_t chw, double nsigma, 
                              int variation, std::string q2_region)
    : T(chw, nsigma),
      m_variation(variation), m_q2_region(q2_region)
  {}


  template<typename T>
  double RPAUniverse<T>::GetRPAWeight( ) const {
    return PlotUtils::GetRPAWeight(*this, T::Getq0True()/1000, T::Getq3True()/1000, m_variation, T::IsProcessingNX());
  }

  template <typename T>
  double RPAUniverse<T>::GetWeightRatioToCV() const {
    return PlotUtils::GetRPAWeight(*this, T::Getq0True()/1000, T::Getq3True()/1000, m_variation, T::IsProcessingNX()) / PlotUtils::GetRPAWeight(*this, T::Getq0True()/1000, T::Getq3True()/1000, 0, T::IsProcessingNX()); //Variation 0 is the CV
  }

  template<typename T>
  void RPAUniverse<T>::SetVariation(int i){ m_variation = i; }


  template<typename T>
  std::string RPAUniverse<T>::ShortName() const { return "RPA_" + m_q2_region; }


  template<typename T>
  std::string RPAUniverse<T>::LatexName() const { return "RPA " + m_q2_region; }


//=================================================================================
// LowQ2Pi
//=================================================================================
  // Constructor
  template<typename T>
  LowQ2PionUniverse<T>::LowQ2PionUniverse(typename T::config_t chw, double nsigma)
    : T(chw, nsigma)
  {
    if(abs(int(nsigma))!=1) 
      throw std::invalid_argument("LowQ2PionUniverse(): nsigma must be +/-1");
  }

  template<typename T>
  double LowQ2PionUniverse<T>::GetLowQ2PiWeight(std::string channel) const { 
    if(!PlotUtils::IsCCRes(*this)) 
      return 1.;
    else
      return PlotUtils::weight_lowq2pi().getWeight(T::GetQ2True()*1e-6 /*GeV^2*/, channel, T::m_nsigma);
  }

  //TODO: Come back to this when I'm ready for Reweighters that provide systematics with a pre-configured channel member.
  /*template <typename T>
  double LowQ2PionUniverse<T>::GetWeightRatioToCV() {
  }*/

  template<typename T>
  std::string LowQ2PionUniverse<T>::ShortName() const { return "LowQ2Pi"; }


  template<typename T>
  std::string LowQ2PionUniverse<T>::LatexName() const { return "LowQ2Pi"; }





#endif // MNVTUNESYSTEMATICS_CXX
