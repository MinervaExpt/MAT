#ifndef MNVTUNESYSTEMATICS_H
#define MNVTUNESYSTEMATICS_H

#include "TreeWrapper.h"

// Helper functions declared/defined in the .cxx
// GetRPASystematicsMap(typename T::config_t chain);
// GetLowQ2PiSystematicsMap(typename T::config_t chain);
// Get2p2hSystematicsMap(typename T::config_t chain);

namespace PlotUtils{
  //=================================================================================
  // 2p2h
  //=================================================================================
    template<class T>
    class Universe2p2h: public T
    {
      public:
        Universe2p2h(typename T::config_t chw, double nsigma, int variation = -1);

        virtual double GetLowRecoil2p2hWeight() const /*override*/;
        double GetWeightRatioToCV() const;

        virtual std::string ShortName() const /*override*/;
        virtual std::string LatexName() const /*override*/;
        virtual bool IsVerticalOnly() const   { return true; }/*override*/;

        void SetVariation(int i);
        int m_variation;
    };


  //=================================================================================
  // RPA
  //=================================================================================
    template<class T>
    class RPAUniverse: public T
    {
      public:
        RPAUniverse(typename T::config_t chw, double nsigma, int variation = -1, 
                    std::string q2_region = "LowQ2");

        virtual double GetRPAWeight() const /*override*/;
        double GetWeightRatioToCV() const;

        virtual std::string ShortName() const /*override*/;
        virtual std::string LatexName() const /*override*/;
        virtual bool IsVerticalOnly() const   { return true; }/*override*/;

        void SetVariation(int i);
        int m_variation;
        std::string m_q2_region;
    };


  //=================================================================================
  // LowQ2Pion
  //=================================================================================
    template <typename T>
    class LowQ2PionUniverse : public T
    {
      public:
        LowQ2PionUniverse(typename T::config_t chw, double nsigma);

        virtual double GetLowQ2PiWeight(std::string channel) const /*override*/;
        //double GetWeightRatioToCV() const; //TODO: Revisit this when I'm ready for Reweighters and Universes that are pre-configured with channel as a member variable.

        virtual std::string ShortName() const /*override*/;
        virtual std::string LatexName() const /*override*/;
        virtual bool IsVerticalOnly() const   { return true; }/*override*/;
    };

  //=================================================================================
  // MnvHadronReweighter (TODO PLAYLIST-DEPENDENT, FIX WITH FLUX WEIGHTS)
  //=================================================================================
  /*
    const double kTrackerFace = 5991.29; // module 27, plane 1
    const double kTrackerBack = 8408.91; // module 80, plane 2
    const double kStandardApothem = 850;
    PlotUtils::MnvHadronReweight& weight_hadron(TTree* truth_tree, TTree* reco_tree,
                                                double front   = PlotUtils::kTrackerFace,
                                                double back    = PlotUtils::kTrackerBack,
                                                double apothem = PlotUtils::kStandardApothem,
                                                std::string cache_dir = "./", 
                                                std::string project_name = "" ) {
      static PlotUtils::MnvHadronReweight* _weight_hadron = 0;
      if(!_weight_hadron){
        _weight_hadron = new PlotUtils::MnvHadronReweight(truth_tree, reco_tree);
        _weight_hadron->useDefaultParticles();
        _weight_hadron->setBasicFiducial(front, back, apothem);
        // Get normalization factors from file, else make the file.
        if (!_weight_hadron->tryLoadingFromFile(cache_dir.c_str() , project_name.c_str())) {
          // And now make the renorm factors and save them
          _weight_hadron->getRenormFactors(cache_dir.c_str(), project_name.c_str());
        } 
      }
      else{
        // Here I would make sure initialized correctly, but HR functions cross
        // that bridge when it comes to it.
      }
      return *_weight_hadron;
    }

    template <typename T>
    class HadronUniverse : public T
    {
      public:
        HadronUniverse(typename T::config_t chw, double nsigma, std::string type)
          : T(chw, nsigma), m_type(type)
        {
          if(abs(int(nsigma))!=1) 
            throw std::invalid_argument("HadronUniverse(): nsigma must be +/-1");
        }

        virtual double GetHadronWeight() const override { 
          return PlotUtils::weight_hadron()
        }

        std::string m_type;
        virtual std::string ShortName() const { return "Hadron_" + m_type; }
        virtual std::string LatexName() const { return "Hadron_" + m_type; }
    };

    template <typename T>
    std::vector<T*> GetHadronSystematics(typename T::config_t chain){ 
      std::vector<T*> ret;
      ret.push_back(new PlotUtils::HadronUniverse<T>(chain, -1), "Pion");
      ret.push_back(new PlotUtils::HadronUniverse<T>(chain, +1), "Pion");
      ret.push_back(new PlotUtils::HadronUniverse<T>(chain, -1), "Proton");
      ret.push_back(new PlotUtils::HadronUniverse<T>(chain, +1), "Proton");
      ret.push_back(new PlotUtils::HadronUniverse<T>(chain, -1), "Neutron");
      ret.push_back(new PlotUtils::HadronUniverse<T>(chain, +1), "Neutron");
      return ret;
    }
  */
}

// Explicit Instantiation, needed for loading into python
// See Ana/PlotUtils/dict/PlotUtilsDict.h

// For template classes, the header needs to know about the definitions
#include "MnvTuneSystematics.cxx"


#endif // MNVTUNESYSTEMATICS_H
