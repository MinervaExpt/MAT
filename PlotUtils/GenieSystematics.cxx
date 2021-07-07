#ifndef GENIESYSTEMATICS_CXX
#define GENIESYSTEMATICS_CXX

#include "NSFDefaults.h"
#include "weightZExp.h"
#include "GenieSystematics.h"
#include "TDecompChol.h"
#include <iostream>

using namespace PlotUtils;

// Helper Functions
namespace PlotUtils{
  //=================================================================================
  // Standard Genie
  //=================================================================================
  
  template <class T>
  std::vector<T*> GetStandardGenieSystematics(typename T::config_t chain, int sigma, bool include_res = true ) {
    std::vector<T*> ret;

    std::vector<std::string> genie_systematics;
      genie_systematics.push_back("AGKYxF1pi"); 
      genie_systematics.push_back("AhtBY");
      genie_systematics.push_back("BhtBY");
      genie_systematics.push_back("CCQEPauliSupViaKF");
      genie_systematics.push_back("CV1uBY");
      genie_systematics.push_back("CV2uBY");
      genie_systematics.push_back("EtaNCEL");
      genie_systematics.push_back("FrAbs_N");
      genie_systematics.push_back("FrAbs_pi");
      genie_systematics.push_back("FrCEx_N");
      genie_systematics.push_back("FrCEx_pi");
      genie_systematics.push_back("FrElas_N");
      genie_systematics.push_back("FrElas_pi");
      genie_systematics.push_back("FrInel_N");
      //genie_systematics.push_back("FrInel_pi"); // Dan says that this knob should not currently be evaluated
                                                  // but that we should revisit this eventually
      genie_systematics.push_back("FrPiProd_N");
      genie_systematics.push_back("FrPiProd_pi");
      genie_systematics.push_back("MFP_N");
      genie_systematics.push_back("MFP_pi");
      genie_systematics.push_back("MaNCEL");
      if( include_res )
      {
        genie_systematics.push_back("MaRES");
        genie_systematics.push_back("MvRES");
        //genie_systematics.push_back("NormCCRES"); // NormCCRES is the normalization component of MaRES, so
                                                    // only one should be used
      }
      genie_systematics.push_back("NormDISCC");
      genie_systematics.push_back("NormNCRES");
      genie_systematics.push_back("RDecBR1gamma");
      genie_systematics.push_back("Rvn2pi");
      genie_systematics.push_back("Rvp2pi");
      genie_systematics.push_back("Theta_Delta2Npi");
      genie_systematics.push_back("VecFFCCQEshape"); 

    for(std::vector<std::string>::const_iterator syst = genie_systematics.begin(); 
            syst != genie_systematics.end(); ++syst){
      ret.push_back( new PlotUtils::GenieUniverse<T>(chain, sigma, *syst) );
    }

    return ret;
  }


  template <class T>
  std::map< std::string, std::vector<T*> > GetStandardGenieSystematicsMap(typename T::config_t chain, bool include_res = true ) {
    std::map< std::string, std::vector<T*> > ret;

    std::vector<double> sigmas;
    sigmas.push_back(-1.);
    sigmas.push_back(+1.);

    std::vector<std::string> genie_systematics;
      genie_systematics.push_back("AGKYxF1pi"); 
      genie_systematics.push_back("AhtBY");
      genie_systematics.push_back("BhtBY");
      genie_systematics.push_back("CCQEPauliSupViaKF");
      genie_systematics.push_back("CV1uBY");
      genie_systematics.push_back("CV2uBY");
      genie_systematics.push_back("EtaNCEL");
      genie_systematics.push_back("FrAbs_N");
      genie_systematics.push_back("FrAbs_pi");
      genie_systematics.push_back("FrCEx_N");
      genie_systematics.push_back("FrCEx_pi");
      genie_systematics.push_back("FrElas_N");
      genie_systematics.push_back("FrElas_pi");
      genie_systematics.push_back("FrInel_N");
      //genie_systematics.push_back("FrInel_pi"); // Dan says that this knob should not currently be evaluated
                                                  // but that we should revisit this eventually
      genie_systematics.push_back("FrPiProd_N");
      genie_systematics.push_back("FrPiProd_pi");
      genie_systematics.push_back("MFP_N");
      genie_systematics.push_back("MFP_pi");
      genie_systematics.push_back("MaNCEL");
      if( include_res )
      {
        genie_systematics.push_back("MaRES");
        genie_systematics.push_back("MvRES");
        //genie_systematics.push_back("NormCCRES"); // NormCCRES is the normalization component of MaRES, so
                                                    // only one should be used
      }
      genie_systematics.push_back("NormDISCC");
      genie_systematics.push_back("NormNCRES");
      genie_systematics.push_back("RDecBR1gamma");
      genie_systematics.push_back("Rvn2pi");
      genie_systematics.push_back("Rvp2pi");
      genie_systematics.push_back("Theta_Delta2Npi");
      genie_systematics.push_back("VecFFCCQEshape"); 

    for(std::vector<double>::const_iterator sigma = sigmas.begin(); 
            sigma != sigmas.end(); ++sigma) {
      for(std::vector<std::string>::const_iterator syst = genie_systematics.begin(); 
              syst != genie_systematics.end(); ++syst){
        ret[*syst].push_back(new PlotUtils::GenieUniverse<T>(chain, *sigma, *syst));
      }
    }
       
    return ret;
  }

  //=================================================================================
  // Genie MaRES and NormCCRes with covariance matrix
  //=================================================================================
  
  template <class T>
  std::vector<T*> GetGenieResPionFitCovSystematics(typename T::config_t chain, int sigma, int row) {
    std::vector<T*> ret;

    ret.push_back( new PlotUtils::GenieMaNormResCovUniverse<T>(chain, sigma, row) );

    return ret;
  }

  template <class T>
  std::map< std::string, std::vector<T*> > GetGenieResPionFitCovSystematicsMap(typename T::config_t chain) {
    std::map< std::string, std::vector<T*> > ret;

    std::vector<double> sigmas;
    sigmas.push_back(-1.);
    sigmas.push_back(+1.);

    for( int iRow = 0; iRow < 2; ++iRow ) {
      for(std::vector<double>::const_iterator sigma = sigmas.begin(); 
              sigma != sigmas.end(); ++sigma) {
        ret[Form("GENIE_Ma_Norm_RES_Cov_Row%d",iRow)].push_back( new PlotUtils::GenieMaNormResCovUniverse<T>(chain, *sigma, iRow) );
      }
    }

    return ret;
  }

  //=================================================================================
  // Genie MaRES and NormCCRES with new systematics
  //=================================================================================
  
  template <class T>
  std::vector<T*> GetGenieResPionFitSystematics(typename T::config_t chain, int sigma) {
    std::vector<T*> ret;

    ret.push_back( new PlotUtils::GenieMaResUniverse<T>(chain, sigma) );
    ret.push_back( new PlotUtils::GenieNormCCResUniverse<T>(chain, sigma) );

    return ret;
  }


  template <class T>
  std::map< std::string, std::vector<T*> > GetGenieResPionFitSystematicsMap(typename T::config_t chain) {
    std::map< std::string, std::vector<T*> > ret;

    std::vector<double> sigmas;
    sigmas.push_back(-1.);
    sigmas.push_back(+1.);

    for(std::vector<double>::const_iterator sigma = sigmas.begin(); 
            sigma != sigmas.end(); ++sigma) {
      //TODO Change
      ret["GENIE_D2_MaRES"].push_back( new PlotUtils::GenieMaResUniverse<T>(chain, *sigma) );
      ret["GENIE_D2_NormCCRES"].push_back( new PlotUtils::GenieNormCCResUniverse<T>(chain, *sigma) );
    }
       
    return ret;
  }

  //=================================================================================
  // Genie MvRES with new systematics
  //=================================================================================
  
  template <class T>
  std::vector<T*> GetGenieEPMvResSystematics(typename T::config_t chain, int sigma) {
    std::vector<T*> ret;

    ret.push_back( new PlotUtils::GenieMvResUniverse<T>(chain, sigma) );

    return ret;
  }


  template <class T>
  std::map< std::string, std::vector<T*> > GetGenieEPMvResSystematicsMap(typename T::config_t chain) {
    std::map< std::string, std::vector<T*> > ret;

    std::vector<double> sigmas;
    sigmas.push_back(-1.);
    sigmas.push_back(+1.);

    for(std::vector<double>::const_iterator sigma = sigmas.begin(); 
            sigma != sigmas.end(); ++sigma) {
      //TODO Change
      ret["GENIE_EP_MvRES"].push_back( new PlotUtils::GenieMvResUniverse<T>(chain, *sigma) );
    }
       
    return ret;
  }

  //=================================================================================
  // Genie Rvx1pi
  //=================================================================================
  
  template <class T>
  std::vector<T*> GetGenieRvx1piSystematics(typename T::config_t chain, int sigma) {
    std::vector<T*> ret;

    std::vector<std::string> genie_systematics;
      genie_systematics.push_back("Rvn1pi");
      genie_systematics.push_back("Rvp1pi");

    for(std::vector<std::string>::const_iterator syst = genie_systematics.begin(); 
            syst != genie_systematics.end(); ++syst){
      ret.push_back( new PlotUtils::GenieRvx1piUniverse<T>(chain, sigma, *syst) );
    }

    return ret;
  }


  template <class T>
  std::map< std::string, std::vector<T*> > GetGenieRvx1piSystematicsMap(typename T::config_t chain) {
    std::map< std::string, std::vector<T*> > ret;

    std::vector<double> sigmas;
    sigmas.push_back(-1.);
    sigmas.push_back(+1.);

    std::vector<std::string> genie_systematics;
      genie_systematics.push_back("Rvn1pi");
      genie_systematics.push_back("Rvp1pi");

    for(std::vector<double>::const_iterator sigma = sigmas.begin(); 
            sigma != sigmas.end(); ++sigma) {
      for(std::vector<std::string>::const_iterator syst = genie_systematics.begin(); 
              syst != genie_systematics.end(); ++syst){
        ret[*syst].push_back(new PlotUtils::GenieRvx1piUniverse<T>(chain, *sigma, *syst));
      }
    }
       
    return ret;
  }

  //=================================================================================
  // Genie FaCCQE
  //=================================================================================
  
  template <class T>
  std::vector<T*> GetGenieFaCCQESystematics(typename T::config_t chain,
                                     int n_universes){
    std::vector<T*> ret;
    // if n_universes == -1, provide the default +/- 1 sigma shifts
    if(n_universes == -1){
      std::vector<std::string> genie_systematics;
        genie_systematics.push_back("MaCCQE");
      
      for(std::vector<std::string>::const_iterator syst = genie_systematics.begin();
              syst != genie_systematics.end(); ++syst){
        // I'm not sure how to pass through a purposeful sigma to give to the below line
        // without confusing z-expansion users
        ret.push_back( new PlotUtils::GenieFaCCQEUniverse<T>(chain, 1.0, -1) );
      }
    }
    // for n_universes > 0, provide the z-expansion weights
    else{
      double nsigma = 1.;
      for(unsigned int i = 0; i < n_universes; ++i){
        ret.push_back( new PlotUtils::GenieFaCCQEUniverse<T>(chain, 1.0, i) );
      }
    }
    return ret;
  }


  template <class T>
  std::map< std::string, std::vector<T*> > GetGenieFaCCQESystematicsMap(typename T::config_t chain,
                                                                 int n_universes) {
    std::map< std::string, std::vector<T*> > ret;
    // if n_universe == -1, provide the default +/- 1 sigma shifts
    if(n_universes == -1){
      std::vector<double> sigmas;
      sigmas.push_back(-1.);
      sigmas.push_back(+1.);

      std::vector<std::string> genie_systematics;
        genie_systematics.push_back("MaCCQE");

      for(std::vector<double>::const_iterator sigma = sigmas.begin(); 
              sigma != sigmas.end(); ++sigma) {
        for(std::vector<std::string>::const_iterator syst = genie_systematics.begin(); 
                syst != genie_systematics.end(); ++syst){
          ret[*syst].push_back(new PlotUtils::GenieFaCCQEUniverse<T>(chain, *sigma, -1));
        }
      }
    }
    // for n_universes > 0, provide the z-expansion weights
    else{
      for (int i = 0; i < n_universes; ++i)
        ret["GENIE_MaCCQE"].push_back(new PlotUtils::GenieFaCCQEUniverse<T>(chain, 1.0, i));
    }
    return ret;
  }


  //=================================================================================
  // All Genie
  //=================================================================================
  
  template <class T>
  std::vector<T*> GetGenieSystematics(typename T::config_t chain, int sigma, bool include_res = true ) {
    std::vector<T*> ret;

    int n_zexp_universes = T::UseZExpansionFaReweight() ? 100 : -1;

    std::vector<T*> genie_systematics_standard = GetStandardGenieSystematics<T>(chain,sigma,include_res);
    std::vector<T*> genie_systematics_Rvx1pi = GetGenieRvx1piSystematics<T>(chain,sigma);
    std::vector<T*> genie_systematics_FaCCQE = GetGenieFaCCQESystematics<T>(chain,n_zexp_universes);

    ret.insert(genie_systematics_standard.begin(),genie_systematics_standard.end());
    ret.insert(genie_systematics_Rvx1pi.begin(),genie_systematics_Rvx1pi.end());
    ret.insert(genie_systematics_FaCCQE.begin(),genie_systematics_FaCCQE.end());
    
    return ret;
  }


  template <class T>
  std::map< std::string, std::vector<T*> > GetGenieSystematicsMap(typename T::config_t chain, bool include_res = true) {
    std::map< std::string, std::vector<T*> > ret;

    int n_zexp_universes = T::UseZExpansionFaReweight() ? 100 : -1;

    std::map< std::string, std::vector<T*> > genie_systematics_map_standard = GetStandardGenieSystematicsMap<T>(chain,include_res);
    std::map< std::string, std::vector<T*> > genie_systematics_map_Rvx1pi = GetGenieRvx1piSystematicsMap<T>(chain);
    std::map< std::string, std::vector<T*> > genie_systematics_map_FaCCQE = GetGenieFaCCQESystematicsMap<T>(chain,n_zexp_universes);

    ret.insert(genie_systematics_map_standard.begin(),genie_systematics_map_standard.end());
    ret.insert(genie_systematics_map_Rvx1pi.begin(),genie_systematics_map_Rvx1pi.end());
    ret.insert(genie_systematics_map_FaCCQE.begin(),genie_systematics_map_FaCCQE.end());

    return ret;
  }


  //=================================================================================
  // Non-resonant pion
  //=================================================================================
  template <typename T>
  bool IsNonResPi(const T& universe) {
    bool is_nonrespi = universe.GetVecElem("truth_genie_wgt_Rvn1pi", 2) < 1.0
                       ||
                       universe.GetVecElem("truth_genie_wgt_Rvp1pi", 2) < 1.0;
    return is_nonrespi;
  }

  //=================================================================================
  // Resonant pion (this is already in MnvTuneSystematics)
  //=================================================================================
  template <typename T>
  bool IsCCResonance(const T& universe) {
    bool is_ccres = universe.GetInt("mc_intType") == 2  // Res
                    &&
                    universe.GetInt("mc_current") == 1; // CC
    return is_ccres;
  }

  //=================================================================================
  // Method to adjust the GENIE weight given a new parameter 
  // Based of CCCohPionUtils GetMa(v)ResWeight
  //=================================================================================
  template <class T>
  double GetGenieParReweight( const T& universe, std::string branch_name, 
                                double newPar, double oldPar, double oldSigma )
  {
    int wgtIdx = (newPar < oldPar) ? 2 : 4;
    // This should return 
    //((difference between parameters)/(1 old sigma))*(shift in weight based off 1 old sigma)
    return 1.0 + ( fabs(newPar - oldPar) / oldSigma )*(universe.GetVecElem(branch_name.c_str(),wgtIdx) - 1.0);
  }
}


//=================================================================================
// Genie
//=================================================================================
// Constructor
template <class T>
GenieUniverse<T>::GenieUniverse(typename T::config_t chw, 
                                          double nsigma, std::string genie_name)
  : T(chw, nsigma), 
    m_name(genie_name), m_branch_name("truth_genie_wgt_" + genie_name)
{
  //std::cout << "Make a PlotUtils Genie Universe" << std::endl;
  
}


// Get Genie Weight
template <class T>
double GenieUniverse<T>::GetGenieWeight() const { 
  double wgt = T::m_nsigma < 0 ? T::GetVecElem(m_branch_name.c_str(), 2):
                                 T::GetVecElem(m_branch_name.c_str(), 4);
  // Protections against unphysical weights due to a potential lack of a
  // constraint in the GENIE FSI reweighter. i.e. bug in GENIE version 
  // used to generate NX sample
  if(m_name.compare("MFP_N")==0&&wgt>10.){
    wgt = 5.0;
  }
  // Any negative weight is unphysical
  if(m_name.compare("MFP_N")==0&&wgt<0.){
    wgt = 0.0;
  }
  // This systematic universe's modified the CV weight,
  // which includes the effect of the nonResPi Reweight
  return wgt*T::GetGenieWeight();
}

template <class T>
double GenieUniverse<T>::GetWeightRatioToCV() const {
  double wgt = T::m_nsigma < 0 ? T::GetVecElem(m_branch_name.c_str(), 2):
                                 T::GetVecElem(m_branch_name.c_str(), 4);
  // Protections against unphysical weights due to a potential lack of a
  // constraint in the GENIE FSI reweighter. i.e. bug in GENIE version 
  // used to generate NX sample
  if(m_name.compare("MFP_N")==0&&wgt>10.){
    wgt = 5.0;
  }
  // Any negative weight is unphysical
  if(m_name.compare("MFP_N")==0&&wgt<0.){
    wgt = 0.0;
  }
  // This systematic universe's modified the CV weight,
  // which includes the effect of the nonResPi Reweight
  return wgt;
}

template <class T>
std::string GenieUniverse<T>::ShortName() const {
  return "GENIE_" + m_name;
};

template <class T>
std::string GenieUniverse<T>::LatexName() const {
  return "GENIE " + m_name;
};

//=================================================================================
// Genie MvRes systematic based on electroproduction data (New: 3%, Old: 10%)
//=================================================================================
// Constructor
template <class T>
GenieMvResUniverse<T>::GenieMvResUniverse(typename T::config_t chw, double nsigma)
  : T(chw, nsigma), 
    m_name("MvRES"), m_branch_name("truth_genie_wgt_MvRES")
{
}

// Get Genie Weight
template <class T>
double GenieMvResUniverse<T>::GetGenieWeight() const { 
  //Should this be CCRES only?
  if( !PlotUtils::IsCCResonance(*this) ) return T::GetGenieWeight();
  double wgt = PlotUtils::GetGenieParReweight( *this, m_branch_name, 
                                               NSFDefaults::GENIE_MvRES + T::m_nsigma*NSFDefaults::ELECTROPROD_MvRES_1Sig,
                                               NSFDefaults::GENIE_MvRES, NSFDefaults::GENIE_MvRES_1Sig );
  return wgt*T::GetGenieWeight();
}

template <class T>
double GenieMvResUniverse<T>::GetWeightRatioToCV() const {
  //Should this be CCRES only?
  if( !PlotUtils::IsCCResonance(*this) ) return T::GetGenieWeight();
  return PlotUtils::GetGenieParReweight( *this, m_branch_name,
                                         NSFDefaults::GENIE_MvRES + T::m_nsigma*NSFDefaults::ELECTROPROD_MvRES_1Sig,
                                         NSFDefaults::GENIE_MvRES, NSFDefaults::GENIE_MvRES_1Sig );
}

template <class T>
std::string GenieMvResUniverse<T>::ShortName() const {
  //TODO: Change
  //return "GENIE_" + m_name;
  return "GENIE_EP_MvRES";
};

template <class T>
std::string GenieMvResUniverse<T>::LatexName() const {
  return "GENIE " + m_name;
};

//=================================================================================
// Genie NormCCRes systematic based on deuterium fit (old 1.00 +/- 20%, new: 1.15 +/- 7%)
//=================================================================================
// Constructor
template <class T>
GenieNormCCResUniverse<T>::GenieNormCCResUniverse(typename T::config_t chw, double nsigma)
  : T(chw, nsigma), 
    m_name("NormCCRES"), m_branch_name("truth_genie_wgt_NormCCRES")
{
}

// Get Genie Weight
template <class T>
double GenieNormCCResUniverse<T>::GetGenieWeight() const { 
  //Should this be CCRES only?
  if( !PlotUtils::IsCCResonance(*this) ) return T::GetGenieWeight();
  //Check if cv reweight is on. If not, I wouldn't recommend changing the uncertainty  
  double wgt = 1.0;
  if( T::UseDeuteriumGeniePiTune() ) {
    // Neeed to divide the CV weight out
    wgt = (NSFDefaults::DEUTERIUM_RES_NORM+T::m_nsigma*NSFDefaults::DEUTERIUM_RES_NORM_1Sig)/NSFDefaults::DEUTERIUM_RES_NORM;
  }
  else { 
    wgt = T::m_nsigma < 0 ? T::GetVecElem(m_branch_name.c_str(), 2):
                            T::GetVecElem(m_branch_name.c_str(), 4);
  }
  
  return wgt*T::GetGenieWeight();
}

template <class T>
double GenieNormCCResUniverse<T>::GetWeightRatioToCV() const {
  //Should this be CCRES only?
  if( !PlotUtils::IsCCResonance(*this) ) return T::GetGenieWeight();
  //Check if cv reweight is on. If not, I wouldn't recommend changing the uncertainty  
  double wgt = 1.0;
  if( T::UseDeuteriumGeniePiTune() ) {
    // Neeed to divide the CV weight out
    wgt = (NSFDefaults::DEUTERIUM_RES_NORM+T::m_nsigma*NSFDefaults::DEUTERIUM_RES_NORM_1Sig)/NSFDefaults::DEUTERIUM_RES_NORM;
  }
  else {
    wgt = T::m_nsigma < 0 ? T::GetVecElem(m_branch_name.c_str(), 2):
                            T::GetVecElem(m_branch_name.c_str(), 4);
  }

  return wgt;
}

template <class T>
std::string GenieNormCCResUniverse<T>::ShortName() const {
//TODO: Change
  //return "GENIE_" + m_name;
  return "GENIE_D2_NormCCRES";
};

template <class T>
std::string GenieNormCCResUniverse<T>::LatexName() const {
  return "GENIE " + m_name;
};

//=================================================================================
// Genie MaRes systematic based on deuterium fit (old 1.12 +/- 20%, new: 0.94 +/- 5%)
//=================================================================================
// Constructor
template <class T>
GenieMaResUniverse<T>::GenieMaResUniverse(typename T::config_t chw, double nsigma)
  : T(chw, nsigma), 
    m_name("MaRES"), m_branch_name("truth_genie_wgt_MaRES")
{
}

// Get Genie Weight
template <class T>
double GenieMaResUniverse<T>::GetGenieWeight() const { 
  //Should this be CCRES only?
  if( !PlotUtils::IsCCResonance(*this) ) return T::GetGenieWeight();
  //Check if cv reweight is on. If not, I wouldn't recommend changing the uncertainty  
  double wgt = 1.0;
  if( T::UseDeuteriumGeniePiTune() ) {
    // Neeed to divide the CV weight out
    wgt = PlotUtils::GetGenieParReweight( *this, m_branch_name, 
                                           NSFDefaults::DEUTERIUM_MaRES + T::m_nsigma*NSFDefaults::DEUTERIUM_MaRES_1Sig,
                                           NSFDefaults::GENIE_MaRES, NSFDefaults::GENIE_MaRES_1Sig ) / 
                                           PlotUtils::GetGenieParReweight(*this, m_branch_name,  
                                                                          NSFDefaults::DEUTERIUM_MaRES,
                                                                          NSFDefaults::GENIE_MaRES, 
                                                                          NSFDefaults::GENIE_MaRES_1Sig );
  }
  else { 
    wgt = T::m_nsigma < 0 ? T::GetVecElem(m_branch_name.c_str(), 2):
                            T::GetVecElem(m_branch_name.c_str(), 4);
  }
  
  
  return wgt*T::GetGenieWeight();
}

template <class T>
double GenieMaResUniverse<T>::GetWeightRatioToCV() const {
  //Should this be CCRES only?
  if( !PlotUtils::IsCCResonance(*this) ) return T::GetGenieWeight();
  //Check if cv reweight is on. If not, I wouldn't recommend changing the uncertainty  
  double wgt = 1.0;
  if( T::UseDeuteriumGeniePiTune() ) {
    // Neeed to divide the CV weight out
    wgt = PlotUtils::GetGenieParReweight( *this, m_branch_name,
                                           NSFDefaults::DEUTERIUM_MaRES + T::m_nsigma*NSFDefaults::DEUTERIUM_MaRES_1Sig,
                                           NSFDefaults::GENIE_MaRES, NSFDefaults::GENIE_MaRES_1Sig ) /
                                           PlotUtils::GetGenieParReweight(*this, m_branch_name,
                                                                          NSFDefaults::DEUTERIUM_MaRES,
                                                                          NSFDefaults::GENIE_MaRES,
                                                                          NSFDefaults::GENIE_MaRES_1Sig );
  }
  else {
    wgt = T::m_nsigma < 0 ? T::GetVecElem(m_branch_name.c_str(), 2):
                            T::GetVecElem(m_branch_name.c_str(), 4);
  }

  return wgt;
}

template <class T>
std::string GenieMaResUniverse<T>::ShortName() const {
//TODO: Change
  //return "GENIE_" + m_name;
  return "GENIE_D2_MaRES";
};

template <class T>
std::string GenieMaResUniverse<T>::LatexName() const {
  return "GENIE " + m_name;
};

//=================================================================================
// Genie Rvx1pi
//=================================================================================
// Constructor
template <class T>
GenieRvx1piUniverse<T>::GenieRvx1piUniverse(typename T::config_t chw,
                                          double nsigma, std::string genie_name)
  : T(chw, nsigma), 
    m_name(genie_name), m_branch_name("truth_genie_wgt_" + genie_name)
{
  
}


// Get Genie Weight
template <class T>
double GenieRvx1piUniverse<T>::GetGenieWeight() const { 
  double wgt = T::m_nsigma < 0 ? T::GetVecElem(m_branch_name.c_str(), 2):
                                 T::GetVecElem(m_branch_name.c_str(), 4);
  // If the nonrespi reweight is being applied, modify the below branches to
  // use alternate weights, provided in GenieSystematics header
  if(T::UseNonResPiReweight()&&PlotUtils::IsNonResPi(*this)){
    double plusMinus = T::m_nsigma < 0 ? -1. : 1.;
    wgt = PlotUtils::kNonResPiWeight + plusMinus * PlotUtils::kNonResPiWeightShift;
    return wgt;
  }
  // This systematic universe's modified the CV weight,
  // which includes the effect of the nonResPi Reweight
  return wgt*T::GetGenieWeight();
}

template <class T>
double GenieRvx1piUniverse<T>::GetWeightRatioToCV() const {
  double wgt = T::m_nsigma < 0 ? T::GetVecElem(m_branch_name.c_str(), 2):
                                 T::GetVecElem(m_branch_name.c_str(), 4);
  // If the nonrespi reweight is being applied, modify the below branches to
  // use alternate weights, provided in GenieSystematics header
  if(T::UseNonResPiReweight()&&PlotUtils::IsNonResPi(*this)){
    double plusMinus = T::m_nsigma < 0 ? -1. : 1.;
    wgt = 1 + plusMinus * PlotUtils::kNonResPiWeightShift / PlotUtils::kNonResPiWeight;
    return wgt;
  }
  // This systematic universe's modified the CV weight,
  // which includes the effect of the nonResPi Reweight
  return wgt;
}

template <class T>
std::string GenieRvx1piUniverse<T>::ShortName() const {
  return "GENIE_" + m_name;
};


template <class T>
std::string GenieRvx1piUniverse<T>::LatexName() const {
  return "GENIE " + m_name;
};

//=================================================================================
// Genie Fitted MaRes+NormRes systematic w/ Cov Matrix
//=================================================================================
// Constructor
template <class T>
GenieMaNormResCovUniverse<T>::GenieMaNormResCovUniverse(typename T::config_t chw, double nsigma, int row)
  : T(chw, nsigma), m_row(row),
    m_ma_branch_name("truth_genie_wgt_MaRES"),m_norm_branch_name("truth_genie_wgt_NormCCRES") 
{
  //Don't know where else to put this for now
  //TMatrixD corr(2,2);
  ////MaRES            NORM RES
  //corr[1][0] = -0.9; corr[1][1] = 1.0;  
  //corr[0][0] = 1.0;  corr[0][1] = -0.9;  

  //TDecompChol decomp_chol(corr);
  //decomp_chol.Decompose();
  //TMatrixD d = decomp_chol.GetU();
  //m_factors.push_back(d[row][0]); m_factors.push_back(d[row][1]);

  //TMatrixD dan(2,2);
  ////MaRES            NORM RES
  //dan[1][0] = -0.9; dan[1][1] = 0.4359;  
  //dan[0][0] = 1.0;  dan[0][1] = 0.0;  
  //m_factors.push_back(dan[row][0]); m_factors.push_back(dan[row][1]);
}

// Get Genie Weight
template <class T>
double GenieMaNormResCovUniverse<T>::GetGenieWeight() const { 
  if( !PlotUtils::IsCCResonance(*this) ) return T::GetGenieWeight();
  //Check if cv reweight is on. If not, I wouldn't recommend changing the uncertainty  
  double wgt = 1.0;
  if( T::UseDeuteriumGeniePiTune() ) {
    double new_ma_wgt = PlotUtils::GetGenieParReweight(*this,"truth_genie_wgt_MaRES", NSFDefaults::DEUTERIUM_MaRES, 
                                                         NSFDefaults::GENIE_MaRES, NSFDefaults::GENIE_MaRES_1Sig );
    //double ma_wgt_sigma = NSFDefaults::DEUTERIUM_RES_NORM * (PlotUtils::GetGenieParReweight( *this, m_ma_branch_name.c_str(), 
    //                                                          NSFDefaults::DEUTERIUM_MaRES + T::m_nsigma*NSFDefaults::DEUTERIUM_MaRES_1Sig,
    //                                                          NSFDefaults::GENIE_MaRES, NSFDefaults::GENIE_MaRES_1Sig )-new_ma_wgt); 
    double ma_wgt_1sigma = NSFDefaults::DEUTERIUM_RES_NORM * (PlotUtils::GetGenieParReweight( *this, m_ma_branch_name.c_str(), 
                                                              T::m_nsigma < 0 ?
                                                              NSFDefaults::DEUTERIUM_MaRES - NSFDefaults::DEUTERIUM_MaRES_1Sig :
                                                              NSFDefaults::DEUTERIUM_MaRES + NSFDefaults::DEUTERIUM_MaRES_1Sig,
                                                              NSFDefaults::GENIE_MaRES, NSFDefaults::GENIE_MaRES_1Sig )-new_ma_wgt); 
    double norm_wgt_1sigma = NSFDefaults::DEUTERIUM_RES_NORM_1Sig*new_ma_wgt;
      
    double wgt_sigma = 0;
    if(m_row == 0)
    {
      wgt_sigma = T::m_nsigma*ma_wgt_1sigma;     
    }
    else if(m_row == 1)
    {
      wgt_sigma = T::m_nsigma*norm_wgt_1sigma*(-0.9+0.4359);
    }

    wgt = NSFDefaults::DEUTERIUM_RES_NORM * new_ma_wgt + wgt_sigma;
    //Need remove the CV reweight for this
    wgt = wgt/( NSFDefaults::DEUTERIUM_RES_NORM * new_ma_wgt );
  }
  else {//Just combine the weights
    wgt = T::m_nsigma < 0 ? T::GetVecElem(m_ma_branch_name.c_str(), 2)*T::GetVecElem(m_norm_branch_name.c_str(), 2) :
                            T::GetVecElem(m_ma_branch_name.c_str(), 4)*T::GetVecElem(m_norm_branch_name.c_str(), 4);
  }
  
  
  return wgt*T::GetGenieWeight();
}

template <class T>
double GenieMaNormResCovUniverse<T>::GetWeightRatioToCV() const {
  if( !PlotUtils::IsCCResonance(*this) ) return 1.;
  //Check if cv reweight is on. If not, I wouldn't recommend changing the uncertainty  
  double wgt = 1.0;
  if( T::UseDeuteriumGeniePiTune() ) {
    double new_ma_wgt = PlotUtils::GetGenieParReweight(*this,"truth_genie_wgt_MaRES", NSFDefaults::DEUTERIUM_MaRES, 
                                                         NSFDefaults::GENIE_MaRES, NSFDefaults::GENIE_MaRES_1Sig );
    //double ma_wgt_sigma = NSFDefaults::DEUTERIUM_RES_NORM * (PlotUtils::GetGenieParReweight( *this, m_ma_branch_name.c_str(), 
    //                                                          NSFDefaults::DEUTERIUM_MaRES + T::m_nsigma*NSFDefaults::DEUTERIUM_MaRES_1Sig,
    //                                                          NSFDefaults::GENIE_MaRES, NSFDefaults::GENIE_MaRES_1Sig )-new_ma_wgt); 
    double ma_wgt_1sigma = NSFDefaults::DEUTERIUM_RES_NORM * (PlotUtils::GetGenieParReweight( *this, m_ma_branch_name.c_str(), 
                                                              T::m_nsigma < 0 ?
                                                              NSFDefaults::DEUTERIUM_MaRES - NSFDefaults::DEUTERIUM_MaRES_1Sig :
                                                              NSFDefaults::DEUTERIUM_MaRES + NSFDefaults::DEUTERIUM_MaRES_1Sig,
                                                              NSFDefaults::GENIE_MaRES, NSFDefaults::GENIE_MaRES_1Sig )-new_ma_wgt); 
    double norm_wgt_1sigma = NSFDefaults::DEUTERIUM_RES_NORM_1Sig*new_ma_wgt;
      
    double wgt_sigma = 0;
    if(m_row == 0)
    {
      wgt_sigma = T::m_nsigma*ma_wgt_1sigma;     
    }
    else if(m_row == 1)
    {
      wgt_sigma = T::m_nsigma*norm_wgt_1sigma*(-0.9+0.4359);
    }

    wgt = NSFDefaults::DEUTERIUM_RES_NORM * new_ma_wgt + wgt_sigma;
    //Need remove the CV reweight for this
    wgt = wgt/( NSFDefaults::DEUTERIUM_RES_NORM * new_ma_wgt );
  }
  else {//Just combine the weights
    wgt = T::m_nsigma < 0 ? T::GetVecElem(m_ma_branch_name.c_str(), 2)*T::GetVecElem(m_norm_branch_name.c_str(), 2) :
                            T::GetVecElem(m_ma_branch_name.c_str(), 4)*T::GetVecElem(m_norm_branch_name.c_str(), 4);
  }
  
  return wgt;
}

template <class T>
std::string GenieMaNormResCovUniverse<T>::ShortName() const {
//TODO: Change
  //return "GENIE_" + m_name;
  return Form("GENIE_Ma_Norm_RES_Cov_Row%d",m_row);
};

template <class T>
std::string GenieMaNormResCovUniverse<T>::LatexName() const {
  return Form("GENIE Ma Norm RES Cov Row %d",m_row);
};

//=================================================================================
// Genie FaCCQE
//=================================================================================
// Constructor
template <class T>
GenieFaCCQEUniverse<T>::GenieFaCCQEUniverse(typename T::config_t chw, 
                                          double nsigma, int universe_number)
  : T(chw, nsigma), 
    m_universe_number(universe_number)
{
  
}


// Get Genie Weight
template <class T>
double GenieFaCCQEUniverse<T>::GetGenieWeight() const { 
  double wgt = T::m_nsigma < 0 ? T::GetVecElem("truth_genie_wgt_MaCCQE", 2):
                                 T::GetVecElem("truth_genie_wgt_MaCCQE", 4);
  // If the nonrespi reweight is being applied, modify the below branches to
  // use alternate weights, provided in GenieSystematics header
  if(m_universe_number == -1){
    // This systematic universe's modified the CV weight,
    // which includes the effect of the nonResPi Reweight
    return wgt*T::GetGenieWeight();
  }
  if(T::GetInt("mc_intType")!=1) return T::GetGenieWeight();  // Only reweight CCQE events
  if(T::GetInt("mc_targetZ")<6) return T::GetGenieWeight(); // Don't reweight hydrogen
  const double q2 = T::GetDouble("mc_Q2")/(1000*1000); // Convert to GeV
  static PlotUtils::weightZExp zExpWeighter = PlotUtils::weightZExp("$MPARAMFILESROOT/data/Reweight/Z_Expansion_Reweight_v2126.root");
  double cvWgt = zExpWeighter.getWeight(q2);
  double univWgt = zExpWeighter.getWeight(q2,m_universe_number);
  return (univWgt/cvWgt)*T::GetGenieWeight();
}

template <class T>
double GenieFaCCQEUniverse<T>::GetWeightRatioToCV() const {
  double wgt = T::m_nsigma < 0 ? T::GetVecElem("truth_genie_wgt_MaCCQE", 2):
                                 T::GetVecElem("truth_genie_wgt_MaCCQE", 4);
  // If the nonrespi reweight is being applied, modify the below branches to
  // use alternate weights, provided in GenieSystematics header
  if(m_universe_number == -1){
    // This systematic universe's modified the CV weight,
    // which includes the effect of the nonResPi Reweight
    return wgt;
  }
  if(T::GetInt("mc_intType")!=1) return 1;  // Only reweight CCQE events
  if(T::GetInt("mc_targetZ")<6) return 1; // Don't reweight hydrogen
  const double q2 = T::GetDouble("mc_Q2")/(1000*1000); // Convert to GeV
  static PlotUtils::weightZExp zExpWeighter = PlotUtils::weightZExp("$MPARAMFILESROOT/data/Reweight/Z_Expansion_Reweight_v2126.root");
  double cvWgt = zExpWeighter.getWeight(q2);
  double univWgt = zExpWeighter.getWeight(q2,m_universe_number);
  return univWgt/cvWgt;
}

template <class T>
std::string GenieFaCCQEUniverse<T>::ShortName() const {
  return "GENIE_MaCCQE";
};


template <class T>
std::string GenieFaCCQEUniverse<T>::LatexName() const {
  return "GENIE MaCCQE";
};

#endif // GENIESYSTEMATICS_CXX
