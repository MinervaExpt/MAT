#ifndef GeanHadronSystematics_cxx
#define GeanHadronSystematics_cxx

#include "GeantHadronSystematics.h"
#include "TreeWrapper.h"
#include <cassert>
#include <iostream>

#include "NSFDefaults.h"

using namespace PlotUtils;

//=================================================================================
// Helper Functions
//=================================================================================
//More information about how to use this (or build a new renormalization file) is in 
//docdb 28556
namespace PlotUtils {
// Static instance of MnvHadronReweight
template <class T>
MnvHadronReweight& weight_hadron(
    T chw = 0, 
    std::string playlist = "ME1A",
    double min_z = NSFDefaults::TrackerFace,
    double max_z = NSFDefaults::TrackerBack,
    double apothem = NSFDefaults::StandardApothem,
    bool useNeutronCVReweight = true,
    bool useElastics = true,
    std::string project_name = "Renorm_Kine_Truth",
    std::string directory = Form("%s/data/mhrwKineRenorm",getenv("PLOTUTILSROOT")) ){
  static MnvHadronReweight* _weight_hadron = 0;
  // First time -- use the input args to make a MHRW
  if (!_weight_hadron) {
    // check chw, project name, FV are all legit
    // TODO
    /*
    if (!chw) {
      std::cerr << "weight_hadron: Static MnvHadronReweighter hasn't been "
                << "initialized yet, but you didn't pass a valid chain to "
                << "this getter. Returning 0.";
      return *_weight_hadron;
    }
    */
    TTree* tree = chw->GetTree();
    assert(tree);
    _weight_hadron = new MnvHadronReweight(NULL, tree);
    _weight_hadron->useDefaultParticles();
    _weight_hadron->setPlaylist(playlist);
    _weight_hadron->setBasicFiducial(min_z, max_z, apothem);
    if( !useElastics ) _weight_hadron->TurnOffElasticReweight();
    if( useNeutronCVReweight ) _weight_hadron->useReweightedNeutronCV();
    _weight_hadron->setProjectName(project_name);
    _weight_hadron->setInputDirectory( directory );
    std::cout << "MnvHadronReweight object initialized \n  Project name "
              << project_name << "\n Playlist " 
              << playlist << "\n Using Elastics "<<useElastics<<"\n Using NeutronCV "<< useNeutronCVReweight
              <<"\n Fiducial vol " << min_z << " < z < "
              << max_z << ", apothem, " << apothem
              << "\n  Default particles (p, n, pi) being considered.\n";
  } 
  //else {  // Not first time -- make sure user isn't trying to change chw,
  //          // project name, or FV. Also, ensure _weight_hadron is well
  //          // constructed.
  //  // check chw is good
  //  // TODO

  //  // check projectname hasn't changed
  //  if (project_name != _weight_hadron->getProjectName()) {
  //    std::cout << "GeantHadronSystematics WARNING: you're getting a "
  //                 "MnvHadronReweight object with a different project name ("
  //              << project_name << ") than the one you initially used.\n";
  //    std::cout << "This isn't an internally consistent thing to do. I won't "
  //                 "allow it. Continuing to use the original project name "
  //              << _weight_hadron->getProjectName() << "\n";
  //  }
  //}
  return *_weight_hadron;
}

// Get standard universes -- set MnvHadronReweight options
template <class T>
std::map<std::string, std::vector<T*> > GetGeantHadronSystematicsMap(
    typename T::config_t chain ) { 

  // Advanced users - write your own GetMap function based off of this one and
  // in it, modify the mhrw here:
    // For neutron vars. The CV is bad.
    // weight_hadron<typename T::config_t>().m_reweightNeutronCV = true;

    // For neutron vars. Elastic collisions just as important as inelastic.
    // weight_hadron<typename T::config_t>().doElasticReweight = true; 

    // Kaons functionality not ready. Lucky you, kaon analyser, you get to
    // finish it.
    // weight_hadron<typename T::config_t>().SetParticle(321, true); 

    // Remove particles.
    // weight_hadron<typename T::config_t>().removeParticle(2212);

  std::map<std::string, std::vector<T*> > ret;
  ret[std::string("GEANT_Proton")].push_back( new GeantHadronUniverse<T>(chain, 1., 2212 )); 
  ret[std::string("GEANT_Proton")].push_back( new GeantHadronUniverse<T>(chain, -1., 2212 ));
  ret[std::string("GEANT_Neutron")].push_back(new GeantHadronUniverse<T>(chain, 1., 2112 ));
  ret[std::string("GEANT_Neutron")].push_back(new GeantHadronUniverse<T>(chain, -1., 2112 ));
  ret[std::string("GEANT_Pion")].push_back(   new GeantHadronUniverse<T>(chain, 1., 211 ));
  ret[std::string("GEANT_Pion")].push_back(   new GeantHadronUniverse<T>(chain, -1., 211 ));
  return ret;
}
}  // namespace PlotUtils

//=================================================================================
// Class
//=================================================================================
// Constructor
template <class T>
GeantHadronUniverse<T>::GeantHadronUniverse(typename T::config_t chw,
                                            double nsigma, int pdg)
    : T(chw, nsigma), m_pdg(pdg) {
  //Want to get the renorm factors once.  Also doing some setup once
  T::SetupMHRWeighter();
  
  // Make a nicer label for the universe
  switch( abs(pdg) ){
    case 2212:
      m_part_name = "Proton";
      break;
    case 2112:
      m_part_name = "Neutron";
      break;
    case 211:
      m_part_name = "Pion";
      break;
    default:
      //I don't think there's a reweight for anything else
      m_part_name = "Unknown";
  }
}

// Get weight
template <class T>
double GeantHadronUniverse<T>::GetGeantHadronWeight() const {
  if( T::IsTruth() ) return 1; //No reweighting for truth events 
  InelXSecWeights weights =
      weight_hadron<PlotUtils::TreeWrapper*>().getWeights(*this);
  double wgt = 1.;
  if (weights.eventHas[m_pdg]) {
    wgt = T::m_nsigma < 0 ? weights.weightDown[m_pdg] : weights.weightUp[m_pdg];
  }
  return wgt * T::GetGeantHadronWeight();
}

template <class T>
double GeantHadronUniverse<T>::GetWeightRatioToCV() const {
  if( T::IsTruth() ) return 1; //No reweighting for truth events 
  InelXSecWeights weights =
      weight_hadron<PlotUtils::TreeWrapper*>().getWeights(*this);
  double wgt = 1.;
  if (weights.eventHas[m_pdg]) {
    wgt = T::m_nsigma < 0 ? weights.weightDown[m_pdg] : weights.weightUp[m_pdg];
  }

  return wgt;
}

#endif  // GeanHadronSystematics_cxx
