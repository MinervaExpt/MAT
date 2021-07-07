#ifndef PLOTUTILSDICT_H 
#define PLOTUTILSDICT_H 1

// Include files for PlotUtils dictionary.

/** @file PlotUtilsDict.h
 *  
 *
 *  @author Jeremy Wolcott <jwolcott@fnal.gov>
 *  @date   2012-11-25
 */
// ============================================================================
// PlotUtils
// ============================================================================

// here we need to include all the header files
// for the classes we want to make dictionaries of

#include <vector>

#include "../PlotUtils/NSFDefaults.h"
#include "../PlotUtils/AnaBinning.h"
#include "../PlotUtils/ArachneUtils.h"
#include "../PlotUtils/Exceptions.h"
#include "../PlotUtils/HistogramUtils.h"
#include "../PlotUtils/MnvFluxConstraint.h" 
#include "../PlotUtils/MnvH1D.h"
#include "../PlotUtils/MnvH2D.h"
#include "../PlotUtils/MnvH3D.h"
#include "../PlotUtils/MnvAnaTuple.h"
#include "../PlotUtils/MnvLatErrorBand.h"
#include "../PlotUtils/MnvLatErrorBand2D.h"
#include "../PlotUtils/MnvLatErrorBand3D.h"
#include "../PlotUtils/MnvPlotter.h"
#include "../PlotUtils/MnvRecoShifter.h"
#include "../PlotUtils/MnvVertErrorBand.h"
#include "../PlotUtils/MnvVertErrorBand2D.h"
#include "../PlotUtils/MnvVertErrorBand3D.h"
#include "../PlotUtils/TargetUtils.h"
#include "../PlotUtils/MnvNormalization.h"
#include "../PlotUtils/POTCounter.h"
#include "../PlotUtils/FluxReweighter.h"
#include "../PlotUtils/FluxReweighterWithWiggleFit.h"
#include "../PlotUtils/ChainWrapper.h"
#include "../PlotUtils/TreeWrapper.h"
#include "../PlotUtils/HyperDimLinearizer.h"
#include "../PlotUtils/PhysicsVariables.h"
#include "../PlotUtils/MnvColors.h"
#include "../PlotUtils/GridCanvas.h"

// PlotUtils weight classes
#include "../PlotUtils/weightRPA.h"
#include "../PlotUtils/weight_2p2h.h"
#include "../PlotUtils/weightLowQ2Pi.h"
#include "../PlotUtils/weightDIS.h"
#include "../PlotUtils/weightZExp.h"

//PlotUtils systematic universes classes (new sys framework)
#include "../PlotUtils/DefaultCVUniverse.h"
#include "../PlotUtils/MinervaUniverse.h"
//#include "../PlotUtils/HistWrapper.h"
//#include "../PlotUtils/Hist2DWrapper.h"
#include "../PlotUtils/FluxSystematics.h"
#include "../PlotUtils/GenieSystematics.h"
#include "../PlotUtils/GeantHadronSystematics.h"
#include "../PlotUtils/MnvTuneSystematics.h"
#include "../PlotUtils/MinosEfficiencySystematics.h"
#include "../PlotUtils/MuonSystematics.h"
#include "../PlotUtils/MuonResolutionSystematics.h"
#include "../PlotUtils/MichelSystematics.h"
#include "../PlotUtils/AngleSystematics.h"
#include "../PlotUtils/ResponseSystematics.h"


//TODO: Do I need this?
#include "../PlotUtils/ErrorHandler.h"

// this garbage is necessary so that gccxml is able to create dictionaries for these custom containers
// (since it otherwise doesn't know which specific version of these templated classes to instantiate)
// see: http://root.cern.ch/root/roottalk/roottalk10/0035.html
// somehow std::map<>s seem to be instantiated somewhere else, so explicit instantiation is not necessary?
#ifdef __GCCXML__
template class std::vector<PlotUtils::MnvEVD::Event>;                                       // the 'Events' typedef
template class std::pair<std::string, std::vector<PlotUtils::MnvEVD::Event> >;              // the 'EventGroup' typedef

template class std::map<std::string, std::vector<std::string> >;
template class std::pair<std::string, std::vector<std::string> >;

template class PlotUtils::BeamAngleXUniverse<PlotUtils::DefaultCVUniverse>;
template class PlotUtils::BeamAngleYUniverse<PlotUtils::DefaultCVUniverse>;
template class PlotUtils::FluxUniverse<PlotUtils::DefaultCVUniverse>;
template class PlotUtils::GenieUniverse<PlotUtils::DefaultCVUniverse>;
template class PlotUtils::GenieRvx1piUniverse<PlotUtils::DefaultCVUniverse>;
template class PlotUtils::GenieFaCCQEUniverse<PlotUtils::DefaultCVUniverse>;
template class PlotUtils::ResponseUniverse<PlotUtils::DefaultCVUniverse>;
template class PlotUtils::MinosEfficiencyUniverse<PlotUtils::DefaultCVUniverse>;
template class PlotUtils::MuonAngleXResolutionUniverse<PlotUtils::DefaultCVUniverse>;
template class PlotUtils::MuonAngleYResolutionUniverse<PlotUtils::DefaultCVUniverse>;
template class PlotUtils::MuonResolutionUniverse<PlotUtils::DefaultCVUniverse>;
template class PlotUtils::MichelEfficiencyUniverse<PlotUtils::DefaultCVUniverse>;
template class PlotUtils::LowQ2PionUniverse<PlotUtils::DefaultCVUniverse>;
template class PlotUtils::MuonUniverseMinerva<PlotUtils::DefaultCVUniverse>;
template class PlotUtils::MuonUniverseMinos<PlotUtils::DefaultCVUniverse>;
template class PlotUtils::RPAUniverse<PlotUtils::DefaultCVUniverse>;
template class PlotUtils::Universe2p2h<PlotUtils::DefaultCVUniverse>;

//template class PlotUtils::HistWrapper<PlotUtils::DefaultCVUniverse>;
//template class PlotUtils::Hist2DWrapper<PlotUtils::DefaultCVUniverse>;

// The std::pair<>s for those std::map<>s don't seem to be generated though.
template class std::pair< std::string, PlotUtils::MnvLatErrorBand* >;
template class std::pair< std::string, PlotUtils::MnvVertErrorBand* >;
template class std::pair< std::string, TMatrixT<double>* >;
template class std::pair< std::string, PlotUtils::MnvVertErrorBand2D* >;
template class std::pair< std::string, PlotUtils::MnvLatErrorBand2D* >;
template class std::pair< std::string, PlotUtils::MnvVertErrorBand3D* >;
template class std::pair< std::string, PlotUtils::MnvLatErrorBand3D* >;

// Use extern keyword because these functions are instantiated in FluxReweighter.cxx already
extern template void PlotUtils::FluxReweighter::AddFluxErrorBand<PlotUtils::MnvH1D>(PlotUtils::MnvH1D*);
extern template void PlotUtils::FluxReweighter::AddFluxErrorBand<PlotUtils::MnvH2D>(PlotUtils::MnvH2D*);
extern  template MnvH1D* FluxReweighter::GetIntegratedFluxReweighted<MnvH1D>( int nuPDG,
                                                                 MnvH1D* template_hist,
                                                                 double min_energy,
                                                                 double max_energy,
                                                                 bool use_muon_correlations);
extern template MnvH2D* FluxReweighter::GetIntegratedFluxReweighted<MnvH2D>( int nuPDG,
                                                                 MnvH2D* template_hist,
                                                                 double min_energy,
                                                                 double max_energy,
                                                                 bool use_muon_correlations);

#endif
#endif // PLOTUTILSDICT_H

