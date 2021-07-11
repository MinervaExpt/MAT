#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ class PlotUtils::MnvLatErrorBand+;
#pragma link C++ class PlotUtils::MnvVertErrorBand+;
#pragma link C++ class PlotUtils::MnvH1D+;
#pragma link C++ class PlotUtils::MnvH2D+;
#pragma link C++ class PlotUtils::TreeWrapper+;
#pragma link C++ class PlotUtils::ChainWrapper+;
#pragma link C++ class PlotUtils::GridCanvas+;
#pragma link C++ class PlotUtils::MnvLatErrorBand2D+;
#pragma link C++ class PlotUtils::MnvVertErrorBand2D+;

//NSF STUFF
#pragma link C++ class PlotUtils::DefaultCVUniverse+;
#pragma link C++ class PlotUtils::HistWrapper<PlotUtils::DefaultCVUniverse>+;
#pragma link C++ class PlotUtils::Hist2DWrapper<PlotUtils::DefaultCVUniverse>+;
#pragma link C++ class PlotUtils::FluxUniverse<PlotUtils::DefaultCVUniverse>+;
#pragma link C++ class PlotUtils::ResponseUniverse<PlotUtils::DefaultCVUniverse>+;
#pragma link C++ class PlotUtils::BeamAngleXUniverse<PlotUtils::DefaultCVUniverse>+;
#pragma link C++ class PlotUtils::BeamAngleYUniverse<PlotUtils::DefaultCVUniverse>+;
#pragma link C++ class PlotUtils::Universe2p2h<PlotUtils::DefaultCVUniverse>+;
#pragma link C++ class PlotUtils::RPAUniverse<PlotUtils::DefaultCVUniverse>+;
#pragma link C++ class PlotUtils::MuonUniverse<PlotUtils::DefaultCVUniverse>+;
#pragma link C++ class PlotUtils::GenieUniverse<PlotUtils::DefaultCVUniverse>+;
#pragma link C++ class PlotUtils::LowQ2PionUniverse<PlotUtils::DefaultCVUniverse>+;
#pragma link C++ class PlotUtils::MinosEfficiencyUniverse<PlotUtils::DefaultCVUniverse>+;
#pragma link C++ class PlotUtils::FluxUniverse<PlotUtils::GenieUniverse>+;
#pragma link C++ class PlotUtils::FluxUniverse<PlotUtils::MuonUniverse>+;
#pragma link C++ class PlotUtils::FluxUniverse<PlotUtils::MuonResolutionUniverse>+;
#pragma link C++ class PlotUtils::MuonResolutionUniverse<PlotUtils::DefaultCVUniverse>+;
#pragma link C++ class PlotUtils::MuonUniverseMinerva<PlotUtils::DefaultCVUniverse>+;
#pragma link C++ class PlotUtils::MuonUniverseMinos<PlotUtils::DefaultCVUniverse>+;
#pragma link C++ class PlotUtils::Universe2p2h<PlotUtils::DefaultCVUniverse>+;
#pragma link C++ class PlotUtils::MinosEfficiencyUniverse<PlotUtils::DefaultCVUniverse>+;
#pragma link C++ class PlotUtils::TrueProtonKECutUniverse<PlotUtils::DefaultCVUniverse>+;
#pragma link C++ class PlotUtils::RecoProtonKECutUniverse<PlotUtils::DefaultCVUniverse>+;
#pragma link C++ class PlotUtils::ProtonScoreUniverse<PlotUtils::DefaultCVUniverse>+;
#pragma link C++ class PlotUtils::GeantHadronUniverse<PlotUtils::DefaultCVUniverse>+;

#pragma link C++ class PlotUtils::MnvFluxConstraint<PlotUtils::MnvH1D,PlotUtils::MnvVertErrorBand>+;
#pragma link C++ class PlotUtils::MnvFluxConstraint<PlotUtils::MnvH2D,PlotUtils::MnvVertErrorBand2D>+;
#endif
