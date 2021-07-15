#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ class MAT::MnvLatErrorBand+;
#pragma link C++ class MAT::MnvVertErrorBand+;
#pragma link C++ class MAT::MnvH1D+;
#pragma link C++ class MAT::MnvH2D+;
#pragma link C++ class MAT::MnvH3D+;
#pragma link C++ class MAT::TreeWrapper+;
#pragma link C++ class MAT::ChainWrapper+;
#pragma link C++ class MAT::GridCanvas+;
#pragma link C++ class MAT::MnvLatErrorBand2D+;
#pragma link C++ class MAT::MnvLatErrorBand3D+;
#pragma link C++ class MAT::MnvVertErrorBand2D+;
#pragma link C++ class MAT::MnvVertErrorBand3D+;

//NSF STUFF
#pragma link C++ class MAT::DefaultCVUniverse+;
#pragma link C++ class MAT::HistWrapper<MAT::DefaultCVUniverse>+;
#pragma link C++ class MAT::Hist2DWrapper<MAT::DefaultCVUniverse>+;
#pragma link C++ class MAT::FluxUniverse<MAT::DefaultCVUniverse>+;
#pragma link C++ class MAT::ResponseUniverse<MAT::DefaultCVUniverse>+;
#pragma link C++ class MAT::BeamAngleXUniverse<MAT::DefaultCVUniverse>+;
#pragma link C++ class MAT::BeamAngleYUniverse<MAT::DefaultCVUniverse>+;
#pragma link C++ class MAT::Universe2p2h<MAT::DefaultCVUniverse>+;
#pragma link C++ class MAT::RPAUniverse<MAT::DefaultCVUniverse>+;
#pragma link C++ class MAT::MuonUniverse<MAT::DefaultCVUniverse>+;
#pragma link C++ class MAT::GenieUniverse<MAT::DefaultCVUniverse>+;
#pragma link C++ class MAT::LowQ2PionUniverse<MAT::DefaultCVUniverse>+;
#pragma link C++ class MAT::MinosEfficiencyUniverse<MAT::DefaultCVUniverse>+;
#pragma link C++ class MAT::FluxUniverse<MAT::GenieUniverse>+;
#pragma link C++ class MAT::FluxUniverse<MAT::MuonUniverse>+;
#pragma link C++ class MAT::FluxUniverse<MAT::MuonResolutionUniverse>+;
#pragma link C++ class MAT::MuonResolutionUniverse<MAT::DefaultCVUniverse>+;
#pragma link C++ class MAT::MuonUniverseMinerva<MAT::DefaultCVUniverse>+;
#pragma link C++ class MAT::MuonUniverseMinos<MAT::DefaultCVUniverse>+;
#pragma link C++ class MAT::Universe2p2h<MAT::DefaultCVUniverse>+;
#pragma link C++ class MAT::MinosEfficiencyUniverse<MAT::DefaultCVUniverse>+;
#pragma link C++ class MAT::TrueProtonKECutUniverse<MAT::DefaultCVUniverse>+;
#pragma link C++ class MAT::RecoProtonKECutUniverse<MAT::DefaultCVUniverse>+;
#pragma link C++ class MAT::ProtonScoreUniverse<MAT::DefaultCVUniverse>+;
#pragma link C++ class MAT::GeantHadronUniverse<MAT::DefaultCVUniverse>+;

#pragma link C++ class MAT::MnvFluxConstraint<MAT::MnvH1D,MAT::MnvVertErrorBand>+;
#pragma link C++ class MAT::MnvFluxConstraint<MAT::MnvH2D,MAT::MnvVertErrorBand2D>+;
#endif
