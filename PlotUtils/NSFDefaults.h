#ifndef MINERVAPHYSICALCONSTANTS_H
#define MINERVAPHYSICALCONSTANTS_H 1


//! Place to put default systematics values for New Systematics Framework

namespace NSFDefaults {
//! Helicities
  static const int AntiNuHelicity  = 2;
  static const int NuHelicity = 1;

//! beam offset and uncertainty
  static const double beamThetaX_Err = 0.001;  // rad (error on beam angleX correction - E. Valencia's - DocDB 11550)
  static const double beamThetaY_Err = 0.0009; // rad (error on beam angleY correction - E. Valencia's - DocDB 11550)
  static const double beamThetaX_Central = 0.0;
  static const double beamThetaY_Central = 0.0;

//! Muon "resolution"
  static const double muonResolution_Err = 0.004; // Fractional

//! muon track angle resolution. B. Zeimer (DocDB 7635), but future talks with Mike K lowered the uncertainty to 2mrad
  //Only use for muon tracks!
  //static const double muon_angle_res = 0.002; //radians

//! Further scrunity of the plots in the above DocDB has come with a fraction uncertainty on the angle of a conservative 2%
  static const double muon_angle_frac_res = 0.02;// Fractional

//! Muon range and curve - from Ana/AnaUtils/src/MuonUtils.cpp
  static const double MinosMuonPRange_Err = 0.00984; // Fractional
  static const double MinosMuonHighPCurvature_Err = 0.006; // Fractional
  static const double MinosMuonLowPCurvature_Err = 0.025; // Fractional

//! Michel tag efficiency (applied to weight) (docdb 27916)
  static const double MichelTagEfficiency_Err = 0.025; // Fractional

//! material assay
  static const double MaterialAssayUncertainty = 17.; // MeV
  static const double MaterialAssayUncertaintyNoNuke = 11.; // MeV

//! Target mass uncertainties (DocDB #6016-v7)
  static const double c_err  = 0.005;// 0.5%
  static const double fe_err = 0.010;// 1%
  static const double pb_err = 0.005;// 0.5%
  static const double ch_err = 0.014;// 1.4%
  static const double h2o_err = 0.02;// 2% (H. Budd, docdb tbd)

//! Pion production parameters (docdb 7451 and 19774)
//! GENIE central parameters
  static const double GENIE_MvRES         = 0.84;
  static const double GENIE_MvRES_1Sig    = 0.1 * GENIE_MvRES;
  static const double GENIE_MaRES         = 1.12;
  static const double GENIE_MaRES_1Sig    = 0.22; //20%
  static const double GENIE_RES_NORM      = 1.00; //100%
  static const double GENIE_RES_NORM_1Sig = 0.20; //20%
//! Reduced MvRES error from pion electroproduction data fit
  static const double ELECTROPROD_MvRES_1Sig     = 0.03*GENIE_MvRES;
//! Pion production parameters and errors from deuterium fit
  static const double DEUTERIUM_MaRES            = 0.94;
  static const double DEUTERIUM_MaRES_1Sig       = 0.05;
  static const double DEUTERIUM_RES_NORM         = 1.15;
  static const double DEUTERIUM_RES_NORM_1Sig    = 0.07;

//! dEdx
//! TODO from where? for what?
  static const double dEdxUncertainty = 40.; // MeV
  static const double dEdxUncertaintyNoNuke = 30.; // MeV

//! binding energies needed for kinematic calculations
  static const double CCQEBindingE_neutron = 34.; //MeV
  static const double CCQEBindingE_proton  = 30.; //MeV
  static const double TrueProtonKECutCentral = 120.; //MeV
  static const double TrueProtonKECutVariation = 20.; // MeV
  static const double RecoProtonKECutCentral = 120.; //MeV
  static const double RecoProtonKECutVariation = 20.; // MeV

//! nuclear target face (the face of the detector)
  static const double nuclearTargetZFace = 4293.04;

//! Boundary between nuclear target and tracker regions
  static const double nuclearTargetZEnd = 5835.0; // mm

//! default fiducial volume used by GeantHadronSystematics/MnvHadronReweighter
  static const double TrackerFace = 5991.29; // mm, module 27, plane 1
  static const double TrackerBack = 8408.91; // mm, module 80, plane 2
  static const double StandardApothem = 850; // mm

}  // end of NSFDefaults namespace

#endif
