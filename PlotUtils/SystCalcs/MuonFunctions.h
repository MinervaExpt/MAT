//==============================================================================
/*! @brief Minerva-standard muon kinematic and systematic functions. Include
           this file inside of your User::CVUniverse class definition.

Public-use functions here defined:

- double GetPmu() const (MeV)
- double GetEmu() const (MeV)
- double GetThetaXmu() const (radians)
- double GetThetaYmu() const (radians)
- double GetThetamu() const (radians w.r.t. incident nu dirn)
- double GetPhimu() const (radians w.r.t. incident nu dirn)
- ROOT::Math::PxPyPzEVector GetMuon4V() const (MeV)
- double GetVertexZ() const
- double GetMinosEfficiencyWeight() const

*/
//==============================================================================

#ifndef MUONFUNCTIONS_H
#define MUONFUNCTIONS_H
//==============================================================================
// Helper functions – not intended for public use
//==============================================================================
// This functionality should *not* be used by standard users. This will allow
// data preservation users in the future to offset the muon energy scale in data
// when it is no longer feasible to reproduce ntuples. So, the warning is
// inverted in that the user will be warned only if they _do_ use the method to
// specify an offset. n.b. 'offset' refers to CV only, and 'shift' refers to
// systematic universes.
static double GetCVMuonMomentumOffset() {
  if (_is_muon_momentum_cv_offset_set && !_has_muon_error_been_thrown) {
    std::cout << "WARNING: YOU HAVE SET A NONZERO CV MUON MOMENTUM SHIFT. "
              << "IF YOU DON'T KNOW WHAT THIS MEANS, "
              << "YOU SHOULDN'T BE DOING THIS! AND PROBABLY ROB IS THE ONLY "
              << "PERSON WHO SHOULD BE DOING THIS." << std::endl;
    _has_muon_error_been_thrown = true;
  }
  return m_muon_momentum_cv_offset;
}

static bool SetCVMuonMomentumOffset(double muon_momentum_cv_offset) {
  if (_is_muon_momentum_cv_offset_set) {
    std::cout << "WARNING: YOU ATTEMPTED SETTING THE MUON MOMENTUM CV SHIFT "
              << "A SECOND TIME. THIS IS NOT ALLOWED FOR CONSISTENCY."
              << std::endl;
    return false;
  } else {
    m_muon_momentum_cv_offset = muon_momentum_cv_offset;
    _is_muon_momentum_cv_offset_set = true;
  }
  return true;
}

// Muon momenta functions used under the hood to calculate GetPmu
// Get total muon momentum w/o CV shift
virtual double GetPmu_nominal() const {                       /* MeV */
  std::string lepton_branch = GetAnaToolName() + "_leptonE";  // AnaTool
  // std::string lepton_branch = "Muon_leptonE";
  double px = GetVecElem(lepton_branch.c_str(), 0);
  double py = GetVecElem(lepton_branch.c_str(), 1);
  double pz = GetVecElem(lepton_branch.c_str(), 2);
  double total_p_nominal = sqrt(px * px + py * py + pz * pz);

  return total_p_nominal;
}

virtual double GetPmuMinos_nominal() const { /* MeV */
  // std::string minos_momentum_branch = "Muon_minos_trk_p";
  std::string minos_momentum_branch = GetAnaToolName() + "_minos_trk_p";
  double minos_p_nominal = GetDouble(minos_momentum_branch.c_str());

  return minos_p_nominal;
}

virtual double GetPmuMinerva() const { /* MeV */
  // this method can be overridden by systematics universes.
  // note that it includes a mechanism for introducing an overall
  // multiplicative factor

  double total_p_nominal = GetPmu_nominal();
  double minos_p_nominal = GetPmuMinos_nominal();
  double minerva_p = total_p_nominal - minos_p_nominal;

  return minerva_p;
}

virtual double GetPmuMinos() const { /* MeV */
  // this method can be overridden by systematics universes.
  // note that it includes a mechanism for introducing an overall
  // multiplicative factor
  double minos_p_nominal = GetPmuMinos_nominal();

  // this CV offset is set by CVUniverse.SetCVMuonMomentumOffset method
  double momentum_offset = GetCVMuonMomentumOffset();
  if (momentum_offset == 0.0) {
    return minos_p_nominal;
  } else {
    double minos_p_shifted = minos_p_nominal * (1. + momentum_offset);
    return minos_p_shifted;
  }
}

//==============================================================================
// User, public functions
//==============================================================================
virtual double GetPmu() const { /* MeV */
  double Pmu_minerva = GetPmuMinerva();
  double Pmu_minos = GetPmuMinos();

  return Pmu_minerva + Pmu_minos;
}

virtual double GetEmu() const { /* MeV */
  double m_mu = MinervaUnits::M_mu;
  double pmu = GetPmu();
  return sqrt(m_mu * m_mu + pmu * pmu);
}

virtual double GetThetaXmu() const { /* radians */
  // Default offset is zero - but allows adding resolution
  std::string thetaBranch = "muon_thetaX";
  return GetDouble(thetaBranch.c_str()) + GetBeamAngleOffsetX();
}

virtual double GetThetaYmu() const { /* radians */
  // Default offset is zero - but allows adding resolution
  std::string thetaBranch = "muon_thetaY";
  return GetDouble(thetaBranch.c_str()) + GetBeamAngleOffsetY();
}

//! Theta is a derived quantity
virtual double GetThetamu() const { /* radians w.r.t. incident nu dirn */
  return theta3D(GetThetaXmu(), GetThetaYmu());
}

virtual double GetPhimu() const { /* radians w.r.t. incident nu dirn */
  return phi3D(GetThetaXmu(), GetThetaYmu());
}

//! 4 vector
virtual ROOT::Math::PxPyPzEVector GetMuon4V() const { /* MeV */
  double p = GetPmu();
  double thetax = GetThetaXmu();
  double thetay = GetThetaYmu();
  double theta = GetThetamu();
  double px = p * std::sin(thetax);
  double py = GetPmu() * std::sin(thetay);
  double pz = GetPmu() * std::cos(theta);
  double E = GetEmu();
  return ROOT::Math::PxPyPzEVector(px, py, pz, E);
}

// Vtx Z position -- required for muon mom systs with recalculation of shifts
virtual double GetVertexZ() const {
  double vertex_Z = GetVecElem("vtx", 2);
  return vertex_Z;
}

virtual double GetMinosEfficiencyWeight() const {
  if (m_is_truth) return 1.;  // No efficiency reweigting for truth events
  if (!isFHC() && !isRHC())
    return 1.;  // LE or nonstandard playlist.  Don't have away of dealing with
                // this yet
  if (IsPlaylistME(GetPlaylist())) {  // The correction factors are different
                                      // between ME and LE
    double pmu = GetPmuMinos() / 1000;  // GetCorrection expects GeV
    return PlotUtils::MinosMuonEfficiencyCorrection::Get(isFHC()).GetCorrection(
        pmu, GetBatchPOT(), isFHC());
  } else {                       // Assume if not ME, then it's LE
    double pmu = GetPmuMinos();  // MnVnormalizer GetCorrection expects MeV
#ifndef __CINT__
    static PlotUtils::MnvNormalizer mnvNormalizer =
        PlotUtils::MnvNormalizer("Eroica", GetPlaylist());
#endif // __CINT__
    return mnvNormalizer.GetCorrection(pmu);
  }
}

#endif  // MUONFUNCTIONS
