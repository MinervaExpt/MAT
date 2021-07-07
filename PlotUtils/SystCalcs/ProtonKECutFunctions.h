#ifndef PROTONKECUTFUNCTIONS_H
#define PROTONKECUTFUNCTIONS_H

/*
static double
    m_true_proton_ke_cut_central;  // = NSFDefaults::TrueProtonKECutCentral;
static double
    m_reco_proton_ke_cut_central;  // = NSFDefaults::RecoProtonKECutCentral;
*/
// Set and get the proton KE cut central value for the BaseUniverse, the
// variation is set when systematics map is set up
static bool SetTrueProtonKECutCentral(
    double cut = NSFDefaults::TrueProtonKECutCentral) {
  m_true_proton_ke_cut_central = cut;
  return true;
}
//! cut for generator level protons for CCQE
virtual double GetTrueProtonKECut() const {
  return m_true_proton_ke_cut_central;
}
// The method you use to get the cut value which will be varied by
// TrueProtonKECutSystematics
static double GetTrueProtonKECutCentral() {
  return m_true_proton_ke_cut_central;
}

// Reco
static bool SetRecoProtonKECutCentral(
    double cut = NSFDefaults::RecoProtonKECutCentral) {
  m_reco_proton_ke_cut_central = cut;
  return true;
}
//! cut for generator level protons for CCQE
virtual double GetRecoProtonKECut() const {
  return m_reco_proton_ke_cut_central;
}
// The method you use to get the cut value which will be varied by
// TrueProtonKECutSystematics
static double GetRecoProtonKECutCentral() {
  return m_reco_proton_ke_cut_central;
}

#endif  // PROTONKECUTFUNCTIONS
