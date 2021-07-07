#ifndef RECOILENERGYFUNCTIONS_H
#define RECOILENERGYFUNCTIONS_H

// Particle response systematic quantity
std::vector<int> m_non_cal_idx;

// Particle response systematics functions
virtual void SetNonCalIndices(std::vector<int>& non_cal_idx) {
  m_non_cal_idx = non_cal_idx;
}
virtual std::vector<int> GetNonCalIndices() { return m_non_cal_idx; }

// Recoil Energy
// Analyzer should specify what their recoil energy is (in MeV) in CVUniverse)
double GetRecoilEnergy() const {
  return GetCalRecoilEnergy() + GetNonCalRecoilEnergy();  //  LOOK AT PlotUtils/SystCalcs/RecoilEnergyFunctions.h
                                                          //  If you're here because GetCalRecoilEnergy or GetNonCalRecoilEnergy doesn't eixt
                                                          //  Please look below
}

////////////////////////////////////////////////////////////////////////
// LOOK HERE
// Need to define these functions in your personal CVUniverse

//virtual double GetCalRecoilEnergy() const {
//  std::cout << "GetCalRecoilEnergy() should be implemented in CVUniverse to"<< std::endl;
//  std::cout << "     return all recoil energy found calorimetrically (MeV)"
//  std::cout << " If all of your recoil energy is find calorimetrically (such as E available
//  std::cout << "   or by spline correction, just put that branch/ calculation here"
//            << std::endl;
//      
//  return calorimetric_recoil_energy;
//}
//
//virtual double GetNonCalRecoilEnergy() const {
//  std::cout << "GetNonCalRecoilEnergy() should be implemented in CVUniverse to"<< std::endl;
//  std::cout << "     return all recoil energy not found calorimetrically like dEdX (MeV)"
//            << std::endl;
//  return non_calorimetric_recoil_energy;
//}

#endif  // RECOILENERGYFUNCTIONS
