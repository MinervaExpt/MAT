#ifndef VARIABLEBASE_H
#define VARIABLEBASE_H

#include <algorithm>
#include <string>
#include <vector>
#include <functional>


//==============================================================================
//
// VARIABLE BASE CLASS
//
//==============================================================================
namespace PlotUtils {

#ifndef __CINT__  // CINT doesn't know about std::function
template <class UNIVERSE>
class VariableBase {
 private:
  //============================================================================
  // Typedef for a pointer to a UNIVERSE::GetValue function,
  // e.g. CVUniverse::GetEnu
  //============================================================================
  // typedef for a function that:
  // * returns a double,
  // * is a member of UNIVERSE, and
  // * acts on a `this` which is a const reference to a UNIVERSE.
  // HMS 8-10-2022 add support for different true and reco binning
  typedef std::function<double(const UNIVERSE&)> PointerToCVUniverseFunction;

 public:
  //============================================================================
  // CTORS
  //============================================================================
  // variable binning
  VariableBase(const std::string name, const std::string xaxis_label,
               const std::vector<double> binning,
               PointerToCVUniverseFunction reco_func = &UNIVERSE::GetDummyVar,
               PointerToCVUniverseFunction true_func = &UNIVERSE::GetDummyVar);
  
  VariableBase(const std::string name, const std::string xaxis_label,
               const std::vector<double> binning,
               const std::vector<double> reco_binning,
               PointerToCVUniverseFunction reco_func = &UNIVERSE::GetDummyVar,
               PointerToCVUniverseFunction true_func = &UNIVERSE::GetDummyVar);

  // uniform binning
  VariableBase(const std::string name, const std::string xaxis_label,
               const int nbins, const double xmin, const double xmax,
               PointerToCVUniverseFunction reco_func = &UNIVERSE::GetDummyVar,
               PointerToCVUniverseFunction true_func = &UNIVERSE::GetDummyVar);
  
  VariableBase(const std::string name, const std::string xaxis_label,
               const int nbins, const double xmin, const double xmax,
               const int nrecobins, const double xminreco, const double xmaxreco,
               PointerToCVUniverseFunction reco_func = &UNIVERSE::GetDummyVar,
               PointerToCVUniverseFunction true_func = &UNIVERSE::GetDummyVar);

  //============================================================================
  // Getters
  //============================================================================
  std::string GetName() const;
  std::string GetAxisLabel() const;
  int GetNBins() const;
  int GetNRecoBins() const;
  bool HasRecoBinning (){ return m_has_reco_binning;} // flag to tell if you need to bother with reco bins
  bool SetRecoBinning (const bool has){m_has_reco_binning = has;}
  std::vector<double> GetBinVec() const;
  std::vector<double> GetRecoBinVec() const;
  void PrintBinning() const;
  void PrintRecoBinning() const;

  //============================================================================
  // GetValue
  //============================================================================
  virtual double GetRecoValue(const UNIVERSE& universe, const int idx1 = -1,
                              const int idx2 = -1) const;
  virtual double GetTrueValue(const UNIVERSE& universe, const int idx1 = -1,
                              const int idx2 = -1) const;
 protected:
  // hms - make protected so can overload in derived classes without introducing dependencies
  std::string m_xaxis_label;
  std::string m_name;
  bool m_has_reco_binning;
  std::vector<double> m_binning;
  std::vector<double> m_reco_binning;
  PointerToCVUniverseFunction m_pointer_to_GetRecoValue;
  PointerToCVUniverseFunction m_pointer_to_GetTrueValue;

  //============================================================================
  // Data Members
  //============================================================================
  
  VariableBase();  // Default is off-limits for the time being

  // Helper functions
  std::vector<double> MakeUniformBinning(const int nbins, const double min,
                                         const double max);
  
  std::vector<double> SplitBinning(const std::vector<double> binning, const int split);
  
  std::vector<double> GetSortedVector(const std::vector<double>& vin);
};
#endif  // __CINT__

}  // namespace PlotUtils

//==============================================================================
//
// EXCLUSIVE VARIABLE (1 ARG) BASE CLASS
//
//==============================================================================
namespace PlotUtils {
#ifndef __CINT__  // CINT doesn't know about std::function
template <class UNIVERSE, class VARIABLE>
class ExclusiveVariable1Arg : public VARIABLE {
 private:
  typedef std::function<double(const UNIVERSE&, int)>
      PointerToCVUniverse1ArgFunction;

 public:
  //============================================================================
  // CTORS
  //============================================================================
  // uniform binning
  ExclusiveVariable1Arg(
      const std::string name, const std::string xaxis_label, const int nbins,
      const double xmin, const double xmax,
      PointerToCVUniverse1ArgFunction reco_func = &UNIVERSE::GetDummyVar1Arg,
      PointerToCVUniverse1ArgFunction true_func = &UNIVERSE::GetDummyVar1Arg);

  // variable binning
  ExclusiveVariable1Arg(
      const std::string name, const std::string xaxis_label,
      const std::vector<double> binning,
      PointerToCVUniverse1ArgFunction reco_func = &UNIVERSE::GetDummyVar1Arg,
      PointerToCVUniverse1ArgFunction true_func = &UNIVERSE::GetDummyVar1Arg);

  //============================================================================
  // GetValue
  //============================================================================
  virtual double GetRecoValue(const UNIVERSE& universe, const int idx1 = -1,
                              const int idx2 = -1) const;
  virtual double GetTrueValue(const UNIVERSE& universe, const int idx1 = -1,
                              const int idx2 = -1) const;

 private:
  //============================================================================
  // Data members
  //============================================================================
  PointerToCVUniverse1ArgFunction pointer_to_GetExclusiveRecoValue;
  PointerToCVUniverse1ArgFunction pointer_to_GetExclusiveTrueValue;
};
#endif  // __CINT__
}  // namespace PlotUtils

//==============================================================================
//
// EXCLUSIVE VARIABLE (2 ARG) BASE CLASS
//
//==============================================================================
namespace PlotUtils {
#ifndef __CINT__  // CINT doesn't know about std::function
template <class UNIVERSE, class VARIABLE>
class ExclusiveVariable2Arg : public VARIABLE {
 private:
  typedef std::function<double(const UNIVERSE&, int, int)>
      PointerToCVUniverse2ArgFunction;

 public:
  //============================================================================
  // CTORS
  //============================================================================
  // uniform binning
  ExclusiveVariable2Arg(
      const std::string name, const std::string xaxis_label, const int nbins,
      const double xmin, const double xmax,
      PointerToCVUniverse2ArgFunction reco_func = &UNIVERSE::GetDummyVar2Arg,
      PointerToCVUniverse2ArgFunction true_func = &UNIVERSE::GetDummyVar2Arg);

  // variable binning
  ExclusiveVariable2Arg(
      const std::string name, const std::string xaxis_label,
      const std::vector<double> binning,
      PointerToCVUniverse2ArgFunction reco_func = &UNIVERSE::GetDummyVar2Arg,
      PointerToCVUniverse2ArgFunction true_func = &UNIVERSE::GetDummyVar2Arg);

  //============================================================================
  // GetValue
  //============================================================================
  virtual double GetRecoValue(const UNIVERSE& universe, const int idx1 = -1,
                              const int idx2 = -1) const;
  virtual double GetTrueValue(const UNIVERSE& universe, const int idx1 = -1,
                              const int idx2 = -1) const;

 private:
  //============================================================================
  // Data members
  //============================================================================
  PointerToCVUniverse2ArgFunction pointer_to_GetExclusiveRecoValue;
  PointerToCVUniverse2ArgFunction pointer_to_GetExclusiveTrueValue;
};
#endif  // __CINT__
}  // namespace PlotUtils

#include "VariableBase.cxx"

#endif  // VARIABLEBASE_H
