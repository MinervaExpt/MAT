#ifndef VARIABLEBASE_CXX
#define VARIABLEBASE_CXX

#include "PlotUtils/VariableBase.h"
#include <iostream>
using namespace PlotUtils;

//==============================================================================
//
// Variable Base
//
//==============================================================================
//==============================================================================
// CTORs
//==============================================================================
// default (off-limits at the moment)
template <class UNIVERSE>
VariableBase<UNIVERSE>::VariableBase()
    : m_name(),
      m_xaxis_label(),
      m_binning(),
      m_pointer_to_GetRecoValue(&UNIVERSE::GetDummyVar),
      m_pointer_to_GetTrueValue(&UNIVERSE::GetDummyVar) {}

// variable binning
template <class UNIVERSE>
VariableBase<UNIVERSE>::VariableBase(const std::string name,
                                     const std::string xaxis_label,
                                     const std::vector<double> binning,
                                     PointerToCVUniverseFunction reco_func,
                                     PointerToCVUniverseFunction true_func)
    : m_name(name),
      m_xaxis_label(xaxis_label),
      m_binning(GetSortedVector(binning)),
      m_pointer_to_GetRecoValue(reco_func),
      m_pointer_to_GetTrueValue(true_func) {}

// uniform binning
template <class UNIVERSE>
VariableBase<UNIVERSE>::VariableBase(const std::string name,
                                     const std::string xaxis_label,
                                     const int nbins, const double xmin,
                                     const double xmax,
                                     PointerToCVUniverseFunction reco_func,
                                     PointerToCVUniverseFunction true_func)
    : m_name(name),
      m_xaxis_label(xaxis_label),
      m_binning(MakeUniformBinning(nbins, xmin, xmax)),
      m_pointer_to_GetRecoValue(reco_func),
      m_pointer_to_GetTrueValue(true_func) {}

//==============================================================================
// Set/Get
//==============================================================================
template <class UNIVERSE>
std::string VariableBase<UNIVERSE>::GetName() const {
  return m_name;
}

template <class UNIVERSE>
std::string VariableBase<UNIVERSE>::GetAxisLabel() const {
  return m_xaxis_label;
}

template <class UNIVERSE>
int VariableBase<UNIVERSE>::GetNBins() const {
  return m_binning.size() - 1;
}

template <class UNIVERSE>
void VariableBase<UNIVERSE>::PrintBinning() const {
  std::cout << GetName() << " binning: ";
  for (const auto b : m_binning) std::cout << b << " ";
  std::cout << "\n";
}

template <class UNIVERSE>
std::vector<double> VariableBase<UNIVERSE>::GetBinVec() const {
  return m_binning;
}

//==============================================================================
// Get Reco and True Values
//==============================================================================
template <class UNIVERSE>
double VariableBase<UNIVERSE>::GetRecoValue(const UNIVERSE& universe,
                                            const int idx1,
                                            const int idx2) const {
  return m_pointer_to_GetRecoValue(universe);
}

template <class UNIVERSE>
double VariableBase<UNIVERSE>::GetTrueValue(const UNIVERSE& universe,
                                            const int idx1,
                                            const int idx2) const {
  return m_pointer_to_GetTrueValue(universe);
}

//==============================================================================
// Binning Helpers
//==============================================================================
template <class UNIVERSE>
std::vector<double> VariableBase<UNIVERSE>::MakeUniformBinning(
    const int nbins, const double min, const double max) {
  double step_size = (max - min) / nbins;
  double arr[nbins + 1];  // +1 because binning arrays include top edge.
  for (int i = 0; i <= nbins; ++i) arr[i] = min + i * step_size;
  const int size = sizeof(arr) / sizeof(*arr);
  std::vector<double> ret(arr, arr + size);
  return ret;
}

template <class UNIVERSE>
std::vector<double> VariableBase<UNIVERSE>::GetSortedVector(
    const std::vector<double>& vin) {
  std::vector<double> vout = vin;
  std::sort(vout.begin(), vout.end());
  return vout;
}

//==============================================================================
//
// ExclusiveVariable1Arg
//
//==============================================================================
// uniform binning
template <class UNIVERSE, class VARIABLE>
ExclusiveVariable1Arg<UNIVERSE, VARIABLE>::ExclusiveVariable1Arg(
    const std::string name, const std::string xaxis_label, const int nbins,
    const double xmin, const double xmax,
    PointerToCVUniverse1ArgFunction reco_func,
    PointerToCVUniverse1ArgFunction true_func)
    : VARIABLE(name, xaxis_label, nbins, xmin, xmax),
      pointer_to_GetExclusiveRecoValue(reco_func),
      pointer_to_GetExclusiveTrueValue(true_func) {}

// variable binning
template <class UNIVERSE, class VARIABLE>
ExclusiveVariable1Arg<UNIVERSE, VARIABLE>::ExclusiveVariable1Arg(
    const std::string name, const std::string xaxis_label,
    const std::vector<double> binning, PointerToCVUniverse1ArgFunction reco_func,
    PointerToCVUniverse1ArgFunction true_func)
    : VARIABLE(name, xaxis_label, binning),
      pointer_to_GetExclusiveRecoValue(reco_func),
      pointer_to_GetExclusiveTrueValue(true_func) {}

// Get the variable's value
template <class UNIVERSE, class VARIABLE>
double ExclusiveVariable1Arg<UNIVERSE, VARIABLE>::GetRecoValue(
    const UNIVERSE& universe, const int idx1, const int idx2) const {
  return pointer_to_GetExclusiveRecoValue(universe, idx1);
}

template <class UNIVERSE, class VARIABLE>
double ExclusiveVariable1Arg<UNIVERSE, VARIABLE>::GetTrueValue(
    const UNIVERSE& universe, const int idx1, const int idx2) const {
  return pointer_to_GetExclusiveTrueValue(universe, idx1);
}

//==============================================================================
//
// ExclusiveVariable2Arg
//
//==============================================================================
// uniform binning
template <class UNIVERSE, class VARIABLE>
ExclusiveVariable2Arg<UNIVERSE, VARIABLE>::ExclusiveVariable2Arg(
    const std::string name, const std::string xaxis_label, const int nbins,
    const double xmin, const double xmax,
    PointerToCVUniverse2ArgFunction reco_func,
    PointerToCVUniverse2ArgFunction true_func)
    : VARIABLE(name, xaxis_label, nbins, xmin, xmax),
      pointer_to_GetExclusiveRecoValue(reco_func),
      pointer_to_GetExclusiveTrueValue(true_func) {}

// variable binning
template <class UNIVERSE, class VARIABLE>
ExclusiveVariable2Arg<UNIVERSE, VARIABLE>::ExclusiveVariable2Arg(
    const std::string name, const std::string xaxis_label,
    const std::vector<double> binning, PointerToCVUniverse2ArgFunction reco_func,
    PointerToCVUniverse2ArgFunction true_func)
    : VARIABLE(name, xaxis_label, binning),
      pointer_to_GetExclusiveRecoValue(reco_func),
      pointer_to_GetExclusiveTrueValue(true_func) {}

// Get the variable's value
template <class UNIVERSE, class VARIABLE>
double ExclusiveVariable2Arg<UNIVERSE, VARIABLE>::GetRecoValue(
    const UNIVERSE& universe, const int idx1, const int idx2) const {
  return pointer_to_GetExclusiveRecoValue(universe, idx1, idx2);
}

template <class UNIVERSE, class VARIABLE>
double ExclusiveVariable2Arg<UNIVERSE, VARIABLE>::GetTrueValue(
    const UNIVERSE& universe, const int idx1, const int idx2) const {
  return pointer_to_GetExclusiveTrueValue(universe, idx1, idx2);
}

#endif  // VARIABLEBASE_CXX
