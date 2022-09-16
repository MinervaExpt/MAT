#ifndef VARIABLE2DBASE_CXX
#define VARIABLE2DBASE_CXX

#include "PlotUtils/Variable2DBase.h"

using namespace PlotUtils;

//==============================================================================
//
// Variable 2D Base
//
//==============================================================================
//==============================================================================
// CTORS
//==============================================================================
template <class UNIVERSE>
Variable2DBase<UNIVERSE>::Variable2DBase(const VariableBase<UNIVERSE>& x,
                                         const VariableBase<UNIVERSE>& y)
    : m_name(x.GetName() + "_" + y.GetName()),
      m_xaxis_label(x.GetAxisLabel()),
      m_yaxis_label(y.GetAxisLabel()),
      m_var_x(new VariableBase<UNIVERSE>(x)),
      m_var_y(new VariableBase<UNIVERSE>(y)) {}

template <class UNIVERSE>
Variable2DBase<UNIVERSE>::Variable2DBase(const std::string name,
                                         const VariableBase<UNIVERSE>& x,
                                         const VariableBase<UNIVERSE>& y)
    : m_name(name),
      m_xaxis_label(x.GetAxisLabel()),
      m_yaxis_label(y.GetAxisLabel()),
      m_var_x(new VariableBase<UNIVERSE>(x)),
      m_var_y(new VariableBase<UNIVERSE>(y)) {}

//==============================================================================
// Set/Get
//==============================================================================
template <class UNIVERSE>
std::string Variable2DBase<UNIVERSE>::SetName(const std::string name) {
  return m_name = name;
}

template <class UNIVERSE>
std::string Variable2DBase<UNIVERSE>::GetName() const {
  return m_name;
}

template <class UNIVERSE>
std::string Variable2DBase<UNIVERSE>::GetAxisLabel() const {
  return m_xaxis_label;
}

template <class UNIVERSE>
std::string Variable2DBase<UNIVERSE>::GetNameX() const {
  return m_var_x->GetName();
}

template <class UNIVERSE>
std::string Variable2DBase<UNIVERSE>::GetNameY() const {
  return m_var_y->GetName();
}

template <class UNIVERSE>
int Variable2DBase<UNIVERSE>::GetNBinsX() const {
  return m_var_x->GetNBins();
}

template <class UNIVERSE>
int Variable2DBase<UNIVERSE>::GetNBinsY() const {
  return m_var_y->GetNBins();
}

template <class UNIVERSE>
std::vector<double> Variable2DBase<UNIVERSE>::GetBinVecX() const {
  return m_var_x->GetBinVec();
}

template <class UNIVERSE>
std::vector<double> Variable2DBase<UNIVERSE>::GetBinVecY() const {
  return m_var_y->GetBinVec();
}

template <class UNIVERSE>
void Variable2DBase<UNIVERSE>::PrintBinningX() const {
  m_var_x->PrintBinning();
}

template <class UNIVERSE>
void Variable2DBase<UNIVERSE>::PrintBinningY() const {
  m_var_y->PrintBinning();
}

template <class UNIVERSE>
int Variable2DBase<UNIVERSE>::GetNRecoBinsX() const {
  return m_var_x->GetNRecoBins();
}

template <class UNIVERSE>
int Variable2DBase<UNIVERSE>::GetNRecoBinsY() const {
  return m_var_y->GetNRecoBins();
}

template <class UNIVERSE>
std::vector<double> Variable2DBase<UNIVERSE>::GetRecoBinVecX() const {
  return m_var_x->GetRecoBinVec();
}

template <class UNIVERSE>
std::vector<double> Variable2DBase<UNIVERSE>::GetRecoBinVecY() const {
  return m_var_y->GetRecoBinVec();
}

template <class UNIVERSE>
void Variable2DBase<UNIVERSE>::PrintRecoBinningX() const {
  m_var_x->PrintRecoBinning();
}

template <class UNIVERSE>
void Variable2DBase<UNIVERSE>::PrintRecoBinningY() const {
  m_var_y->PrintRecoBinning();
}

//==============================================================================
// GetValues
//==============================================================================
template <class UNIVERSE>
double Variable2DBase<UNIVERSE>::GetRecoValueX(const UNIVERSE& universe,
                                               const int idx1,
                                               const int idx2) const {
  return m_var_x->GetRecoValue(universe, idx1, idx2);
}

template <class UNIVERSE>
double Variable2DBase<UNIVERSE>::GetRecoValueY(const UNIVERSE& universe,
                                               const int idx1,
                                               const int idx2) const {
  return m_var_y->GetRecoValue(universe, idx1, idx2);
}

template <class UNIVERSE>
double Variable2DBase<UNIVERSE>::GetTrueValueX(const UNIVERSE& universe,
                                               const int idx1,
                                               const int idx2) const {
  return m_var_x->GetTrueValue(universe, idx1, idx2);
}

template <class UNIVERSE>
double Variable2DBase<UNIVERSE>::GetTrueValueY(const UNIVERSE& universe,
                                               const int idx1,
                                               const int idx2) const {
  return m_var_y->GetTrueValue(universe, idx1, idx2);
}

#endif  // VARIABLE2DBASE_CXX
