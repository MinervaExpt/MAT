#ifndef VARIABLE2DBASE_H
#define VARIABLE2DBASE_H

#include <memory>

#include "PlotUtils/VariableBase.h"

namespace PlotUtils {

#ifndef __CINT__
template <class UNIVERSE>
class Variable2DBase {
 public:
  //============================================================================
  // CTORS
  //============================================================================
  Variable2DBase(const VariableBase<UNIVERSE>& x,
                 const VariableBase<UNIVERSE>& y);
  Variable2DBase(const std::string name,
                 const VariableBase<UNIVERSE>& x,
                 const VariableBase<UNIVERSE>& y);

 public:
  //============================================================================
  // Setters/Getters
  //============================================================================
  std::string SetName(const std::string name);
  std::string GetAxisLabel() const;
  std::string GetName() const;
  std::string GetNameX() const;
  std::string GetNameY() const;
  int GetNBinsX() const;
  int GetNBinsY() const;
  std::vector<double> GetBinVecX() const;
  std::vector<double> GetBinVecY() const;
  void PrintBinningX() const;
  void PrintBinningY() const;

  //============================================================================
  // Get Value
  //============================================================================
  double GetRecoValueX(const UNIVERSE& universe, const int idx1 = -1,
                              const int idx2 = -1) const;
  double GetRecoValueY(const UNIVERSE& universe, const int idx1 = -1,
                              const int idx2 = -1) const;
  double GetTrueValueX(const UNIVERSE& universe, const int idx1 = -1,
                              const int idx2 = -1) const;
  double GetTrueValueY(const UNIVERSE& universe, const int idx1 = -1,
                              const int idx2 = -1) const;

 protected:
  std::string m_xaxis_label;
  std::string m_yaxis_label;

 private:
  //============================================================================
  // Member Variables
  //============================================================================
  std::string m_name;

  std::unique_ptr<VariableBase<UNIVERSE>> m_var_x;
  std::unique_ptr<VariableBase<UNIVERSE>> m_var_y;

  Variable2DBase();  // off-limits
};
#endif

}  // namespace PlotUtils

#include "Variable2DBase.cxx"

#endif  // VARIABLE2DBASE_H
