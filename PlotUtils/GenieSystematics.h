#ifndef GENIESYSTEMATICS_H
#define GENIESYSTEMATICS_H

#include <string>
#include "TMatrixD.h"
// Helper functions declared/defined in the .cxx
// GetStandardGenieSystematicsMap(typename T::config_t chain);
// GetGenieRvx1piSystematicsMap(typename T::config_t chain);
// GetGenieFaCCQESystematicsMap(typename T::config_t chain, int n_universes);


namespace PlotUtils{
  const double kNonResPiWeight = 0.43;
  const double kNonResPiWeightShift = 0.04;

  //=================================================================================
  // Standard Genie
  //=================================================================================
  template<class T>
  class GenieUniverse : public T
  {
    public:
      GenieUniverse(typename T::config_t chw, double nsigma, std::string genie_name);

      // I can't get things to compile with "override"
      virtual double GetGenieWeight() const /* override */;
      double GetWeightRatioToCV() const;

      virtual std::string ShortName() const /* override */;
      virtual std::string LatexName() const /* override */;
      virtual bool IsVerticalOnly() const   { return true; }/* override */;

      const std::string m_name;
      const std::string m_branch_name;
  };

  //=================================================================================
  // Genie MvRes systematic based on electroproduction data (New: 3%, Old: 10%)
  //=================================================================================
  template<class T>
  class GenieMvResUniverse : public T
  {
    public:
      GenieMvResUniverse(typename T::config_t chw, double nsigma);

      virtual double GetGenieWeight() const /* override */;
      double GetWeightRatioToCV() const;

      virtual std::string ShortName() const /* override */;
      virtual std::string LatexName() const /* override */;
      virtual bool IsVerticalOnly() const   { return true; }/* override */;

      const std::string m_name;
      const std::string m_branch_name;
  };

  //=================================================================================
  // Genie MaRes systematic based on deuterium fit (old 1.12 +/- 20%, new: 0.94 +/- 5%)
  //=================================================================================
  template<class T>
  class GenieNormCCResUniverse : public T
  {
    public:
      GenieNormCCResUniverse(typename T::config_t chw, double nsigma);

      virtual double GetGenieWeight() const /* override */;
      double GetWeightRatioToCV() const;

      virtual std::string ShortName() const /* override */;
      virtual std::string LatexName() const /* override */;
      virtual bool IsVerticalOnly() const   { return true; }/* override */;

      const std::string m_name;
      const std::string m_branch_name;
  };

  //=================================================================================
  // Genie MaRes systematic based on deuterium fit (old 1.12 +/- 20%, new: 0.94 +/- 5%)
  //=================================================================================
  template<class T>
  class GenieMaResUniverse : public T
  {
    public:
      GenieMaResUniverse(typename T::config_t chw, double nsigma);

      virtual double GetGenieWeight() const /* override */;
      double GetWeightRatioToCV() const;

      virtual std::string ShortName() const /* override */;
      virtual std::string LatexName() const /* override */;
      virtual bool IsVerticalOnly() const   { return true; }/* override */;

      const std::string m_name;
      const std::string m_branch_name;
  };

  //=================================================================================
  // Genie Rv{n,p}1pi systematics; nonrespi reweight
  //=================================================================================
  template<class T>
  class GenieRvx1piUniverse : public T
  {
    public:
      GenieRvx1piUniverse(typename T::config_t chw, double nsigma, std::string genie_name);

      virtual double GetGenieWeight() const /* override */;
      double GetWeightRatioToCV() const;

      virtual std::string ShortName() const /* override */;
      virtual std::string LatexName() const /* override */;
      virtual bool IsVerticalOnly() const   { return true; }/* override */;

      const std::string m_name;
      const std::string m_branch_name;
  };

  //=================================================================================
  // Genie Fitted MaRes+NormRes systematic w/ Cov Matrix
  //=================================================================================
  template<class T>
  class GenieMaNormResCovUniverse : public T
  {
    public:
      GenieMaNormResCovUniverse(typename T::config_t chw, double nsigma, int row);

      virtual double GetGenieWeight() const /* override */;
      double GetWeightRatioToCV() const;

      virtual std::string ShortName() const /* override */;
      virtual std::string LatexName() const /* override */;
      virtual bool IsVerticalOnly() const   { return true; }/* override */;

      int m_row;
      std::vector<double> m_factors; 
      const std::string m_ma_branch_name;  
      const std::string m_norm_branch_name;
  };

  //=================================================================================
  // Genie FaCCQE (MaCCQE) systematics; z-expansion reweight
  //=================================================================================
  template<class T>
  class GenieFaCCQEUniverse : public T
  {
    public:
      GenieFaCCQEUniverse(typename T::config_t ch, double nsigma, int universe_number = -1);

      virtual double GetGenieWeight() const /* override */;
      double GetWeightRatioToCV() const;

      virtual std::string ShortName() const /* override */;
      virtual std::string LatexName() const /* override */;
      virtual bool IsVerticalOnly() const   { return true; }/* override */;

      int m_universe_number;
  };

}

// Explicit Instantiation, needed for loading into python
// See Ana/PlotUtils/dict/PlotUtilsDict.h

// For template classes, the header needs to know about the definitions
#include "GenieSystematics.cxx"

#endif // GENIESYSTEMATICS_H
