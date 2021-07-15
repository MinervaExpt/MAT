//File: Cut.h
//Brief: A toolkit for defining selection criteria that can be shared
//       across analyses.  A MAT::Cut takes a CVUniverse and optionally
//       an EVENT and returns a boolean about whether this combination
//       passes the Cut.  A MAT::Cut can also be used to generate a
//       table summarizing cut performance.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef PLOTUTILS_CUT_H
#define PLOTUTILS_CUT_H

//Local includes
//#include "PlotUtils/Table.h"

//c++ includes
#include <memory> //std::unqiue_ptr<>
#include <algorithm> //std::all_of<>()
#include <vector>
#include <string> 

#ifdef __GCCXML__
#define override
#endif

namespace MAT
{
  namespace detail
  {
    struct empty {};
  }

  //Intended for reconstructed quantities
  template <class UNIVERSE, class EVENT = detail::empty>
  class Cut
  {
    public:
    Cut(const std::string& name): m_Name(name), m_SignalPassed(0), m_TotalPassed(0) {}
#ifndef __GCCXML__
      virtual ~Cut() = default;
#endif
      //Check whether the pair (univ, evt) passes this Cut.  Keep
      //some handy statistics about events that passed.
      bool passesCut(const UNIVERSE& univ, EVENT& evt, const double weight = 1, const bool isSignal = true)
      {
        const bool result = checkCut(univ, evt);
        m_SignalPassed += isSignal * result * weight;
        m_TotalPassed += result * weight;
        return result;
      }

      inline const std::string& getName() const { return m_Name; }
      inline double getSignalPassed() const { return m_SignalPassed; }
      inline double getTotalPassed() const { return m_TotalPassed; }

      void resetStats()
      {
        m_SignalPassed = 0;
        m_TotalPassed = 0;
      }      

    private:
      //Override this function from a base class to use Cut.
    virtual bool checkCut(const UNIVERSE& univ, EVENT& evt) const =0;

      //Metrics for summarizing Cut performance
      const std::string m_Name; //Name of this Cut for the summary
      double m_SignalPassed; //Sum of event weights for which passesCut() returned true
      double m_TotalPassed; //Number of times a UNIVERSE, EVENT pair passed this Cut
  };

#ifndef __GCCXML__
  //Alias for a container of Cut<> because that's too long!
  template <class UNIVERSE, class EVENT = detail::empty>
  using cuts_t = std::vector<std::unique_ptr<Cut<UNIVERSE, EVENT> > >;
#endif
  //Intended for truth quantities.  These should always work on the Truth tree.
  template <class UNIVERSE>
  class SignalConstraint
  {
    public:
#ifndef __GCCXML__
      virtual ~SignalConstraint() = default;
#endif
      SignalConstraint(const std::string& name): m_Name(name),m_TotalPassed(0) {}

      virtual bool passes(const UNIVERSE& univ, const double weight = 0)
      {
        const bool result = checkConstraint(univ);
        //If we ever need to keep statistics here, we don't need to change the user inte
        //Default value of weight is 0 to avoid double counting 
        m_TotalPassed += result * weight;
        return result;
      }

      inline const std::string& getName() const { return m_Name; }
      inline double getTotalPassed() const { return m_TotalPassed; }

      void resetStats()
      {
        m_TotalPassed = 0;
      }

    private:
      //Override this function from a base class to use SignalConstraint
      virtual bool checkConstraint(const UNIVERSE& univ) const = 0;
      const std::string m_Name; //Name of this SignalConstraint for the summary
      double m_TotalPassed; //Number of times a UNIVERSE, EVENT pair passed this Constraint  
  };

#ifndef __GCCXML__
  //Alias for a container of SignalConstraint<> because that's too long!
  template <class UNIVERSE>
  using constraints_t = std::vector<std::unique_ptr<SignalConstraint<UNIVERSE> > >;
#endif

  //Classes are complicated, and you probably only want to cut on functions from
  //your CVUniverse anyway, right?  Let me wrap over a member function pointer
  //for you.
  template <class UNIVERSE, class var_t, var_t(UNIVERSE::*var)() const, class EVENT = detail::empty>
  class Maximum: public Cut<UNIVERSE, EVENT>
  {
    public:
      Maximum(const var_t max, const std::string& name): Cut<UNIVERSE, EVENT>(name), m_Max(max)
      {
      }
#ifndef __GCCXML__
      virtual ~Maximum() = default;
#endif
    private:
      bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
      {
        return (univ.*var)() <= m_Max;
      }

      const var_t m_Max; //Maximum value of UNIVERSE::var that passes this Cut
  };

  template <class UNIVERSE, class var_t, var_t(UNIVERSE::*var)() const, class EVENT = detail::empty>
  class Minimum: public Cut<UNIVERSE, EVENT>
  {
    public:
      Minimum(const var_t min, const std::string& name): Cut<UNIVERSE, EVENT>(name), m_Min(min)
      {
      }
#ifndef __GCCXML__
      virtual ~Minimum() = default;
#endif
    private:
      bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
      {
        return (univ.*var)() >= m_Min;
      }

      const var_t m_Min; //Minimum value of UNIVERSE::var that passes this Cut
  };

  template <class UNIVERSE, class var_t, var_t(UNIVERSE::*var)() const, class EVENT = detail::empty>
  class IsSame: public Cut<UNIVERSE, EVENT>
  {
    public:
      IsSame(const var_t matches, const std::string& name): Cut<UNIVERSE, EVENT>(name), m_Matches(matches)
      {
      }
#ifndef __GCCXML__
      virtual ~IsSame() = default;
#endif
    private:
      bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
      {
        return (univ.*var)() == m_Matches;
      }
                                                                                
      const var_t m_Matches; //Maximum value of UNIVERSE::var that passes this Cut
  };
}
#ifdef __GCCXML__
#undef override
#endif

#endif //PLOTUTILS_CUT_H
