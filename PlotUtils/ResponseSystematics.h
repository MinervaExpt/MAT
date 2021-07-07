#ifndef RESPONSESYSTEMATICS_H
#define RESPONSESYSTEMATICS_H

#include "TreeWrapper.h"


/*
This implementation doesn't include the Neutron Systematics. 
Currently Neutron Systematics is supposed to handled by Hadron Reweight systematics. 
If Hadron Reweight systematics is not used, neutron systematics should be turned on. 
Further decision on this case pending (based on Conversation with Trung Jan 17, 2020)
*/
namespace PlotUtils
{
  //=================================================================================
  // Response 
  //=================================================================================
  template<class T>
  class ResponseUniverse : public T
  {
    public:
      ResponseUniverse(typename T::config_t chw, double nsigma, std::string response_name, bool use_new_part_resp = false);

      double GetRecoilShift() const;
      virtual double GetCalRecoilEnergy() const; /* override */
      //virtual double GetRecoilEnergy() const; /* override */

      virtual std::string ShortName() const /*override*/;
      virtual std::string LatexName() const /*override*/;
      virtual bool IsVerticalOnly()  const  { return false; }/*override*/

      const std::string m_name;
      const std::string m_branch_name;

      //For the new particle response systematics
      double m_frac_unc; //Wish this could be a const, but oh well
      const bool m_use_new_part_resp;
      const std::string m_particle_passive_name;
      const std::string m_container_particle_passive_name;
  };

}

// For template classes, the header needs to know about the definitions
#include "ResponseSystematics.cxx"

#endif // RESPONSESYSTEMATICS_H
