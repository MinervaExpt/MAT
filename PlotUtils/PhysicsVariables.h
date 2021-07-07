

#ifndef PHYSICSVARIABLES_H
#define PHYSICSVARIABLES_H 1

//#ifndef PLOTUTILS_STANDALONE
//#include "Kernel/MinervaPhysicalConstants.h"
//#else
#include "PlotUtilsPhysicalConstants.h"
//#endif

#include "NSFDefaults.h"  // Get binding energies here



namespace PlotUtils{

    // Assume CCQE Kinematics, formula differs slightly for nu vs anti-nu, all energies are MeV.
    // Q^2 calc also computes E for a CCQE hypothesis and uses that value.
    double nuEnergyCCQE( double lep_energy, double lep_p, double lep_theta, int charge );
    double qSquaredCCQE( double lep_energy, double lep_p, double lep_theta, int charge );
    double nuEnergyCCQE( double lep_energy, double lep_p, double lep_theta, int charge, double bindingE, double lep_mass = MinervaUnits::M_mu );
    double qSquaredCCQE( double lep_energy, double lep_p, double lep_theta, int charge, double bindingE, double lep_mass = MinervaUnits::M_mu );

    /// @f$Q^{2} = E_{\nu}E_{\mu}sin(\frac{\theta_{\mu}}{2})^{2} @f$
    double qSquared( double nu_energy, double lep_energy, double lep_theta,
                     double lep_mass = MinervaUnits::M_mu );
  /// Return the Mass of the struck nucleon based on observed lepton charge (int)
    double struckNucleonMass( int leptonCharge );

 
    /// Return the Mass of the struck nucleon based on Neutrino Helicity
    double struckNucleonMass( int  helicity );

    /// Get W given nucleon mass, recoil energy, and Q^2
    double W( double mass, double recoilE, double qsquared );
    /// Get W^2 given nucleon mass, recoil energy, and Q^2
    double WSquared( double mass, double recoilE, double qsquared );
    /// Get Bjorken-x given Q^2, nucleon mass, and recoil energy.
    double xBjorken( double qsquared, double mass, double recoilE );

    /// Return a NeutrinoInt::NeutrinoHelicity based on lepton charge (int)
    int getHelicity( int charge );
   
    
    /// Return the default binding energy in use.  Depends on neutrino helicity...
    double getDefaultBindingE( int  hel );

    //const double m_CCQEBindingE_proton = 34.;
    //const double m_CCQEBindingE_neutron = 30.;

    /// Return an integer lepton charge based on NeutrinoInt::NeutrinoHelicity.
    inline int getCharge( int helicity )
    {

      if( helicity == 2   ) return  1;
      if( helicity == 1) return -1;
     
    return 0;
    }

   


};


#endif
