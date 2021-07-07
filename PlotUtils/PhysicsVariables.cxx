#include <iostream>
#include <cmath>

#include "PhysicsVariables.h"

namespace PlotUtils{

//===============================================================================================
// Assume CCQE Kinematics, formula differs slightly for nu vs anti-nu, return units are GeV-type.
// Q^2 calc also computes E for a CCQE hypothesis and uses that value.
//===============================================================================================
double qSquaredCCQE( double lep_energy, double lep_p, double lep_theta, int charge )
{
  return PlotUtils::qSquaredCCQE( lep_energy, lep_p, lep_theta, charge, PlotUtils::getDefaultBindingE(getHelicity(charge)) );
}
//=============================================================================
double qSquaredCCQE( double lep_energy, double lep_p, double lep_theta, int charge, double bindingE, double lep_mass )
{
  // verbose() << "qSquaredCCQE w/ binding E = " << bindingE << " and charge = " << charge << std::endl;
  double nu_energy = PlotUtils::nuEnergyCCQE( lep_energy, lep_p, lep_theta, charge, bindingE, lep_mass );
  if( nu_energy < 0 ) return -1000000;
  double Qsquared  = 2.0 * nu_energy * (lep_energy - lep_p * std::cos(lep_theta)) - std::pow(lep_mass,2);
  // verbose() << " Qsquared = " << Qsquared << " (MeV/c)^2" << std::endl;
  return Qsquared;
}
//=============================================================================
double nuEnergyCCQE( double lep_energy, double lep_p, double lep_theta, int charge )
{
  return PlotUtils::nuEnergyCCQE( lep_energy, lep_p, lep_theta, charge, PlotUtils::getDefaultBindingE(PlotUtils::getHelicity(charge)) );
}
//=============================================================================
double nuEnergyCCQE( double lep_energy, double lep_p, double lep_theta, int charge, double bindingE, double lep_mass )
{
  // verbose() << "nuEnergyCCQE w/ binding E = " << bindingE << " and charge = " << charge << std::endl;
  // verbose() << "  lep_energy = " << lep_energy << "; lep_p = " << lep_p << "; lep_theta = " << lep_theta << std::endl;
  double nu_energy = -1000;
  if( charge > 0 ) {
    double nu_energy_num = std::pow(MinervaUnits::M_n,2) - std::pow(MinervaUnits::M_p - bindingE,2)
      - std::pow(lep_mass,2) + 2.0*(MinervaUnits::M_p - bindingE)*lep_energy;
    double nu_energy_den = 2.0*(MinervaUnits::M_p - bindingE - lep_energy + lep_p*std::cos(lep_theta));
    // verbose() << " Nu E numerator   = " << nu_energy_num << std::endl;
    // verbose() << " Nu E denominator = " << nu_energy_den << std::endl;
    if( nu_energy_den ) nu_energy = nu_energy_num / nu_energy_den;
  } else if ( charge < 0 ) {
    double nu_energy_num = std::pow(MinervaUnits::M_p,2) - std::pow(MinervaUnits::M_n - bindingE,2)
      - std::pow(lep_mass,2) + 2.0*(MinervaUnits::M_n - bindingE)*lep_energy;
    double nu_energy_den = 2.0*(MinervaUnits::M_n - bindingE - lep_energy + lep_p*std::cos(lep_theta));
    // verbose() << " Nu E numerator   = " << nu_energy_num << std::endl;
    // verbose() << " Nu E denominator = " << nu_energy_den << std::endl;
    if( nu_energy_den ) nu_energy = nu_energy_num / nu_energy_den;
  }
  // verbose() << " Nu Energy = " << nu_energy << " MeV" << std::endl;
  return nu_energy;
}

//=============================================================================
// Return Q-Squared using the general formula
//=============================================================================
double qSquared( double nu_energy, double lep_energy, double lep_theta, double lep_mass )
{
  double lep_p = lep_energy*lep_energy - lep_mass*lep_mass;
  if(lep_p >= 0){
    lep_p = sqrt(lep_p);
  }
  else{
    std::cout << "Got a negative lepton momentum in qSquared " << lep_energy << " " << lep_mass << std::endl;
    lep_p = 0.;
  }

  return  2.0 * nu_energy * (lep_energy - lep_p * std::cos(lep_theta)) - (lep_mass*lep_mass);
}

//=============================================================================
// Return the Mass of the struck nucleon based on observed muon charge (int)
//=============================================================================
double struckNucleonMass( int leptonCharge )
{
  if( -1 == leptonCharge ) { //struck nucleon is neutron
    return MinervaUnits::M_neutron;
  } else if( 1 == leptonCharge ) { //struck nucleon is proton
    return MinervaUnits::M_proton;
  } else if( 0 == leptonCharge ) {
    return (MinervaUnits::M_neutron + MinervaUnits::M_proton)/2.0;
  } else {
    std::cerr  << "Unrecognized charge for the Lepton! Setting mass of struck nucleon to -1! " << std::endl;
  }
  return -1.0;
}


////=============================================================================
//// Return the Mass of the struck nucleon based on observed muon charge (Particle class enum)
////=============================================================================
//double struckNucleonMass( Minerva::Particle::Q leptonCharge )
//{
//  return struckNucleonMass( getCharge( leptonCharge ) );
//}


////=============================================================================
//// Return the Mass of the struck nucleon based on Neutrino Helicity
////=============================================================================
//double struckNucleonMass( Minerva::NeutrinoInt::NeutrinoHelicity helicity )
//{
//  return struckNucleonMass( getCharge( helicity ) );
//}


//=============================================================================
// Get W given nucleon mass, recoil energy, and Q^2
//=============================================================================
double W( double mass, double recoilE, double qsquared )
{
  double w2 = PlotUtils::WSquared( mass, recoilE, qsquared );
  double w  = (w2<0) ? -1.0 : sqrt(w2);
  return w;
}


//=============================================================================
// Get W^2 given nucleon mass, recoil energy, and Q^2
//=============================================================================
double WSquared( double mass, double recoilE, double qsquared )
{
  return ( std::pow(mass, 2) + (2 * mass * recoilE) - qsquared );
}


//=============================================================================
// Get Bjorken-x given Q^2, nucleon mass, and recoil energy.
//=============================================================================
double xBjorken( double qsquared, double mass, double recoilE )
{
  if( 0 == mass    ) { std::cerr << "Zero Mass in xBjorken!"   << std::endl; return -1.0; }
  if( 0 == recoilE ) { std::cerr << "Zero Recoil in xBjorken!" << std::endl; return -1.0; }
  return 0.5 * qsquared / ( mass * recoilE );
}

double getDefaultBindingE( int hel )
{
  
	// this is a bit unkind to anyone in the future
	// who tries to use this method for a neutral current analysis...
	// but I regret nothing
	switch (hel)
	{
    //static const double m_CCQEBindingE_neutron = 34.;
    //static const double m_CCQEBindingE_proton  = 30.;
    case 1:
			return NSFDefaults::CCQEBindingE_neutron;

    case 2:
			return NSFDefaults::CCQEBindingE_proton;

		default:
			std::cerr << "Binding energy for unknown helicity requested. Assuming neutrino" << std::endl;
			return PlotUtils::getDefaultBindingE(hel);
	}
}

//=============================================================================
// Return a NeutrinoInt::NeutrinoHelicity based on lepton charge.
//=============================================================================
int getHelicity ( int leptonCharge )
{
  const int AntiNuHelicity = 2;
  const int NuHelicity = 1;
 if( leptonCharge > 0 ) return AntiNuHelicity;
 if( leptonCharge < 0 ) return NuHelicity;
  return  -999;
}
//int  getHelicity( Minerva::Particle::Q leptonCharge )
//{
//  return this->getHelicity( getCharge( leptonCharge ) );
//}
}
