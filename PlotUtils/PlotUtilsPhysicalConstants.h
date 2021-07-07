#ifndef MINERVAUNITS_H
#define MINERVAUNITS_H

// modified to run standalone

//#include "CLHEP/Units/PhysicalConstants.h"

//#####################################################################
//
// MinervaPhysicalConstants
// 
// One-stop-shopping for a variety of unit conventions and physical
// constants used in MINERvA software.  
//
// Please, ALWAYS use this header and NEVER define constants in local 
// code.
//
// Important note on changing units in CLHEP style:
//        1) MULTIPLY by the CURRENT unit
//   and  
//        2) DIVIDE by the DESIRED unit
//
// e.g. M_pi0 is given in MeV.  To change to GeV, one does:
//      double mass_in_GeV = M_pi0 * CLHEP::MeV / CLHEP::GeV;
//
//######################################################################

namespace MinervaUnits {
 
  // ================================================================== 
  // Units
  // ==================================================================
  
  static const double inch = 0.0254;
  static const double ft   = 12.0*inch;

  //static const double electron_charge_fc =   1e15;


  // ================================================================== 
  // MINERvA/NuMI Specific
  // ==================================================================

  static const double numi_beam_angle_rad = -0.05887;


  // ================================================================== 
  // Derived Units/Conversion Factors
  // ==================================================================


  // ================================================================== 
  // Mathematical Constants
  // ==================================================================


  // ================================================================== 
  // Physical Constants, speeds, times, distances etc.
  // ==================================================================

  //static const double c_light = CLHEP::c_light; // mm/ns


  // ================================================================== 
  // Particle Properties
  // ==================================================================
 
  //-- Bosons  -- hm these are in GeV 
  static const double M_W = 80.398  ;     // W Boson
  static const double M_Z = 91.1876 ;     // Z Boson
  
  //-- Baryons
  static const double M_p = 938.272013 ;    // Proton
  static const double M_proton = M_p;
  static const double M_n = 939.56536;   // Neutron
  static const double M_neutron = M_n;

  //-- Mesons
  static const double M_pi0  = 134.9766 ; // Neutral Pion
  static const double M_pion = 139.5701 ; // Charged Pion
  static const double M_k0   = 497.648  ; // Neutral Kaon
  static const double M_kaon = 493.677  ; // Charged Kaon

  //-- Leptons
  static const double M_e = 0.510998910;  // Electron
  static const double M_electron = M_e;
  static const double M_mu = 105.6583 ;   // Muon
  static const double M_muon = M_mu;
  static const double M_tau = 1776.99 ;   // Tau


  static const int AntiNuHelicity  = 2;
  static const int NuHelicity = 1;

  static const double m_CCQEBindingE_neutron = 34.;
  static const double m_CCQEBindingE_proton  = 30.;

}



#endif
