#ifndef PARTICLERESPONSECONSTANTS_H
#define PARTICLERESPONSECONSTANTS_H 1
#include <iostream>

//Placing the standard particle response fractional uncertainties here
// $ANAUTILSROOT/src/ParticleResponseUtils.cpp has the standard responses used by the Anatool
// https://cdcvs.fnal.gov/redmine/projects/minerva/wiki/NX_systematics has references on where these comes from
namespace PlotUtils
{
  namespace PartRespDefaults
  {
    const double protonResponse          = 0.035 ; // %
    const double protonResponseLowEkin   = 0.04  ; // %  
    const double protonResponseMidEkin   = 0.035 ; // %
    const double protonResponseHighEkin  = 0.03  ; // %
    const double neutronResponseLowEkin  = 0.25  ; // %
    const double neutronResponseMidEkin  = 0.10  ; // %
    const double neutronResponseHighEkin = 0.20  ; // %
    const double mesonResponse           = 0.05  ; // %
    const double electromagneticResponse = 0.03  ; // %
    const double otherParticleResponse   = 0.20  ; // %
    const double muonResponse            = 0.024 ; // %
    const double xtalkResponse           = 0.20  ; // %

    double GetDefaultPartRespFracUnc( std::string particle );
  }
}

#endif
