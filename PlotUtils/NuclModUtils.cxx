#include "NuclModUtils.h"
#include <TMath.h>

using namespace MAT;

//constants for me
namespace
{
  //handy constants
  const double Mn = (.9396 + .9383 ) / 2.;
  const double frac_c_in_ch = .88;

  //constants used for BY
/*  const double BY_A = .419;
  const double BY_B = .223;
*/
  const double BY_A = 0.538;
  const double BY_B = 0.305;
  //fit parameters for E139 fit
  const double E139_c[3] = {1.69029097E-2, 1.80889367E-2, 5.04268396E-3};
  const double E139_a[9] = {-6.98871401E-2, 2.18888887E0, -2.46673765E1, 1.45290967E2, -4.97236711E2, 1.01312929E3, -1.20839250E3, 7.75766802E2, -2.05872410E2 };

}

NuclModUtils&  NuclModUtils::Get()
{
  static NuclModUtils singleton;
  return singleton;
}


//===============
// Bodek-Yang scaling variables
//===============
double NuclModUtils::BY_ScalingVar_Xi_TM( double x, double q2 ) const
{
  const double nu = q2/(2.*x*Mn);
  const double xi_tm = q2 / ( Mn*nu*(1. + sqrt(1. + q2/TMath::Power(nu,2) ) ) );
  return xi_tm;
}

double NuclModUtils::BY_ScalingVar_Xi_w( double x, double q2 ) const
{
  double num = 2*x*(q2 + BY_B);
  double den = q2 * ( 1. + TMath::Sqrt( 1. + TMath::Power(2.*x*Mn,2)/q2 ) ) + 2*BY_A*x;
  if( den > 0. )
    return num / den;
  else
    return 0.;
}

//================
// Bodek-Yang 2003
//================
double NuclModUtils::BY_Free_to_D( double x ) const
{
  const double x2 = x * x;
  const double x3 = x2* x;
  const double x4 = x3* x;
  const double x5 = x4* x;
  return .985*(1. + .422*x - 2.745*x2 + 7.57*x3 - 10.335*x4 + 5.422*x5);
}

double NuclModUtils::BY_D_to_Fe( double x ) const
{
  return 1.096 - .364*x - .278*TMath::Exp(-21.94*x) + 2.772*TMath::Power(x,14.417);
}

//=======================
// Bodek-Yang 2013
//=======================
double NuclModUtils::BY13_D_to_Fe( double xi_tm) const
{
  //xi_tm is the scaling variable
  const double x = xi_tm;

  return 1.096 - 0.38*x - 0.3*TMath::Exp(-23*x) + 8.*TMath::Power(x,15);
}

double NuclModUtils::BY13_C_to_Fe( double xi_tm) const
{
  //xi_tm is the scaling variable
  const double x = xi_tm;
  const double x2 = x * x;
  const double x3 = x2* x;
  const double x4 = x3* x;
  const double x5 = x4* x;
  return .919 + 1.844*x - 12.73*x2 + 36.89*x3 - 46.77*x4 + 21.22*x5;
}

double NuclModUtils::BY13_Fe_to_Pb( double xi_tm) const
{
  //xi_tm is the scaling variable
  const double x = xi_tm;
  const double x2 = x * x;
  const double x3 = x2* x;
  const double x4 = x3* x;
  const double x5 = x4* x;
  const double x6 = x5* x;
  return .932 + 2.461*x - 24.23*x2 + 101.03*x3 - 203.47*x4 + 193.85*x5 - 69.82*x6;
}

//=====================
// E139's fit to CL data
//=====================
double NuclModUtils::E139( double x, double A ) const
{
  //caller doesn't care about goodfit.  it's a dummy
  bool goodfit = false;
  return E139( x, A, goodfit );
}

double NuclModUtils::E139( double x, double A, bool &goodfit ) const
{
  goodfit = true;

  // do nothing for D or free nuc
  if( A < 2.5 )
    return 1.;

  //check x limits
  if( x < .0085 )
  {
    x = .0085;
    goodfit = false;
  }
  if( .7 < x )
  {
    x = .7;
    goodfit = false;
  }

  double c = 0.;
  for( int i = 0; i < 3; ++i )
    c += E139_c[i]*TMath::Power( TMath::Log(x), i );
  c = TMath::Exp(c);

  double a = 0.;
  for( int i = 0; i < 9; ++i )
    a += E139_a[i]*TMath::Power( x, i );

  return c * TMath::Power( A, a );
}  


//===============================
// GENIE style event weights
//===============================
double NuclModUtils::GENIE_BY( double x, int A ) const
{
  //special case A=tracker
  if( 0 == A )
    return GENIE_BY(x,12)*frac_c_in_ch + (1.-frac_c_in_ch);

  double f = 1.;

  //scale to deuterium for all nuclei
  if( 2 <= A )
    f *= BY_Free_to_D(x);

  //scale to isoscalar iron for all nuclei larger than deuterium
  if( 2 < A )
    f *= BY_D_to_Fe(x);

  return f;
}

double NuclModUtils::GENIE_BY_bug( double x, double q2, int A ) const
{
  //special case A=tracker
  if( 0 == A )
    return GENIE_BY_bug(x, q2, 12)*frac_c_in_ch + (1.-frac_c_in_ch);

  //get the scaling variable genie uses
  double xi_w = BY_ScalingVar_Xi_w( x, q2 );

  //calculate BY modification with xi_w
  return GENIE_BY( xi_w, A );
  
}


double NuclModUtils::GENIE_BY13( double x, double q2, int A ) const
{
  //special case A=tracker
  if( 0 == A )
    return GENIE_BY13(x, q2, 12)*frac_c_in_ch + (1.-frac_c_in_ch);

  double f = 1.;

  //scale to deuterium for all nuclei
  if( 2 <= A )
    f *= BY_Free_to_D(x);

  if( 2 < A )
  {
    //get BY13 scaling variable
    const double xi_tm = BY_ScalingVar_Xi_TM( x, q2 );
    const double d_to_fe = BY13_D_to_Fe( xi_tm );

    if( 12 == A )
    {
      //we have a specific prediction for carbon
      const double c_to_fe = BY13_C_to_Fe( xi_tm );
      f *= ( d_to_fe / c_to_fe );
    }
    else if( 56 == A )
    {
      //we have a specific prediction for iron
      f *= d_to_fe;
    }
    else if( 207 == A )
    {
      //we have a specific prediction for lead
      const double fe_to_pb = BY13_Fe_to_Pb( xi_tm );
      f *= ( d_to_fe * fe_to_pb );
    }
    else
    {
      //use iron when we have no specific prediction for a nucleus
      //perhaps use E193's A^alpha to scale between the nearest nuclei?
      f *= d_to_fe;
    }

  }//end if nucleus is greater than deuterium

  return f;
}

double NuclModUtils::GENIE_E139( double x, int A ) const
{
  //special case A=tracker
  if( 0 == A )
    return GENIE_E139(x,12)*frac_c_in_ch + (1.-frac_c_in_ch);

  double f = 1.;

  //scale to deuterium for all nuclei
  if( 2 <= A )
    f *= BY_Free_to_D(x);

  //scale using E139 fit for all nuclei (does nothing to D)
  if( 2 < A )
    f *= E139( x, A );

  return f;
}
