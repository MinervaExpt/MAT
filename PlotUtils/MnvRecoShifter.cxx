#ifndef MNV_MnvRecoShifter_cxx
#define MNV_MnvRecoShifter_cxx 1

#include "MnvRecoShifter.h"

#include "TError.h"
#include "HistogramUtils.h" //for NotPhysicsShift

// Heidi says include from here instead of the big build

/*#ifndef PLOTUTILS_STANDALONE
#include "Kernel/MinervaPhysicalConstants.h"
#else
#include "PlotUtilsPhysicalConstants.h"
#endif
*/

#include <math.h>
#include <iostream>

#define DEBUGLEVEL   0
#define DEBUGQ2SHIFT 5  // If DEBUGLEVEL > DEBUGQ2SHIFT, print

using namespace PlotUtils;


//=================================
// BEGIN - Constructors/Destructors
//=================================

// Default constructor - default to neutrino.
MnvRecoShifter::MnvRecoShifter() :
  m_isAntiNu( false ),
  m_nUniverses( 1000 )
{ 
  SetDefaultEnergyLimits();
  GetRandomShiftVectors( m_nUniverses );
}

// Constructor with nu/antinu declaration.
MnvRecoShifter::MnvRecoShifter( bool isAnti) :
  m_isAntiNu( isAnti ),
  m_nUniverses( 1000 )
{
  SetDefaultEnergyLimits();
  GetRandomShiftVectors( m_nUniverses );
}

// Constructor with nu/antinu and nUniverses declaration
MnvRecoShifter::MnvRecoShifter( bool isAnti, int nUniverses) :
  m_isAntiNu( isAnti ),
  m_nUniverses( nUniverses )
{
  SetDefaultEnergyLimits();
  GetRandomShiftVectors( m_nUniverses );
}

// Default destructor
MnvRecoShifter::~MnvRecoShifter()
{}

//=================================
// END - Constructors/Destructors
//=================================


//=================================
// BEGIN - Calculators
//=================================

//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::Calc_q2( double& q2, const double muonE, const double hadronicE, const double muonTheta ) const
{
  //! Sanity checks on the variables
  if( muonE < 0. || hadronicE < 0. )
    return false;

  //! q2 = 4. * muonE * neutrinoE * pow( sin( muonTheta / 2. ), 2. );
  q2 = 4. * muonE * (muonE+hadronicE) * pow( sin( muonTheta / 2. ), 2. );

  return true;
}

//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::Calc_w( double& w, const double muonE, const double hadronicE, const double muonTheta ) const
{
  //! Sanity checks on the variables
  if( muonE < 0. || hadronicE < 0. )
    return false;

  //! What about mass?  Split difference between proton and neutron for now.
  const double mass = ( MinervaUnits::M_p + MinervaUnits::M_n ) / 2.;

  //! First calculate Q^2.
  double q2 = 0.0;
  bool didCalculateQ2 = Calc_q2( q2, muonE, hadronicE, muonTheta );

  //! Next, calculate W^2.
  double w2 = pow( mass, 2 ) + (2 * mass * hadronicE) - q2;

  //! If W^2 is positive, calculate W and return true (anded with didCalculateQ2).
  if( w2 > 0 ) { 
    w = sqrt( w2 );
    return didCalculateQ2;
  }

  //! Else, if W^2 is not positive, set W to 0 and return false to signal an unphysical W.
  w = 0.0;
  return false;
}

//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::Calc_x( double& x, const double muonE, const double hadronicE, const double muonTheta ) const
{
  //! Sanity checks on the variables
  if( muonE < 0. ||  hadronicE < 0. )
    return false;

  //! What about mass?  Split difference between proton and neutron for now.
  const double mass = ( MinervaUnits::M_p + MinervaUnits::M_n ) / 2.;

  //! get q2 from the calculator
  double q2 = 0.;
  if( ! Calc_q2( q2, muonE, hadronicE, muonTheta ) )
    return false;

  //! x = 0.5 * q2 / ( mass * hadronicE );
  x = 0.5 * q2 / ( mass * hadronicE );

  return true;

}

//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::Calc_y( double& y, const double muonE, const double hadronicE ) const
{
  //! Sanity checks on the variables
  if( muonE < 0. || hadronicE < 0. )
    return false;

  //! neutrino E = muonE + hadronicE
  const double E = muonE + hadronicE;

  //! If E is positive, then y=hadronicE/E.
  //! If E is not positive ( E=0 ), return false because we failed
  if( E > 0. )
    y = hadronicE / E;
  else
    return false;

  return true;
}
//=================================
// END - Calculators
//=================================
//==========================================
// BEGIN - CCQE Calculators and Scale Shifts
//==========================================

bool MnvRecoShifter::Calc_ccqe_nuE( double& nuE, const double lepP, const double lepTheta, const double bindingE, double lepMass ) const
{

  const double lepE = sqrt( pow( lepP, 2 ) + pow( lepMass, 2 ) );

  if( m_isAntiNu ) 
  {
    double nu_energy_num = pow( MinervaUnits::M_n, 2 ) - pow( MinervaUnits::M_p - bindingE, 2 ) - pow( lepMass, 2 ) + 2.0 * ( MinervaUnits::M_p - bindingE ) * lepE;
    double nu_energy_den = 2.0 * ( MinervaUnits::M_p - bindingE - lepE + lepP * cos( lepTheta ) );
    if ( nu_energy_den ) 
      nuE = nu_energy_num / nu_energy_den;
    else 
      return false;
  }
  else
  {
    double nu_energy_num = pow( MinervaUnits::M_p, 2 ) - pow( MinervaUnits::M_n - bindingE, 2 ) - pow( lepMass, 2 ) + 2.0 * ( MinervaUnits::M_n - bindingE ) * lepE;
    double nu_energy_den = 2.0 * ( MinervaUnits::M_n - bindingE - lepE + lepP * cos( lepTheta ) );
    if ( nu_energy_den ) 
      nuE = nu_energy_num / nu_energy_den;
    else 
      return false;
  }

  return true;
}

bool MnvRecoShifter::Calc_ccqe_q2( double& q2, const double lepP, const double lepTheta, const double bindingE, double lepMass ) const
{ 

  const double lepE = sqrt( pow( lepP, 2 ) + pow( lepMass, 2 ) );

  double nuE = 0.;
  if ( ! Calc_ccqe_nuE( nuE, lepP, lepTheta, bindingE ) )
    return false;
  if( nuE < 0 ) 
    return false;

  q2 = 2.0 * nuE * ( lepE - lepP * cos( lepTheta ) ) - pow( lepMass, 2 );

  return true;
}

bool MnvRecoShifter::Calc_tki_vars( dVec& vars, const TVector3 &nu3P,const TLorentzVector &lepton, const TLorentzVector &nucleon, const double NucleusMass, const double ISNucleonMass, const double bindingE) const
{
  XYZVector Nu(nu3P.X(), nu3P.Y(), nu3P.Z() );
  XYZTVector Lepton(lepton.X(), lepton.Y(), lepton.Z(), lepton.E() );
  XYZTVector Nucleon(nucleon.X(), nucleon.Y(), nucleon.Z(), nucleon.E() );

  return Calc_tki_vars( vars, Nu, Lepton, Nucleon, NucleusMass, ISNucleonMass, bindingE );
}

bool MnvRecoShifter::Calc_tki_vars( dVec& vars, const XYZVector &nu3P,const XYZTVector &lepton, const XYZTVector &nucleon, const double NucleusMass, const double ISNucleonMass, const double bindingE) const
{
  double RemnantNucleusMass = NucleusMass - ISNucleonMass + bindingE;
  XYZVector lepton3P = lepton.Vect();
  XYZVector nucleon3P = nucleon.Vect();
  double Elepton = lepton.E(), Enucleon = nucleon.E(), Enu = -999.;

  XYZVector coordz = nu3P.Unit();
  XYZVector coordx = lepton3P.Cross(coordz).Unit();
  XYZVector coordy = coordz.Cross(coordx).Unit();

  //prime coordinate of coord(x,y,z)
  XYZVector pNucleon3P( nucleon3P.Dot(coordx),nucleon3P.Dot(coordy),nucleon3P.Dot(coordz));
  XYZVector pLepton3P( lepton3P.Dot(coordx),lepton3P.Dot(coordy),lepton3P.Dot(coordz));

  XYZVector pNucleonPT( pNucleon3P.X(), pNucleon3P.Y(), 0 );
  XYZVector pLeptonPT( pLepton3P.X(), pLepton3P.Y(), 0 );

  XYZVector dpTvec = pNucleonPT+pLeptonPT;

  double R = NucleusMass + pLepton3P.Z() + pNucleon3P.Z() - Elepton - Enucleon;
  double dPL = 0.5*( R- ( RemnantNucleusMass*RemnantNucleusMass + dpTvec.Perp2() )/R  );
  double Pn = TMath::Power(dpTvec.Perp2() + dPL*dPL, 0.5 );
  Enu = pLepton3P.Z() + pNucleon3P.Z()-dPL;

  double dalphaT = TMath::ACos( -pLeptonPT.Unit().Dot( dpTvec.Unit() )) * 180/TMath::Pi();
  double dpT = dpTvec.R();
  double dpTx = dpTvec.X();
  double dpTy = dpTvec.Y();
  double dphiT = TMath::ACos( -pLeptonPT.Unit().Dot( pNucleonPT.Unit() ) ) *180/TMath::Pi();
  double sign = (dpTx>0)? +1 : -1;

  vars= dVec({ Enu, Pn, dpT,dPL, dpTx, dpTy, dalphaT, dphiT, sign });  
  return true;
}

bool MnvRecoShifter::Calc_expected_nucleon( TLorentzVector& expected, TVector3 &nu3P, TLorentzVector &lepton, double ISMass, double FSMass, double BindingE )
{

  //rotate to nu3P = (0,0,1)
  TVector3 nuAxis(0,0,1);
  TRotation R1;//( ACos( nuAxis.Dot( nu3P.Unit() )), nu3P.Cross(nuAxis) );
  double rotationAngle = ACos(nuAxis.Dot( nu3P.Unit() ));
  TVector3 rotationAxis = nu3P.Cross(nuAxis) ;
  R1.AngleAxis( rotationAngle, rotationAxis );
  //AxisAngle R1( nu3P.Cross( nuAxis) , ACos( nuAxis.Dot( nu3P.Unit() )) );
  //LorentzVector lep = (LorentzVector) lepton;
  
  TLorentzVector leptonR1 = (lepton);
  leptonR1.Transform(R1);
  //TLorentzVector leptonR1(lepton3P.X(), lepton3P.Y(), lepton3P.Z(), lepton.E());
  
  //(v00v) + (M000) = (Em pxpy pz ) + (EM -px-py pmz)
  // v-EM = Em-M, v- pmz = pz
  //EM - pmz = pz - Em + M0
  // A = pz - Em + M0
  // M1^2 + px^2 + py^2 + pmz^2 = A^2 + pmz^2 + 2A*pmz, 
  // A = pz - Em + Mi, Mi = M0 - Eb
  // B^2 = M1^2 + px^2 + py^2
  // B^2 = A^2 + 2A*pmz ===> pmz = (B^2 - A^2 )/(2A)
  double B2 = FSMass*FSMass + leptonR1.X()*leptonR1.X() + leptonR1.Y()*leptonR1.Y();
  double A = leptonR1.Z() - leptonR1.E() + (ISMass-BindingE);
  double pMz = (B2/A-A)/2.;
  double pMx = -leptonR1.X(), pMy = -leptonR1.Y();
  TLorentzVector nucleonR1( pMx, pMy, pMz, sqrt( pMx*pMx + pMy*pMy + pMz*pMz + FSMass*FSMass ) );
  
  //rotate back
  //expected = TVector3(R1.Inverse()*nucleonR1) ;  
  expected = TLorentzVector(nucleonR1.Transform(R1.Inverse())) ;  

  return true;


}


bool MnvRecoShifter::Calc_expected_nucleon( XYZTVector& expected, XYZVector &nu3P, XYZTVector &lepton, double ISMass, double FSMass, double BindingE )
{
  TLorentzVector expected1;
  TLorentzVector lepton1( lepton.X(), lepton.Y(), lepton.Z(), lepton.E() );
  TVector3 nu3P1( nu3P.X(), nu3P.Y(), nu3P.Z());
  bool ret = Calc_expected_nucleon( expected1, nu3P1, lepton1,  ISMass,  FSMass,  BindingE ) ;
  expected.SetXYZT( expected1.X(), expected1.Y(), expected1.Z(), expected1.E() );

  return ret;
}



bool MnvRecoShifter::ComputeNeutronAngularVars( dVec & res, TVector3 &nu, TLorentzVector& lepton, TVector3 &targetVec,double ISMass, double FSMass, double BindingE  )
{ 
  
  TVector3 tgt = targetVec;
  TLorentzVector expVec;
  bool ret = Calc_expected_nucleon( expVec, nu, lepton, ISMass, FSMass, BindingE );

  TVector3 tgt_expected = tgt.Unit()*expVec.Vect().Mag(); 
  // this is the inferred neutron momentum, it has same magnitude as the expected momentum, but with the measured direction as the targetVec
  TVector3 coordz = expVec.Vect().Unit();
  TVector3 coordx = coordz.Cross( nu ).Unit();
  TVector3 coordy = coordz.Cross( coordx );

  double dx = coordx.Dot( tgt ), dy = coordy.Dot(tgt), dz = coordz.Dot(tgt);
  double dxi = coordx.Dot( tgt_expected ), dyi = coordy.Dot(tgt_expected);//, dzi = coordz.Dot(tgt_expected);
  double dPPerp = dx, dPReact = dy;
  double dPPerpI = dxi, dPReactI = dyi;
  
  double ReactPlaneAngle = ATan2( dy, dz );
  //double PerpPlaneAngle = ATan2( dy, Power( (dz*dz+dx*dx), 0.5 ) );
  double PerpPlaneAngle = ATan2( dx, dz );
  double Theta = ACos( expVec.Vect().Unit().Dot( tgt ) );
  res = dVec({PerpPlaneAngle, ReactPlaneAngle, Theta, dPPerp, dPReact, dPPerpI, dPReactI});
  //ret.push_back(PerpPlaneAngle);
  //ret.push_back(ReactPlaneAngle);
  return ret;
}

bool MnvRecoShifter::ComputeNeutronAngularVars( dVec & res, XYZVector &nu, XYZTVector& lepton, XYZVector &targetVec,double ISMass, double FSMass, double BindingE  )
{ 
  TVector3 Nu( nu.X(), nu.Y(), nu.Z() );
  TVector3 TargetVec( targetVec.X(), targetVec.Y(), targetVec.Z() );
  TLorentzVector Lep( lepton.X(), lepton.Y(), lepton.Z(), lepton.E() );
  return ComputeNeutronAngularVars( res, Nu, Lep, TargetVec, ISMass, FSMass, BindingE );
 }


bool MnvRecoShifter::GetShifts_muon( TVector3& shift, const double error, const int nUniverse, const double nSigma ) const
{
  //valid universe numbers: [ 1, m_nUniverses ]
  if ( nUniverse < 1 || m_nUniverses < nUniverse )
  {
    Error("MnvRecoShifter::GetShifts_beam_tki", "You need to specify a valid Universe number.");
    return false;
  }
  double a = error*nSigma;
  double sx = m_rShiftMuonDirX[nUniverse-1];
  double sy = m_rShiftMuonDirY[nUniverse-1];
  double sz = m_rShiftMuonDirZ[nUniverse-1];
  shift.SetXYZ( a*sx, a*sy, a*sz );

  return true;
}


bool MnvRecoShifter::GetShifts_hadron( TVector3& shift, const double error, const int nUniverse, const double nSigma ) const
{
  //valid universe numbers: [ 1, m_nUniverses ]
  if ( nUniverse < 1 || m_nUniverses < nUniverse )
  {
    Error("MnvRecoShifter::GetShifts_beam_tki", "You need to specify a valid Universe number.");
    return false;
  }
  double a = error*nSigma;
  double sx = m_rShiftHadronDirX[nUniverse-1];
  double sy = m_rShiftHadronDirY[nUniverse-1];
  double sz = m_rShiftHadronDirZ[nUniverse-1];
  shift.SetXYZ( a*sx, a*sy, a*sz );

  return true;
}


bool MnvRecoShifter::GetShifts_beamTheta( TVector3& shifted_beam, const double error, const std::string axis, const int nUniverse, const double nSigma ) const
{
  //valid universe numbers: [ 1, m_nUniverses ]
  
  if ( nUniverse < 1 || m_nUniverses < nUniverse )
  {
    Error("MnvRecoShifter::GetShifts_beamTheta", "You need to specify a valid Universe number.");
    return false;
  }
  if( axis != "x" || axis != "y" )
  {
    Error("MnvRecoShifter::GetShifts_beamTheta", "You need to specify a valid axis (x or y)");
    return false;
  }
  //double sx = m_rShiftBeamDirX[nUniverse-1];
  //double sy = m_rShiftBeamDirY[nUniverse-1];
  //double sz = m_rShiftBeamDirZ[nUniverse-1];

  //TVector3 beam(0,sin(numi_beam_angle_rad),cos(numi_beam_angle_rad));
  //1. rotate beam to beam' by theta x or y  in the beam frame
  //2. rotate beam' into the detector frame
  TVector3 beam(0,0,1);
  if (axis == "x")
  {
    double shift_theta = error*nSigma*m_rShiftBeamThetaX[nUniverse-1];
    beam.RotateX(shift_theta);
    beam.RotateX(-MinervaUnits::numi_beam_angle_rad);
  }
  if (axis == "y")
  {
    double shift_theta = error*nSigma*m_rShiftBeamThetaY[nUniverse-1];
    beam.RotateY(shift_theta);
    beam.RotateX(-MinervaUnits::numi_beam_angle_rad);
  }
  shifted_beam = beam;
  return true;
}

//-----------------------------------
// Muon momentum shifts
//-----------------------------------
bool MnvRecoShifter::GetShift_ccqe_muonP( double& shift, const double muonPErr,
    const int nUniverse /* = -1 */, const double nSigma /*= 1.*/ ) const
{

  //valid universe numbers: [ 1, m_nUniverses ]
  if ( nUniverse < 1 || m_nUniverses < nUniverse )
  {
    Error("MnvRecoShifter::GetShift_ccqe_muonP", "You need to specify a valid Universe number.");
    return false;
  }

  shift = muonPErr * m_rShiftMuon[nUniverse-1] * nSigma ;

  return true;
}

bool MnvRecoShifter::GetShift_ccqe_muonE_muonP( double& shift, 
    const double muonP, const double muonPErr,
    const int nUniverse /* = -1 */, const double nSigma /*= 1.*/ ) const
{

  bool ok = true;
  double muonP_shift = 0., muonE_prime = 0., muonE = 0.;

  ok = ok && GetShift_ccqe_muonP( muonP_shift, muonPErr, nUniverse, nSigma); 
  muonE_prime = sqrt( pow( muonP + muonP_shift, 2 ) + pow( MinervaUnits::M_mu, 2 ) );
  muonE = sqrt( pow( muonP, 2 ) + pow( MinervaUnits::M_mu, 2 ) );

  shift = muonE_prime - muonE;

  return ok;
}

bool MnvRecoShifter::GetShift_ccqe_muonT_muonP( double& shift,
    const double muonP, const double muonPErr,
    const int nUniverse /* = -1 */, const double nSigma /*=1.*/) const
{

  bool ok = true;
  double muonP_shift = 0., muonT_prime = 0., muonT = 0.;

  ok = ok && GetShift_ccqe_muonP( muonP_shift, muonPErr, nUniverse, nSigma); 
  muonT_prime = sqrt( pow( muonP + muonP_shift, 2 ) + pow( MinervaUnits::M_mu, 2 ) ) - MinervaUnits::M_mu;
  muonT = sqrt( pow( muonP, 2 ) + pow( MinervaUnits::M_mu, 2 ) ) - MinervaUnits::M_mu;

  shift = muonT_prime - muonT;

  return ok;
}

bool MnvRecoShifter::GetShift_ccqe_nuE_muonP( double &shift, const double muonP, const double muonTheta, const double bindingE, const double muonPErr, const int nUniverse /* = -1 */, const double nSigma /*= 1.*/) const
{

  bool ok = true;
  double nuE_prime = 0., nuE = 0., muonP_shift = 0.;

  ok = ok && GetShift_ccqe_muonP( muonP_shift, muonPErr, nUniverse, nSigma);
  ok = ok && Calc_ccqe_nuE( nuE_prime, muonP + muonP_shift, muonTheta, bindingE );
  ok = ok && Calc_ccqe_nuE( nuE, muonP, muonTheta, bindingE );

  shift = nuE_prime - nuE;

  return ok;
}

bool MnvRecoShifter::GetShift_ccqe_q2_muonP(double &shift, const double muonP, const double muonTheta, const double bindingE, const double muonPErr, const int nUniverse /* = -1 */, const double nSigma /*= 1.*/) const
{

  bool ok = true;
  double q2_prime = 0., q2 = 0., muonP_shift = 0.;

  ok = ok && GetShift_ccqe_muonP( muonP_shift, muonPErr, nUniverse, nSigma); 
  ok = ok && Calc_ccqe_q2( q2_prime, muonP + muonP_shift, muonTheta, bindingE);
  ok = ok && Calc_ccqe_q2( q2, muonP, muonTheta, bindingE);

  shift = q2_prime - q2 ;

  return ok;
}

bool MnvRecoShifter::GetShift_ccqe_muonPT_muonP(double &shift, const double muonP, const double muonTheta, const double muonPErr, const int nUniverse, const double nSigma) const 
{ 

  bool ok = true; 
  double muonP_shift = 0.0, muonPT_prime = 0.0, muonPT = 0.0; 

  ok = ok && GetShift_ccqe_muonP(muonP_shift, muonPErr, nUniverse, nSigma); 
  muonPT_prime = (muonP + muonP_shift)*sin(muonTheta); 
  muonPT       = muonP*sin(muonTheta); 

  shift = muonPT_prime - muonPT; 

  return ok; 
} 

bool MnvRecoShifter::GetShift_ccqe_muonPZ_muonP(double &shift, const double muonP, const double muonTheta, const double muonPErr, const int nUniverse, const double nSigma) const 
{ 

  bool ok = true; 
  double muonP_shift = 0.0, muonPZ_prime = 0.0, muonPZ = 0.0; 

  ok = ok && GetShift_ccqe_muonP(muonP_shift, muonPErr, nUniverse, nSigma); 
  muonPZ_prime = (muonP + muonP_shift)*cos(muonTheta); 
  muonPZ       = muonP*cos(muonTheta); 

  shift = muonPZ_prime - muonPZ; 

  return ok; 
} 

//-----------------------------------
// Muon angle shifts
//-----------------------------------
bool MnvRecoShifter::GetShift_ccqe_muonTheta( double& shift, const double muonThetaErr,
    const int nUniverse /* = -1 */, const double nSigma /*= 1.*/ ) const
{

  //valid universe numbers: [ 1, m_nUniverses ]
  if ( nUniverse < 1 || m_nUniverses < nUniverse )
  {
    Error("MnvRecoShifter::GetShift_ccqe_muonTheta", "You need to specify a valid Universe number.");
    return false;
  }

  shift = muonThetaErr * m_rShiftMuon[nUniverse-1] * nSigma ;

  return true;
}

bool MnvRecoShifter::GetShift_ccqe_muonThetaXY( double& shift, const double muonTheta, const double muonThetaX, const double muonThetaX_Err, const double muonThetaY, const double muonThetaY_Err, const int nUniverse /* = -1 */, const double nSigma /*= 1.*/ ) const
{

  bool ok = true; 
  double muonTheta_prime = 0.;
  double muonThetaX_shift = 0.0, muonThetaY_shift = 0.0; 
  double cos_newThetaX = 0.0, cos_newThetaY = 0.0, cos_muonTheta_shifted = 0.0;  

  ok = ok && GetShift_ccqe_muonTheta(muonThetaX_shift, muonThetaX_Err, nUniverse, nSigma); 

  ok = ok && GetShift_ccqe_muonTheta(muonThetaY_shift, muonThetaY_Err, nUniverse, nSigma); 
                                                                                    
  cos_newThetaX = cos(muonThetaX + muonThetaX_shift); 
  cos_newThetaY = cos(muonThetaY + muonThetaY_shift); 

  //This is merely: sec^2(theta_X) + sec^2(theta_Y) = sec^2(theta) +1 
  cos_muonTheta_shifted = 1/(sqrt(1.0/pow(cos_newThetaX, 2) + 1.0/pow(cos_newThetaY, 2) -1) ); 
  muonTheta_prime = acos(cos_muonTheta_shifted); 

  shift = muonTheta_prime - muonTheta;

  return ok;
}

bool MnvRecoShifter::GetShift_ccqe_nuE_muonTheta( double &shift, const double muonP, const double muonTheta, const double bindingE, const double muonThetaErr, const int nUniverse /* = -1 */, const double nSigma /*= 1.*/) const
{

  bool ok = true;
  double nuE_prime = 0., nuE = 0., muonTheta_shift = 0.;

  ok = ok && GetShift_ccqe_muonP( muonTheta_shift, muonThetaErr, nUniverse, nSigma);
  ok = ok && Calc_ccqe_nuE( nuE_prime, muonP, muonTheta + muonTheta_shift, bindingE );
  ok = ok && Calc_ccqe_nuE( nuE, muonP, muonTheta, bindingE );

  shift = nuE_prime - nuE;

  return ok;
}

bool MnvRecoShifter::GetShift_ccqe_q2_muonTheta(double &shift, const double muonP, const double muonTheta, const double bindingE, const double muonThetaErr, const int nUniverse /* = -1 */, const double nSigma /*= 1.*/) const
{

  bool ok = true;
  double q2_prime = 0., q2 = 0., muonTheta_shift = 0.;

  ok = ok && GetShift_ccqe_muonTheta( muonTheta_shift, muonThetaErr, nUniverse, nSigma); 
  ok = ok && Calc_ccqe_q2( q2_prime, muonP, muonTheta + muonTheta_shift, bindingE);
  ok = ok && Calc_ccqe_q2( q2, muonP, muonTheta, bindingE);

  shift = q2_prime - q2 ;

  return ok;
}

bool MnvRecoShifter::GetShift_ccqe_q2_muonThetaXY( double &shift, const double muonP, const double muonTheta, const double bindingE, const double muonThetaX, const double muonThetaX_Err, const double muonThetaY, const double muonThetaY_Err, const int nUniverse, const double nSigma ) const 
{

  bool ok = true; 
  double q2_prime = 0., q2 = 0.;
  double muonThetaX_shift = 0.0, muonThetaY_shift = 0.0; 
  double cos_newThetaX = 0.0, cos_newThetaY = 0.0, cos_muonTheta_shifted = 0.0;  

  ok = ok && GetShift_ccqe_muonTheta(muonThetaX_shift, muonThetaX_Err, nUniverse, nSigma);

  ok = ok && GetShift_ccqe_muonTheta(muonThetaY_shift, muonThetaY_Err, nUniverse, nSigma); 
                                                                                    
  cos_newThetaX = cos(muonThetaX + muonThetaX_shift); 
  cos_newThetaY = cos(muonThetaY + muonThetaY_shift); 

  //This is merely: sec^2(theta_X) + sec^2(theta_Y) = sec^2(theta) +1 
  cos_muonTheta_shifted = 1/(sqrt(1.0/pow(cos_newThetaX, 2) + 1.0/pow(cos_newThetaY, 2) -1) ); 
  double muonTheta_shifted = acos(cos_muonTheta_shifted); 

  ok = ok && Calc_ccqe_q2( q2_prime, muonP, muonTheta_shifted, bindingE);
  ok = ok && Calc_ccqe_q2( q2, muonP, muonTheta, bindingE);

  shift = q2_prime - q2 ;

  return ok;
}

bool MnvRecoShifter::GetShift_ccqe_muonPT_muonTheta(double &shift, const double muonP, const double muonTheta, const double muonThetaErr, const int nUniverse, const double nSigma) const 
{ 

  bool ok = true; 
  double muonTheta_shift = 0.0, muonPT_prime = 0.0, muonPT = 0.0; 

  ok = ok && GetShift_ccqe_muonTheta(muonTheta_shift, muonThetaErr, nUniverse, nSigma); 
  muonPT_prime = muonP*sin(muonTheta+muonTheta_shift); 
  muonPT       = muonP*sin(muonTheta); 

  shift = muonPT_prime - muonPT; 

  return ok; 
} 

bool MnvRecoShifter::GetShift_ccqe_muonPT_muonThetaXY(double &shift, const double muonP, const double muonTheta, const double muonThetaX, const double muonThetaX_Err, const double muonThetaY, const double muonThetaY_Err, const int nUniverse, const double nSigma) const 
{

  bool ok = true; 
  double muonThetaX_shift = 0.0, muonThetaY_shift = 0.0, muonPT_prime = 0.0, muonPT = 0.0; 
  double cos_newThetaX = 0.0, cos_newThetaY = 0.0, cos_muonTheta_shifted = 0.0;  

  ok = ok && GetShift_ccqe_muonTheta(muonThetaX_shift, muonThetaX_Err, nUniverse, nSigma); 

  ok = ok && GetShift_ccqe_muonTheta(muonThetaY_shift, muonThetaY_Err, nUniverse, nSigma); 

  cos_newThetaX = cos(muonThetaX + muonThetaX_shift); 
  cos_newThetaY = cos(muonThetaY + muonThetaY_shift); 

  //This is merely: sec^2(theta_X) + sec^2(theta_Y) = sec^2(theta) +1 
  cos_muonTheta_shifted = 1/(sqrt(1.0/pow(cos_newThetaX, 2) + 1.0/pow(cos_newThetaY, 2) -1) ); 

  muonPT_prime = muonP*sin(acos(cos_muonTheta_shifted)); 
  muonPT       = muonP*sin(muonTheta); 

  shift = muonPT_prime - muonPT; 

  return ok; 
}

bool MnvRecoShifter::GetShift_ccqe_muonPZ_muonTheta(double &shift, const double muonP, const double muonTheta, const double muonThetaErr, const int nUniverse, const double nSigma) const 
{ 

  bool ok = true; 
  double muonTheta_shift = 0.0, muonPZ_prime = 0.0, muonPZ = 0.0; 

  ok = ok && GetShift_ccqe_muonTheta(muonTheta_shift, muonThetaErr, nUniverse, nSigma); 
  muonPZ_prime = muonP*cos(muonTheta+muonTheta_shift); 
  muonPZ       = muonP*cos(muonTheta); 

  shift = muonPZ_prime - muonPZ; 

  return ok; 
} 

bool MnvRecoShifter::GetShift_ccqe_muonPZ_muonThetaXY(double &shift, const double muonP, const double muonTheta, const double muonThetaX, const double muonThetaX_Err, const double muonThetaY, const double muonThetaY_Err, const int nUniverse, const double nSigma) const 
{ 
 
  bool ok = true; 
  double muonThetaX_shift = 0.0, muonThetaY_shift = 0.0, muonPZ_prime = 0.0, muonPZ = 0.0; 
  double cos_newThetaX = 0.0, cos_newThetaY = 0.0, cos_muonTheta_shifted = 0.0;  

  ok = ok && GetShift_ccqe_muonTheta(muonThetaX_shift, muonThetaX_Err, nUniverse, nSigma); 

  ok = ok && GetShift_ccqe_muonTheta(muonThetaY_shift, muonThetaY_Err, nUniverse, nSigma); 

  cos_newThetaX = cos(muonThetaX + muonThetaX_shift); 
  cos_newThetaY = cos(muonThetaY + muonThetaY_shift); 

  //This is merely: sec^2(theta_X) + sec^2(theta_Y) = sec^2(theta) +1 
  cos_muonTheta_shifted = 1/(sqrt(1.0/pow(cos_newThetaX, 2) + 1.0/pow(cos_newThetaY, 2) -1) ); 

  muonPZ_prime = muonP*cos_muonTheta_shifted; 
  muonPZ       = muonP*cos(muonTheta); 

  shift = muonPZ_prime - muonPZ; 

  return ok; 
}

//-----------------------------------
// Binding energy shifts 
//-----------------------------------
bool MnvRecoShifter::GetShift_ccqe_q2_bindingEnergy(double &shift, const double muonP, const double muonTheta, const double bindingE, const double bindingE_Err, const int nUniverse /* = -1 */, const double nSigma /*= 1.*/) const
{

  bool ok = true;
  double q2_prime = 0., q2 = 0., bindingE_shift = 0.;

  ok = ok && GetShift_ccqe_muonP( bindingE_shift, bindingE_Err, nUniverse, nSigma);
 
  ok = ok && Calc_ccqe_q2( q2_prime, muonP, muonTheta, bindingE+bindingE_shift);
  ok = ok && Calc_ccqe_q2( q2, muonP, muonTheta, bindingE);

  shift = q2_prime - q2 ;

  return ok;
}

//-----------------------------------
// Recoil energy shifts 
//-----------------------------------
bool MnvRecoShifter::GetShift_ccqe_recoilE( double& shift, const double muonRecoilerr,
    const int nUniverse /* = -1 */, const double nSigma /*= 1.*/ ) const
{

  //valid universe numbers: [ 1, m_nUniverses ]
  if ( nUniverse < 1 || m_nUniverses < nUniverse )
  {
    Error("MnvRecoShifter::GetShift_ccqe_muonTheta", "You need to specify a valid Universe number.");
    return false;
  }

  shift = muonRecoilerr * m_rShiftHadron[nUniverse-1] * nSigma ;

  return true;
}

//==========================================
// END - CCQE Calculators and Scale Shifts
//==========================================
//================================= 
// BEGIN - Muon Angle Scale Shifts
//================================= 
bool MnvRecoShifter::GetShifts_muonTheta( dVec& shifts,
    const double muonTheta, const double muonE, const double hadronicE,
    const double muonThetaErr) const
{
  shifts.clear();
  double shift_muonTheta = 0.;
  bool allOK = true;
  //! Loop over all universes of deviation from CV and find the shift amount
  for( dVecCItr i = m_rShiftMuon.begin(); i != m_rShiftMuon.end(); ++i )
  {
    bool ok = GetShift_muonTheta( shift_muonTheta, muonTheta, muonE, hadronicE, muonThetaErr, *i );
    allOK = allOK && ok;
    shifts.push_back( shift_muonTheta ); //record shift in muonTheta
  }
  return allOK;
}


//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShift_muonTheta( double& shift,
    const double muonTheta, const double muonE, const double hadronicE,
    const double muonThetaErr,
    const double nSigma /*= 1.*/ )const
{
  bool ok = true;
  double muonEerr = 0.;

  double muonE_prime     = muonE     + muonEerr*nSigma;
  double hadronicE_prime = hadronicE;
  double E_prime         = muonE_prime + hadronicE_prime;

  double muonTheta_prime = muonTheta + muonThetaErr * nSigma;

  //! If any energy is out of bounds, then set the shift to the NotPhysicalShiftNumber
  if( !ok || OutOfBounds_E( E_prime ) || OutOfBounds_muonE( muonE_prime ) || OutOfBounds_hadronicE( hadronicE_prime ) )
    shift = MnvHist::NotPhysicalShiftNumber;
  else
    shift = muonTheta_prime - muonTheta;

  return ok;
}

//~~~~~~~~~~~~~~~~~~~

//================================= 
// BEGIN - Muon Energy Scale Shifts
//================================= 

bool MnvRecoShifter::GetShifts_muonE( dVec& shifts,
    const double muonE, const double hadronicE,
    const double muonEerr) const
{
  shifts.clear();
  double shift_muonE = 0.;
  bool allOK = true;
  //! Loop over all universes of deviation from CV and find the shift amount
  for( dVecCItr i = m_rShiftMuon.begin(); i != m_rShiftMuon.end(); ++i )
  {
    bool ok = GetShift_muonE( shift_muonE, muonE, hadronicE, muonEerr, *i );
    allOK = allOK && ok;
    shifts.push_back( shift_muonE ); //record shift in muonE
  }
  return allOK;
}


//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShift_muonE( double& shift,
    const double muonE, const double hadronicE,
    const double muonEerr,
    const double nSigma /*= 1.*/ )const
{
  bool ok = true;

  double muonE_prime     = muonE     + muonEerr*nSigma;
  double hadronicE_prime = hadronicE;
  double E_prime         = muonE_prime + hadronicE_prime;

  //! If any energy is out of bounds, then set the shift to the NotPhysicalShiftNumber
  if( !ok || OutOfBounds_E( E_prime ) || OutOfBounds_muonE( muonE_prime ) || OutOfBounds_hadronicE( hadronicE_prime ) )
    shift = MnvHist::NotPhysicalShiftNumber;
  else
    shift = muonE_prime - muonE;

  return ok;
}

//~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShifts_E_muonE( dVec& shifts,
    const double muonE, const double hadronicE,
    const double muonEerr) const
{
  shifts.clear();
  double shift_E = 0.;
  bool allOK = true;
  //! Loop over all universes of deviation from CV and find the shift amount
  for( dVecCItr i = m_rShiftMuon.begin(); i != m_rShiftMuon.end(); ++i )
  {
    bool ok = GetShift_E_muonE( shift_E, muonE, hadronicE, muonEerr, *i );
    allOK = allOK && ok;
    shifts.push_back( shift_E ); //record shift in E
  }
  return allOK;
}

//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShift_E_muonE( double& shift,
    const double muonE, const double hadronicE, const double muonEerr,
    const double nSigma /*= 1.*/) const
{
  bool ok = true;

  //! Calculate the shifted quantities
  double muonE_shift     = 0.;

  ok = ok && GetShift_muonE( muonE_shift, muonE, hadronicE, muonEerr, nSigma );

  //! E = hadronicE + muonE
  double E = hadronicE + muonE;
  double muonE_prime = muonE + muonE_shift;
  double hadronicE_prime = hadronicE;
  double E_prime = hadronicE_prime + muonE_prime;

  //! If any energy is out of bounds, then set the shift to the NotPhysicalShiftNumber
  if( !ok || OutOfBounds_E( E_prime ) || OutOfBounds_muonE( muonE_prime ) || OutOfBounds_hadronicE( hadronicE_prime ) )
    shift = MnvHist::NotPhysicalShiftNumber;
  else
    shift = E_prime - E;

  return true;
}

//~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShifts_q2_muonE( dVec& shifts,
    const double muonE, const double hadronicE, const double muonTheta,
    const double muonEerr) const
{
  shifts.clear();
  double shift_q2 = 0.;
  bool allOK = true;
  //! Loop over all universes of deviation from CV and find the shift amount
  for( dVecCItr i = m_rShiftMuon.begin(); i != m_rShiftMuon.end(); ++i )
  {
    bool ok = GetShift_q2_muonE( shift_q2, muonE, hadronicE, muonTheta, muonEerr, *i );
    allOK = allOK && ok;
    shifts.push_back( shift_q2 ); //record shift in q2
  }
  return allOK;
}

//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShift_q2_muonE( double& shift,
    const double muonE, const double hadronicE, const double muonTheta, const double muonEerr,
    const double nSigma /*= 1.*/) const
{
  bool ok = true;
#if DEBUGLEVEL > DEBUGQ2SHIFT
  std::cout << "GetShift_q2_muonE" << std::endl;
  std::cout << "  shift = " << shift << "; muonE = " << muonE << "; hadronicE = " << hadronicE 
    << "; muonTheta = " << muonTheta << "; muonEerr = " << muonEerr << "; nSigma = " << nSigma << std::endl;
#endif

  //! Calculate the shifted quantities
  //! how do we handle correlations?
  double muonE_shift     = 0.;
  double hadronicE_shift = 0.;
  double muonTheta_shift = 0.;

  ok = ok && GetShift_muonE( muonE_shift, muonE, hadronicE, muonEerr, nSigma);
  // ok = ok && GetShift_muonTheta_muonE( ok, muonTheta, muonE, muonTerr, nSigma );// <--- someday?
#if DEBUGLEVEL > DEBUGQ2SHIFT
  std::cout << " After GetShift_muonE..." << std::endl;
  std::cout << " muonE_shift = " << muonE_shift << std::endl;
#endif

  //! Calculate Q2 with CV and shifted values
  double q2 = 0, q2_prime = 0.;
  ok = ok && Calc_q2( q2,       muonE, hadronicE, muonTheta );
  ok = ok && Calc_q2( q2_prime, muonE + muonE_shift, hadronicE + hadronicE_shift, muonTheta + muonTheta_shift );
#if DEBUGLEVEL > DEBUGQ2SHIFT
  std::cout << " q2 = " << q2 << "; q2_prime = " << q2_prime << "; shift = " << (q2_prime - q2) << std::endl;
#endif


  double muonE_prime = muonE + muonE_shift;
  double hadronicE_prime = hadronicE + hadronicE_shift;
  double E_prime = hadronicE_prime + muonE_prime;

  //! If any energy is out of bounds, then set the shift to the NotPhysicalShiftNumber
  if( !ok || OutOfBounds_E( E_prime ) || OutOfBounds_muonE( muonE_prime ) || OutOfBounds_hadronicE( hadronicE_prime ) )
    shift = MnvHist::NotPhysicalShiftNumber;
  else
    shift = q2_prime - q2;

  return ok;
}

//~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShifts_w_muonE( dVec& shifts,
    const double muonE, const double hadronicE, const double muonTheta,
    const double muonEerr) const
{
  shifts.clear();
  double shift_w = 0.;
  bool allOK = true;
  //! Loop over all universes of deviation from CV and find the shift amount
  for( dVecCItr i = m_rShiftMuon.begin(); i != m_rShiftMuon.end(); ++i )
  {
    bool ok = GetShift_w_muonE( shift_w, muonE, hadronicE, muonTheta, muonEerr, *i );
    allOK = allOK && ok;
    shifts.push_back( shift_w ); //record shift in W
  }
  return allOK;
}

//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShift_w_muonE( double& shift,
    const double muonE, const double hadronicE, const double muonTheta, const double muonEerr,
    const double nSigma /*= 1.*/)const
{
  bool ok = true;

  //! Calculate the shifted quantities
  //! how do we handle correlations?
  double muonE_shift     = 0.;
  double hadronicE_shift = 0.;
  double muonTheta_shift = 0.;

  ok = ok && GetShift_muonE( muonE_shift, muonE, hadronicE, muonEerr, nSigma );
  // ok = ok && GetShift_muonTheta_muonE( ok, muonTheta, muonE, muonTerr, nSigma );// <--- someday?

  //! Calculate W with CV and shifted values
  double w = 0, w_prime = 0.;
  ok = ok && Calc_w( w,       muonE, hadronicE, muonTheta );
  ok = ok && Calc_w( w_prime, muonE + muonE_shift, hadronicE + hadronicE_shift, muonTheta + muonTheta_shift );

  double muonE_prime = muonE + muonE_shift;
  double hadronicE_prime = hadronicE + hadronicE_shift;
  double E_prime = hadronicE_prime + muonE_prime;

  //! If any energy is out of bounds (or W is unphysical), then set the shift to the NotPhysicalShiftNumber
  if( !ok || OutOfBounds_E( E_prime ) || OutOfBounds_muonE( muonE_prime ) || OutOfBounds_hadronicE( hadronicE_prime ) )
    shift = MnvHist::NotPhysicalShiftNumber;
  else
    shift = w_prime - w;

  return ok;
}

//~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShifts_x_muonE( dVec& shifts,
    const double muonE, const double hadronicE, const double muonTheta,
    const double muonEerr) const
{
  shifts.clear();
  double shift_x = 0.;
  bool allOK = true;
  //! Loop over all universes of deviation from CV and find the shift amount
  for( dVecCItr i = m_rShiftMuon.begin(); i != m_rShiftMuon.end(); ++i )
  {
    bool ok = GetShift_x_muonE( shift_x, muonE, hadronicE, muonTheta, muonEerr, *i );
    allOK = allOK && ok;
    shifts.push_back( shift_x ); //record shift in x
  }
  return allOK;
}

//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShift_x_muonE( double& shift,
    const double muonE, const double hadronicE, const double muonTheta, const double muonEerr,
    const double nSigma /*= 1.*/) const
{
  bool ok = true;

  //! Calculate the shifted quantities
  double muonE_shift     = 0.;
  double hadronicE_shift = 0.;
  double muonTheta_shift = 0.; 

  ok = ok && GetShift_muonE( muonE_shift, muonE, hadronicE, muonEerr, nSigma );
  // ok = ok && GetShift_muonTheta_muonE( ok, muonTheta, muonE, muonTerr, nSigma );// <--- someday?

  //! Calculate x with CV and shifted values
  double x = 0., x_prime = 0.;
  ok = ok && Calc_x( x,       muonE, hadronicE, muonTheta );
  ok = ok && Calc_x( x_prime, muonE + muonE_shift, hadronicE + hadronicE_shift, muonTheta + muonTheta_shift );

  double muonE_prime = muonE + muonE_shift;
  double hadronicE_prime = hadronicE + hadronicE_shift;
  double E_prime = hadronicE_prime + muonE_prime;

  //! If any energy is out of bounds (or W is unphysical), then set the shift to the NotPhysicalShiftNumber
  if( !ok || OutOfBounds_E( E_prime ) || OutOfBounds_muonE( muonE_prime ) || OutOfBounds_hadronicE( hadronicE_prime ) )
    shift = MnvHist::NotPhysicalShiftNumber;
  else
    shift = x_prime - x;

  return ok;
}

//~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShifts_y_muonE( dVec& shifts,
    const double muonE, const double hadronicE,
    const double muonEerr) const
{
  shifts.clear();
  double shift_y = 0.;
  bool allOK = true;
  //! Loop over all universes of deviation from CV and find the shift amount
  for( dVecCItr i = m_rShiftMuon.begin(); i != m_rShiftMuon.end(); ++i )
  {
    bool ok = GetShift_y_muonE( shift_y, muonE, hadronicE, muonEerr, *i );
    allOK = allOK && ok;
    shifts.push_back( shift_y ); //record shift in y
  }
  return allOK;
}

//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShift_y_muonE( double& shift,
    const double muonE, const double hadronicE, const double muonEerr,
    const double nSigma /*= 1.*/) const
{
  bool ok = true;

  //! Calculate the shifted quantities
  double muonE_shift     = 0.;
  double hadronicE_shift = 0.;

  ok = ok && GetShift_muonE( muonE_shift, muonE, hadronicE, muonEerr, nSigma );

  //! Calculate y with CV and shifted values
  double y = 0., y_prime = 0.;
  ok = ok && Calc_y( y, muonE, hadronicE );
  ok = ok && Calc_y( y_prime, muonE + muonE_shift, hadronicE + hadronicE_shift );

  double muonE_prime = muonE + muonE_shift;
  double hadronicE_prime = hadronicE + hadronicE_shift;
  double E_prime = hadronicE_prime + muonE_prime;

  //! If any energy is out of bounds (or W is unphysical), then set the shift to the NotPhysicalShiftNumber
  if( !ok || OutOfBounds_E( E_prime ) || OutOfBounds_muonE( muonE_prime ) || OutOfBounds_hadronicE( hadronicE_prime ) )
    shift = MnvHist::NotPhysicalShiftNumber;
  else
    shift = y_prime - y;

  return ok;
}

//================================= 
// END - Muon Energy Shifts
//================================= 

//================================= 
// BEGIN - Hadronic Energy Shifts
//=================================

//~~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShifts_hadronicE( dVec& shifts,
    const double muonE, const double hadronicE ) const
{
  shifts.clear();
  double shift_hadronicE = 0.;
  bool allOK = true;
  //! Loop over all universes of deviation from CV and find the shift amount
  for( dVecCItr i = m_rShiftHadron.begin(); i != m_rShiftHadron.end(); ++i )
  {
    bool ok = GetShift_hadronicE( shift_hadronicE, muonE, hadronicE, *i );
    //cout << "Ehad shift: " << shift_hadronicE << endl;
    allOK = allOK && ok;
    shifts.push_back( shift_hadronicE ); //record shift in hadronicE 
  }
  return allOK;
}

//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShift_hadronicE( double& shift,
    const double muonE, const double hadronicE,
    const double nSigma) const
{
  //! Old method - Use 10% for all : shift = hadronicE * 0.10;
  //! New method - J. Devan's function. Careful about units!
  bool ok = true;

  //! Change hadronic energy input from MeV to GeV.
  double nml = hadronicE * (1./1000.);
  //! Calculate the shift.
  shift      = CCnumu_calorimetric_err( nml, m_isAntiNu );
  //! Shift arrives as a percentage, not an absolute. Change it to an absolute.
  shift     *= nml;
  //! Change the shift units back from GeV to MeV.
  shift     *= 1000.;
  //! Now, multiply the result by the number (and sign!) of "sigma" requested.
  shift     *= nSigma;

  //! Calculate the shifted energies
  double muonE_prime     = muonE; //silly but consistent
  double hadronicE_prime = hadronicE + shift;
  double E_prime = hadronicE_prime + muonE_prime;

  //! If any shifted energy is out of bounds, then set the shift to the NotPhysicalShiftNumber
  if( !ok || OutOfBounds_E( E_prime ) || OutOfBounds_muonE( muonE_prime ) || OutOfBounds_hadronicE( hadronicE_prime ) )
    shift = MnvHist::NotPhysicalShiftNumber;

  return ok;
}

//~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShifts_E_hadronicE( dVec& shifts, 
    const double muonE, const double hadronicE ) const
{
  shifts.clear();
  double shift_E = 0.;
  bool allOK = true;
  //! Loop over all universes of deviation from CV and find the shift amount
  for( dVecCItr i = m_rShiftHadron.begin(); i != m_rShiftHadron.end(); ++i )
  {
    bool ok = GetShift_E_hadronicE( shift_E, muonE, hadronicE, *i );
    allOK = allOK && ok;
    shifts.push_back( shift_E ); //record shift in E
  }
  return allOK;
}

//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShift_E_hadronicE( double& shift,
    const double muonE, const double hadronicE, 
    const double nSigma)const
{
  bool ok = true;

  //! Calculate the shifted quantities
  double hadronicE_shift = 0.;

  ok = ok && GetShift_hadronicE( hadronicE_shift, muonE, hadronicE, nSigma );

  double muonE_prime     = muonE;
  double hadronicE_prime = hadronicE + hadronicE_shift;

  //! E = hadronicE + muonE
  double E       = hadronicE + muonE;
  double E_prime = hadronicE_prime + muonE_prime;

  //! If any energy is out of bounds, then set the shift to the NotPhysicalShiftNumber
  if( !ok || OutOfBounds_E( E_prime ) || OutOfBounds_muonE( muonE_prime ) || OutOfBounds_hadronicE( hadronicE_prime ) )
    shift = MnvHist::NotPhysicalShiftNumber;
  else
    shift = E_prime - E;

  return ok;
}

//~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShifts_q2_hadronicE( dVec& shifts,
    const double muonE, const double hadronicE, const double muonTheta ) const
{
  shifts.clear();
  double shift_q2 = 0.;
  bool allOK = true;
  //! Loop over all universes of deviation from CV and find the shift amount
  for( dVecCItr i = m_rShiftHadron.begin(); i != m_rShiftHadron.end(); ++i )
  {
    bool ok = GetShift_q2_hadronicE( shift_q2, muonE, hadronicE, muonTheta, *i );
    allOK = allOK && ok;
    shifts.push_back( shift_q2 );
  }
  return allOK;
}


//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShift_q2_hadronicE( double& shift,
    const double muonE, const double hadronicE, const double muonTheta,
    const double nSigma) const
{
  bool ok = true;
#if DEBUGLEVEL > DEBUGQ2SHIFT
  std::cout << "GetShift_q2_hadronicE" << std::endl;
  std::cout << "  shift = " << shift << "; muonE = " << muonE << "; hadronicE = " << hadronicE 
    << "; muonTheta = " << muonTheta << "; nSigma = " << nSigma << std::endl;
#endif

  //! Calculate the shifted quantities
  //! how do we handle correlations?
  double muonE_shift     = 0.;
  double hadronicE_shift = 0.;
  double muonTheta_shift = 0.;

  ok = ok && GetShift_hadronicE( hadronicE_shift, muonE, hadronicE, nSigma );
#if DEBUGLEVEL > DEBUGQ2SHIFT
  std::cout << " After GetShift_hadronicE..." << std::endl;
  std::cout << " hadronicE_shift = " << hadronicE_shift << std::endl;
#endif

  //! Calculate Q2 with CV and shifted values
  double q2 = 0, q2_prime = 0.;
  ok = ok && Calc_q2( q2, muonE, hadronicE, muonTheta );
  ok = ok && Calc_q2( q2_prime, muonE + muonE_shift, hadronicE + hadronicE_shift, muonTheta + muonTheta_shift );
#if DEBUGLEVEL > DEBUGQ2SHIFT
  std::cout << " q2 = " << q2 << "; q2_prime = " << q2_prime << "; shift = " << (q2_prime - q2) << std::endl;
#endif

  double muonE_prime     = muonE     + muonE_shift;
  double hadronicE_prime = hadronicE + hadronicE_shift;
  double E_prime         = muonE_prime + hadronicE_prime;

  //! If any energy is out of bounds, then set the shift to the NotPhysicalShiftNumber
  if( !ok || OutOfBounds_E( E_prime ) || OutOfBounds_muonE( muonE_prime ) || OutOfBounds_hadronicE( hadronicE_prime ) )
    shift = MnvHist::NotPhysicalShiftNumber;
  else
    shift = q2_prime - q2;

  return ok;
}

//~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShifts_w_hadronicE( dVec& shifts,
    const double muonE, const double hadronicE, const double muonTheta ) const
{
  shifts.clear();
  double shift_w = 0.;
  bool allOK = true;
  //! Loop over all universes of deviation from CV and find the shift amount
  for( dVecCItr i = m_rShiftHadron.begin(); i != m_rShiftHadron.end(); ++i )
  {
    bool ok = GetShift_w_hadronicE( shift_w, muonE, hadronicE, muonTheta, *i );
    allOK = allOK && ok;
    shifts.push_back( shift_w );
  }
  return allOK;
}

//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShift_w_hadronicE( double& shift,
    const double muonE, const double hadronicE, const double muonTheta,
    const double nSigma) const
{
  bool ok = true;

  //! Calculate the shifted quantities
  //! how do we handle correlations?
  double muonE_shift     = 0.;
  double hadronicE_shift = 0.;
  double muonTheta_shift = 0.;

  ok = ok && GetShift_hadronicE( hadronicE_shift, muonE, hadronicE, nSigma );

  //! Calculate W with CV and shifted values
  double w = 0, w_prime = 0.;
  ok = ok && Calc_w( w, muonE, hadronicE, muonTheta );
  ok = ok && Calc_w( w_prime, muonE + muonE_shift, hadronicE + hadronicE_shift, muonTheta + muonTheta_shift );

  double muonE_prime     = muonE     + muonE_shift;
  double hadronicE_prime = hadronicE + hadronicE_shift;
  double E_prime         = muonE_prime + hadronicE_prime;

  //! If any energy is out of bounds, then set the shift to the NotPhysicalShiftNumber
  if( !ok || OutOfBounds_E( E_prime ) || OutOfBounds_muonE( muonE_prime ) || OutOfBounds_hadronicE( hadronicE_prime ) )
    shift = MnvHist::NotPhysicalShiftNumber;
  else
    shift = w_prime - w;

  return ok;
}

//~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShifts_x_hadronicE( dVec& shifts,
    const double muonE, const double hadronicE, const double muonTheta ) const
{
  shifts.clear();
  double shift_x = 0.;
  bool allOK = true;
  //! Loop over all universes of deviation from CV and find the shift amount
  for( dVecCItr i = m_rShiftHadron.begin(); i != m_rShiftHadron.end(); ++i )
  {
    bool ok = GetShift_x_hadronicE( shift_x, muonE, hadronicE, muonTheta, *i );
    allOK = allOK && ok;
    shifts.push_back( shift_x );
  }
  return allOK;
}

//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShift_x_hadronicE( double& shift,
    const double muonE, const double hadronicE, const double muonTheta,
    const double nSigma)const
{
  bool ok = true;

  //! Calculate the shifted quantities
  double muonE_shift     = 0.;
  double hadronicE_shift = 0.;
  double muonTheta_shift = 0.;

  ok = ok && GetShift_hadronicE( hadronicE_shift, muonE, hadronicE, nSigma);

  //! Calculate x with CV and shifted values
  double x = 0., x_prime = 0.;
  ok = ok && Calc_x( x, muonE, hadronicE, muonTheta );
  ok = ok && Calc_x( x_prime, muonE + muonE_shift, hadronicE + hadronicE_shift, muonTheta + muonTheta_shift );

  double muonE_prime     = muonE     + muonE_shift;
  double hadronicE_prime = hadronicE + hadronicE_shift;
  double E_prime         = muonE_prime + hadronicE_prime;

  //! If any energy is out of bounds, then set the shift to the NotPhysicalShiftNumber
  if( !ok || OutOfBounds_E( E_prime ) || OutOfBounds_muonE( muonE_prime ) || OutOfBounds_hadronicE( hadronicE_prime ) )
    shift = MnvHist::NotPhysicalShiftNumber;
  else
    shift = x_prime - x;

  return ok;
}

//~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShifts_y_hadronicE( dVec& shifts,
    const double muonE, const double hadronicE ) const
{
  shifts.clear();
  double shift_y = 0.;
  bool allOK = true;
  //! Loop over all universes of deviation from CV and find the shift amount
  for( dVecCItr i = m_rShiftHadron.begin(); i != m_rShiftHadron.end(); ++i )
  {
    bool ok = GetShift_y_hadronicE( shift_y, muonE, hadronicE, *i );
    allOK = allOK && ok;
    shifts.push_back( shift_y );
  }
  return allOK;
}


//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShift_y_hadronicE( double& shift,
    const double muonE, const double hadronicE, 
    const double nSigma)const
{
  bool ok = true;

  //! Calculate the shifted quantities
  double muonE_shift     = 0.;
  double hadronicE_shift = 0.;

  ok = ok && GetShift_hadronicE( hadronicE_shift, muonE, hadronicE, nSigma );

  //! Calculate y with CV and shifted values
  double y = 0., y_prime = 0.;
  ok = ok && Calc_y( y, muonE, hadronicE );
  ok = ok && Calc_y( y_prime, muonE + muonE_shift, hadronicE + hadronicE_shift );

  double muonE_prime     = muonE     + muonE_shift;
  double hadronicE_prime = hadronicE + hadronicE_shift;
  double E_prime         = muonE_prime + hadronicE_prime;

  //! If any energy is out of bounds, then set the shift to the NotPhysicalShiftNumber
  if( !ok || OutOfBounds_E( E_prime ) || OutOfBounds_muonE( muonE_prime ) || OutOfBounds_hadronicE( hadronicE_prime ) )
    shift = MnvHist::NotPhysicalShiftNumber;
  else
    shift = y_prime - y;

  return ok;
}

//~~~~~~~~~~~~~~~~~~~
// The attached function returns the calorimetric systematic error from the histograms in J. Devan's 
// presentation (docdb 7605: http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=7605). 
// The first argument, nml, is the neutrino minus lepton energy in GeV. Set the second argument, 
// anti, to true if you want the anti-neutrino result.
double MnvRecoShifter::CCnumu_calorimetric_err( double nml, 
    bool anti /*= false */ ) const
{
  const int n_bins = 31;
  const double x[] = {0.15,0.252356,0.359659,0.47241,0.591194,0.716692,0.849712,0.991212,1.14235,1.30454,1.47952,1.66949,1.87728,2.10657,2.36236,2.65161,2.98448,3.37659,3.85402,4.4654,5.32045,6.79052,8.75278,10.7528,12.7528,14.7528,16.7528,18.7528,20.7528,22.7528,24.7528};
  const double anti_y[] = {0.141761,0.133744,0.112329,0.111782,0.109596,0.105645,0.101442,0.0989989,0.0964361,0.0949933,0.0923315,0.0901037,0.0886676,0.0855434,0.0828689,0.0814455,0.0786826,0.0769916,0.0758833,0.0754399,0.0737412,0.0723754,0.0699579,0.0676747,0.0667669,0.0642251,0.0625776,0.0652039,0.0613044,0.0615586,0.0585469};
  const double numu_y[] = {0.107806,0.102049,0.0928338,0.0923239,0.0924461,0.0929755,0.0929152,0.0923734,0.0921905,0.0913156,0.0905763,0.0891249,0.0872907,0.0850411,0.0834782,0.0807297,0.078877,0.0770799,0.076058,0.074665,0.0738296,0.0725464,0.0704331,0.0671367,0.0658373,0.0643413,0.0627155,0.0627347,0.0613374,0.0596863,0.0603769};
  const double *y = (anti)? anti_y : numu_y;

  if(nml < x[0]) return y[0];
  for(int i = 1; i < n_bins; i++)
  {
    if(nml < x[i]) return y[i-1]+(y[i]-y[i-1])/(x[i]-x[i-1])*(nml-x[i-1]); // linear interpolation
  }
  return y[n_bins-1];
}


//================================= 
// END - Hadronic Energy Shifts
//================================= 




//================================= 
// BEGIN - Muon Theta Shifts
//=================================

double MnvRecoShifter::GetMuonTheta_err( ) const 
{
  // Muon angle residual fits for CC inclusive show
  // ThetaX and ThetaY each have sigma < .7mrad 
  // and mean consistent enough with 0.
  // See TN 25 for resolution plots
  // Be generous and round up to 1mrad.
  return .001; 
}

//~~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShifts_muonTheta( dVec& shifts,
    const double muonTheta) const
{
  shifts.clear();
  double shift_muonTheta = 0.;
  bool allOK = true;
  // Loop over all universes of deviation from CV and find the shift amount
  for( dVecCItr i = m_rShiftMuonTheta.begin(); i != m_rShiftMuonTheta.end(); ++i )
  {
    bool ok = GetShift_muonTheta( shift_muonTheta, muonTheta, *i );
    allOK = allOK && ok;
    shifts.push_back( shift_muonTheta ); //record shift in muonTheta
  }
  return allOK;
}

//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShift_muonTheta( double& shift,
    const double muonTheta,
    const double nSigma) const
{
  bool ok = true;

  // Calculate the shift.
  shift      = GetMuonTheta_err();
  // Now, multiply the result by the number (and sign!) of "sigma" requested.
  shift     *= nSigma;

  // Calculate the shifted energies
  double muonTheta_prime = muonTheta + shift;

  //need to check for negative angles?
  //NotPhysicalShift's name suggests we use that, but I think it's a different situation because the theta shift shouldn't really go negative.
  if( muonTheta_prime < 0. )
  {
    muonTheta_prime = 0.;
    shift = -muonTheta;
  }

  // If any shifted angle is out of bounds, then set the shift to the NotPhysicalShiftNumber
  if( !ok || OutOfBounds_muonTheta( muonTheta_prime ) )
    shift = MnvHist::NotPhysicalShiftNumber;

  return ok;
}


//~~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShifts_q2_muonTheta( dVec& shifts,
    const double muonE, const double hadronicE, const double muonTheta) const
{
  shifts.clear();
  double shift_q2 = 0.;
  bool allOK = true;
  // Loop over all universes of deviation from CV and find the shift amount
  for( dVecCItr i = m_rShiftMuonTheta.begin(); i != m_rShiftMuonTheta.end(); ++i )
  {
    bool ok = GetShift_q2_muonTheta( shift_q2, muonE, hadronicE, muonTheta, *i );
    allOK = allOK && ok;
    shifts.push_back( shift_q2 ); //record shift in q2
  }
  return allOK;
}


//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShift_q2_muonTheta( double& shift,
    const double muonE, const double hadronicE, const double muonTheta,
    const double nSigma) const
{
  bool ok = true;

  // Calculate the shifted quantities
  double muonE_shift     = 0.;
  double hadronicE_shift = 0.;
  double muonTheta_shift = 0.;

  ok = ok && GetShift_muonTheta( muonTheta_shift, muonTheta, nSigma );

  // Calculate Q2 with CV and shifted values
  double q2 = 0, q2_prime = 0.;
  ok = ok && Calc_q2( q2, muonE, hadronicE, muonTheta );
  ok = ok && Calc_q2( q2_prime, muonE + muonE_shift, hadronicE + hadronicE_shift, muonTheta + muonTheta_shift );

  double muonTheta_prime = muonTheta + muonTheta_shift;

  // If any energy is out of bounds, then set the shift to the NotPhysicalShiftNumber
  if( !ok || OutOfBounds_muonTheta( muonTheta_prime ) )
    shift = MnvHist::NotPhysicalShiftNumber;
  else
    shift = q2_prime - q2;

  return ok;
}

//~~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShifts_x_muonTheta( dVec& shifts,
    const double muonE, const double hadronicE, const double muonTheta) const
{
  shifts.clear();
  double shift_x = 0.;
  bool allOK = true;
  // Loop over all universes of deviation from CV and find the shift amount
  for( dVecCItr i = m_rShiftMuonTheta.begin(); i != m_rShiftMuonTheta.end(); ++i )
  {
    bool ok = GetShift_x_muonTheta( shift_x, muonE, hadronicE, muonTheta, *i );
    allOK = allOK && ok;
    shifts.push_back( shift_x ); //record shift in x
  }
  return allOK;
}

//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShift_x_muonTheta( double& shift,
    const double muonE, const double hadronicE, const double muonTheta,
    const double nSigma) const
{
  bool ok = true;

  // Calculate the shifted quantities
  double muonE_shift     = 0.;
  double hadronicE_shift = 0.;
  double muonTheta_shift = 0.;

  ok = ok && GetShift_muonTheta( muonTheta_shift, muonTheta, nSigma );

  // Calculate x with CV and shifted values
  double x = 0, x_prime = 0.;
  ok = ok && Calc_x( x, muonE, hadronicE, muonTheta );
  ok = ok && Calc_x( x_prime, muonE + muonE_shift, hadronicE + hadronicE_shift, muonTheta + muonTheta_shift );

  double muonTheta_prime = muonTheta + muonTheta_shift;

  // If any energy is out of bounds, then set the shift to the NotPhysicalShiftNumber
  if( !ok || OutOfBounds_muonTheta( muonTheta_prime ) )
    shift = MnvHist::NotPhysicalShiftNumber;
  else
    shift = x_prime - x;

  return ok;
}


bool MnvRecoShifter::GetShifts_w_muonTheta( dVec& shifts,
    const double muonE, const double hadronicE, const double muonTheta) const
{
  shifts.clear();
  double shift_w = 0.;
  bool allOK = true;
  // Loop over all universes of deviation from CV and find the shift amount
  for( dVecCItr i = m_rShiftMuonTheta.begin(); i != m_rShiftMuonTheta.end(); ++i )
  {
    bool ok = GetShift_w_muonTheta( shift_w, muonE, hadronicE, muonTheta, *i );
    allOK = allOK && ok;
    shifts.push_back( shift_w ); //record shift in w
  }
  return allOK;
}

//~~~~~~~~~~~~~~~~~~~
bool MnvRecoShifter::GetShift_w_muonTheta( double& shift,
    const double muonE, const double hadronicE, const double muonTheta,
    const double nSigma) const
{
  bool ok = true;

  // Calculate the shifted quantities
  double muonE_shift     = 0.;
  double hadronicE_shift = 0.;
  double muonTheta_shift = 0.;

  ok = ok && GetShift_muonTheta( muonTheta_shift, muonTheta, nSigma );

  // Calculate w with CV and shifted values
  double w = 0, w_prime = 0.;
  ok = ok && Calc_w( w, muonE, hadronicE, muonTheta );
  ok = ok && Calc_w( w_prime, muonE + muonE_shift, hadronicE + hadronicE_shift, muonTheta + muonTheta_shift );

  double muonTheta_prime = muonTheta + muonTheta_shift;

  // If any energy is out of bounds, then set the shift to the NotPhysicalShiftNumber
  if( !ok || OutOfBounds_muonTheta( muonTheta_prime ) )
    shift = MnvHist::NotPhysicalShiftNumber;
  else
    shift = w_prime - w;

  return ok;
}


//===========================================================
// BEGIN - Get a Vector of Gaus Random Initialized variables
//         for Lateral multiuniverses filling purposes
//===========================================================
bool MnvRecoShifter::GetRandomShiftVector(std::vector<double>& rShift, const double seed, const int nUniverses /* = 1000 */ )
{
  if (nUniverses<=0)
  {
    Error("MnvRecoShifter::GetRandomShiftVector", "Cannot set number of Universes less or equal than 0");
    return false;
  }
  rShift.clear();

  //! If 2 universes are specified, then make the -1 and +1 sigma
  if( 2 == nUniverses )
  {
    rShift.push_back(-1);
    rShift.push_back(+1);
  }
  else
  {
    if( nUniverses < 10 )
      Warning( "MnvRecoShifter::GetRandomShiftVector", Form( "Using a Gaussian distribution of errors in only %d universes.  Results may not make sense.",nUniverses) );

    gRandom->SetSeed(seed); //Setting a fixed seed
    for (int i=0; i < nUniverses ; ++i)
    {
      rShift.push_back( gRandom->Gaus(0,1) );
    }

  }
  return true;
}

bool MnvRecoShifter::GetRandomShiftVectors(const int nUniverses /* = 1000 */ )
{
  bool ok = true;
  // We set the Seeds for every observable variable 
  SetSeeds();
  // Clear the current random shift vectors
  m_rShiftMuon.clear();
  m_rShiftHadron.clear();
  m_rShiftMuonTheta.clear();
  // Fill Muon Shift Vector
  ok = ok && GetRandomShiftVector(m_rShiftMuon, m_seedMuon, nUniverses);
  // Fill Hadron Shift Vector
  ok = ok && GetRandomShiftVector(m_rShiftHadron, m_seedHadron, nUniverses);
  // Fill MuonTheta Shift Vector
  ok = ok && GetRandomShiftVector(m_rShiftMuonTheta, m_seedMuonTheta, nUniverses);
  // Fill Beam XYZ Shift Vector
  //ok = ok && GetRandomShiftVector(m_rShiftBeamDirX, m_seedBeamX, nUniverses);
  //ok = ok && GetRandomShiftVector(m_rShiftBeamDirY, m_seedBeamY, nUniverses);
  //ok = ok && GetRandomShiftVector(m_rShiftBeamDirZ, m_seedBeamZ, nUniverses);
  ok = ok && GetRandomShiftVector(m_rShiftBeamThetaX, m_seedBeamThetaX, nUniverses);
  ok = ok && GetRandomShiftVector(m_rShiftBeamThetaY, m_seedBeamThetaY, nUniverses);

  // Fill Muon XYZ Shift Vector
  ok = ok && GetRandomShiftVector(m_rShiftMuonDirX, m_seedMuonX, nUniverses);
  ok = ok && GetRandomShiftVector(m_rShiftMuonDirY, m_seedMuonY, nUniverses);
  ok = ok && GetRandomShiftVector(m_rShiftMuonDirZ, m_seedMuonZ, nUniverses);

  // Fill Hadron XYZ Shift Vector
  ok = ok && GetRandomShiftVector(m_rShiftHadronDirX, m_seedHadronX, nUniverses);
  ok = ok && GetRandomShiftVector(m_rShiftHadronDirY, m_seedHadronY, nUniverses);
  ok = ok && GetRandomShiftVector(m_rShiftHadronDirZ, m_seedHadronZ, nUniverses);

  return ok;
}

void MnvRecoShifter::SetSeeds()
{
  m_seedMuon = 123456;
  m_seedHadron = 654321;
  m_seedMuonTheta = 13579;

  //m_seedBeamX = 12409;
  //m_seedBeamY = 28751;
  //m_seedBeamZ = 1042;
  m_seedBeamThetaX = 12409;
  m_seedBeamThetaY = 28751;
                 
  m_seedMuonX = 25464;
  m_seedMuonY = 12822;
  m_seedMuonZ = 30426;
                 
  m_seedHadronX = 16985;
  m_seedHadronY = 6999;
  m_seedHadronZ = 5422;

}

void MnvRecoShifter::SetDefaultEnergyLimits()
{
  //! By default there are no limits, make them extreme
  m_lowerLimit_E = -99999999.;
  m_upperLimit_E = +99999999.;

  m_lowerLimit_hadronicE = -99999999.;
  m_upperLimit_hadronicE = +99999999.;

  m_lowerLimit_muonE = -99999999.;
  m_upperLimit_muonE = +99999999.;

  m_lowerLimit_muonTheta = -99999999.;
  m_upperLimit_muonTheta = +99999999.;
}

//============================================================
// END - Get a Vector of Gaus Random Initialized variables
//============================================================

#endif
