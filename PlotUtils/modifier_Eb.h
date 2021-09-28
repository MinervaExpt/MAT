#ifndef modifier_eb_h
#define modifier_eb_h

#include <cmath>
#include <functional>
#include <algorithm>

#include <TMath.h>
#include <TVector3.h>

#include "PlotUtils/ChainWrapper.h"


#include <iostream>
#include <map>
#include <vector>

using namespace std;
using namespace TMath;

namespace PlotUtils{
  class ChainWrapper;

  class KinematicsModifier
  {
    public:
    KinematicsModifier();
    ~KinematicsModifier(){};


    template<class V>
    void CorrectQE(  V &muon, V &nucleon, V& mc_incomingPartVec,V& mc_initNucVec, const int A, const int Z,  int mode = 0, double  Eb=0.025, double Ex = 0.010, double Coulomb = 0.0031);

    template<class T, class V>
    void CorrectQE( T* evt, V& muon, V& nucleon, int mode = 0 );

    template<class V>
    void CorrectQE( PlotUtils::ChainWrapper* chw, Long64_t entry, V& muon, V &nucleon, int mode = 0 );

    template<class T, class V>
    V GetCorrectedMuon( T* evt, V muon, int mode = 0 );
    template<class V>
    V GetCorrectedMuon( PlotUtils::ChainWrapper* chw, Long64_t entry, V muon, int mode = 0 );

    template<class T, class V>
    V GetCorrectedNucleon( T* evt, V nucleon, int mode = 0 );
    template<class V>
    V GetCorrectedNucleon( PlotUtils::ChainWrapper* chw, Long64_t entry, V nucleon, int mode = 0 );

    //energy should be in GeV
    template<class V>
    void ModifyEnergy( V& part, double deltaE=0 );

    double GetGenieBE( const int A, const int Z );
    double GetCoulomb( const int A, const int Z );
    pair<double,double> GetUoptPars( const int A, const int Z );
    double GetUopt( const int A, const int Z, double q3Ksq );//(q3+K)^2 K-->Fermi motion
    double GetExcitationEnergy( const int A, const int Z, bool isNeutrino );


    private:

    vector<int> Ex_exist;
    vector<double> Ex_nu, Ex_antinu;

  };

KinematicsModifier::KinematicsModifier()
{
  Ex_exist=vector<int>({1,3,6,8,13,14,18,20,23,26,28,79,82});
  Ex_nu=vector<double>({0,12.2,10.0, 10.2,21.6,12.4,21.8,19.8,17.0,19.0,16.8,19.5,16.9});
  Ex_antinu=vector<double>({0,12.2,10.1, 10.9,21.6,12.4,17.8,19.4,17.0,19.0,16.8,19.5,14.7});
}


template<class V> 
void KinematicsModifier::ModifyEnergy( V& part, double deltaE )
{
  double E = part.E();
  TVector3 p3( part.X(), part.Y(), part.Z() );
  double p = TMath::Sqrt( p3.X()*p3.X() + p3.Y()*p3.Y() + p3.Z()*p3.Z() );
  if (p==0) return;
  double M = part.M();
  double Eprime = E+deltaE;
  double pprime = TMath::Sqrt(Eprime*Eprime - M*M);
  p3 = p3*(pprime/p);
  part.SetXYZT( p3.X(), p3.Y(), p3.Z(), Eprime );
  return;
}
double KinematicsModifier::GetCoulomb( const int A, const int Z )
{
  double columb = 0;
  if (A == 6)  columb =  1.4; // lithium
  if (A == 12) columb =  3.1; // carbon
  if (A == 16) columb =  3.5; // oxygen
  if (A == 24) columb =  5.1; // magnesium --> use Al27
  if (A == 40) columb =  6.3; // argon
  if (A == 48) columb =  8.1; // Ti48 --> use V50
  if (A == 56) columb =  8.9; // 56 iron
  if (A == 58) columb =  9.8; // 58 nickel
  if (A == 197) columb =  18.5; // 197 gold
  if (A >= 206 && A <= 208) columb =  18.9; // 208 lead

  // Specialty in MINERvA,picked off numbers by hand.
  if (A == 28) columb =  5.5;  // silicon
  if (A == 27) columb =  5.1;  // aluminum
  if (A == 14) columb =  3.3;  // nitrogen --> average between carbon and oxygen
  if (A == 55) columb =  8.9;  // iron55 or manganese55
  if (A == 35) columb =  6.3;  // chlorine  --> use argon
  if (A == 4) columb =  1.4;  // helium? --> use lithium

  return columb/1000.;   // GENIE defaults to something like this, yep no coulomb.
}

double KinematicsModifier::GetGenieBE( const int A, const int Z )
{

  double Eb = 8.0;
  // From UserPhysicsOptions.xml file
  // these are hard coded, so we can get an exact match
  // but beware if you try to port this code to any new GENIE.
  
  if (A == 1) Eb = 0; // hydrogen
  if (A == 6) Eb = 17.0; // lithium
  if (A == 12) Eb = 25.0; // carbon
  if (A == 16) Eb = 27.0; // oxygen
  if (A == 24) Eb = 32.0; // magnesium
  if (A == 40) Eb = 29.5; // argon
  if (A == 48) Eb = 30.0; // Ti48
  if (A == 56) Eb = 36.0; // 56 iron
  if (A == 58) Eb = 36.0; // 58 nickel
  if (A >= 206 && A <= 208) Eb = 44.0; // 208 lead

  // Specialty in MINERvA,picked off numbers by hand.
  if (A == 28) Eb = 8.219751;  // silicon
  if (A == 27) Eb = 8.115287;  // aluminum
  if (A == 14) Eb = 7.185166;  // nitrogen
  if (A == 55) Eb = 8.653063;  // iron55 or manganese55
  if (A == 35) Eb = 8.347164;  // chlorine  

  if (A == 4) Eb = 5.0;  // this is rough, 

  //else
  // this is a problem, all other numbers come from a semi-empirical binding energy formula
  // need to back off the QE precision for those.
  return Eb/1000.;   // GENIE defaults to something like this.  Needs to be exact.  Check it.
}

pair<double,double> KinematicsModifier::GetUoptPars( const int A, const int Z )
{
  if( A <6  ) return  make_pair(0.,0.); //hydrogen
  if( A <12  ) return make_pair(0.0281,-0.0046);//lithium
  if( A ==12 ) return make_pair(0.0421,-0.0296);//carbon and oxygen
  if( A ==16 ) return make_pair(0.0421,-0.0296);
  if( A ==27 ) return make_pair(0.0396,-0.0286);
  if( A ==28 ) return make_pair(0.0396,-0.0286);
  if( A ==40 ) return make_pair(0.0525,-0.0378);
  if( A ==40 ) return make_pair(0.0525,-0.0378);
  if( A ==50 ) return make_pair(0.0544,-0.0332);
  if( A ==56 ) return make_pair(0.0544,-0.0332);
  if( A ==197) return make_pair(0.0873,-0.0352);
  if( A >=206 && A<=208   ) return make_pair(0.0873,-0.0352);
  return make_pair(0.,0.);
}
double KinematicsModifier::GetUopt( const int A, const int Z, double q3Ksq )//(q3+K)^2 K-->Fermi motion
{
  pair<double, double> par = this->GetUoptPars(A,Z);
  return min( par.first*q3Ksq+par.second, 0. );
}


double KinematicsModifier::GetExcitationEnergy( const int A, const int Z, bool isNeutrino )
{
  double Ex = 0;
  auto it = find( this->Ex_exist.begin(), this->Ex_exist.end(), A );
  if( it == this->Ex_exist.end() ) return 0;
  int index = it - this->Ex_exist.begin();
  if( isNeutrino ) Ex = this->Ex_nu[index];
  else Ex = this->Ex_antinu[index];
  return Ex/1000.;
}

template<class V> 
void KinematicsModifier::CorrectQE(  V &muon, V &nucleon, V& mc_incomingPartVec,V& mc_initNucVec, const int A, const int Z,  int mode, double  Eb, double Ex, double Coulomb)
{
  //#Mode 0: no correction
  //#Mode 1: Just change binding energy, GENIE style
  //#Mode 2: UFSI Correction for Default GENIE without Veff
  //#Mode 3: UFSI Correction for Default GENIE with Veff
  //#Mode 4: UFSI Correction for Default GENIE with Veff and Veff from lepton/proton


  if( mode == 0 ) return;

  this->ModifyEnergy( muon, -Ex );
  this->ModifyEnergy( nucleon, Eb );

  if( mode == 1 ) return;

  V q = mc_incomingPartVec - muon;
  V k = mc_initNucVec;
  double q0 = q.E();
  TVector3 q3(q.X(), q.Y(), q.Z() );
  double q3ksq = TMath::Power((k+q).P(),2);
  double Uopt = this->GetUopt( A, Z, q3ksq );

  if( mode == 2 )
  {
    this->ModifyEnergy(nucleon, Uopt);
    this->ModifyEnergy(muon, -Uopt);
  }
  else if( mode == 3 )
  {
    this->ModifyEnergy(nucleon, Uopt+Coulomb);
    this->ModifyEnergy(muon, -Uopt-Coulomb);
  }
  return;
}

template<class T, class V> 
void KinematicsModifier::CorrectQE( T* evt, V& muon, V &nucleon, int mode  )
{
  if( evt->mc_intType != 1 && evt->mc_current !=1 ) return;
  V mc_incomingPartVec(  evt->mc_incomingPartVec[0]/1000., 
                                      evt->mc_incomingPartVec[1]/1000., 
                                      evt->mc_incomingPartVec[2]/1000., 
                                      evt->mc_incomingPartVec[3]/1000. );
  V mc_initNucVec( evt->mc_initNucVec[0]/1000.,
                                evt->mc_initNucVec[1]/1000.,
                                evt->mc_initNucVec[2]/1000.,
                                evt->mc_initNucVec[3]/1000. );
  const int A = evt->mc_targetA;
  const int Z = evt->mc_targetZ;
  bool isNeutrino = ( evt->mc_incoming > 0 );
  double Eb = GetGenieBE(A,Z);
  double Coulomb = GetCoulomb(A,Z);
  if (!isNeutrino) Coulomb*=-1;
  double Ex = GetExcitationEnergy(A,Z,isNeutrino);

  return this->CorrectQE(   muon,  nucleon,  mc_incomingPartVec, mc_initNucVec, A, Z,  mode, Eb, Ex, Coulomb);
}


template<class V> 
void KinematicsModifier::CorrectQE( PlotUtils::ChainWrapper* chw, Long64_t entry, V& muon, V &nucleon, int mode  )
{
  if( chw->GetInt("mc_intType", entry) != 1 && chw->GetInt("mc_current", entry)!=1 ) return;
  vector<double> ipV = chw->GetValueVector<double>("mc_incomingPartVec", entry);
  V mc_incomingPartVec(  ipV[0]/1000., 
                                      ipV[1]/1000., 
                                      ipV[2]/1000., 
                                      ipV[3]/1000. );
  vector<double> inV = chw->GetValueVector<double>("mc_initNucVec", entry);
  V mc_initNucVec( inV[0]/1000.,
                                inV[1]/1000.,
                                inV[2]/1000.,
                                inV[3]/1000. );
  const int A = chw->GetInt("mc_targetA",entry);
  const int Z = chw->GetInt("mc_targetZ",entry);
  bool isNeutrino = (  chw->GetInt("mc_incoming",entry)> 0 );
  double Eb = GetGenieBE(A,Z);
  double Coulomb = GetCoulomb(A,Z);
  if (!isNeutrino) Coulomb*=-1;
  double Ex = GetExcitationEnergy(A,Z,isNeutrino);

  return this->CorrectQE(   muon,  nucleon,  mc_incomingPartVec, mc_initNucVec, A, Z,  mode, Eb, Ex, Coulomb);
}


template<class T, class V> 
V KinematicsModifier::GetCorrectedMuon( T* evt, V muon, int mode  )
{
  if( evt->mc_intType != 1 && evt->mc_current !=1 ) return muon;
  V mc_incomingPartVec(  evt->mc_incomingPartVec[0]/1000., 
                                      evt->mc_incomingPartVec[1]/1000., 
                                      evt->mc_incomingPartVec[2]/1000., 
                                      evt->mc_incomingPartVec[3]/1000. );
  V mc_initNucVec( evt->mc_initNucVec[0]/1000.,
                                evt->mc_initNucVec[1]/1000.,
                                evt->mc_initNucVec[2]/1000.,
                                evt->mc_initNucVec[3]/1000. );
  const int A = evt->mc_targetA;
  const int Z = evt->mc_targetZ;
  bool isNeutrino = ( evt->mc_incoming > 0 );
  double Eb = GetGenieBE(A,Z);
  double Coulomb = GetCoulomb(A,Z);
  if (!isNeutrino) Coulomb*=-1;
  double Ex = GetExcitationEnergy(A,Z,isNeutrino);

  V nucleon(0,0,0,1);
  this->CorrectQE(   muon,  nucleon,  mc_incomingPartVec, mc_initNucVec, A, Z,  mode, Eb, Ex, Coulomb);
  return muon;
}


template<class V> 
V KinematicsModifier::GetCorrectedMuon( PlotUtils::ChainWrapper* chw, Long64_t entry, V muon, int mode  )
{
  if( chw->GetInt("mc_intType", entry) != 1 && chw->GetInt("mc_current", entry)!=1 ) return muon;
  vector<double> ipV = chw->GetValueVector<double>("mc_incomingPartVec", entry);
  V mc_incomingPartVec(  ipV[0]/1000., 
                                      ipV[1]/1000., 
                                      ipV[2]/1000., 
                                      ipV[3]/1000. );
  vector<double> inV = chw->GetValueVector<double>("mc_initNucVec", entry);
  V mc_initNucVec( inV[0]/1000.,
                                inV[1]/1000.,
                                inV[2]/1000.,
                                inV[3]/1000. );
  const int A = chw->GetInt("mc_targetA",entry);
  const int Z = chw->GetInt("mc_targetZ",entry);
  bool isNeutrino = (  chw->GetInt("mc_incoming",entry)> 0 );
  double Eb = GetGenieBE(A,Z);
  double Coulomb = GetCoulomb(A,Z);
  if (!isNeutrino) Coulomb*=-1;
  double Ex = GetExcitationEnergy(A,Z,isNeutrino);

  V nucleon(0,0,0,1);

  this->CorrectQE(   muon,  nucleon,  mc_incomingPartVec, mc_initNucVec, A, Z,  mode, Eb, Ex, Coulomb);
  return muon;
}


template<class T, class V> 
V KinematicsModifier::GetCorrectedNucleon( T* evt, V nucleon, int mode  )
{
  if( evt->mc_intType != 1 && evt->mc_current !=1 ) return nucleon;
  V mc_incomingPartVec(  evt->mc_incomingPartVec[0]/1000., 
                                      evt->mc_incomingPartVec[1]/1000., 
                                      evt->mc_incomingPartVec[2]/1000., 
                                      evt->mc_incomingPartVec[3]/1000. );
  V mc_initNucVec( evt->mc_initNucVec[0]/1000.,
                                evt->mc_initNucVec[1]/1000.,
                                evt->mc_initNucVec[2]/1000.,
                                evt->mc_initNucVec[3]/1000. );
  const int A = evt->mc_targetA;
  const int Z = evt->mc_targetZ;
  bool isNeutrino = ( evt->mc_incoming > 0 );
  double Eb = GetGenieBE(A,Z);
  double Coulomb = GetCoulomb(A,Z);
  if (!isNeutrino) Coulomb*=-1;
  double Ex = GetExcitationEnergy(A,Z,isNeutrino);

  V muon(0,0,0,1);
  this->CorrectQE(   muon,  nucleon,  mc_incomingPartVec, mc_initNucVec, A, Z,  mode, Eb, Ex, Coulomb);
  return nucleon;
}


template<class V>
V KinematicsModifier::GetCorrectedNucleon( PlotUtils::ChainWrapper* chw, Long64_t entry, V nucleon, int mode  )
{
  if( chw->GetInt("mc_intType", entry) != 1 && chw->GetInt("mc_current", entry)!=1 ) return nucleon;
  vector<double> ipV = chw->GetValueVector<double>("mc_incomingPartVec", entry);
  V mc_incomingPartVec(  ipV[0]/1000., 
                                      ipV[1]/1000., 
                                      ipV[2]/1000., 
                                      ipV[3]/1000. );
  vector<double> inV = chw->GetValueVector<double>("mc_initNucVec", entry);
  V mc_initNucVec( inV[0]/1000.,
                                inV[1]/1000.,
                                inV[2]/1000.,
                                inV[3]/1000. );
  const int A = chw->GetInt("mc_targetA",entry);
  const int Z = chw->GetInt("mc_targetZ",entry);
  bool isNeutrino = (  chw->GetInt("mc_incoming",entry)> 0 );
  double Eb = GetGenieBE(A,Z);
  double Coulomb = GetCoulomb(A,Z);
  if (!isNeutrino) Coulomb*=-1;
  double Ex = GetExcitationEnergy(A,Z,isNeutrino);

  V muon(0,0,0,1);

  this->CorrectQE(   muon,  nucleon,  mc_incomingPartVec, mc_initNucVec, A, Z,  mode, Eb, Ex, Coulomb);
  return nucleon;
}







}

















#endif
