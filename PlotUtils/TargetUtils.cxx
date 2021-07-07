#include "TargetUtils.h"

#include <TError.h>
#include <TString.h>
#include <TMath.h>
#include "TSystem.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace PlotUtils;


// singleton
TargetUtils& TargetUtils::Get()
{
  static TargetUtils singleton;
  return singleton;
}


//! Get the area of a hexagon (mm^2)
double TargetUtils::GetHexArea( double apothem /* = 850. */ ) const 
{ 
  return pow(apothem,2) * ( 6. / sqrt(3.) ); 
}

//! Is (x,y) coordinate inside the hexagon
bool TargetUtils::IsInHexagon( double x, double y, double apothem /*= 850. */ ) const
{
   double lenOfSide = apothem*(2/sqrt(3)); 
   double slope     = (lenOfSide/2.0)/apothem;
   double xp        = fabs(x);
   double yp        = fabs(y);
   
   if( (xp*xp + yp*yp) < apothem*apothem )             return true;
   else if( xp <= apothem && yp*yp < lenOfSide/2.0 )   return true; 
   else if( xp <= apothem && yp < lenOfSide-xp*slope ) return true;

   return false;
}

//! Get the pass of a number of planes
double TargetUtils::GetTrackerMass( double nPlanes, bool isMC, double apothem /* = 850. */ ) const
{
  const double area = GetHexArea( apothem );
  const double dens = isMC ? TargetProp::Scint::areal_density_MC : TargetProp::Scint::areal_density;
  const double massOfOnePlane = area * dens;
  return nPlanes * massOfOnePlane;
}

double TargetUtils::GetTrackerMass( double minZ, double maxZ, bool isMC, double apothem /* = 850. */ ) const
{
    double nPlanes = GetNPlanes(minZ, maxZ);
    return GetTrackerMass(nPlanes, isMC, apothem);
}

double TargetUtils::GetTrackerNAtoms( double nPlanes, bool isMC, double apothem /* = 850. */ ) const
{
  const double mass = GetTrackerMass( nPlanes, isMC, apothem );
  const double nAtoms_H = mass * GetTrackerElementMassFraction( 1, isMC ) * TargetProp::AtomsPerGram::H;
  const double nAtoms_C = mass * GetTrackerElementMassFraction( 6, isMC ) * TargetProp::AtomsPerGram::C;
  const double nAtoms_O = mass * GetTrackerElementMassFraction( 8, isMC ) * TargetProp::AtomsPerGram::O;
  const double nAtoms_Al = mass * GetTrackerElementMassFraction( 13, isMC ) * TargetProp::AtomsPerGram::Al;
  const double nAtoms_Si = mass * GetTrackerElementMassFraction( 14, isMC ) * TargetProp::AtomsPerGram::Si;
  const double nAtoms_Cl = mass * GetTrackerElementMassFraction( 17, isMC ) * TargetProp::AtomsPerGram::Cl;
  const double nAtoms_Ti = mass * GetTrackerElementMassFraction( 22, isMC ) * TargetProp::AtomsPerGram::Ti;
  const double nAtoms_N = mass * GetTrackerElementMassFraction( 7, isMC ) * TargetProp::AtomsPerGram::N;

  return nAtoms_C + nAtoms_H + nAtoms_O + nAtoms_Al + nAtoms_Ti + nAtoms_Si + nAtoms_Cl + nAtoms_N;
}

double TargetUtils::GetTrackerNAtoms( double minZ, double maxZ, bool isMC, double apothem /* = 850. */ ) const
{
    int nPlanes = GetNPlanes(minZ, maxZ);
    return GetTrackerNAtoms(nPlanes, isMC, apothem);
}

// Get the number of protons in some planes of tracker
double TargetUtils::GetTrackerNProtons( double nPlanes, bool isMC, double apothem /* = 850. */ )const
{
  const double mass = GetTrackerMass( nPlanes, isMC, apothem );
  const double nProtons_H  = mass * GetTrackerElementMassFraction( 1, isMC ) * TargetProp::AtomsPerGram::H  * TargetProp::ProtonsPerAtom::H;
  const double nProtons_C  = mass * GetTrackerElementMassFraction( 6, isMC ) * TargetProp::AtomsPerGram::C  * TargetProp::ProtonsPerAtom::C;
  const double nProtons_O  = mass * GetTrackerElementMassFraction( 8, isMC ) * TargetProp::AtomsPerGram::O  * TargetProp::ProtonsPerAtom::O;
  const double nProtons_Al  = mass * GetTrackerElementMassFraction( 13, isMC ) * TargetProp::AtomsPerGram::Al  * TargetProp::ProtonsPerAtom::Al;
  const double nProtons_Si  = mass * GetTrackerElementMassFraction( 14, isMC ) * TargetProp::AtomsPerGram::Si  * TargetProp::ProtonsPerAtom::Si;
  const double nProtons_Cl  = mass * GetTrackerElementMassFraction( 17, isMC ) * TargetProp::AtomsPerGram::Cl  * TargetProp::ProtonsPerAtom::Cl;
  const double nProtons_Ti  = mass * GetTrackerElementMassFraction( 22, isMC ) * TargetProp::AtomsPerGram::Ti  * TargetProp::ProtonsPerAtom::Ti;
  const double nProtons_N  = mass * GetTrackerElementMassFraction( 7, isMC ) * TargetProp::AtomsPerGram::N  * TargetProp::ProtonsPerAtom::N;

  return nProtons_C + nProtons_H + nProtons_O + nProtons_Al + nProtons_Ti + nProtons_Si + nProtons_Cl + nProtons_N;
}

double TargetUtils::GetTrackerNProtons( double minZ, double maxZ, bool isMC, double apothem /* = 850. */ )const
{
  int nPlanes = GetNPlanes(minZ, maxZ);
  return GetTrackerNProtons(nPlanes, isMC, apothem);
}

// Get the number of neutrons in some planes of tracker
double TargetUtils::GetTrackerNNeutrons( double nPlanes, bool isMC, double apothem /* = 850. */ )const
{
const double mass = GetTrackerMass( nPlanes, isMC, apothem );
const double nNeutrons_H  = mass * GetTrackerElementMassFraction( 1, isMC ) * TargetProp::AtomsPerGram::H  * TargetProp::NeutronsPerAtom::H;
const double nNeutrons_C  = mass * GetTrackerElementMassFraction( 6, isMC ) * TargetProp::AtomsPerGram::C  * TargetProp::NeutronsPerAtom::C;
const double nNeutrons_O  = mass * GetTrackerElementMassFraction( 8, isMC ) * TargetProp::AtomsPerGram::O  * TargetProp::NeutronsPerAtom::O;
const double nNeutrons_Al  = mass * GetTrackerElementMassFraction( 13, isMC ) * TargetProp::AtomsPerGram::Al  * TargetProp::NeutronsPerAtom::Al;
const double nNeutrons_Si  = mass * GetTrackerElementMassFraction( 14, isMC ) * TargetProp::AtomsPerGram::Si  * TargetProp::NeutronsPerAtom::Si;
const double nNeutrons_Cl  = mass * GetTrackerElementMassFraction( 17, isMC ) * TargetProp::AtomsPerGram::Cl  * TargetProp::NeutronsPerAtom::Cl;
const double nNeutrons_Ti  = mass * GetTrackerElementMassFraction( 22, isMC ) * TargetProp::AtomsPerGram::Ti  * TargetProp::NeutronsPerAtom::Ti;
const double nNeutrons_N  = mass * GetTrackerElementMassFraction( 7, isMC ) * TargetProp::AtomsPerGram::N  * TargetProp::NeutronsPerAtom::N;

return nNeutrons_C + nNeutrons_H + nNeutrons_O + nNeutrons_Al + nNeutrons_Ti + nNeutrons_Si + nNeutrons_Cl + nNeutrons_N;
}

double TargetUtils::GetTrackerNNeutrons( double minZ, double maxZ, bool isMC, double apothem /* = 850. */ )const
{
    int nPlanes = GetNPlanes(minZ, maxZ);
    return GetTrackerNNeutrons(nPlanes, isMC, apothem);
}

// Get the number of nucleons in some planes of tracker
double TargetUtils::GetTrackerNNucleons( double nPlanes, bool isMC, double apothem /* = 850. */ ) const
{
  return GetTrackerNProtons( nPlanes, isMC, apothem ) + GetTrackerNNeutrons( nPlanes, isMC, apothem );
}

double TargetUtils::GetTrackerNNucleons( double minZ, double maxZ, bool isMC, double apothem /* = 850. */ ) const
{
    int nPlanes = GetNPlanes(minZ, maxZ);
    return GetTrackerNNucleons(nPlanes, isMC, apothem);
}

// Get the number of atoms in some planes of tracker for element with Z protons
double TargetUtils::GetTrackerElementNAtoms( int elementZ, double nPlanes, bool isMC, double apothem /* = 850. */ ) const
{
  const double mass = GetTrackerMass( nPlanes, isMC, apothem );
  const double massFrac = GetTrackerElementMassFraction( elementZ, isMC );
  const double atomsPerGram = GetTrackerElementAtomsPerGram( elementZ );

  return mass * massFrac * atomsPerGram;
}

double TargetUtils::GetTrackerElementNAtoms( int elementZ, double minZ, double maxZ, bool isMC, double apothem /* = 850. */ ) const
{
    int nPlanes = GetNPlanes(minZ, maxZ);
    return GetTrackerElementNAtoms(elementZ, nPlanes, isMC, apothem);
}

double TargetUtils::GetTrackerElementMassFraction( int elementZ, bool isMC ) const
{
  switch(elementZ)
  {
    case 1:
      return isMC ? TargetProp::Scint::MassFractionMC::H : TargetProp::Scint::MassFraction::H;
    case 6:
      return isMC ? TargetProp::Scint::MassFractionMC::C : TargetProp::Scint::MassFraction::C;
    case 7:
      return isMC ? TargetProp::Scint::MassFractionMC::N : TargetProp::Scint::MassFraction::N;
    case 8:
      return isMC ? TargetProp::Scint::MassFractionMC::O : TargetProp::Scint::MassFraction::O;
    case 13:
      return isMC ? TargetProp::Scint::MassFractionMC::Al : TargetProp::Scint::MassFraction::Al;
    case 14:
      return isMC ? TargetProp::Scint::MassFractionMC::Si : TargetProp::Scint::MassFraction::Si;
    case 17:
      return isMC ? TargetProp::Scint::MassFractionMC::Cl : TargetProp::Scint::MassFraction::Cl;
    case 22:
      return isMC ? TargetProp::Scint::MassFractionMC::Ti : TargetProp::Scint::MassFraction::Ti;
  }

  Error( "TargetUtils::GetTrackerElementMassFraction", Form( "I don't have any elements with Z = %d in the tracker.", elementZ) );
  throw 1;
}

double TargetUtils::GetTrackerElementA( int elementZ ) const
{
  switch(elementZ)
  {
    case 1:
      return TargetProp::AtomicMass::H - elementZ;
    case 6:
      return TargetProp::AtomicMass::C - elementZ;
    case 7:
      return TargetProp::AtomicMass::N - elementZ;
    case 8:
      return TargetProp::AtomicMass::O - elementZ;
    case 13:
      return TargetProp::AtomicMass::Al - elementZ;
    case 14:
      return TargetProp::AtomicMass::Si - elementZ;
    case 17:
      return TargetProp::AtomicMass::Cl - elementZ;
    case 22:
      return TargetProp::AtomicMass::Ti - elementZ;
  }

  Error( "TargetUtils::GetTrackerElementA", Form( "I don't have any elements with Z = %d in the tracker.", elementZ) );
  throw 1;
}

double TargetUtils::GetTrackerElementAtomsPerGram( int elementZ ) const
{
  switch(elementZ)
  {
    case 1:
      return TargetProp::AtomsPerGram::H;
    case 6:
      return TargetProp::AtomsPerGram::C;
    case 7:
      return TargetProp::AtomsPerGram::N;
    case 8:
      return TargetProp::AtomsPerGram::O;
    case 13:
      return TargetProp::AtomsPerGram::Al;
    case 14:
      return TargetProp::AtomsPerGram::Si;
    case 17:
      return TargetProp::AtomsPerGram::Cl;
    case 22:
      return TargetProp::AtomsPerGram::Ti;
  }

  Error( "TargetUtils::GetTrackerElementAtomsPerGram", Form( "I don't have any elements with Z = %d in the tracker.", elementZ) );
  throw 1;
}


// Get the number of C atoms in some planes of tracker
double TargetUtils::GetTrackerNCarbonAtoms( double nPlanes, bool isMC, double apothem /* = 850. */ ) const
{
  return GetTrackerElementNAtoms( 6, nPlanes, isMC, apothem );
}

double TargetUtils::GetTrackerNCarbonAtoms( double minZ, double maxZ, bool isMC, double apothem /* = 850. */ ) const
{
    int nPlanes = GetNPlanes(minZ, maxZ);
    return GetTrackerNCarbonAtoms(nPlanes, isMC, apothem);
}

//! In tracker?
bool TargetUtils::InTracker( double vtx_x, double vtx_y, double vtx_z, double apothem/* = 850. */ )
{
  if( !InTrackerZ( vtx_z ) ) return false;

  // Is the vertex in the hexagon?  If not, boot out
  if( !IsInHexagon( vtx_x, vtx_y, apothem ) ) return false;
   
  return true; 
}
bool TargetUtils::InTrackerZ( double vtx_z )
{
  //Is it in the tracker volume?
  if( vtx_z < PlotUtils::TargetProp::Tracker::Face ) return false;
  if( PlotUtils::TargetProp::Tracker::Back < vtx_z ) return false;

  return true; 
}

//! In nuclear target region?
bool TargetUtils::InNukeRegion( double vtx_x, double vtx_y, double vtx_z, double apothem/* = 850. */ )
{
  if( !InNukeRegionZ( vtx_z ) ) return false;

  // Is the vertex in the hexagon?  If not, boot out
  if( !IsInHexagon( vtx_x, vtx_y, apothem ) ) return false;

  return true;
}
bool TargetUtils::InNukeRegionZ( double vtx_z )
{
  //Is it in the nuclear target region volume?
  if( vtx_z < PlotUtils::TargetProp::NukeRegion::Face ) return false;
  if( PlotUtils::TargetProp::NukeRegion::Back < vtx_z ) return false;

  return true; 
}
  

//! Get the area of a passive target (mm^2)
double TargetUtils::GetPassiveTargetArea( int targetID, int targetZ, double apothem /* = 850. */ ) const
{
  using namespace std;
  const double hexArea = GetHexArea( apothem );

  if( targetZ < 0 )
    return hexArea;

  //plastic reference targets have no division of materials, so do not subtract the excluded buffer zone
  double distToDivCut = ( targetID < 10 ) ? distToDivCut_ : 0.;

  const int refTarg = targetID % 10;

  if( refTarg == 6) return hexArea;

  if( refTarg == 4 )
  {
    if( targetZ == 82 )
      return hexArea;
  }

  if( refTarg == 1 || refTarg == 2 || refTarg == 5 )
  {
    if( apothem / sqrt(3.) < distToDivCut )
    {
      //note: length of side = apothem * ( 2 / sqrt(3) )
      //      When the distance to division cut is larger than half the length of a side the buffer region is no longer a rectangle.
      //      A more complicated calculation is needed, but probably is not necessary.  If you saw the error and exception you probably read this, so pound out some geometry.

      Error( "TargetUtils::GetPassiveTargetArea", "Buffer area calculation for target 1,2,3 not valid when the distance to division cut is greater than one half of the length of a hexagon side.");
      throw(1);
    }

    const double dmzArea = TargetProp::offset_pb_fe * ( 2 * apothem ); //rectange that goes from center to division line
    const double bufferArea = distToDivCut * ( 2 * apothem );          //rectange cut out of fe/pb because it is too close to division

    if( targetZ == 26 )
      return hexArea / 2. + dmzArea - bufferArea;
    if( targetZ == 82 )
      return hexArea / 2. - dmzArea - bufferArea;
  }

  if( refTarg == 3 )
  {
    const double centerToPoint = 2*apothem/sqrt(3.);

    if( targetZ == 6 )
    {
      // excluded area is a trapezoid, which can be described as a rectangle of two triangles cut out (one for either edge of hexagon)
      // height=distToDivCut, width=2*apothem, triangle={height=distToDivCut, angle opposite height=60deg}
      const double areaOfTrangle   = 0.5 * distToDivCut * ( distToDivCut * tan( TMath::Pi() / 6. ) );
      const double areaOfRectangle = distToDivCut*2*centerToPoint;
      const double bufferArea = areaOfRectangle - 2*areaOfTrangle;

      return hexArea / 2. - bufferArea;
    }
    else if( 26 == targetZ )
    {
      //the excluded area is 2 rhombi
      const double areaOfRhombusA = centerToPoint * distToDivCut; 
      const double areaOfRhombusB = (centerToPoint-distToDivCut) * distToDivCut;
      const double bufferArea     = areaOfRhombusA + areaOfRhombusB;

      return hexArea / 3. - bufferArea;
    }
    else if( targetZ == 82 )
    {
      //the excluded area is 2 trapezoids (A and B)
      const double length1A = centerToPoint;
      const double length2A = centerToPoint - distToDivCut;
      const double length1B = centerToPoint - distToDivCut;
      const double length2B = centerToPoint - distToDivCut - 2*distToDivCut*tan( TMath::Pi() / 6. );
      const double areaOfTrapA = .5 * distToDivCut * ( length1A + length2A );
      const double areaOfTrapB = .5 * distToDivCut * ( length1B + length2B );
      const double bufferArea = areaOfTrapA + areaOfTrapB;

      return hexArea / 6. - bufferArea;
    }
  }

  Error( "GetArea", Form("Unable to GetArea for targetID %d and targetZ %d", targetID, targetZ ) );

  throw 1;
}

//! Get the areal mass of a passive target ( g/mm^2 )
double TargetUtils::GetPassiveTargetArealDensity( int targetID, int targetZ, bool isMC ) const
{
  //! targetID above 10 means faux target, so return the tracker areal mass
  if( targetID > 10 )
    return TargetProp::Scint::areal_density * TargetProp::Faux::n_planes;

  if( isMC )
  {
    if( 1 == targetID && 26 == targetZ ) return TargetProp::DensityMC::Fe * TargetProp::ThicknessMC::Tgt1::Fe;
    if( 1 == targetID && 82 == targetZ ) return TargetProp::DensityMC::Pb * TargetProp::ThicknessMC::Tgt1::Pb;

    if( 2 == targetID && 26 == targetZ ) return TargetProp::DensityMC::Fe * TargetProp::ThicknessMC::Tgt2::Fe;
    if( 2 == targetID && 82 == targetZ ) return TargetProp::DensityMC::Pb * TargetProp::ThicknessMC::Tgt2::Pb;

    if( 3 == targetID &&  6 == targetZ ) return TargetProp::DensityMC::C  * TargetProp::ThicknessMC::Tgt3::C;
    if( 3 == targetID && 26 == targetZ ) return TargetProp::DensityMC::Fe * TargetProp::ThicknessMC::Tgt3::Fe;
    if( 3 == targetID && 82 == targetZ ) return TargetProp::DensityMC::Pb * TargetProp::ThicknessMC::Tgt3::Pb;

    if( 4 == targetID && 82 == targetZ ) return TargetProp::DensityMC::Pb * TargetProp::ThicknessMC::Tgt4::Pb;

    if( 5 == targetID && 26 == targetZ ) return TargetProp::DensityMC::Fe * TargetProp::ThicknessMC::Tgt5::Fe;
    if( 5 == targetID && 82 == targetZ ) return TargetProp::DensityMC::Pb * TargetProp::ThicknessMC::Tgt5::Pb;

    //Special code for water target
    if( 6 == targetID && ( 8 == targetZ || 1 == targetZ ) ) return TargetProp::DensityMC::H2O * TargetProp::WaterTarget::ThicknessMC; 
  }
  else
  {
    if( 1 == targetID && 26 == targetZ ) return TargetProp::Density::Fe * TargetProp::Thickness::Tgt1::Fe;
    if( 1 == targetID && 82 == targetZ ) return TargetProp::Density::Pb * TargetProp::Thickness::Tgt1::Pb;

    if( 2 == targetID && 26 == targetZ ) return TargetProp::Density::Fe * TargetProp::Thickness::Tgt2::Fe;
    if( 2 == targetID && 82 == targetZ ) return TargetProp::Density::Pb * TargetProp::Thickness::Tgt2::Pb;

    if( 3 == targetID &&  6 == targetZ ) return TargetProp::Density::C  * TargetProp::Thickness::Tgt3::C;
    if( 3 == targetID && 26 == targetZ ) return TargetProp::Density::Fe * TargetProp::Thickness::Tgt3::Fe;
    if( 3 == targetID && 82 == targetZ ) return TargetProp::Density::Pb * TargetProp::Thickness::Tgt3::Pb;

    if( 4 == targetID && 82 == targetZ ) return TargetProp::Density::Pb * TargetProp::Thickness::Tgt4::Pb;

    if( 5 == targetID && 26 == targetZ ) return TargetProp::Density::Fe * TargetProp::Thickness::Tgt5::Fe;
    if( 5 == targetID && 82 == targetZ ) return TargetProp::Density::Pb * TargetProp::Thickness::Tgt5::Pb;

    //Special code for water target
    if( 6 == targetID && ( 8 == targetZ || 1 == targetZ ) ) return TargetProp::Density::H2O * TargetProp::WaterTarget::Thickness; 
  }

  Error( "TargetUtils::GetPassiveTargetArealDensity", Form("Could not find areal density for target with ID %d, Z %d, MC? %d", targetID, targetZ, isMC) );
  return 0.;
}

//! Get the mass of a passive target (g)
double TargetUtils::GetPassiveTargetMass( int targetID, int targetZ, bool isMC, double apothem /* = 850. */ ) const
{
  const double arealDens = GetPassiveTargetArealDensity( targetID, targetZ, isMC );
  const double area      = GetPassiveTargetArea( targetID, targetZ, apothem );

  return arealDens*area;
}

// Get the number of protons in a target section
double TargetUtils::GetPassiveTargetNProtons( int targetID, int targetZ, bool isMC, double apothem /* = 850 */ ) const
{
  const double mass = GetPassiveTargetMass( targetID, targetZ, isMC, apothem );

  if( targetID > 10 )
  {
    //! For faux targets, use the tracker calculcators then scale because we may not be using the whole plane's area
    const double fullTrackerMass  = GetTrackerMass( TargetProp::Faux::n_planes, isMC, apothem );
    const double trackerNProt     = GetTrackerNProtons( TargetProp::Faux::n_planes, isMC, apothem );
    const double massFrac = ( fabs(mass) > 1E-6 ) ? mass / fullTrackerMass : 0.;
    return massFrac * trackerNProt;
  }

  if( targetID == 6 ) return mass * TargetProp::AtomsPerGram::H2O * TargetProp::ProtonsPerAtom::H2O; 
  if(  6 == targetZ ) return mass * TargetProp::AtomsPerGram::C   * TargetProp::ProtonsPerAtom::C;
  if( 26 == targetZ ) return mass * TargetProp::AtomsPerGram::Fe  * TargetProp::ProtonsPerAtom::Fe;
  if( 82 == targetZ ) return mass * TargetProp::AtomsPerGram::Pb  * TargetProp::ProtonsPerAtom::Pb;

  Error( "TargetUtils::GetPassiveTargetNProtons", Form("Could not get answer for Z = %d", targetZ) );

  return 0.;
}

// Get the number of neutrons in a target section
double TargetUtils::GetPassiveTargetNNeutrons( int targetID, int targetZ, bool isMC, double apothem /* = 850 */ ) const
{
  const double mass = GetPassiveTargetMass( targetID, targetZ, isMC, apothem );

  if( targetID > 10 )
  {
    //! For faux targets, use the tracker calculcators then scale because we may not be using the whole plane's area
    const double fullTrackerMass  = GetTrackerMass( TargetProp::Faux::n_planes, isMC, apothem );
    const double trackerNNeut     = GetTrackerNNeutrons( TargetProp::Faux::n_planes, isMC, apothem );
    const double massFrac = ( fabs(mass) > 1E-6 ) ? mass / fullTrackerMass : 0.;
    return massFrac * trackerNNeut;
  }

  if( targetID == 6 ) return mass * TargetProp::AtomsPerGram::H2O * TargetProp::NeutronsPerAtom::H2O; 
  if(  6 == targetZ ) return mass * TargetProp::AtomsPerGram::C   * TargetProp::NeutronsPerAtom::C;
  if( 26 == targetZ ) return mass * TargetProp::AtomsPerGram::Fe  * TargetProp::NeutronsPerAtom::Fe;
  if( 82 == targetZ ) return mass * TargetProp::AtomsPerGram::Pb  * TargetProp::NeutronsPerAtom::Pb;

  Error( "TargetUtils::GetPassiveTargetNNeutrons", Form("Could not get answer for Z = %d", targetZ) );

  return 0.;
}

// Get the number of nucleons (p+n) in a passive target section
double TargetUtils::GetPassiveTargetNNucleons( int targetID, int targetZ, bool isMC, double apothem /* = 850. */ ) const
{
  return GetPassiveTargetNProtons( targetID, targetZ, isMC, apothem ) + GetPassiveTargetNNeutrons( targetID, targetZ, isMC, apothem );
}

//! In MC target 1?
bool TargetUtils::InTarget1MC( double vtx_x, double vtx_y, double vtx_z, int mc_nucleiZ, double apothem /* = 850. */, bool excludeBuffer /* = false */ )
{
  if( mc_nucleiZ != 26 && mc_nucleiZ != 82 ) return false;
  if( !InTarget1VolMC( vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return false;
  return true;
}

//! In MC target 2?
bool TargetUtils::InTarget2MC( double vtx_x, double vtx_y, double vtx_z, int mc_nucleiZ, double apothem /* = 850. */, bool excludeBuffer /* = false */ )
{
  if( mc_nucleiZ != 26 && mc_nucleiZ != 82 ) return false;
  if( !InTarget2VolMC( vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return false;
  return true;
}

//! In MC target 3?
bool TargetUtils::InTarget3MC( double vtx_x, double vtx_y, double vtx_z, int mc_nucleiZ, double apothem /* = 850. */, bool excludeBuffer /* = false */ )
{
  if( mc_nucleiZ != 6 && mc_nucleiZ != 26 && mc_nucleiZ != 82 ) return false;
  if( !InTarget3VolMC( vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return false;
  return true;
}

//! In MC target 4?
bool TargetUtils::InTarget4MC( double vtx_x, double vtx_y, double vtx_z, int mc_nucleiZ, double apothem /* = 850. */ )
{
  if( mc_nucleiZ != 82 ) return false;
  if( !InTarget4VolMC( vtx_x, vtx_y, vtx_z, apothem ) ) return false;
  return true;
}

//! In MC target 5?
bool TargetUtils::InTarget5MC( double vtx_x, double vtx_y, double vtx_z, int mc_nucleiZ, double apothem /* = 850. */, bool excludeBuffer /* = false */ )
{
  if( mc_nucleiZ != 26 && mc_nucleiZ != 82 ) return false;
  if( !InTarget5VolMC( vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return false;
  return true;
}

//! In MC water target?
bool TargetUtils::InWaterTargetMC( double vtx_x, double vtx_y, double vtx_z, int mc_nucleiZ, double apothem /* = 850. */ )
{
  if( mc_nucleiZ != 8 && mc_nucleiZ != 1 ) return false;
  if( !InWaterTargetVolMC( vtx_x, vtx_y, vtx_z, apothem ) ) return false;
  return true;
}

bool TargetUtils::InPassiveTargetMC( double vtx_x, double vtx_y, double vtx_z, int mc_nucleiZ, double apothem /* = 850. */, bool excludeBuffer /* = false */ )
{
  if( InTarget1MC( vtx_x, vtx_y, vtx_z, mc_nucleiZ, apothem, excludeBuffer ) ) return true;
  if( InTarget2MC( vtx_x, vtx_y, vtx_z, mc_nucleiZ, apothem, excludeBuffer ) ) return true;
  if( InTarget3MC( vtx_x, vtx_y, vtx_z, mc_nucleiZ, apothem, excludeBuffer ) ) return true;
  if( InTarget4MC( vtx_x, vtx_y, vtx_z, mc_nucleiZ, apothem )                ) return true;
  if( InTarget5MC( vtx_x, vtx_y, vtx_z, mc_nucleiZ, apothem, excludeBuffer ) ) return true;
  if( InWaterTargetMC( vtx_x, vtx_y, vtx_z, mc_nucleiZ, apothem )            ) return true;
  return false;
}

bool TargetUtils::InCarbonTargetVolMC( double vtx_x, double vtx_y, double vtx_z, double apothem, bool excludeBuffer )
{
  if( InCarbon3VolMC( vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return true;
  return false;
}

bool TargetUtils::InIronTargetVolMC( double vtx_x, double vtx_y, double vtx_z, double apothem, bool excludeBuffer )
{
  if( InIron1VolMC( vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return true;
  if( InIron2VolMC( vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return true;
  if( InIron3VolMC( vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return true;
  if( InIron5VolMC( vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return true;
  return false;
}

bool TargetUtils::InLeadTargetVolMC( double vtx_x, double vtx_y, double vtx_z, double apothem, bool excludeBuffer )
{
  if( InLead1VolMC( vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return true;
  if( InLead2VolMC( vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return true;
  if( InLead3VolMC( vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return true;
  if( InLead4VolMC( vtx_x, vtx_y, vtx_z, apothem ) ) return true;
  if( InLead5VolMC( vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return true;
  return false;
}

//! In the correct MC volume?
bool TargetUtils::InTarget1VolMC( double vtx_x, double vtx_y, double vtx_z, double apothem, bool excludeBuffer )
{
  if( InIron1VolMC(   vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return true;
  if( InLead1VolMC(   vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return true;
  return false;
}

bool TargetUtils::InTarget2VolMC( double vtx_x, double vtx_y, double vtx_z, double apothem, bool excludeBuffer )
{
  if( InIron2VolMC(   vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return true;
  if( InLead2VolMC(   vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return true;
  return false;
}

bool TargetUtils::InTarget3VolMC( double vtx_x, double vtx_y, double vtx_z, double apothem, bool excludeBuffer )
{
  if( InCarbon3VolMC( vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return true;
  if( InIron3VolMC(   vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return true;
  if( InLead3VolMC(   vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return true;
  return false;
}

bool TargetUtils::InTarget4VolMC( double vtx_x, double vtx_y, double vtx_z, double apothem )
{
  if( InLead4VolMC(   vtx_x, vtx_y, vtx_z, apothem ) ) return true;
  return false;
}

bool TargetUtils::InTarget5VolMC( double vtx_x, double vtx_y, double vtx_z, double apothem, bool excludeBuffer )
{
  if( InIron5VolMC(   vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return true;
  if( InLead5VolMC(   vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return true;
  return false;
}

bool TargetUtils::InWaterTargetVolMC( double vtx_x, double vtx_y, double vtx_z, double apothem )
{
  //Is this in the correct z range?
  if( !InWaterTargetZMC( vtx_z ) ) return false;

  //Is this in the hexagon?
  if ( !IsInHexagon(vtx_x, vtx_y, apothem) ) return false;

  return true; 
}

bool TargetUtils::InPassiveTargetVolMC( double vtx_x, double vtx_y, double vtx_z, double apothem /* = 850. */, bool excludeBuffer /* = false */ )
{
  if( InTarget1VolMC( vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return true;
  if( InTarget2VolMC( vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return true;
  if( InTarget3VolMC( vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return true;
  if( InTarget4VolMC( vtx_x, vtx_y, vtx_z, apothem )                ) return true;
  if( InTarget5VolMC( vtx_x, vtx_y, vtx_z, apothem, excludeBuffer ) ) return true;
  if( InWaterTargetVolMC( vtx_x, vtx_y, vtx_z, apothem )            ) return true;
  return false;
}

//! In the correct z range?
bool TargetUtils::InTarget1ZMC( double vtx_z, int nucleiZ )
{
  // Need to know which target section you're looking at
  // Iron (26) or lead (82)
  if( !OnIron( nucleiZ ) && !OnLead( nucleiZ ) ) return false;

  double width = 0;  
  //Get the width of target  
  if( OnIron( nucleiZ ) ) width = PlotUtils::TargetProp::ThicknessMC::Tgt1::Fe/2;
  if( OnLead( nucleiZ ) ) width = PlotUtils::TargetProp::ThicknessMC::Tgt1::Pb/2;

  // Is the vertex z in the target?
  if( fabs( vtx_z - GetTarget1CenterZMC() ) <= width ) return true;
  return false; 
}

bool TargetUtils::InTarget2ZMC( double vtx_z, int nucleiZ )
{
  // Need to know which target section you're looking at
  // Iron (26) or lead (82)
  if( !OnIron( nucleiZ ) && !OnLead( nucleiZ ) ) return false;

  double width = 0;  
  //Get the width of target  
  if( OnIron( nucleiZ ) ) width = PlotUtils::TargetProp::ThicknessMC::Tgt2::Fe/2;
  if( OnLead( nucleiZ ) ) width = PlotUtils::TargetProp::ThicknessMC::Tgt2::Pb/2;

  // Is the vertex z in the target?
  if( fabs( vtx_z - GetTarget2CenterZMC() ) <= width ) return true;
  return false; 
}

bool TargetUtils::InTarget3ZMC( double vtx_z, int nucleiZ )
{
  // Need to know which target section you're looking at
  // Carbon (6), iron (26) or lead (82)
  if( !OnCarbon( nucleiZ ) && !OnIron( nucleiZ ) && !OnLead( nucleiZ ) ) return false;

  double width = 0;  
  //Get the width of target  
  if( OnCarbon( nucleiZ) )  width = PlotUtils::TargetProp::ThicknessMC::Tgt3::C/2;
  if( OnIron( nucleiZ ) )   width = PlotUtils::TargetProp::ThicknessMC::Tgt3::Fe/2;
  if( OnLead( nucleiZ ) )   width = PlotUtils::TargetProp::ThicknessMC::Tgt3::Pb/2;

  // Is the vertex z in the target?
  if( ( nucleiZ != 6 ) && fabs( vtx_z - GetTarget3CenterZMC() ) <= width ) return true;
  else if( ( nucleiZ == 6 ) && fabs( vtx_z - GetTarget3CarbonCenterZMC() ) <= width ) return true;
  return false; 
}

bool TargetUtils::InTarget4ZMC( double vtx_z, int nucleiZ )
{
  // Need to know which target section you're looking at
  // Lead (82)
  if( !OnLead( nucleiZ ) ) return false;

  double width = 0;  
  //Get the width of target  
  if( OnLead( nucleiZ ) ) width = PlotUtils::TargetProp::ThicknessMC::Tgt4::Pb/2;

  // Is the vertex z in the target?
  if( fabs( vtx_z - GetTarget4CenterZMC() ) <= width ) return true;
  return false; 
}

bool TargetUtils::InTarget5ZMC( double vtx_z, int nucleiZ )
{
  // Need to know which target section you're looking at
  // Iron (26) or lead (82)
  if( !OnIron( nucleiZ ) && !OnLead( nucleiZ ) ) return false;

  double width = 0;  
  //Get the width of target  
  if( OnIron( nucleiZ ) ) width = PlotUtils::TargetProp::ThicknessMC::Tgt5::Fe/2;
  if( OnLead( nucleiZ ) ) width = PlotUtils::TargetProp::ThicknessMC::Tgt5::Pb/2;

  // Is the vertex z in the target?
  if( fabs( vtx_z - GetTarget5CenterZMC() ) <= width ) return true;
  return false; 
}

bool TargetUtils::InWaterTargetZMC( double vtx_z )
{
  if( vtx_z < PlotUtils::TargetProp::WaterTarget::Face ) return false;
  if( PlotUtils::TargetProp::WaterTarget::Back < vtx_z ) return false;

  return true;
}

bool TargetUtils::InIron1VolMC(double vtx_x, double vtx_y, double vtx_z, double apothem, bool excludeBuffer/*default=false*/)
{
   //Is this in the correct z range?
   if ( !InTarget1ZMC( vtx_z, 26 ) ) return false;

   //Is this in the hexagon?
   if ( !IsInHexagon(vtx_x, vtx_y, apothem) ) return false;

   double u = GetCoordU(vtx_x, vtx_y);
   double udist = u - PlotUtils::TargetProp::offset_pb_fe;

   //Is this on the right side on the divide?
   //  If we want to exclude the buffer, is it outside the buffer?  
   return ( udist < 0.0 && ( excludeBuffer ? fabs(udist) > distToDivCut_ : true ) );
}

bool TargetUtils::InLead1VolMC(double vtx_x, double vtx_y, double vtx_z, double apothem/*default=850.0mm*/, bool excludeBuffer/*default=false*/)
{
   //Is this in the correct z range?
   if ( !InTarget1ZMC( vtx_z, 82 ) ) return false;

   //Is this in the hexagon?
   if ( !IsInHexagon(vtx_x, vtx_y, apothem) ) return false;

   double u = GetCoordU(vtx_x, vtx_y);
   double udist = u - PlotUtils::TargetProp::offset_pb_fe;
   //Is this on the right side on the divide?
   //  If we want to exclude the buffer, is it outside the buffer?  
   return ( udist >= 0.0 && ( excludeBuffer ? fabs(udist) > distToDivCut_ : true ) );
}

bool TargetUtils::InIron2VolMC(double vtx_x, double vtx_y, double vtx_z, double apothem, bool excludeBuffer)
{
   //Is this in the correct z range?
   if ( !InTarget2ZMC( vtx_z, 26 ) ) return false;

   //Is this in the hexagon?
   if ( !IsInHexagon(vtx_x, vtx_y, apothem) ) return false;

   double d = GetCoordD(vtx_x, vtx_y);
   double ddist = d - PlotUtils::TargetProp::offset_pb_fe;

   //Is this on the right side on the divide?
   //  If we want to exclude the buffer, is it outside the buffer?  
   return ( ddist < 0.0 && ( excludeBuffer ? fabs(ddist) > distToDivCut_ : true ) );
}

bool TargetUtils::InLead2VolMC(double vtx_x, double vtx_y, double vtx_z, double apothem, bool excludeBuffer)
{
   //Is this in the correct z range?
   if ( !InTarget2ZMC( vtx_z, 82 ) ) return false;

   //Is this in the hexagon?
   if ( !IsInHexagon(vtx_x, vtx_y, apothem) ) return false;

   double d = GetCoordD(vtx_x, vtx_y);
   double ddist = d - PlotUtils::TargetProp::offset_pb_fe;

   //Is this on the right side on the divide?
   //  If we want to exclude the buffer, is it outside the buffer?  
   return ( ddist >= 0.0 && ( excludeBuffer ? fabs(ddist) > distToDivCut_ : true ) );
}

bool TargetUtils::InCarbon3VolMC(double vtx_x, double vtx_y, double vtx_z, double apothem, bool excludeBuffer)
{
   //Is this in the correct z range?
   if ( !InTarget3ZMC( vtx_z, 6 ) ) return false;

   //Is this in the hexagon?
   if ( !IsInHexagon(vtx_x, vtx_y, apothem) ) return false;

   double c = GetCoordC(vtx_x, vtx_y);
   double cdist = c;

   //Is this on the right side on the divide?
   //  If we want to exclude the buffer, is it outside the buffer?  
   return ( cdist >= 0.0 && ( excludeBuffer ? fabs(cdist) > distToDivCut_ : true ) );
}

bool TargetUtils::InIron3VolMC(double vtx_x, double vtx_y, double vtx_z, double apothem, bool excludeBuffer)
{
   //Is this in the correct z range?
   if ( !InTarget3ZMC( vtx_z, 26 ) ) return false;

   //Is this in the hexagon?
   if ( !IsInHexagon(vtx_x, vtx_y, apothem) ) return false;

   double c = GetCoordC(vtx_x, vtx_y);
   double cdist = c;

   //Is this on the right side on the divide?
   //  If we want to exclude the buffer, is it outside the buffer?  
   return ( cdist < 0.0 && vtx_x < 0.0 && ( excludeBuffer ? fabs(cdist) > distToDivCut_ && fabs(vtx_x) > distToDivCut_ : true ) );
}

bool TargetUtils::InLead3VolMC(double vtx_x, double vtx_y, double vtx_z, double apothem, bool excludeBuffer)
{
   //Is this in the correct z range?
   if ( !InTarget3ZMC( vtx_z, 82 ) ) return false;

   //Is this in the hexagon?
   if ( !IsInHexagon(vtx_x, vtx_y, apothem) ) return false;

   double c = GetCoordC(vtx_x, vtx_y);
   double cdist = c;

   //Is this on the right side on the divide?
   //  If we want to exclude the buffer, is it outside the buffer?  
   return ( cdist < 0.0 && vtx_x >= 0.0 && ( excludeBuffer ? fabs(cdist) > distToDivCut_ && fabs(vtx_x) > distToDivCut_ : true ) );
}

bool TargetUtils::InLead4VolMC(double vtx_x, double vtx_y, double vtx_z, double apothem)
{
   //Is this in the correct z range?
   if ( !InTarget4ZMC( vtx_z, 82 ) ) return false;

   //Is this in the hexagon?
   if ( !IsInHexagon(vtx_x, vtx_y, apothem) ) return false;

   return true;
}

bool TargetUtils::InIron5VolMC(double vtx_x, double vtx_y, double vtx_z, double apothem, bool excludeBuffer)
{
   //Is this in the correct z range?
   if ( !InTarget5ZMC( vtx_z, 26 ) ) return false;

   //Is this in the hexagon?
   if ( !IsInHexagon(vtx_x, vtx_y, apothem) ) return false;

   double u = GetCoordU(vtx_x, vtx_y);
   double udist = u - PlotUtils::TargetProp::offset_pb_fe;

   //Is this on the right side on the divide?
   //  If we want to exclude the buffer, is it outside the buffer?  
   return ( udist < 0.0 && ( excludeBuffer ? fabs(udist) > distToDivCut_ : true ) );
}

bool TargetUtils::InLead5VolMC(double vtx_x, double vtx_y, double vtx_z, double apothem, bool excludeBuffer)
{
   //Is this in the correct z range?
   if ( !InTarget5ZMC( vtx_z, 82 ) ) return false;

   //Is this in the hexagon?
   if ( !IsInHexagon(vtx_x, vtx_y, apothem) ) return false;

   double u = GetCoordU(vtx_x, vtx_y);
   double udist = u - PlotUtils::TargetProp::offset_pb_fe;

   //Is this on the right side on the divide?
   //  If we want to exclude the buffer, is it outside the buffer?  
   return ( udist >= 0.0 && ( excludeBuffer ? fabs(udist) > distToDivCut_ : true ) );
}

double TargetUtils::GetCoordU(double vtx_x, double vtx_y)
{//see equation 5 of docdb 9069
  return -vtx_x*cos(TMath::Pi()/6) + vtx_y*sin(TMath::Pi()/6);
} 

double TargetUtils::GetCoordD(double vtx_x, double vtx_y)
{
  return vtx_x*cos(TMath::Pi()/6) + vtx_y*sin(TMath::Pi()/6);
} 

double TargetUtils::GetCoordC(double vtx_x, double vtx_y)
{      
  return vtx_x*sin(TMath::Pi()/6) + vtx_y*cos(TMath::Pi()/6);
}

double TargetUtils::GetNPlanes(double minZ, double maxZ) const
{
    const char* plotutils = gSystem->Getenv("PLOTUTILSROOT");
    if(!plotutils || !strlen(plotutils)){
        std::cout<<"$PLOTUTILSROOT is not set. Can't find flux histograms"<<std::endl;
        std::exit(1);
    }

    std::string plotutils_dir(plotutils);
    std::string file_dir = plotutils_dir + "/data/Minerva_Planes.txt"; 
    std::ifstream file_planes;
    file_planes.open(file_dir.c_str());

    if (!file_planes.is_open()){
        std::cerr<<"WARNING! File cannot be opened"<<std::endl;
        exit(1);
    }

    std::string line;
    double dummy;
    double temp_plane_z;
    std::vector<double> plane_z;
    while(!file_planes.eof()){
        getline(file_planes,line);
        std::stringstream line_ss(line);
        line_ss>>dummy>>dummy>>dummy>>dummy>>dummy>>temp_plane_z;  
       
        plane_z.push_back(temp_plane_z);
    }

    double nplanes = 0;
    for (unsigned int i = 0; i < plane_z.size(); ++i){
        if (plane_z[i] < minZ) continue;
        else if (plane_z[i] > maxZ) break;

        nplanes++;
    }
    
    file_planes.close();

    return nplanes;
}
