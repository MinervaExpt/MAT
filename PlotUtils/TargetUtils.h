#ifndef MNV_TargetUtils
#define MNV_TargetUtils 1

namespace PlotUtils
{
  /*! Define target properties
    See doc 6016 for measurements

    For MC numbers, see tagged version v10r4p1 of MinervaDDDB, specifically
  http://cdcvs.fnal.gov/cgi-bin/public-cvs/cvsweb-public.cgi/AnalysisFramework/Det/MinervaDDDB/DDDB/materials/MinervaMaterials.xml?rev=1.11&content-type=text/x-cvsweb-markup&cvsroot=mnvsoft
  http://cdcvs.fnal.gov/cgi-bin/public-cvs/cvsweb-public.cgi/AnalysisFramework/Det/MinervaDDDB/DDDB/MINERVA/geometry.xml?rev=1.18&content-type=text/x-cvsweb-markup&cvsroot=mnvsoft
  http://cdcvs.fnal.gov/cgi-bin/public-cvs/cvsweb-public.cgi/AnalysisFramework/Det/MinervaDDDB/DDDB/Plane/plane_geometry.xml?rev=1.8&content-type=text/x-cvsweb-markup&cvsroot=mnvsoft

    Densities are all areal densities in g/mm^2
    Thicknesses are in mm
    Masses are in g
   */
  namespace TargetProp
  {
    //! How far from center is the divide between lead and iron in passive targets 1,2,5?
    const double offset_pb_fe = 205.0; //mm
    const double offset_c_pbfe= 0.0; //mm

    const double gcm2_to_gmm2 = 1. / ( 10.*10. ); //convert g/cm^2 to g/mm^2
    const double gcm3_to_gmm3 = 1. / ( 10.*10.*10. ); //convert g/cm^3 to g/mm^3
    const double avagadro = 6.0221412927E23;

    //! What are some scintillator properties?
    namespace Scint
    {
      const double areal_density = 1.9872 * gcm2_to_gmm2;
      const double areal_density_err = 0.028 * gcm2_to_gmm2;

      const double areal_density_MC = 1.9872 * gcm2_to_gmm2; //From NX MC 

      //composition by mass from docdb 11686
      // Ron's estimate (doc 6016, updated after destructive plane measurement)
      namespace MassFraction
      {
        const double H  = 0.0818;
        const double C  = 0.8851;
        const double O  = 0.0250;
        const double Ti = 0.0047;
        //const double O  = 0.0268;
        //const double Ti = 0.0070;
        // Correcting Ti and O fraction, see
        // DocDB 28673
        const double Al = 0.0007;
        const double Si = 0.0007;
        const double Cl = 0.0020;
        const double N = 0.0000;
      }
      //Eroica MC fractions found with MinervaMaterialScanAlg
      namespace EroicaMassFractionMC
      {
        const double C  = 0.8784;
        const double H  = 0.0819;
        const double O  = 0.0286;
        const double Ti = 0.0046;
        const double Al = 0.0013;
        const double Si = 0.0014;
        const double Cl = 0.0038;
        const double N  = 0.0000;
      }
      // NX MC fraction (note they are the same as data
      // Correction by Tejin -- MC fraction different from data
      // Refer to DocDB 28673
      namespace NXMassFractionMC
      {
        const double C  = 0.8896;
        const double H  = 0.07533; 
        const double O  = 0.02432;
        const double Ti = 0.00697;
        const double Al = 0.001613;
        const double Si = 0.001613;
        const double Cl = 0.000062;
        const double N  = 0.000651;
      }
      namespace MassFractionMC
      {
        const double C  = NXMassFractionMC::C ;
        const double H  = NXMassFractionMC::H ; 
        const double O  = NXMassFractionMC::O ;
        const double Ti = NXMassFractionMC::Ti;
        const double Al = NXMassFractionMC::Al;
        const double Si = NXMassFractionMC::Si;
        const double Cl = NXMassFractionMC::Cl;
        const double N  = NXMassFractionMC::N ;
      }
      //Composition by total number of atoms 
      //Derived from quantities above
      // number fraction = Number of atoms/ total number of atoms
      // number of atoms per gram scint = ( Mass fraction/Atomic Mass )*avogadro
      namespace NumAtomFraction
      {
        const double C  = 0.470530;
        const double H  = 0.518181;
        const double O  = 0.009977;
        const double Ti = 0.000627;
        const double Al = 0.000166;
        const double Si = 0.000159;
        const double Cl = 0.000360;
        const double N  = 0.000000;
      }
      namespace EroicaNumAtomFractionMC
      {
        const double C  = 0.467379;
        const double H  = 0.519272;
        const double O  = 0.011424;
        const double Ti = 0.000614;
        const double Al = 0.000308;
        const double Si = 0.000319;
        const double Cl = 0.000685;
        const double N  = 0.000000;
      }
      namespace NXNumAtomFractionMC
      {
        const double C  = 0.491700;
        const double H  = 0.496144;
        const double O  = 0.010091;
        const double Ti = 0.000967;
        const double Al = 0.000397;
        const double Si = 0.000381;
        const double Cl = 0.000012;
        const double N  = 0.000309;
      }
      namespace NumAtomFractionMC
      {
        const double C  = NXNumAtomFractionMC::C ;
        const double H  = NXNumAtomFractionMC::H ; 
        const double O  = NXNumAtomFractionMC::O ;
        const double Ti = NXNumAtomFractionMC::Ti;
        const double Al = NXNumAtomFractionMC::Al;
        const double Si = NXNumAtomFractionMC::Si;
        const double Cl = NXNumAtomFractionMC::Cl;
        const double N  = NXNumAtomFractionMC::N ;
      }
    }

    namespace AtomicMass
    {
      const double C  = 12.0107;
      const double H  = 1.00794;
      const double O  = 15.9994;
      const double Ti = 47.867;
      const double Al = 26.982;
      const double Si = 28.0855;
      const double Cl = 35.453;

      const double Fe = (.99*55.845 + .01*54.938045); //99% Fe, 1% Mn
      const double Pb = 207.2;
      const double N = 14.007;
    }

    //How many atoms in on gram
    namespace AtomsPerGram
    {
      const double C  = avagadro / AtomicMass::C;
      const double H  = avagadro / AtomicMass::H;
      const double O  = avagadro / AtomicMass::O;
      const double Ti = avagadro / AtomicMass::Ti;
      const double Al = avagadro / AtomicMass::Al;
      const double Si = avagadro / AtomicMass::Si;
      const double Cl = avagadro / AtomicMass::Cl;
      const double N  = avagadro / AtomicMass::N;

      const double Fe = avagadro / AtomicMass::Fe;
      const double Pb = avagadro / AtomicMass::Pb;
      const double H2O = avagadro / ( AtomicMass::O + 2*AtomicMass::H );
    }

    namespace ProtonsPerAtom
    {
      const double C  = 6;
      const double H  = 1;
      const double O  = 8;
      const double Ti = 22;
      const double Al = 13;
      const double Si = 14;
      const double Cl = 17;

      const double N = 7;

      const double Fe = (.99*26 + .01*25); //99% Fe, 1%Mn
      const double Pb = 82;
      const double H2O = 2*H + O;
    }

    //derived as atomic mass - nProtons
    namespace NeutronsPerAtom
    {
      const double C  = AtomicMass::C - ProtonsPerAtom::C;
      const double H  = AtomicMass::H - ProtonsPerAtom::H;
      const double O  = AtomicMass::O - ProtonsPerAtom::O;
      const double Al = AtomicMass::Al - ProtonsPerAtom::Al;
      const double Ti = AtomicMass::Ti - ProtonsPerAtom::Ti;
      const double Si = AtomicMass::Si - ProtonsPerAtom::Si;
      const double Cl = AtomicMass::Cl - ProtonsPerAtom::Cl;
      const double N = AtomicMass::N - ProtonsPerAtom::N;

      const double Fe = AtomicMass::Fe - ProtonsPerAtom::Fe;
      const double Pb = AtomicMass::Pb - ProtonsPerAtom::Pb;
      const double H2O = ( AtomicMass::O + 2*AtomicMass::H ) - ProtonsPerAtom::H2O;
    }


    //! What is the measured density of passive targets?
    namespace Density
    {
      const double C  = 1.739 * gcm3_to_gmm3;
      const double Fe = 7.834 * gcm3_to_gmm3;
      const double Pb = 11.29 * gcm3_to_gmm3;
      const double C_err  = 0.01 * gcm3_to_gmm3;
      const double Fe_err = .034 * gcm3_to_gmm3;
      const double Pb_err = 0.03 * gcm3_to_gmm3;

      //Mass for the water in the water target for a fid vol with transverse radius of 90 cm is
      //  115.1945 gal +/- 2% -> 436.1 kg +/- 8.7 (H. Budd)
      //Mass for the water in the water target with a hexagon apothem 
      //  of 85 cm is 113.65 gal +/- 2% -> 430.2 kg +/- 8.6 (H. Budd)
      //The volume of the water target is difficult to quantify (lenticular shape)
      //Assuming a hexagonal prism  of 22 cm thickness and 85 cm apothem
      const double H2O = 0.7813 * gcm3_to_gmm3;
      const double H2O_err = 0.016 * gcm3_to_gmm3;
    }

    //! Idealized densities used in MC for passive targets
    namespace DensityMC
    {
      const double C   = 1.739 * gcm3_to_gmm3;
      const double Fe  = 7.87  * gcm3_to_gmm3;
      const double Pb  = 11.35 * gcm3_to_gmm3;

      //Mass of the water in the MC is 459.12 kg for a transverse radius of 90 cm (J. Kleykamp)
      //Using the ratio of hexagon apothem mass to water volume in 90 cm radius leads to 
      //  452.96 kg of water in a volume with hexagon apothem of 85 cm
      //The volume of the water target is difficult to quantify (lenticular shape)
      //Assuming a hexagonal prism of 22 cm thickness and 85 cm apothem
      const double H2O = 0.822  * gcm3_to_gmm3;
    }

    //! What is the measured thicknes of passive targets?
    namespace Thickness
    {
      namespace Tgt1
      {
        const double Fe = 25.67;
        const double Pb = 25.78;
        const double Fe_err = .06;
        const double Pb_err = 0.12;
      }
      namespace Tgt2
      {
        const double Fe = 25.63;
        const double Pb = 25.81;
        const double Fe_err = 0.06;
        const double Pb_err = 0.16;
      }
      namespace Tgt3
      {
        const double C  = 76.20;
        const double Fe = 25.73;
        const double Pb = 25.63;
        const double C_err  = 0.05;
        const double Fe_err = 0.04;
        const double Pb_err = 0.04;
      }
      namespace Tgt4
      {
        const double Pb = 7.95;
        const double Pb_err = .005;
      }
      namespace Tgt5
      {
        const double Fe = 12.89;
        const double Pb = 13.17;
        const double Fe_err = .06;
        const double Pb_err = .07;
      }
    } //end of measured thicknesses

    //! Idealized thicknesses used in the MC for passive targets
    namespace ThicknessMC
    {
      namespace Tgt1
      {
        const double Fe = 25.75;
        const double Pb = 25.75;
      }
      namespace Tgt2
      {
        const double Fe = 25.75;
        const double Pb = 25.75;
      }
      namespace Tgt3
      {
        const double C  = 76.3;
        const double Fe = 25.75;
        const double Pb = 25.75;
      }
      namespace Tgt4
      {
        const double Pb = 8.;
      }
      namespace Tgt5
      {
        const double Fe = 13.;
        const double Pb = 13.;
      }
    }

    //! Center of targets used in the MC for the passive targets (mm)
    //! Offset is due to studies showing MC targets are not where we originally thought (docdb 11497)
    //! Also see Det/MinervaDDDB/DDDB/MINERVA/geometry.xml
    namespace CenterZMC
    {
      namespace Tgt1
      {
        const double center_nooffset = 4478.016;
        const double offset          = 4.3; 
      }
      namespace Tgt2
      {
        const double center_nooffset = 4698.824;
        const double offset          = 3.5;
      }
      namespace Tgt3
      {
        const double center_nooffset = 4915.595;
        const double offset          = 5.1;
      }
      namespace Tgt3C
      {
        const double center_nooffset = 4940.82;
        const double offset          = 5.1;
      }
      namespace Tgt4
      {
        const double center_nooffset = 5641.68;
        const double offset          = 0.0;
      }
      namespace Tgt5
      {
        const double center_nooffset = 5774.32;
        const double offset          = 2.0;
      }
      
    }

    //! Faux target properties for nuclear target analyses (6 modules)
    namespace Faux
    {
      const double n_planes = 12.0;
    }

    namespace Tracker
    {
      const double Face = 5991.29;
      const double Back = 8408.91;
    }
    namespace WaterTarget
    {
      const double Face = 5200.0;
      const double Back = 5420.0;
      const double Thickness   = Back - Face; 
      const double ThicknessMC = Back - Face; 
    }
    namespace NukeRegion
    {
      const double Face = 4293.04;
      const double Back = 5835.0;
    }
  } //end TargetProp namespace


  //==========================================
  // TargetUtils class
  //==========================================
  //! Class to do mass and target number calcuclations
  class TargetUtils
  {
    public:
      //! Default constructor
      TargetUtils() :
        distToDivCut_(25.),
        useCenterZMCOffset_(true)
    {};

      //! Default destructor
      virtual ~TargetUtils() {};

      //! singleton gettor
      static TargetUtils& Get();

      //! Get the area of a hexagon (mm^2)
      double GetHexArea( double apothem = 850. ) const;
      //! Is (x,y) coordinate inside the hexagon
      bool IsInHexagon( double x, double y, double apothem = 850. ) const;

      //==========================================
      // Functions for Tracker
      //==========================================
      //! Get the pass of a number of planes of tracker (g)
      double GetTrackerMass( double nPlanes, bool isMC, double apothem = 850. ) const;
      double GetTrackerMass( double minZ, double maxZ, bool isMC, double apothem = 850. ) const;
      //! Get the number of atoms in tracker
      double GetTrackerNAtoms( double nPlanes, bool isMC, double apothem = 850. ) const;
      double GetTrackerNAtoms( double minZ, double maxZ, bool isMC, double apothem = 850. ) const;
      //! Get the number of protons in tracker
      double GetTrackerNProtons( double nPlanes, bool isMC, double apothem = 850. ) const;
      double GetTrackerNProtons( double minZ, double maxZ, bool isMC, double apothem = 850. ) const;
      //! Get the number of neutrons in tracker
      double GetTrackerNNeutrons( double nPlanes, bool isMC, double apothem = 850. ) const;
      double GetTrackerNNeutrons( double minZ, double maxZ, bool isMC, double apothem = 850. ) const;
      //! Get the number of nucleons (n+p) in tracker
      double GetTrackerNNucleons( double nPlanes, bool isMC, double apothem = 850. ) const;
      double GetTrackerNNucleons( double minZ, double maxZ, bool isMC, double apothem = 850. ) const;

      //get information about elements of the tracker
      double GetTrackerElementMassFraction( int targetZ, bool isMC ) const;
      double GetTrackerElementA( int targetZ ) const;
      double GetTrackerElementAtomsPerGram( int targetZ ) const;

      //! Get the number of EL atoms in tracker
      double GetTrackerElementNAtoms( int elementZ, double nPlanes, bool isMC, double apothem = 850. ) const;
      double GetTrackerElementNAtoms( int elementZ, double minZ, double maxZ, bool isMC, double apothem = 850. ) const;
      //! Get the number of C atoms in tracker (used frequently so is gets a special function)
      double GetTrackerNCarbonAtoms( double nPlanes, bool isMC, double apothem = 850. ) const;
      double GetTrackerNCarbonAtoms( double minZ, double maxZ, bool isMC, double apothem = 850. ) const;

      //! In tracker?
      bool InTrackerZ( double vtx_z );
      bool InTracker( double vtx_x, double vtx_y, double vtx_z, double apothem = 850. );

      //! In nuclear target region?
      bool InNukeRegionZ( double vtx_z );
      bool InNukeRegion( double vtx_x, double vtx_y, double vtx_z, double apothem = 850. );

      //==========================================
      // Functions for Passive Targets
      //==========================================
      //! Get the area of a passive target (mm^2)
      double GetPassiveTargetArea( int targetID, int targetZ, double apothem = 850. ) const;
      //! Get the areal density of a passive target ( g/mm^2 )
      double GetPassiveTargetArealDensity( int targetID, int targetZ, bool isMC ) const;
      //! Get the mass of a passive target (g)
      double GetPassiveTargetMass( int targetID, int targetZ, bool isMC, double apothem = 850. ) const;
      //! Get the number of protons in a passive target
      double GetPassiveTargetNProtons( int targetID, int targetZ, bool isMC, double apothem = 850. ) const;
      //! Get the number of protons in a passive target
      double GetPassiveTargetNNeutrons( int targetID, int targetZ, bool isMC, double apothem = 850. ) const;
      //! Get the number of nucleons( n+p ) in a passive target
      double GetPassiveTargetNNucleons( int targetID, int targetZ, bool isMC, double apothem = 850. ) const;

      //! Get Target 1 Center Z MC
      double GetTarget1CenterZMC() const       { return TargetProp::CenterZMC::Tgt1::center_nooffset +  ( useCenterZMCOffset_ ? TargetProp::CenterZMC::Tgt1::offset : 0 ); } ;
      //! Get Target 2 Center Z MC
      double GetTarget2CenterZMC() const       { return TargetProp::CenterZMC::Tgt2::center_nooffset +  ( useCenterZMCOffset_ ? TargetProp::CenterZMC::Tgt2::offset : 0 ); } ;
      //! Get Target 3 Center Z MC
      double GetTarget3CenterZMC() const       { return TargetProp::CenterZMC::Tgt3::center_nooffset +  ( useCenterZMCOffset_ ? TargetProp::CenterZMC::Tgt3::offset : 0 ); } ;
      //! Get Target 3 Carbon Center Z MC
      double GetTarget3CarbonCenterZMC() const { return TargetProp::CenterZMC::Tgt3C::center_nooffset + ( useCenterZMCOffset_ ? TargetProp::CenterZMC::Tgt3C::offset : 0 ); } ;
      //! Get Target 4 Center Z MC
      double GetTarget4CenterZMC() const       { return TargetProp::CenterZMC::Tgt4::center_nooffset +  ( useCenterZMCOffset_ ? TargetProp::CenterZMC::Tgt4::offset : 0 ); } ;
      //! Get Target 5 Center Z MC
      double GetTarget5CenterZMC() const       { return TargetProp::CenterZMC::Tgt5::center_nooffset +  ( useCenterZMCOffset_ ? TargetProp::CenterZMC::Tgt5::offset : 0 ); } ;

      // Is this a particular nuclei?
      bool OnCarbon( int nucleiZ ) { return nucleiZ == 6;  }
      bool OnIron(   int nucleiZ ) { return nucleiZ == 26; }
      bool OnLead(   int nucleiZ ) { return nucleiZ == 82; }
      bool OnWater(  int nucleiZ ) { return nucleiZ == 8 || nucleiZ == 1; }

      //! Note, the water target geometry is far more complicated than shown here.  This is a first pass approximation
      //! In the MC target volume and has the correct nuclei? Vtx inputs are mm
      bool InTarget1MC( double vtx_x, double vtx_y, double vtx_z, int mc_nucleiZ, double apothem = 850., bool excludeBuffer = false );
      bool InTarget2MC( double vtx_x, double vtx_y, double vtx_z, int mc_nucleiZ, double apothem = 850., bool excludeBuffer = false );
      bool InTarget3MC( double vtx_x, double vtx_y, double vtx_z, int mc_nucleiZ, double apothem = 850., bool excludeBuffer = false );
      bool InTarget4MC( double vtx_x, double vtx_y, double vtx_z, int mc_nucleiZ, double apothem = 850. );
      bool InTarget5MC( double vtx_x, double vtx_y, double vtx_z, int mc_nucleiZ, double apothem = 850., bool excludeBuffer = false );
      bool InWaterTargetMC( double vtx_x, double vtx_y, double vtx_z, int mc_nucleiZ, double apothem = 850. );
      bool InPassiveTargetMC( double vtx_x, double vtx_y, double vtx_z, int mc_nucleiZ, double apothem = 850., bool excludeBuffer = false );     
 
      //Using MC volumes, are you in the volume of a target of the given material?  inputs are mm
      bool InCarbonTargetVolMC( double vtx_x, double vtx_y, double vtx_z, double apothem = 850., bool excludeBuffer= false );
      bool InIronTargetVolMC( double vtx_x, double vtx_y, double vtx_z, double apothem = 850., bool excludeBuffer= false );
      bool InLeadTargetVolMC( double vtx_x, double vtx_y, double vtx_z, double apothem = 850., bool excludeBuffer= false );

      //Using MC volumes, are you in the volume in a given target?  inputs are mm
      bool InTarget1VolMC( double vtx_x, double vtx_y, double vtx_z, double apothem = 850., bool excludeBuffer= false );
      bool InTarget2VolMC( double vtx_x, double vtx_y, double vtx_z, double apothem = 850., bool excludeBuffer= false );
      bool InTarget3VolMC( double vtx_x, double vtx_y, double vtx_z, double apothem = 850., bool excludeBuffer= false );
      bool InTarget4VolMC( double vtx_x, double vtx_y, double vtx_z, double apothem = 850. );
      bool InTarget5VolMC( double vtx_x, double vtx_y, double vtx_z, double apothem = 850., bool excludeBuffer= false );
      bool InWaterTargetVolMC( double vtx_x, double vtx_y, double vtx_z, double apothem = 850. );
      bool InPassiveTargetVolMC( double vtx_x, double vtx_y, double vtx_z, double apothem = 850., bool excludeBuffer = false );     

      //Using MC volumes, are you in the z range for a given target?  inputs are mm
      bool InTarget1ZMC( double vtx_z, int nucleiZ );
      bool InTarget2ZMC( double vtx_z, int nucleiZ );
      bool InTarget3ZMC( double vtx_z, int nucleiZ );
      bool InTarget4ZMC( double vtx_z, int nucleiZ );
      bool InTarget5ZMC( double vtx_z, int nucleiZ );
      bool InWaterTargetZMC( double vtx_z );

      //Using MC volumes, which target/material are you in? inputs are mm
      bool InIron1VolMC( double vtx_x, double vtx_y, double vtx_z, double apothem = 850., bool excludeBuffer = false );
      bool InLead1VolMC( double vtx_x, double vtx_y, double vtx_z, double apothem = 850., bool excludeBuffer = false );
      bool InIron2VolMC( double vtx_x, double vtx_y, double vtx_z, double apothem = 850., bool excludeBuffer = false );
      bool InLead2VolMC( double vtx_x, double vtx_y, double vtx_z, double apothem = 850., bool excludeBuffer = false );
      bool InCarbon3VolMC( double vtx_x, double vtx_y, double vtx_z, double apothem = 850., bool excludeBuffer = false );
      bool InIron3VolMC( double vtx_x, double vtx_y, double vtx_z, double apothem = 850., bool excludeBuffer = false );
      bool InLead3VolMC( double vtx_x, double vtx_y, double vtx_z, double apothem = 850., bool excludeBuffer = false );
      bool InLead4VolMC( double vtx_x, double vtx_y, double vtx_z, double apothem = 850. );
      bool InIron5VolMC( double vtx_x, double vtx_y, double vtx_z, double apothem = 850., bool excludeBuffer = false );
      bool InLead5VolMC( double vtx_x, double vtx_y, double vtx_z, double apothem = 850., bool excludeBuffer = false );

      double GetCoordC( double vtx_x, double vtx_y);
      double GetCoordD( double vtx_x, double vtx_y);
      double GetCoordU( double vtx_x, double vtx_y);

      double GetDistToDivCut( ) const  { return distToDivCut_; };
      void SetDistToDivCut( double x ) { distToDivCut_ = x; };

      bool GetCenterZMCOffset( ) const                   { return useCenterZMCOffset_; };
      void SetCenterZMCOffset( bool useCenterZMCOffset )  { useCenterZMCOffset_ = useCenterZMCOffset; };

      double GetNPlanes(double minZ, double maxZ) const;
    private:
      double distToDivCut_; ///< Points this close to the division between passive target sections are not allowed (default = 25mm)
      bool useCenterZMCOffset_; // Do you offset the center of the targets by the offsets? Should be on

  }; // end class TargetUtils

}//end namespace PlotUtils

#endif //MNV_TargetUtils
