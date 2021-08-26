#ifndef MNV_MnvNuclearModelWeight
#define MNV_MnvNuclearModelWeight 1

#include <string>
#include <TF1.h>
#include <TF2.h>
#include <TH2F.h>
#include <TMath.h>
#include <TObject.h>
#include <TFile.h>
#include <vector>
#include "PlotUtils/NuclModUtils.h"



namespace PlotUtils
{


  //==========================================
  // MnvNuclearModelWeight class
  //==========================================
  //! Class to reweight Resonacne and DIS events to different nuclear models
  //! Autor Joel Mousseau
  //! Send all complaints to /dev/null, brownies to WH10XO
  class MnvNuclearModelWeight : public NuclModUtils
  {
    public:
      //! Default constructor
      MnvNuclearModelWeight();

      //! Specify playlist and possibly analysis
      explicit MnvNuclearModelWeight( const std::string& model );

      //! Default destructor
      virtual ~MnvNuclearModelWeight();


      //=========================================================
      // Each model can have it's nuclear effects parametrized one of two ways.
      // One: A flat f(x) for each given nuclei. This assumes F2 == xF3 Call these function f_C, f_Fe and f_Pb. isFlat = true
      // Two: A ratio of F2_A / F2_C and xF3_A / xF3_C. This assumes F2 != xF3 Call these F2_C, F2_Fe, F2_Pb, xF3_C, xF3_Fe and xF3_Pb. isFlat = false
      // todo: make them read-only to outsiders. Why not
      //=========================================================

      //!Is this model or fit an Overall scaling of F2 and xF3?
      bool isOverallScaling;
       
      
      //!Overall Scalings
      TF1 f_C;
      TF1 f_Fe;
      TF1 f_Pb;
      
      //!If you need a function to go from free nucleon to deuterium
      TF1 f_D;
      
      
      //!Structure Functions.
      TH2F *F1_P;
      TH2F *F1_N;
      TH2F *F1_CH;
      TH2F *F1_C;
      TH2F *F1_Fe;
      TH2F *F1_Pb;
      
      TH2F *F2_P;
      TH2F *F2_N;
      TH2F *F2_CH;
      TH2F *F2_C;
      TH2F *F2_Fe;
      TH2F *F2_Pb;
      
      TH2F *xF3_P;
      TH2F *xF3_N;
      TH2F *xF3_CH;
      TH2F *xF3_C;
      TH2F *xF3_Fe;
      TH2F *xF3_Pb;
      
      //!Genie's Structure Functions. 
      
      //! These are for bound Fe
      
      //! Bound proton
      TH2F *GF1FeP;
      TH2F *GF2FeP;
      TH2F *GxF3FeP;
      
      //! Bound neutron
      TH2F *GF1FeN;
      TH2F *GF2FeN;
      TH2F *GxF3FeN;
      
      //! Free proton
      TH2F *GF1P;
      TH2F *GF2P;
      TH2F *GxF3P;
      
      //! Free neutron
      TH2F *GF1N;
      TH2F *GF2N;
      TH2F *GxF3N;
      
      //!Actually Calculate the Weight
      double CalculateWeight(double Q2, double x, double y, double Enu, int tgtNucleon, int Target_A);
      
      //!Calculate the Differential Cross Section
      double CalculateDifferentialCrossSection(double Q2, double x, double y, double Enu, int Target_A);
      double CalculateDifferentialCrossSection(int q2Bin, int xBin, double x, double y, double Enu, int Target_A);
      
      //!Calculate GENIE'S Cross Section, which is Z independent under our ansantz 
      double CalculateGENIEDifferentialCrossSection(double Q2, double x, double y, double Enu, int Target_A);
      double CalculateGENIEDifferentialCrossSection(int q2Bin, int xBin, double x, double y, double Enu, int Target_A);      

    private:
      //This code does all of the math
      
      
      //!Calculators for R
      double CalcR(double Q2, double x, int fit) const;
      double CalcR(double Q2, double x) const;
      
      //! Not an agnle. Used for alterante R calculation
      double Theta(double Q2, double x) const;      
      
      //! What model you are re-weighting to
      std::string model_ ;
      
      //! Location of Genie Structure Functions
      std::string GenieSFLoc;
      
      //! Location of the structure functions for a parciular model.
      std::string ModelSFLoc;
      
      //! Kinematic Range the Input fit is valid over
      double minQ2;
      double maxQ2;
      double minX;
      double maxX;
      

      
      

  }; // end class MnvNuclearModelWeight

}//end namespace PlotUtils

#endif //MNV_MnvNuclearModelWeight
