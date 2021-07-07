#ifndef MNV_MnvHadronReweight
#define MNV_MnvHadronReweight 1

#include "PlotUtilsPhysicalConstants.h"
#include "TargetUtils.h"
#include <TH2D.h>
#include <TTree.h>
#include <map>
#include <iostream>
#include <cstdlib>
#include <TVector3.h>
#include <algorithm>

class TF1;
/*class TTree;
class TBranch;
struct Int_t;
struct Double_t;*/
                //Tar          //nuc          //pdg          //shift
typedef std::map< int, std::map< int, std::map< int, std::map< int, TH2D* > > > > TruthKineRenormWeights;
namespace PlotUtils
{
  const double pion_mass    = 139.570; // MeV
  const double neutron_mass = 939.565; // MeV
  const double proton_mass  = 938.272; // MeV
  const double muon_mass    = 105.658; // Mev
  const double kaon_mass    = 493.677; // MeV
  
  const int PDGPION = 211;
  const int PDGPROTON = 2212;
  const int PDGKAON = 321;
  const int PDGNEUTRON = 2112;
  const int PDGNEUTRONCV = 2112 + 10000;
  
  const int UP = 0;
  const int DOWN = 1;
  const int CV   = 2;

  const int minKineEntries = 3;

  const char* const mhrw_shift_names[3] = {"up","down","cv"};

  const int mhrw_nTars       = 12;                                     //tracker, water 
  const int mhrw_targets[12] = {  1,  1,  2,  2, 3,  3,  3,  4,  5,  5, 6,  7 };
  const int mhrw_nuclei[12]  = { 26, 82, 26, 82, 6, 26, 82, 82, 26, 82, 6,  8 };

  const int mhrw_nNTR_Scint  = 7;

  struct LeadingFSParticle
  {
    double p;
    double theta;
  };

  struct InelXSecWeights
  {
    // 
    std::map<int, double> weightUp;
    std::map<int, double> weightDown;
    // Whether event has this type of particle
    std::map<int, bool> eventHas;
  };
 
  struct HadronReweightInfo {

    int nPaths;
    int nPoints;

    std::vector<int> pDG; // the pdg of the particle
    std::vector<int> intCode; // The interaction type
    
    // Tracks the track ID
    std::vector<int> trackID;
    std::vector<int> nuke; // The type of nucleus that the particle is in
    
    std::vector<double> initialE; // Energy
    std::vector<double> finalE;
    std::vector<double> initialSigmaE; // Error in energy (for dedx estimation)
    std::vector<double> finalSigmaE;
    std::vector<double> columnarDensity; // the integrated areal density using path length tool
    std::vector<double> posX; // The initial position of the particle
    std::vector<double> posY;
    std::vector<double> posZ;
  };
  
  static const int AL = 4096; // need enough array size or else it will crash
  
  class HadronBranchContainer
  {
    public:
      
      //! Default constructor
      HadronBranchContainer(TTree* fchain);
      
      //! Default destructor
      virtual ~HadronBranchContainer();
      
      Int_t GetEntry(Long64_t entry = 0, Int_t getall = 0);

      template <class T>
      void FillUnivEntry( const T& univ )
      {
        truth_hadronReweightNPaths  = univ.GetInt("truth_hadronReweightNPaths" );
        truth_hadronReweightNPoints = univ.GetInt("truth_hadronReweightNPoints");

        truth_hadronReweightIntCode_sz = univ.GetInt("truth_hadronReweightIntCode_sz");
        for( int i = 0; i < truth_hadronReweightIntCode_sz; ++i ) truth_hadronReweightIntCode[i] = univ.GetVecElemInt("truth_hadronReweightIntCode",i);   

        truth_hadronReweightNuke_sz = univ.GetInt("truth_hadronReweightNuke_sz");
        for( int i = 0; i < truth_hadronReweightNuke_sz; ++i ) truth_hadronReweightNuke[i] = univ.GetVecElemInt("truth_hadronReweightNuke",i);   

        truth_hadronReweightPDG_sz = univ.GetInt("truth_hadronReweightPDG_sz");
        for( int i = 0; i < truth_hadronReweightPDG_sz; ++i ) truth_hadronReweightPDG[i] = univ.GetVecElemInt("truth_hadronReweightPDG",i);   

        truth_hadronReweightTrackID_sz = univ.GetInt("truth_hadronReweightTrackID_sz");
        for( int i = 0; i < truth_hadronReweightTrackID_sz; ++i ) truth_hadronReweightTrackID[i] = univ.GetVecElemInt("truth_hadronReweightTrackID",i);   

        truth_hadronReweightColumnarDensity_sz = univ.GetInt("truth_hadronReweightColumnarDensity_sz");
        for( int i = 0; i < truth_hadronReweightColumnarDensity_sz; ++i ) truth_hadronReweightColumnarDensity[i] = univ.GetVecElem("truth_hadronReweightColumnarDensity",i);   

        truth_hadronReweightFinalE_sz = univ.GetInt("truth_hadronReweightFinalE_sz");
        for( int i = 0; i < truth_hadronReweightFinalE_sz; ++i ) truth_hadronReweightFinalE[i] = univ.GetVecElem("truth_hadronReweightFinalE",i);   
        truth_hadronReweightFinalSigmaE_sz = univ.GetInt("truth_hadronReweightFinalSigmaE_sz");
        for( int i = 0; i < truth_hadronReweightFinalSigmaE_sz; ++i ) truth_hadronReweightFinalSigmaE[i] = univ.GetVecElem("truth_hadronReweightFinalSigmaE",i);   

        truth_hadronReweightInitialE_sz = univ.GetInt("truth_hadronReweightInitialE_sz");
        for( int i = 0; i < truth_hadronReweightInitialE_sz; ++i ) truth_hadronReweightInitialE[i] = univ.GetVecElem("truth_hadronReweightInitialE",i);   
        truth_hadronReweightInitialSigmaE_sz = univ.GetInt("truth_hadronReweightInitialSigmaE_sz");
        for( int i = 0; i < truth_hadronReweightInitialSigmaE_sz; ++i ) truth_hadronReweightInitialSigmaE[i] = univ.GetVecElem("truth_hadronReweightInitialSigmaE",i);   

        truth_hadronReweightPosX_sz = univ.GetInt("truth_hadronReweightPosX_sz");
        for( int i = 0; i < truth_hadronReweightPosX_sz; ++i ) truth_hadronReweightPosX[i] = univ.GetVecElem("truth_hadronReweightPosX",i);   
        truth_hadronReweightPosY_sz = univ.GetInt("truth_hadronReweightPosY_sz");
        for( int i = 0; i < truth_hadronReweightPosY_sz; ++i ) truth_hadronReweightPosY[i] = univ.GetVecElem("truth_hadronReweightPosY",i);   
        truth_hadronReweightPosZ_sz = univ.GetInt("truth_hadronReweightPosZ_sz");
        for( int i = 0; i < truth_hadronReweightPosZ_sz; ++i ) truth_hadronReweightPosZ[i] = univ.GetVecElem("truth_hadronReweightPosZ",i);     

        truth_hadronReweightIntCodePerSegment_sz = univ.GetInt("truth_hadronReweightIntCodePerSegment_sz");
        for( int i = 0; i < truth_hadronReweightIntCodePerSegment_sz; ++i ) truth_hadronReweightIntCodePerSegment[i] = univ.GetVecElemInt("truth_hadronReweightIntCodePerSegment",i);     

        mc_targetZ = univ.GetInt("mc_targetZ");
        for( int i = 0; i < 3; ++i ) mc_vtx[i] = univ.GetVecElem("mc_vtx",i);

        //Info on FS particles
        mc_nFSPart = univ.GetInt("mc_nFSPart");
        if (mc_nFSPart >= AL || mc_nFSPart <0){
           std::cout <<  "Warning: MnvHadronReweight -  mc_nFSPart out of bounds" << mc_nFSPart << std::endl;
           mc_nFSPart = 0;
        }
      
        for( int i = 0; i < mc_nFSPart; ++i )
        {
          mc_FSPartPDG[i] = univ.GetVecElemInt("mc_FSPartPDG",i);
          mc_FSPartE[i]   = univ.GetVecElem("mc_FSPartE",i);
          mc_FSPartPx[i]  = univ.GetVecElem("mc_FSPartPx",i);
          mc_FSPartPy[i]  = univ.GetVecElem("mc_FSPartPy",i);
          mc_FSPartPz[i]  = univ.GetVecElem("mc_FSPartPz",i);
        }

      }
      void setOSF(bool isOSF){fOSF=isOSF;};
      bool fOSF;
      bool defaultmode;
    
      Int_t           truth_hadronReweightNPaths;
      Int_t           truth_hadronReweightNPoints;
      Int_t           truth_hadronReweightIntCode_sz;
      Int_t           truth_hadronReweightIntCode[AL];   //[truth_hadronReweightIntCode_sz]
      Int_t           truth_hadronReweightNuke_sz;
      Int_t           truth_hadronReweightNuke[AL];   //[truth_hadronReweightNuke_sz]
      Int_t           truth_hadronReweightPDG_sz;
      Int_t           truth_hadronReweightPDG[AL];   //[truth_hadronReweightPDG_sz]
      Int_t           truth_hadronReweightTrackID_sz;
      Int_t           truth_hadronReweightTrackID[AL];   //[truth_hadronReweightTrackID_sz]
      Int_t           truth_hadronReweightColumnarDensity_sz;
      Double_t        truth_hadronReweightColumnarDensity[AL];   //[truth_hadronReweightColumnarDensity_sz]
      Int_t           truth_hadronReweightFinalE_sz;
      Double_t        truth_hadronReweightFinalE[AL];   //[truth_hadronReweightFinalE_sz]
      Int_t           truth_hadronReweightFinalSigmaE_sz;
      Double_t        truth_hadronReweightFinalSigmaE[AL];   //[truth_hadronReweightFinalSigmaE_sz]
      Int_t           truth_hadronReweightInitialE_sz;
      Double_t        truth_hadronReweightInitialE[AL];   //[truth_hadronReweightInitialE_sz]
      Int_t           truth_hadronReweightInitialSigmaE_sz;
      Double_t        truth_hadronReweightInitialSigmaE[AL];   //[truth_hadronReweightInitialSigmaE_sz]
      Int_t           truth_hadronReweightPosX_sz;
      Double_t        truth_hadronReweightPosX[AL];   //[truth_hadronReweightPosX_sz]
      Int_t           truth_hadronReweightPosY_sz;
      Double_t        truth_hadronReweightPosY[AL];   //[truth_hadronReweightPosY_sz]
      Int_t           truth_hadronReweightPosZ_sz;
      Double_t        truth_hadronReweightPosZ[AL];   //[truth_hadronReweightPosZ_sz]  
      Int_t           truth_hadronReweightIntCodePerSegment_sz;
      Int_t           truth_hadronReweightIntCodePerSegment[AL];   //[truth_hadronReweightPosZ_sz]  

      Int_t           mc_targetZ;
      Double_t        mc_vtx[3];

      //Info on FS particles
      Int_t           mc_nFSPart;
      Int_t           mc_FSPartPDG[AL];
      Double_t        mc_FSPartE[AL];
      Double_t        mc_FSPartPx[AL];
      Double_t        mc_FSPartPy[AL];
      Double_t        mc_FSPartPz[AL];

      TBranch        *b_truth_hadronReweightNPaths;   //!
      TBranch        *b_truth_hadronReweightNPoints;   //!
      TBranch        *b_truth_hadronReweightIntCode_sz;   //!
      TBranch        *b_truth_hadronReweightIntCode;   //!
      TBranch        *b_truth_hadronReweightNuke_sz;   //!
      TBranch        *b_truth_hadronReweightNuke;   //!
      TBranch        *b_truth_hadronReweightPDG_sz;   //!
      TBranch        *b_truth_hadronReweightPDG;   //!
      TBranch        *b_truth_hadronReweightTrackID_sz;   //!
      TBranch        *b_truth_hadronReweightTrackID;   //!
      TBranch        *b_truth_hadronReweightColumnarDensity_sz;   //!
      TBranch        *b_truth_hadronReweightColumnarDensity;   //!
      TBranch        *b_truth_hadronReweightFinalE_sz;   //!
      TBranch        *b_truth_hadronReweightFinalE;   //!
      TBranch        *b_truth_hadronReweightFinalSigmaE_sz;   //!
      TBranch        *b_truth_hadronReweightFinalSigmaE;   //!
      TBranch        *b_truth_hadronReweightInitialE_sz;   //!
      TBranch        *b_truth_hadronReweightInitialE;   //!
      TBranch        *b_truth_hadronReweightInitialSigmaE_sz;   //!
      TBranch        *b_truth_hadronReweightInitialSigmaE;   //!
      TBranch        *b_truth_hadronReweightPosX_sz;   //!
      TBranch        *b_truth_hadronReweightPosX;   //!
      TBranch        *b_truth_hadronReweightPosY_sz;   //!
      TBranch        *b_truth_hadronReweightPosY;   //!
      TBranch        *b_truth_hadronReweightPosZ_sz;   //!
      TBranch        *b_truth_hadronReweightPosZ;   //!
      TBranch        *b_truth_hadronReweightIntCodePerSegment_sz;   //!
      TBranch        *b_truth_hadronReweightIntCodePerSegment;   //!

      TBranch        *b_mc_targetZ;
      TBranch        *b_mc_vtx;

      TBranch        *b_mc_nFSPart;
      TBranch        *b_mc_FSPartPDG;
      TBranch        *b_mc_FSPartE;
      TBranch        *b_mc_FSPartPx;
      TBranch        *b_mc_FSPartPy;
      TBranch        *b_mc_FSPartPz;
      
      TTree          *fChain;
  };
  
  class MnvHadronReweight
  {
 public:
      //static MnvH1D* getPlotFromVertError(MnvH1D* hist, const std::string &name, int number);
      
      //static MnvHadronReweight* get(const char* projectName); 
      static MnvHadronReweight* get(); // TODO remove singleton because there's a lot of initializing todo, like fiducial volume and project name
      
      //! Default constructor
      MnvHadronReweight(TTree* truth=0, TTree* data=0, bool hdXSec = false);
      
      //! Default destructor
      virtual ~MnvHadronReweight();
      
      void setTruthTree(TTree* truth);

      void setDataTree(TTree* data);
      
      //Original Renormalization (single factor)
      void getRenormFactors(const char* filename, const char* projectname, TTree* truth=0);
      
      bool tryLoadingFromFile(const char* filename, const char* projectname);
      
      bool saveRenormFactorsToFile(const char* filename, const char* projectname);
      
      std::string makefullfilename(const std::string filename, const std::string projectname);
      
      void setRenorm(bool nRenorm) { renorm = nRenorm; };
      
      bool getRenorm() { return renorm; };
      
      //Kinematic renormalization
      void getTruthKineRenorm();

      bool tryLoadingKineFile();
      
      bool saveKineRenormFactorsToFile();
      
      std::string makekinefilename();

      std::string getTargetNucleiName( int tar, int nuc );

      /*  
      template <class T>
      std::map<int, LeadingFSParticle> getLeadingFSParticles( const T& univ )
      {
               //pdg            idx  energy
        std::map<int, std::pair<int, double > > leadingPartIdxE;

        for( int iPart = 0; iPart < univ.GetInt("mc_nFSPart"); ++iPart )
        {
          const int pdg = univ.GetVecElemInt("mc_FSPartPDG",iPart);
          bool bPart     = particles.find(pdg) != particles.end();
          bool bAntiPart = particles.find( getAntiPDG(pdg) ) != particles.end();

          // Is the particle or antiparticle in particles vector?
          if( !bPart && !bAntiPart ) continue;

          // If it is the antiparticle, do we care?
          if( bAntiPart && !particles[getAntiPDG(pdg)] ) continue;

          //Make the pdg we store this match those in particles vector
          const int part_pdg = bPart ? pdg : getAntiPDG(pdg);

          if( leadingPartIdxE.find(part_pdg) == leadingPartIdxE.end() )
          {
            leadingPartIdxE[part_pdg].first  = iPart;
            leadingPartIdxE[part_pdg].second = univ.GetVecElem("mc_FSPartE",iPart);
          }
          else
          {
            if( univ.GetVecElem("mc_FSPartE",iPart) > leadingPartIdxE[part_pdg].second )
            {
              leadingPartIdxE[part_pdg].first  = iPart;
              leadingPartIdxE[part_pdg].second = univ.GetVecElem("mc_FSPartE",iPart);
            }
          } 
        }

        std::map<int, LeadingFSParticle> FSParts;

        //Loop through map
        //Make TVector3, make magnitude and angle, import into FSParts
        for( auto &lPart : leadingPartIdxE )
        {
          const int pdg = lPart.first;
          const int fs_idx = lPart.second.second;    

          TVector3 p3FS( univ.GetVecElem("mc_FSPartPx",fs_idx), univ.GetVecElem("mc_FSPartPy",fs_idx), univ.GetVecElem("mc_FSPartPz",fs_idx) );
          p3FS.RotateX(MinervaUnits::numi_beam_angle_rad); //Do I actually need this if it's internal?

          FSParts[pdg].p     = p3FS.Mag();
          FSParts[pdg].theta = p3FS.Theta()*180/3.14159;
        }
      
        return FSParts;
      }
      */ 

      void setRenormKine(bool nRenorm) { renormKine = nRenorm; };
      
      bool getRenormKine() { return renormKine; };
      
      void initKineRenormWeights( TruthKineRenormWeights& tkrw, int initialValue = 0.0 );

      void initializeWeights(InelXSecWeights& weights, double initialValue = 1.0);

      InelXSecWeights getWeights( Long64_t entry, bool usetruth = false );

      template <class T>
      InelXSecWeights getWeights( const T& univ,  bool usetruth = false )
      {
        //Check if weights have been calculated for this entry
        if( univ.GetEntry() == m_current_entry_weight ) return m_current_weights;
        else m_current_entry_weight = univ.GetEntry();

        //Need to have a HadronBranchContainer
        HadronBranchContainer* b;
        if (usetruth) b = fTruth;
        else b = fData;

        b->FillUnivEntry( univ );

        return getWeights( b );
         
      };



      TTree* makeTree(const char* filename, const char* treename = "reweight_tree");
      TTree* makeTree(const char* filename, bool usetruth, const char* treename = "reweight_tree");
      
      std::vector<double> getGeant4Weight(int pdg, std::string xsec, double p_i, double p_f, double path_mm, double delta);

      //bool getWeight(const int nuke, const int pdg, const double density, const double Ei, const double Ef, const double delta, double& up, double& down, double& val);
      bool getWeight(const int nuke, const int pdg, const int intCode, const double density, const double Ei, const double Ef, const double delta, double& up, double& down);
      bool getWeight(InelXSecWeights& weights, const int intCode, const double valA, const double valB, const int pdg, const bool treatAnti);
      bool inFiducial(HadronBranchContainer* b, int iPrimTraj);
      bool inFiducial(const double x, const double y, const double z);
      double getDelta(HadronBranchContainer* b, int iPrimTraj);
      double getDelta(const int pdg, const int target, const double ke);
      double getDeltaElastic(const int pdg, const int nuke, const double ke);
      double NeutronReweightAmount(double neutronKE);
      double getMass(const int pdg);
      inline int getAntiPDG(const int pdg)
      {
        // Leave open the possibility that there are other ways to do antiparticles
        return -pdg;
      }
      inline int getNormalPDG(const int pdg)
      {
        // Leave open the possibility that there are other ways to do antiparticles
        return abs(pdg);
      }

      void setReadoutVolume( std::string volname );
      void setBasicFiducial(double minZ = 5800.0, double maxZ = 8600.0, double apothem = 894.45);
      void getBasicFiducial(double& minZ, double& maxZ, double& apothem);
      void setFiducialVolumeType(int type);
      void setdefaultmode(bool mode = false);
      void setParticle(int pdg, bool treatAntiParticleSame);
      void removeParticle(int pdg);
      void useDefaultParticles();
      double getCorrectedDensity(const int nuke, const double density);
      void setDeltaScale(double newscale);
      double ezrand();
      void changeLowEElasticXSec(double multiplier);
      void useHDXSec(bool useHDXSec = true);
      enum XSecType
      {
        xsectype_inelastic,
        xsectype_elastic,
        xsectype_total
      };
      bool getXSecFunc(int pdg, int nuke, XSecType xsectype, TF1& function);
      bool getWeightWithElastic(const int nuke, const int pdg, const int intCode, const double density, const double Ei, const double Ef, const double delta_inelastic, const double delta_elastic, double& up, double& down);
      void TurnOffElasticReweight();
      void TurnOnFakeElastics();
      double doWeightCalcInteracted(double rhoX, double sigma_partial_data, double sigma_partial_geant, double sigma_total_data, double sigma_total_geant);
      double doWeightCalcNoninteracted(double rhoX, double sigma_total_data, double sigma_total_geant);
      double doWeightCalcInelasticOnly(double rhoX, double sigma_inelastic_data, double sigma_inelastic_geant);
      double getAverageCrossSection2(int pdg, int nuke, double T_i, double T_f, XSecType xsectype, bool isCVCalc = false);
      double getAverageNewNeutronCrossSection(int nuke, double Ti, double Tf, XSecType xsectype);
      int GetNFakeElastics();
      void ResetFakeElasticsSeed();
      void TurnOnElasticReweightForParticlesOtherThanNeutron();
      void useReweightedNeutronCV();
      void setReweightedNeutronCV( bool reweightNeutronCV )
      {
        m_reweightNeutronCV = reweightNeutronCV;
        if( m_reweightNeutronCV ) useReweightedNeutronCV();
      }
      bool getWeightToNewNeutronXSec(int nuke, int intCode, double density, double Ei, double Ef, double& wgt);

      double reweightNeutronCV(Long64_t entry, bool usetruth = false );

      template <class T>
      double reweightNeutronCV( const T& univ, bool usetruth = false ) 
      {
        if (!m_reweightNeutronCV) return 1.0;
        if( univ.GetEntry() == m_current_entry_neutronCV ) return m_current_rwgt_neutronCV;
        else m_current_entry_neutronCV = univ.GetEntry();
                
        // Get the branch containing all event data
        HadronBranchContainer* b;
        if (usetruth) b = fTruth;
        else b = fData;

        b->FillUnivEntry( univ );

        return reweightNeutronCV( b );
      }


      void setInputDirectory( std::string directory ) 
      { 
        if( !m_directory.empty() && directory != m_directory ) std::cout<<"Changing directory "<<m_directory<<" to "<<directory<<std::endl;
        m_directory = directory; 
      }
      std::string getInputDirectory() { return m_directory; }

      void setProjectName( std::string projectname ) 
      { 
        if( !m_projectname.empty() && projectname != m_projectname ) std::cout<<"Changing projectname "<<m_projectname<<" to "<<projectname<<std::endl;
        m_projectname = projectname; 
      }
      std::string getProjectName() { return m_projectname; }

      void setPlaylist( std::string playlist ) 
      { 
        std::string tmp_playlist = playlist;
        //String formatting
        if( playlist.find("minerva") != std::string::npos  || playlist.find("Minerva") != std::string::npos )
        {
          tmp_playlist = playlist.substr( 7, playlist.size()-7 );
        }

        std::transform( tmp_playlist.begin(), tmp_playlist.end(), tmp_playlist.begin(), ::toupper);
        
        if( !m_playlist.empty() && tmp_playlist != m_playlist ) std::cout<<"Changing playlist "<<m_playlist<<" to "<<tmp_playlist<<std::endl;

        m_playlist = tmp_playlist; 
      }
      std::string getPlaylist() { return m_playlist; }

 private:

   static TargetUtils* m_TargetUtils;

   static MnvHadronReweight* instance;
   // TODO, have 3 options for reweight amount
   // 1 flat value,
   // 2 TFormula as function of single True variable
   // 3 TH1 as function of single true variable (with option to do lerp from center of each bin)
   bool renorm;
   // 
   InelXSecWeights renormWeights;
   // Scales delta by a constant
   double deltascale;
   bool doElasticReweight;

   bool renormKine;   
   TruthKineRenormWeights kineRenormWeights;
   TObjArray m_kineRenormHists;

   // (PDG, Treat Antiparticle As Same)
   std::map<int, bool> particles;
   
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   TTree          *fTruthChain;   //!pointer to the analyzed TTree or TChain


   
   HadronBranchContainer *fData;
   HadronBranchContainer *fTruth;
   
   std::map<int, double> deltas;
   std::map<int, int> specialfunctiondelta;
   bool bNeutrinoMode;
   std::string m_directory;
   std::string m_projectname;
   std::string m_playlist;  
   std::string m_kinefilename;
 
   int fiducialType;
   double fiducialMinZ;
   double fiducialMaxZ;
   double fiducialApothem;
   bool defaultmode;
   double minimum_fake_elastic_threshold;
   bool look_for_fake_elastics;
   int nfakeelastics;
   bool elasticReweightForAllParticles;
   bool m_reweightNeutronCV;

   int m_current_entry_neutronCV;
   int m_current_entry_weight;
   InelXSecWeights m_current_weights;
   double m_current_rwgt_neutronCV;

   InelXSecWeights getWeights( HadronBranchContainer* b );
   double reweightNeutronCV( HadronBranchContainer *b );

   std::map<int, LeadingFSParticle> getLeadingFSParticles( HadronBranchContainer* b );
   std::pair<int,int> getTargetNuclei( HadronBranchContainer* b ); 
  
   bool isInHexagon(double x, double y, double apothem);

  };

} // end namespace HadronReweight


#endif //MNV_MnvHadronReweight
