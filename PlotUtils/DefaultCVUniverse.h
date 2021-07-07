#ifndef DEFAULTCVUNIVERSE_H
#define DEFAULTCVUNIVERSE_H

#include <iostream>
#include "ChainWrapper.h"
#include "FluxReweighter.h"
#include "Math/Vector4D.h"
#include "NSFDefaults.h"
#include "PlotUtilsPhysicalConstants.h"

namespace PlotUtils {
  /*! @brief Base class for all user and systematic universes
  
    @author Ben Messerly
  */
  class DefaultCVUniverse {
    protected:
    //! helper functions for calculating physics angles from measured angles
      double phi3D(double thetaX, double thetaY) const;
      double theta3D(double thetaX, double thetaY) const;
      double calcq0(const double Enu, const double Elep) const;
      double calcq3(const double Q2, const double Enu, const double Elep) const;
  
    // member variables
      PlotUtils::TreeWrapper* m_chw; // chain of events
      double m_nsigma;               // n sigma that universe event shifts
      Long64_t m_entry;              // current chain entry

    // Particle response systematic quantity
      std::vector<int> m_non_cal_idx;
  
    public:
      DefaultCVUniverse() {}  // Not to be used...yet.  Need for python
  
      // User hook to add more constructor parameters!  To define your
      // own config_t, just add this line in your CVUnvierse:
      //
      // typdef SomeStruct config_t;
      //
      // If you require a c++11 compiler, you can use "using" instead.
      typedef PlotUtils::ChainWrapper* config_t;
  
      // I'm going to take TreeWrapper* and not config_t as constructor
      // argument to hack this so user-specific unvierses don't have to
      // change to keep using ChainWrapper*.
  
      //! Constructor
      DefaultCVUniverse(PlotUtils::TreeWrapper* config,
                        double nsigma = 0);  //, bool is_vertical_only = false);
  
      // Constructors we may wish to add in the future
      // DefaultCVUniverse(const std::string& anatool_name,
      //                  const std::string& playlist,
      //                  double nsigma=0)
      //  : m_chw(makeChainWrapperPtr(anatool_name, playlist)),
      //    m_nsigma(nsigma), m_entry(-1)
      //  {}
  
      // DefaultCVUniverse(const std::string& anatool_name, double nsigma=0)
      //  : m_chw(new PlotUtils::ChainWrapper(anatool_name.c_str())),
      //    m_nsigma(nsigma), m_entry(-1)
      //  {}
  
      //! Destructor
      virtual ~DefaultCVUniverse() {}
  
      //! Get Weights
      virtual double GetGenieWeight() const;
      virtual double GetRPAWeight() const;
      virtual double GetLowRecoil2p2hWeight() const;
      virtual double GetLowQ2PiWeight(std::string channel) const;
      virtual double GetCoherentPiWeight(double thpi_true, double tpi_true) const ;
      virtual double GetMKWeight() const;
      virtual double GetMinosEfficiencyWeight() const;
      virtual double GetTargetMassWeight() const;
      virtual double GetFluxAndCVWeight(double Enu = -99. /*GeV*/,
                                        int nu_pdg = -99) const;
      virtual double GetFSIWeight( int iWeight ) const;
      virtual double GetGeantHadronWeight() const;
      virtual double GetMichelEfficiencyWeight() const;

      //! Get the batch pot.  Needed for the Minos Efficiency weight, reused
      //! when finding the error
      virtual double GetBatchPOT() const;
  
      //! Get true event quantities -- these are needed to calculate several
      //! event weights
      virtual double GetEnuTrue() const;      /* MeV */
      virtual double GetElepTrue() const;     /* MeV */
      virtual double GetPlepTrue() const;     /* MeV */
      virtual double GetThetalepTrue() const; /* radians w.r.t. incident nu dirn */
      virtual double GetPhilepTrue() const;/* radians w.r.t. incident nu dirn */
  
      virtual double GetQ2True() const; /* MeV^2 */
      virtual double Getq0True() const; /* MeV */
      virtual double Getq3True() const; /* MeV */
      virtual int GetTargetZTrue() const; /* atomic number (Z) of struck nucleus */
      virtual double GetVertexZTrue() const; /* longitudinal coord of interaction vertex; mm */
  
      // Get Reco Muon Variables
      // Everything shall be derived from Pmu, theta_x and theta_y,
      // all of which can be overridden
      virtual ROOT::Math::PxPyPzEVector GetMuon4V() const; /* MeV */
      virtual double GetEmu() const;      /* MeV */
      virtual double GetPmu() const;      /* MeV */
      virtual double GetThetamu() const;  /* radians w.r.t. incident nu dirn */
      virtual double GetPhimu() const;    /* radians w.r.t. incident nu dirn */
      virtual double GetThetaXmu() const; /* radians */
      virtual double GetThetaYmu() const; /* radians */
  
      // Muon momenta functions used under the hood to calculate GetPmu
      virtual double GetPmu_nominal() const;      /* MeV */
      virtual double GetPmuMinerva() const;       /* MeV */
      virtual double GetPmuMinos() const;         /* MeV */
      virtual double GetPmuMinos_nominal() const; /* MeV */
  
      //! these allow you to implement beam angle offsets and universes for all
      //! particles. helper functions phi3D and theta3D let you calculate real phi
      //! and theta from thetax, thetay + these offsets.
  
      virtual double GetBeamAngleOffsetX() const;
      virtual double GetBeamAngleOffsetY() const;
  
      //! Particle response systematic functions
      virtual double GetRecoilEnergy() const;
      virtual double GetCalRecoilEnergy() const;
      virtual double GetNonCalRecoilEnergy() const;
      virtual double GetVertexZ() const;
      //! cut for generator level protons for CCQE
      virtual double GetTrueProtonKECut() const;
      virtual double GetRecoProtonKECut() const;
	
      //! Get Dummy Values
      double GetDummyVar() const { return -999.; }
      double GetDummyVar1Arg(const int idx1) const {
        (void)idx1;
        return -999.;
      }
      double GetDummyVar2Arg(const int idx1, const int idx2) const {
        (void)idx1;
        (void)idx2;
        return -999.;
      }
  
      //! Generic Branch Getters
      int GetInt(const char* name) const;
      bool GetBool(const char* name) const;
      double GetDouble(const char* name) const;
      int GetVecElemInt(const char* name, const int i) const;
      double GetVecElem(const char* name, const int i) const;
      double GetVecElem(const char* name, const int i, const int j) const;
      std::vector<int> GetVecInt(const char* name) const;
      std::vector<double> GetVecDouble(const char* name) const;
      std::vector<std::vector<double> > GetVecOfVecDouble(const char* name) const;
      template <typename T>
      std::vector<T> GetVec(const char* name) const {
        return m_chw->GetValueVector<T>(name, m_entry);
      }
  
      //! Class Member Setters/Getters
      void SetEntry(Long64_t entry) { m_entry = entry; }
      int GetEntry() const { return m_entry; }
      double GetSigma() const { return m_nsigma; }
  
      //! Names
      virtual std::string ShortName() const { return "cv"; }
      virtual std::string LatexName() const { return "Central value"; }
      virtual std::string GetAnaToolName() const {
        return m_chw->GetTree()->GetName();
      }
  
      //! "Vertical" = only affects weights, doesn't shift variables
      //! Use this function to optimize your cuts
      //! Override me in child classes!
      virtual bool IsVerticalOnly() const { return false; }
  
      // Some systematics (e.g. muon efficiency) should only be applied to reco,
      // hence the need for m_is_truth
      static void SetTruth(bool is_truth) { m_is_truth = is_truth; }
      bool IsTruth() const { return m_is_truth; }
  
      // Particle response systematics functions
      virtual void SetNonCalIndices(std::vector<int>& non_cal_idx) { m_non_cal_idx = non_cal_idx; }
      virtual std::vector<int> GetNonCalIndices() { return m_non_cal_idx; }
  
      // When this is on, NRP weights are applied to the CV weight
      static bool UseNonResPiReweight();
      static bool SetNonResPiReweight(bool use_constraint);
  
      static bool UseZExpansionFaReweight();
      static bool SetZExpansionFaReweight(bool use_constraint);
  
      //! Playlists
      static std::string GetPlaylist();
      static bool IsPlaylistME(std::string playlist);
      static bool SetPlaylist(std::string playlist);
  
      //! MINERvA Software Processing
      static bool SetProcessingEroica();
      static bool IsProcessingNX();
  
      //! Horn current
      static bool isFHC() { return m_is_FHC; }
      static bool isRHC() { return m_is_RHC; }
      static bool SetHornCurrent(std::string& playlist);
  
      //! Analysis neutrino identity -- used for flux weights
      static int GetAnalysisNuPDG();
      static bool SetAnalysisNuPDG(int nu_pdg);
  
      static bool UseNuEConstraint();
      static bool SetNuEConstraint(bool use_constraint);
  
      static int GetNFluxUniverses();
      static bool SetNFluxUniverses(int n_flux_universes);
  
      static double GetCVMuonMomentumOffset();
      static bool SetCVMuonMomentumOffset(double muon_momentum_cv_offset);
    
      static double GetTrueProtonKECutCentral();
      static bool SetTrueProtonKECutCentral(double m_true_proton_ke_cut_central);
    
      static double GetRecoProtonKECutCentral();
      static bool SetRecoProtonKECutCentral(double m_reco_proton_ke_cut_central);
    private:
      static bool m_use_nuE_constraint;
      static bool m_use_nonResPi_reweight;
      static bool m_use_zExpansionFa_reweight;
      static bool m_is_truth;
      static std::string m_playlist;
      static std::string m_processing;
      static bool m_is_RHC;
      static bool m_is_FHC;
      static int m_analysis_nu_pdg;
      static int m_n_flux_universes;
      static double m_muon_momentum_cv_offset;
      static double m_true_proton_ke_cut_central;
      static double m_reco_proton_ke_cut_central;
  };

}

#endif  // DEFAULTCVUNIVERSE
