
#ifndef weight_fsi_h
#define weight_fsi_h

#include <vector>

#include <Rtypes.h>


/*!
  This is a class for taking a Tproton histogram and apply reweighing on the no FSI to mimic the effects of a corrected elastic FSI. The current parameter for Transverse Imbalance Study, Phys. Rev. Lett. 121, 022504 (2018)/arXiv:1805.05486 is 1.5

  11/28/19: major revision by Rik and Trung to calculate both elastic and absorption weights

   How to use the new class:

   1) outside and before the event loop:
        weight_fsi* ptr = new weight_fsi;
   2) inside the event loop:
        ptr->CalcWeights(arguments); // if use the non-chw interface, where arguments are listed below in the code OR
        ptr->CalcWeights(chw, entry);// if use the chw interface
        double elastic_weight = ptr->GetElasticWeight(config); // config = 1 or 2
        double absorption_weight = ptr->GetAbsorptionWeight();
   3) outside and after the event loop
        delete ptr;

   config == 1: elastic --> no fsi
   config == 2: elastic --> other fsi
  
*/

namespace PlotUtils{
    
    class ChainWrapper;


    class intranuke_particle {
      public:
      intranuke_particle(int pdg, int fate, double ke) :
        __pdg(pdg), __fate(fate), __ke(ke), __weight1(-1.0), __weight2(-1.0) {}

        int __pdg;
        int __fate;
        double __ke;
        double __weight1;
        double __weight2; // reserved
    };

    class fate_code {
      public:

        fate_code();

        void Reset();
        void calcFates(PlotUtils::ChainWrapper* chw,
                       Long64_t entry,
                       bool verbose = false);

        void calcFates(int mc_incoming,
                       int mc_primaryLepton,
                       int mc_charm,
                       int mc_intType,
                       int mc_targetA,
                       int mc_targetZ,
                       int mc_er_nPart,
                       const std::vector<int>& mc_er_ID,
                       const std::vector<int>& mc_er_status,
                       const std::vector<int>& mc_er_FD,
                       const std::vector<int>& mc_er_LD,
                       const std::vector<int>& mc_er_mother,
                       const std::vector<double>& mc_er_Px,
                       const std::vector<double>& mc_er_Py,
                       const std::vector<double>& mc_er_Pz,
                       const std::vector<double>& mc_er_E,
                       Long64_t entry,                      // to check if it's the same as last entry
                       bool verbose = false);
        
        std::vector<intranuke_particle> GetIntranukeParticles();

      private:

        double getGenieBEinMeV( const int A, const int Z = 0 );
        
        std::vector<intranuke_particle> fIntranukeParticles;
        Long64_t fLastEntry;

      public:
        
        int __count_proton;
        int __count_neutron;
        int __count_piplus;
        int __count_pizero;
        int __count_piminus;
        
    };




    
    class weight_fsi
    {
      public:
        
        weight_fsi();

        void Reset();
        double GetElasticWeight(int config);
        double GetAbsorptionWeight();
        void UseTrackingThreshold();
                
        int calcWeights(PlotUtils::ChainWrapper* chw,
                        Long64_t entry,
                        bool verbose = false);

        int calcWeights(int mc_incoming,
                        int mc_primaryLepton,
                        int mc_charm,
                        int mc_intType,
                        int mc_targetA,
                        int mc_targetZ,
                        int mc_resID,
                        int mc_er_nPart,
                        const int* mc_er_ID_arr,
                        const int* mc_er_status_arr,
                        const int* mc_er_FD_arr,
                        const int* mc_er_LD_arr,
                        const int* mc_er_mother_arr,
                        const double* mc_er_Px_arr,
                        const double* mc_er_Py_arr,
                        const double* mc_er_Pz_arr,
                        const double* mc_er_E_arr,
                        Long64_t entry,             // to check if it's the same as last entry
                        bool verbose = false);
        
        int calcWeights(int mc_incoming,
                        int mc_primaryLepton,
                        int mc_charm,
                        int mc_intType,
                        int mc_targetA,
                        int mc_targetZ,
                        int mc_resID,
                        int mc_er_nPart,
                        const std::vector<int>& mc_er_ID,
                        const std::vector<int>& mc_er_status,
                        const std::vector<int>& mc_er_FD,
                        const std::vector<int>& mc_er_LD,
                        const std::vector<int>& mc_er_mother,
                        const std::vector<double>& mc_er_Px,
                        const std::vector<double>& mc_er_Py,
                        const std::vector<double>& mc_er_Pz,
                        const std::vector<double>& mc_er_E,
                        Long64_t entry,             // to check if it's the same as last entry
                        bool verbose = false);

        int __count_elastic_nucleon;
        int __count_elastic_pion;
        double __elastic_weight_nucleon;
        double __elastic_weight_pion;
        bool __use_tracking_threshold;
        
      private:

        double getGenieBEinMeV( const int A, const int Z = 0 );

        int GetQEElasticWeight(double &weightElastic1, double &weightElastic2, int fatecodePN, double inputKEgev, int inputA, int pdg);

        int GetQEElasticFSIWeightFromParam(double &weight1, double &weight2, double inputKEgev, int inputA, int pdg);


        double m_breakpointlow; // These are in GeV
        double m_breakpointhigh; 
        double m_highcutoff;
        
        double fAbsorptionWeight;
        double fElasticWeight1;
        double fElasticWeight2;

        Long64_t fLastEntry;

      

  };

}

#endif

