
#ifndef weight_fsi_cai_h
#define weight_fsi_cai_h

#include <fstream>  //ifstream
#include <iostream> //cout
#include <map>
#include <vector>

#include <TString.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TFile.h>
#include <TMath.h>
#include "genie_particle.h"
#include "math.h"
#include "assert.h"
//For Compatibility with ROOT compiler
//uncomment the following:
//#define ROOT

/*!
  This is a class for taking a Tproton histogram and apply reweighing on the no FSI to mimic the effects of a corrected elastic FSI. The current parameter for Transverse Imbalance Study, Phys. Rev. Lett. 121, 022504 (2018)/arXiv:1805.05486 is 1.5
*/

namespace PlotUtils{

  class weight_fsi_cai
  {
    public:
      //Constructor: Read in params from a filename
      weight_fsi_cai(const TString  f, int finalstate_pdg = 2212, bool reweighQE = true, bool neutrinoMode = true, int config=2):theta_b(3.3402*TMath::DegToRad())
      { 
        read(f); //Read in params from file
        m_FSPDG = finalstate_pdg;
        m_reweighQE = reweighQE;
        m_config = config;
        m_neutrinoMode = neutrinoMode;
        m_breakpointlow = 0.110;  // These are in GeV
        m_breakpointhigh = 1.01;
        m_highcutoff = 3.8;
 
      }    
      
      TString  filename;
      TFile* fFSIWeight;
      TH1D*  hElaFSIWeight;
      TH1D*  hNoFSIWeight;


      //Initializer
      void read(const TString  f);

      void setConfig( int config ) { m_config = config; }

      int getQEFSIMode( int mc_intType, int mc_targetA, int mc_targetZ, 
          int mc_er_nPart, const int* mc_er_ID, const int* mc_er_status, 
          const int* mc_er_FD, const int* mc_er_LD, const int* mc_er_mother, 
          const double* mc_er_Px, const double* mc_er_Py, const double* mc_er_Pz,
          const double* mc_er_E, double &Ti);
      int getResFSIMode( int mc_intType, int mc_targetA, int mc_targetZ, 
          int mc_er_nPart, const int* mc_er_ID, const int* mc_er_status, 
          const int* mc_er_FD, const int* mc_er_LD, const int* mc_er_mother, 
          const double* mc_er_Px, const double* mc_er_Py, const double* mc_er_Pz,
          const double* mc_er_E, double &Ti);

      double getWeight( int mc_intType, int mc_targetA, int mc_targetZ, 
          int mc_er_nPart, const int* mc_er_ID, const int* mc_er_status, 
          const int* mc_er_FD, const int* mc_er_LD, const int* mc_er_mother, 
          const double* mc_er_Px, const double* mc_er_Py, const double* mc_er_Pz,
          const double* mc_er_E);

      double getWeightPoly( int mc_intType, int mc_targetA, int mc_targetZ, 
          int mc_er_nPart, int* mc_er_ID, int* mc_er_status, 
          int* mc_er_FD, int* mc_er_LD, int* mc_er_mother, 
          double* mc_er_Px, double* mc_er_Py, double* mc_er_Pz,
          double* mc_er_E);



      void setThetaBeam( double theta ) { theta_b = theta; };
      void setNeutrinoMode( bool neutrino = true ){ m_neutrinoMode = neutrino ; };

    private:
      double getGenieBE( const int A, const int Z = 0 );
      double getWeightInternal( double Ti, int mode ); //Kinetic Energy GeV, mode = {1:NoFsi,3:Ela,4:Inela}
      double* polySwitcher(double inputKEgev,  double* low, double *mid, double* high );
      //Kinetic Energy GeV, mode = {1:QE,3:Ela,4:Inela},config = {0: do nothing, 1:Lauren 1.5,2:ela->nofsi, 3: do nothing, 4:ela->otherFSI
      double getWeightInternalRik( double Ti, int A, int mode, int config );
      int m_FSPDG;
      bool m_reweighQE;
      bool m_neutrinoMode;
      int theta_b;
      int m_config;

      double m_breakpointlow, m_breakpointhigh, m_highcutoff;

      //=========Rik's weights ======================
      //config == 2

      double lowEpolyHe4[10] = { 39.32856672, -3930.858433, 198879.6178, -5992910.025, 114791528.5, -1431038143, 1.155973558e+10, -5.830506854e+10, 1.667805605e+11, -2.064812952e+11};
      double midEpolyHe4[10] = { 2.448600795, -21.2815073, 175.9270834, -837.1818884, 2452.148013, -4549.790084, 5362.000843, -3893.423638, 1589.436743, -279.1563134};
      double highEpolyHe4[10] = { 1.156044571, -0.03764619617, 0.005627980527, 0., 0., 0., 0., 0., 0., 0.};
 

      double lowEpolyC12[10] = { 17.82282668, 1551.604933, -170931.3066, 6480983.701, -134101728.6, 1696097963, -1.350258893e+10, 6.620396447e+10, -1.829733656e+11, 2.184157212e+11};
      double midEpolyC12[10] = { 4.174891296, -45.54402543, 345.7329055, -1484.999656, 3875.448189, -6356.522435, 6588.903725, -4194.602797, 1499.087099, -230.4838552};
      double highEpolyC12[10] = { 1.231471752, -0.05270919298, 0.008157195446, 0., 0., 0., 0., 0., 0., 0.};


      double lowEpolyO16[10] = {-11.4283, 4987.76, -340194., 1.1199e7, -2.16573e8, 2.63567e9, -2.04893e10, 9.89449e10, -2.70812e11, 3.21319e11};
      double midEpolyO16[10] = {4.85878, -55.0283, 418.162, -1798.99, 4700.67, -7716.98, 8001.46, -5089.56, 1814.09, -277.456};
      double highEpolyO16[10] = {1.26765, -0.0583833, 0.00853185, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0};


      double lowEpolySi28[10] = { 323.0786729, -3756.15957, -732485.8957, 37802444.09, -882957457.5, 1.192714772e+10, -9.873018114e+10, 4.952931632e+11, -1.385578606e+12, 1.661480791e+12};
      double midEpolySi28[10] = { 24.61820776, -403.5263055, 3189.486319, -14362.55838, 40075.08404, -71634.16657, 82127.95804, -58397.4893, 23440.19114, -4058.274728};
      double highEpolySi28[10] = { 1.40181183, -0.09856004905, 0.01487881423, 0., 0., 0., 0., 0., 0., 0.};


      double lowEpolyAr40[10] = {-12.4288, 46744.8, -2.47696e6, 5.40691e7, -5.84637e8, 2.5746e9, 7.87892e9, -1.43636e11, 6.19253e11, -9.44823e11};
      double midEpolyAr40[10] = {31.0973, -517.199, 4089.58, -18404.7, 51320.3, -91701.5, 105134., -74779.8, 3032.9, -5203.64};
      double highEpolyAr40[10] = {1.47462, -0.111806, 0.0166705, 0.0, 0.0,   0.0, 0.0, 0.0, 0.0, 0.0};

      double lowEpolyFe56[10] = {2533.41, -140121., 3.27464e6, -4.09799e7, 2.88262e8, -1.07774e9, 1.67053e9, 0.0, 0.0, 0.0};
      double midEpolyFe56[10] = {31.9074, -525.4, 4131.33, -18520.4, 51485.8, -91755.4, 104940., -74465.9, 29837.6, -5157.84};
      double highEpolyFe56[10] = {1.52685, -0.127437, 0.0198632, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0};


      double lowEpolyPb208[10] = {6301.34, -347755., 8.07485e6, -1.00173e8, 6.97538e8, -2.57953e9, 3.95315e9, 0.0, 0.0, 0.0};
      double midEpolyPb208[10] = {61.6161, -1022.96, 7948.2, -35212.3, 96827.2, -170874., 193711., -136368., 54248.1, -9316.09};
      double highEpolyPb208[10] = {1.90105, -0.195386, 0.0286129, 0.0, 0.0,   0.0, 0.0, 0.0, 0.0, 0.0};

      //======================= antinu
      double lowEpolyHe4Anu[10] = { 76.21450444, -7284.080281, 350409.6144, -10144233.29, 188295757, -2289715752, 1.813028138e+10, -8.996810885e+10, 2.539078252e+11, -3.108221822e+11};
      double midEpolyHe4Anu[10] = { 3.819189505, -41.52923432, 328.7530687, -1528.585143, 4443.389389, -8283.145704, 9883.531893, -7291.837115, 3026.956436, -540.240589};
      double highEpolyHe4Anu[10] = { 1.1301232, -0.02361999224, 0.003182730455, 0., 0., 0., 0., 0., 0., 0.};
  

      double lowEpolyC12Anu[10] = { 91.8661909, -3841.369929, 17803.80266, 2466426.335, -77777993.43, 1160008348, -1.006039598e+10, 5.189316382e+10, -1.481598665e+11, 1.807900697e+11};
      double midEpolyC12Anu[10] = { 6.87059066, -87.0000531, 658.0906236, -2865.934219, 7713.872241, -13257.88765, 14592.16131, -9963.08823, 3845.829084, -641.7424001};
      double highEpolyC12Anu[10] = { 1.191527119, -0.02739756607, 0.003201116244, 0., 0., 0., 0., 0., 0., 0.};

      // this looks radically different.  Huh?
      double lowEpolyO16Anu[10] = { -2.15186349, 6545.892859, -480405.8834, 16120949.72, -313111105.5, 3803439779, -2.941740429e+10, 1.410667677e+11, -3.829205886e+11, 4.502057086e+11};
      double midEpolyO16Anu[10] = { 7.669365414, -96.49518275, 722.5640045, -3118.543517, 8308.536397, -14117.09555, 15347.33815, -10345.77359, 3942.650035, -649.6469222};
      double highEpolyO16Anu[10] = { 1.2325431, -0.04131824065, 0.005846054663, 0., 0., 0., 0., 0., 0., 0.};

      double lowEpolySi28Anu[10] = { 1942.270031, -145015.1637, 4619879.234, -74613366.29, 537963557.4, 938586048, -4.898302179e+10, 3.844432317e+11, -1.353278306e+12, 1.872616614e+12};
      double midEpolySi28Anu[10] = { 45.15947021, -752.4355742, 5833.867552, -25789.37417, 70848.97048, -125102.6277, 142078.2779, -100275.8391, 40006.37836, -6891.096095};
      double highEpolySi28Anu[10] = { 1.342152187, -0.06350811765, 0.008383670794, 0., 0., 0., 0., 0., 0., 0.};
      
      double lowEpolyAr40Anu[10] = { 1675.010584, -72524.21889, 1329865.841, -13064326.4, 72102012.51, -211141714.4, 255512775, 0., 0., 0.};
      double midEpolyAr40Anu[10] = { 53.50240758, -895.3341982, 6935.402331, -30619.68803, 84012.55091, -148171.6227, 168093.1796, -118513.9921, 47235.85732, -8128.517134};
      double highEpolyAr40Anu[10] = { 1.419983761, -0.08493487878, 0.01244785896, 0., 0., 0., 0., 0., 0., 0.};
      
      double lowEpolyFe56Anu[10] = { 3241.107116, -165590.9076, 3582114.609, -41549932.77, 271131806.6, -940872093, 1354028645, 0., 0., 0.};
      double midEpolyFe56Anu[10] = { 54.64574278, -907.6702308, 7006.764962, -30871.98759, 84591.8376, -149056.734, 168992.8464, -119102.8865, 47461.72133, -8167.162769};
      double highEpolyFe56Anu[10] = { 1.443437845, -0.07152361183, 0.00858548564, 0., 0., 0., 0., 0., 0., 0.};

      double lowEpolyPb208Anu[10] = { 6635.553396, -337195.0226, 7238352.564, -83226870.84, 538080687.7, -1849719333, 2637236834, 0., 0., 0.};
      double midEpolyPb208Anu[10] = { 96.25604035, -1597.334355, 12200.22151, -53178.48721, 144204.542, -251602.1285, 282621.0265, -197471.025, 78062.22019, -13333.6168};
      double highEpolyPb208Anu[10] = { 1.783187249, -0.121317479, 0.01488332601, 0., 0., 0., 0., 0., 0., 0.};

      



      //config == 4
      //updated weights from Rik

      double lowEpolyHe4w2[10] = { 0.8525499509, 106.5613564, -1349.864879, -129938.4005, 5706909.418, -106576730.4, 1102103401, -6559583996, 2.108872155e+10, -2.843225144e+10};
      double midEpolyHe4w2[10] ={ 1.401592311, 0.5340768238, 35.9410343, -351.4320019, 1465.967144, -3424.247666, 4784.59239, -3969.621213, 1803.437515, -345.4102984};

      double highEpolyHe4w2[10] = { 1.158515008, 0.0001009351013, -0.0001265995844, 0., 0., 0., 0., 0., 0., 0.};

      
      double lowEpolyC12w2[10] = { -1.09785639, 353.8971638, -14438.15969, 249803.2065, -961138.7354, -33007787.3, 592230359.6, -4425023560, 1.622582081e+10, -2.391266333e+10};
      double midEpolyC12w2[10] = { 1.405735843, 0.5803084824, 33.33315803, -323.6899774, 1330.142548, -3060.151485, 4217.656957, -3458.103359, 1555.362503, -295.3732419};

      double highEpolyC12w2[10] = { 1.16047678, -0.002525807226, 0.0006168334479, 0., 0., 0., 0., 0., 0., 0.};

      
      double lowEpolyO16w2[10] = { -1.600476996, 435.9526759, -19971.39691, 452148.723, -5407872.958, 28338610.09, 57010988.68, -1558249122, 7621399671, -1.283154539e+10};
      double midEpolyO16w2[10] = { 1.448614827, -0.4064702399, 42.84475735, -374.3659837, 1495.641126, -3402.083242, 4662.460835, -3809.164928, 1708.297187, -323.5100567};
      double highEpolyO16w2[10] = { 1.154381145, 0.004104892683, -0.0008410639001, 0., 0., 0., 0., 0., 0., 0.};

      double lowEpolySi28w2[10] = { -11.14679259, 1619.43455, -82842.44543, 2334771.427, -40502433.64, 451579395.5, -3251861347, 1.464296207e+10, -3.753948905e+10, 4.186264037e+10};
      double midEpolySi28w2[10] = { 1.402604068, 0.6760161318, 32.45166596, -320.1710152, 1323.602905, -3055.847889, 4221.183186, -3465.956677, 1560.301404, -296.4792756};
      double highEpolySi28w2[10] = { 1.15626498, 0.001690036606, -0.0002929508535, 0., 0., 0., 0., 0., 0., 0.};
      
      double lowEpolyAr40w2[10] = { -15.91880593, 2230.471033, -116235.7268, 3358158.198, -59926068.7, 688940043.9, -5124126962, 2.385706678e+10, -6.326926353e+10, 7.298762608e+10};
      double midEpolyAr40w2[10] =  { 1.370127894, 1.510992183, 23.88831308, -274.2049331, 1180.476078, -2784.735678, 3904.85482, -3244.201353, 1474.709714, -282.5053889};
      double highEpolyAr40w2[10] = { 1.157107574, 0.001164884679, -0.0003047637304, 0., 0., 0., 0., 0., 0., 0.};

      double lowEpolyFe56w2[10] = { -17.55855486, 2416.340033, -125628.1801, 3636120.214, -65207995.22, 755278024.5, -5671137415, 2.669951459e+10, -7.169405909e+10, 8.382943525e+10};
      double midEpolyFe56w2[10] = { 1.364364871, 1.617991401, 22.78343352, -266.8142579, 1149.241143, -2702.282161, 3770.631447, -3114.17639, 1406.301429, -267.5046458};
      double highEpolyFe56w2[10] = { 1.158295471, -0.0006771545028, 0.0002835267823, 0., 0., 0., 0., 0., 0., 0.};

      double lowEpolyPb208w2[10] = { -22.77635755, 3056.796849, -158957.7289, 4605296.37, -82631008.84, 956843795.2, -7177117819, 3.372961824e+10, -9.035554312e+10, 1.053493431e+11};
      double midEpolyPb208w2[10] = { 1.380837506, 1.340129355, 24.57882428, -272.4572318, 1158.247826, -2708.591884, 3770.535156, -3112.007007, 1405.812079, -267.6767955};
      double highEpolyPb208w2[10] = { 1.156387959, 0.001577782994, -0.0003703094127, 0., 0., 0., 0., 0., 0., 0.};
      

  };

}

#endif

