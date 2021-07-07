#ifndef weight_fsi_absorption_h
#define weight_fsi_absorption_h

#include <Rtypes.h>


namespace PlotUtils {

    class ChainWrapper;
    
    double gethAFSIAbsorptionWeight(int mc_incoming,
                                    int mc_primaryLepton,
                                    int mc_charm,
                                    int mc_intType,
                                    int mc_targetA,
                                    int mc_targetZ,
                                    int mc_er_nPart,
                                    const int* mc_er_ID,
                                    const int* mc_er_status,
                                    const int* mc_er_FD,
                                    const int* mc_er_LD,
                                    const int* mc_er_mother,
                                    const double* mc_er_Px,
                                    const double* mc_er_Py,
                                    const double* mc_er_Pz,
                                    const double* mc_er_E,
                                    bool verbose = false);



    
    double gethAFSIAbsorptionWeight(PlotUtils::ChainWrapper* chw,
                                    Long64_t entry,
                                    bool verbose = false);
    
}


#endif
