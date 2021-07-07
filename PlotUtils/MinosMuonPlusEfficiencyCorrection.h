#ifndef MinosMuonPlusEfficiencyCorrection_h
#define MinosMuonPlusEfficiencyCorrection_h


class TF1;

namespace PlotUtils {
    
        /// how to use:
        /// call anywhere inside your event loop
        /// for example p_mu (MINOS) = 3.0 GeV/c and numi_pot = 32e12 pot
        /// MinosMuonEfficiencyCorrection::Get().GetCorrection(3.0, 32.0);
        /// returns the correction to be applied to MC
    class MinosMuonPlusEfficiencyCorrection {
      public:
        
        static MinosMuonPlusEfficiencyCorrection& Get();

            // p_mu in GeV
            // numi_pot in unit of 1.e12
        double GetCorrection(double p_mu, double batch_pot);

            // muon theta in degree
            // see Anushree's DocDB 20760, v1, p. 32
            // actually it's only theta_mu dependent at the time this function is implemented
        double GetCorrectionErr(double p_mu, double theta_mu, double batch_pot);
        
      private:
        MinosMuonPlusEfficiencyCorrection();
        MinosMuonPlusEfficiencyCorrection(const MinosMuonPlusEfficiencyCorrection&) {}
        ~MinosMuonPlusEfficiencyCorrection();

            // solve for A and B, correction = 1 + A * intensity + B * intensity^2
        void solve_eq(double a1, double b1,double c1,
                      double a2, double b2,double c2,
                      double& A, double& B);
        
        TF1* __tf1_correction_curve_lo;
        TF1* __tf1_correction_curve_hi;
        TF1* __tf1_correction_curve_err; // muon theta dependent
        
        double __pmin;
        double __pmax;
        double __theta_min;
        double __theta_max;
    };

};

#endif
