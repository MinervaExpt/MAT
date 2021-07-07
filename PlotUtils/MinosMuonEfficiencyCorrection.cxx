#include <cmath>
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <limits>

#include <TF1.h>

#include "MinosMuonEfficiencyCorrection.h"


using std::sqrt;
using std::pow;

PlotUtils::MinosMuonEfficiencyCorrection::MinosMuonEfficiencyCorrection(bool isNeutrino /*=true*/)
{
    __pmin = 1.0;
    __pmax = 4.0;


  if(isNeutrino){//neutrino

        // low-intensity correction curve
    double p0  =  0.705549;//   +/-   0.0292251;
    double p1  = 0.499846;//   +/-   0.0684638;
    double p2  = -0.364875;//   +/-   0.0600353;
    double p3 = 0.131063;//   +/-   0.0248497;
    double p4  = -0.0228281;//   +/-   0.00489447;
    double p5  = 0.00154069;//   +/-   0.000369448;

    __tf1_correction_curve_lo = new TF1("correction_low",
                                        "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x",
                                        __pmin, __pmax);
    
    __tf1_correction_curve_lo->SetParameter(0, p0);
    __tf1_correction_curve_lo->SetParameter(1, p1);
    __tf1_correction_curve_lo->SetParameter(2, p2);
    __tf1_correction_curve_lo->SetParameter(3, p3);
    __tf1_correction_curve_lo->SetParameter(4, p4);
    __tf1_correction_curve_lo->SetParameter(5, p5);

        // high-intensity correction curve
   
    p0  =     0.959697;  //   +/-   0.0948765   
    p1  =    -0.139194;  // +/-   0.219377    
    p2  =     0.175986;  // +/-   0.189988    
    p3  =   -0.0802815;  // +/-   0.0777433   
    p4  =    0.0162377;  // +/-   0.0151462   
    p5  =  -0.00122259; //   +/-   0.00113086   
    __tf1_correction_curve_hi = new TF1("correction_high",
                                        "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x",
                                        __pmin, __pmax);
    
    __tf1_correction_curve_hi->SetParameter(0, p0);
    __tf1_correction_curve_hi->SetParameter(1, p1);
    __tf1_correction_curve_hi->SetParameter(2, p2);
    __tf1_correction_curve_hi->SetParameter(3, p3);
    __tf1_correction_curve_hi->SetParameter(4, p4);
    __tf1_correction_curve_hi->SetParameter(5, p5);
    
  }
  else{//anti-neutrino

        // low-intensity correction curve
        double p0  = 1.06348;
	double p1  = 0.49754;
	double p2  = -1.07352;
	double p3  = 0.810622;

	    __tf1_correction_curve_lo = new TF1("correction_low",
                                      "[0] - [1]/x - [2]/(x*x) - [3]/(x*x*x) ",
                                      __pmin, __pmax);
   
    __tf1_correction_curve_lo->SetParameter(0, p0);
    __tf1_correction_curve_lo->SetParameter(1, p1);
    __tf1_correction_curve_lo->SetParameter(2, p2);
    __tf1_correction_curve_lo->SetParameter(3, p3);

        // high-intensity correction curve
 	p0 = 1.04643;
	p1 = 0.37246;
	p2 = -0.75515;
	p3 = 0.571880;
  
          __tf1_correction_curve_hi = new TF1("correction_high",
                                        "[0] - [1]/x - [2]/(x*x) - [3]/(x*x*x)",
                                        __pmin, __pmax);
 
    __tf1_correction_curve_hi->SetParameter(0, p0);
    __tf1_correction_curve_hi->SetParameter(1, p1);
    __tf1_correction_curve_hi->SetParameter(2, p2);
    __tf1_correction_curve_hi->SetParameter(3, p3);

  }

      // systematic uncertainty as function of muon theta
    __theta_min = 0.0;
    __theta_max = 40.0;

    double a =  0.99079;
    double b = -4.24729e-4;
    double c =  7.20874e-5;
    
    __tf1_correction_curve_err = new TF1("correction_err",
                                         "[0] - [1] * x - [2] * x * x",
                                         __theta_min, __theta_max);

    __tf1_correction_curve_err->SetParameter(0, a);
    __tf1_correction_curve_err->SetParameter(1, b);
    __tf1_correction_curve_err->SetParameter(2, c);
    
}

PlotUtils::MinosMuonEfficiencyCorrection::~MinosMuonEfficiencyCorrection() {}

PlotUtils::MinosMuonEfficiencyCorrection& PlotUtils::MinosMuonEfficiencyCorrection::Get(bool isNeutrino /* true*/)
{
    static MinosMuonEfficiencyCorrection singleton(isNeutrino);
    
    return singleton;
}

double PlotUtils::MinosMuonEfficiencyCorrection::GetCorrection(double p_mu /*GeV*/, double pot, bool isNeutrino /* true*/)
{

static double pot_lo;
static double pot_hi;

  if(isNeutrino){//neutrino
     pot_lo  = 3.9385;
     pot_hi  = 8.0311;
  }
  else{//antineutrino
     pot_lo  = 6.74048;  //5A
     pot_hi  = 8.26951;  //6B
  }

    static const double pot_lo2 = pot_lo * pot_lo;
    static const double pot_hi2 = pot_hi * pot_hi;

    p_mu = std::min(p_mu, __pmax);
    p_mu = std::max(p_mu, __pmin);

    double corr_lo = __tf1_correction_curve_lo->Eval(p_mu);
    double corr_hi = __tf1_correction_curve_hi->Eval(p_mu);



    double A = 0.0;
    double B = 0.0;
    solve_eq(pot_hi, pot_hi2, corr_hi, pot_lo, pot_lo2, corr_lo, A, B);
    double corr = 1.0 + A * pot + B * pot * pot;


    return corr;

}

double PlotUtils::MinosMuonEfficiencyCorrection::GetCorrectionErr(double, double theta_mu /*degrees*/, double, bool isNeutrino)
{
    theta_mu = std::min(theta_mu, __theta_max);
    theta_mu = std::max(theta_mu, __theta_min);
    
    double overall = 0.01;
    
    double theta_dependent = 1.0 - __tf1_correction_curve_err->Eval(theta_mu);
    double total_err = pow(overall,2.0);
    if (isNeutrino) total_err += pow(theta_dependent,2.0);

    total_err = sqrt(total_err);
    
    return total_err;
}



void PlotUtils::MinosMuonEfficiencyCorrection::solve_eq(double a1, double b1,double c1,
                                                        double a2, double b2,double c2,
                                                        double& A, double& B)
{
    c1 += -1.0;
    c2 += -1.0;
    
    if (std::abs((a1 * b2) - (b1 * a2)) < std::numeric_limits<double>::epsilon()) {
        
        std::cerr << "The system has no solution." << std::endl;
        return;
        
    } else {
        
        A = ((c1*b2) - (b1*c2))/((a1*b2) - (b1*a2));
        B = ((a1*c2) - (c1*a2))/((a1*b2) - (b1*a2));
    }
}



