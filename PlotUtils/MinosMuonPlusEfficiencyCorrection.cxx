//This code is uploaded on April 8 to cvs.
//I have used ME5A and ME6B antineutrino playlists.
//I will do it again once a low intensity antinuetrino playlist is available. 

#include <cmath>
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <limits>

#include <TF1.h>

#include "MinosMuonPlusEfficiencyCorrection.h"


using std::sqrt;
using std::pow;

PlotUtils::MinosMuonPlusEfficiencyCorrection::MinosMuonPlusEfficiencyCorrection()
{
    __pmin = 1.0;
    __pmax = 4.0;
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

PlotUtils::MinosMuonPlusEfficiencyCorrection::~MinosMuonPlusEfficiencyCorrection() {}

PlotUtils::MinosMuonPlusEfficiencyCorrection& PlotUtils::MinosMuonPlusEfficiencyCorrection::Get()
{
    static MinosMuonPlusEfficiencyCorrection singleton;
    
    return singleton;
}

double PlotUtils::MinosMuonPlusEfficiencyCorrection::GetCorrection(double p_mu, double pot)
{

    static const double pot_lo  = 6.74048;  //5A
    static const double pot_lo2 = pot_lo * pot_lo;
    static const double pot_hi  = 8.26951;  //6B
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

double PlotUtils::MinosMuonPlusEfficiencyCorrection::GetCorrectionErr(double, double theta_mu, double)
{
    theta_mu = std::min(theta_mu, __theta_max);
    theta_mu = std::max(theta_mu, __theta_min);
    
    double overall = 0.01;
    double theta_dependent = 1.0 - __tf1_correction_curve_err->Eval(theta_mu);
    double total_err = sqrt(pow(overall,2.0) + pow(theta_dependent,2.0));
    
    return total_err;
}



void PlotUtils::MinosMuonPlusEfficiencyCorrection::solve_eq(double a1, double b1,double c1,
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



