#include "weightLowQ2Pi.h"

using namespace PlotUtils;

weightLowQ2Pi::weightLowQ2Pi() 
  : weights_file(nullptr),
    h_cvweight_JOINT(nullptr),
    h_cvweight_shift_JOINT(nullptr),
    h_cvweight_NU1PI(nullptr),
    h_cvweight_shift_NU1PI(nullptr),
    h_cvweight_NUNPI(nullptr),
    h_cvweight_shift_NUNPI(nullptr),
    h_cvweight_NUPI0(nullptr),
    h_cvweight_shift_NUPI0(nullptr),
    h_cvweight_NUBARPI0(nullptr),
    h_cvweight_shift_NUBARPI0(nullptr)
{
  char *mparalocation = std::getenv("MPARAMFILESROOT");
  std::string f = std::string(mparalocation) + "/data/Reweight/lowQ2pi_weights.root";
  read(f);
}


void weightLowQ2Pi::read(const std::string f)
{
  weights_file = TFile::Open(f.c_str(),"READONLY");

  if (weights_file){
    h_cvweight_JOINT          = (TGraph*)weights_file->Get("JOINT_weights");
    h_cvweight_shift_JOINT    = (TGraph*)weights_file->Get("JOINT_weight_shifts");
    h_cvweight_NU1PI          = (TGraph*)weights_file->Get("NU1PI_weights");
    h_cvweight_shift_NU1PI    = (TGraph*)weights_file->Get("NU1PI_weight_shifts");
    h_cvweight_NUNPI          = (TGraph*)weights_file->Get("NUNPI_weights");
    h_cvweight_shift_NUNPI    = (TGraph*)weights_file->Get("NUNPI_weight_shifts");
    h_cvweight_NUPI0          = (TGraph*)weights_file->Get("NUPI0_weights");
    h_cvweight_shift_NUPI0    = (TGraph*)weights_file->Get("NUPI0_weight_shifts");
    h_cvweight_NUBARPI0       = (TGraph*)weights_file->Get("NUBARPI0_weights");
    h_cvweight_shift_NUBARPI0 = (TGraph*)weights_file->Get("NUBARPI0_weight_shifts");
    std::cout << "Done readings weights and weight shifts from file " 
              << f << std::endl;
  }
  else
    throw std::runtime_error(std::string("weightLowQ2Pi::read: Could not read file ") + f);
}


double weightLowQ2Pi::getWeight(const double Q2 /*(GeV/c)^2*/, std::string channel, int shift)
{
  // shift
  //  0 : CV weight
  // -1 : weaker weighting universe
  // +1 : stronger weighting universe
  if (shift != -1 && shift != 0 && shift != 1)
    throw std::invalid_argument("weightLowQ2Pi::getWeight: Invalid shift. "
                                "Options are -1, 0, +1");

  double cvweight = 1.0;
  double cvweight_shift = 0.0;

  if (Q2 < 0.0 || Q2 >= 0.7) {
    cvweight = 1.0;
    cvweight_shift = 0.0;
  }
  else {
    if (channel == "MINOS"){
      cvweight = getMinosWeight(Q2);
    }
    else if (channel == "JOINT"){
      assert(h_cvweight_JOINT);
      assert(h_cvweight_shift_JOINT);
      cvweight = h_cvweight_JOINT->Eval(Q2);
      if (shift != 0) cvweight_shift = h_cvweight_shift_JOINT->Eval(Q2);
    }
    else if (channel == "NU1PI"){
      assert(h_cvweight_NU1PI);
      assert(h_cvweight_shift_NU1PI);
      cvweight = h_cvweight_NU1PI->Eval(Q2);
      if (shift != 0) cvweight_shift = h_cvweight_shift_NU1PI->Eval(Q2);
    }
    else if (channel == "NUNPI"){
      assert(h_cvweight_NUNPI);
      assert(h_cvweight_shift_NUNPI);
      cvweight = h_cvweight_NUNPI->Eval(Q2);
      if (shift != 0) cvweight_shift = h_cvweight_shift_NUNPI->Eval(Q2);
    }
    else if (channel == "NUPI0"){
      assert(h_cvweight_NUPI0);
      assert(h_cvweight_shift_NUPI0);
      cvweight = h_cvweight_NUPI0->Eval(Q2);
      if (shift != 0) cvweight_shift = h_cvweight_shift_NUPI0->Eval(Q2);
    }
    else if (channel == "NUBARPI0"){
      assert(h_cvweight_NUBARPI0);
      assert(h_cvweight_shift_NUBARPI0);
      cvweight = h_cvweight_NUBARPI0->Eval(Q2);
      if (shift != 0) cvweight_shift = h_cvweight_shift_NUBARPI0->Eval(Q2);
    }
    else
      throw std::invalid_argument( 
        std::string("weightLowQ2Pi:getWeight: Invalid channel ")+ channel +
        "Options are: MINOS, JOINT, NU1PI, NUNPI, NUPI0, NUBARPI0");
  }

  return cvweight + shift*cvweight_shift;
} 


double weightLowQ2Pi::getMinosWeight(const double Q2/*(GeV/c)^2*/ /*, int variation*/)
{
  const double Aminos         = 1.01;
  //const double Aminos_err     = 0.0;
  const double Q0minos        = 0.156;
  //const double Q0minos_err    = 0.0;
    
   if (Q2 < 0.0 || Q2 >= 0.7) return 1.0; //paranoia

   double denom = 1 + exp(1 - sqrt(Q2)/Q0minos);

   return Aminos/denom;
}

PlotUtils::weightLowQ2Pi& PlotUtils::weight_lowq2pi() {
  static PlotUtils::weightLowQ2Pi* _weight_lowq2pi = new PlotUtils::weightLowQ2Pi();
  return *_weight_lowq2pi;
}    
