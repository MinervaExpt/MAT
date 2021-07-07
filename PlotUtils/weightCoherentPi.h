#ifndef weightCoherentPi_h
#define weightCoherentPi_h

#include <TFile.h>
#include <TH1D.h>
#include <TString.h>

#include <cassert>
#include <iostream>
#include <stdexcept>

namespace PlotUtils {

class weightCoherentPi {
 public:
  weightCoherentPi();
  // explicit coherent_reweight(const std::string& coherent_reweight_filename);

  double get_pion_energy_weight(const double epi_in_GeV);
  double get_pion_theta_weight(const double theta_in_degree);
  double get_combined_weight(const double theta_in_degree,
                             const double epi_in_GeV);

  double get_pion_energy_weight_error(const double epi);
  double get_pion_theta_weight_error(const double theta_in_degree);

  ~weightCoherentPi();

 private:
  TFile* weights_file;

  TH1D* get_bin_error_as_histogram(const TH1D* input);

  TH1D* __h1d_energy_cv;
  TH1D* __h1d_theta_cv;

  TH1D* __h1d_energy_err;
  TH1D* __h1d_theta_err;

  double __epi_max;
  double __theta_max;

  void read(const std::string f);
};

PlotUtils::weightCoherentPi& weight_coherent();
}  // namespace PlotUtils

#endif
