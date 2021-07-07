#include "weightCoherentPi.h"

using namespace PlotUtils;

weightCoherentPi::weightCoherentPi()
    : weights_file(nullptr),
      __h1d_energy_cv(NULL),
      __h1d_theta_cv(NULL),
      __h1d_energy_err(NULL),
      __h1d_theta_err(NULL),
      __epi_max(1e20),
      __theta_max(1e20) {
  char* mparalocation = std::getenv("MPARAMFILESROOT");
  std::string f = std::string(mparalocation) +
                  "/data/Reweight/cohPi_weights_tracker_ME.root";
  read(f);
}

void weightCoherentPi::read(const std::string f) {
  weights_file = TFile::Open(f.c_str(), "READONLY");
  //    TFile* f = new TFile(coherent_reweight_filename.c_str(), "read");
  assert(weights_file);

  __h1d_energy_cv = static_cast<TH1D*>(weights_file->Get("h1d_epi_xsec_ratio"));
  __h1d_theta_cv =
      static_cast<TH1D*>(weights_file->Get("h1d_theta_xsec_ratio"));

  __h1d_energy_err = get_bin_error_as_histogram(__h1d_energy_cv);
  __h1d_theta_err = get_bin_error_as_histogram(__h1d_theta_cv);

  __epi_max = __h1d_energy_cv->GetXaxis()->GetXmax();
  __theta_max = __h1d_theta_cv->GetXaxis()->GetXmax();

  // printf("q2_max epi_max theta_max %10.4f %10.4f %10.4f\n", __q2_max,
  // __epi_max, __theta_max);
}

weightCoherentPi::~weightCoherentPi() {
  delete __h1d_energy_cv;
  delete __h1d_theta_cv;

  delete __h1d_energy_err;
  delete __h1d_theta_err;
}

TH1D* weightCoherentPi::get_bin_error_as_histogram(const TH1D* input) {
  TH1D* output =
      static_cast<TH1D*>(input->Clone(Form("%s_error", input->GetName())));
  output->Reset();
  for (int bin = 1; bin <= input->GetNbinsX(); ++bin) {
    double weight_error = input->GetBinError(bin);
    output->SetBinContent(bin, weight_error);
  }

  return output;
}

// coherent pion weights
double weightCoherentPi::get_pion_energy_weight(double epi) {
  if (epi > __epi_max) return 1.0;

  return __h1d_energy_cv->GetBinContent(__h1d_energy_cv->FindBin(epi));
}

double weightCoherentPi::get_pion_theta_weight(double theta) {
  if (theta > __theta_max) return 1.0;

  return __h1d_theta_cv->GetBinContent(__h1d_theta_cv->FindBin(theta));
}

double weightCoherentPi::get_combined_weight(double theta, double epi) {
  if (theta > __theta_max && epi > __epi_max) return 1.0;
  if (theta > __theta_max)
    return __h1d_energy_cv->GetBinContent(__h1d_energy_cv->FindBin(epi));
  if (epi > __epi_max)
    return __h1d_theta_cv->GetBinContent(__h1d_theta_cv->FindBin(theta));
  return __h1d_theta_cv->GetBinContent(__h1d_theta_cv->FindBin(theta)) *
         __h1d_energy_cv->GetBinContent(__h1d_energy_cv->FindBin(epi));
}

// coherent pion weight errors
double weightCoherentPi::get_pion_energy_weight_error(double epi) {
  if (epi > __epi_max) return 0.0;

  return __h1d_energy_err->GetBinContent(__h1d_energy_cv->FindBin(epi));
}

double weightCoherentPi::get_pion_theta_weight_error(double theta) {
  if (theta > __theta_max) return 0.0;

  return __h1d_theta_err->GetBinContent(__h1d_theta_err->FindBin(theta));
}

PlotUtils::weightCoherentPi& PlotUtils::weight_coherent() {
  static PlotUtils::weightCoherentPi* _weight_coherent =
      new PlotUtils::weightCoherentPi();
  return *_weight_coherent;
}
