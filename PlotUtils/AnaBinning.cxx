#ifndef MNV_ANABINNING_cxx
#define MNV_ANABINNING_cxx 1

#include "AnaBinning.h"
#include <iostream>

using namespace PlotUtils;

AnaBinning& AnaBinning::Get(){
  static AnaBinning singleton;
  return singleton;
}

//==================================================
axis_binning AnaBinning::GetProtonEnergyBinsGeV( )const{
  
  axis_binning tmp;

  std::vector<double> bins;
  int nbins = 50;
  double min = 0, max = 2.0;  // GeV
  double binw = (max-min)/(double)nbins;
  for( int i = 0; i <= nbins; i++ ){
    bins.push_back( min + i*binw );
  }
  tmp.bin_edges = bins;
  tmp.nbins = nbins;
  tmp.min = min;
  tmp.max = max;
  tmp.uniform = true;
  tmp.default_width = binw;

  return tmp;
}
//==================================================
axis_binning AnaBinning::GetBlobEnergyBinsGeV( )const{
  
  axis_binning tmp;

  std::vector<double> bins;
  int nbins = 100;
  double min = 0, max = 2.0;  // GeV
  double binw = (max-min)/(double)nbins;
  for( int i = 0; i <= nbins; i++ ){
    bins.push_back( min + i*binw );
  }
  tmp.bin_edges = bins;
  tmp.nbins = nbins;
  tmp.min = min;
  tmp.max = max;
  tmp.uniform = true;
  tmp.default_width = binw;

  return tmp;
}
//==================================================
axis_binning AnaBinning::GetBlobEnergyBinsMeV( )const{

  axis_binning tmp;

  std::vector<double> bins;
  int nbins = 100;
  double min = 0, max = 500; // MeV
  double binw = (max-min)/(double)nbins;

  for( int i = 0; i <= nbins; i++ ){
    bins.push_back( min + i*binw );
  }
  tmp.bin_edges = bins;
  tmp.nbins = nbins;
  tmp.min = min;
  tmp.max = max;
  tmp.uniform = true;
  tmp.default_width = binw;

  return tmp;
}
//==================================================
axis_binning AnaBinning::GetQ2BinsGeV( )const{

  axis_binning tmp;
  tmp.default_width = 0.025;
  tmp.uniform = false;

  std::vector<double> bins;

  bins.push_back( 0.000 );
  bins.push_back( 0.025 );
  bins.push_back( 0.050 );
  bins.push_back( 0.1 );
  bins.push_back( 0.2 );
  bins.push_back( 0.4 );
  bins.push_back( 0.8 );
  bins.push_back( 1.2 );
  bins.push_back( 2.0 );

  tmp.bin_edges = bins;
  tmp.nbins = bins.size()-1;
  tmp.min = bins.front();
  tmp.max = bins.back();

  return tmp;
}
//==================================================
axis_binning AnaBinning::GetQ2UniformBinsGeV( )const{

  axis_binning tmp;

  std::vector<double> bins;
  int nbins = 80;
  double min = 0.0, max = 2.0;  // GeV^2
  double binw = (max-min)/(double)nbins;

  for( int i = 0; i <= nbins; i++ ){
    bins.push_back( min + i*binw );
  }
  tmp.bin_edges = bins;
  tmp.nbins = nbins;
  tmp.min = min;
  tmp.max = max;
  tmp.uniform = true;
  tmp.default_width = binw;

  return tmp;
}
//==================================================
axis_binning AnaBinning::GetQ2FineBinsGeV( )const{

  axis_binning tmp;
  tmp.default_width = 0.025;
  tmp.uniform = false;

  std::vector<double> bins;

  bins.push_back(0.0);
  for( int i = 1; i <= 8; i++ ){
    bins.push_back( bins.back() + tmp.default_width );
  }
  for( int i = 1; i <= 4; i++ ){
    bins.push_back( bins.back() + 2*tmp.default_width );
  }
  for( int i = 1; i <= 6; i++ ){
    bins.push_back( bins.back() + 4*tmp.default_width );
  }
  tmp.bin_edges = bins;
  tmp.nbins = bins.size()-1;
  tmp.min = bins.front();
  tmp.max = bins.back();

  return tmp;
}
//==================================================
axis_binning AnaBinning::GetQ2FineBinsGeV_05( )const{

  axis_binning tmp;
  tmp.default_width = 0.05;
  tmp.uniform = false;

  std::vector<double> bins;

  bins.push_back(0.0);
  for( int i = 1; i <= 12; i++ ){
    bins.push_back( bins.back() + tmp.default_width );
  }
  for( int i = 1; i <= 4; i++ ){
    bins.push_back( bins.back() + 2*tmp.default_width );
  }
  tmp.bin_edges = bins;
  tmp.nbins = bins.size()-1;
  tmp.min = bins.front();
  tmp.max = bins.back();

  return tmp;
}
//==================================================
axis_binning AnaBinning::GetNeutrinoEnergyBinsGeV( )const{

  return this->GetMuonEnergyBinsGeV();
}
//==================================================
axis_binning AnaBinning::GetNeutrinoEnergyUniformBinsGeV( )const{

  return this->GetMuonEnergyUniformBinsGeV();
}
//==================================================
axis_binning AnaBinning::GetMuonEnergyBinsGeV( )const{

  axis_binning tmp;
  tmp.default_width = 0.500;
  tmp.uniform = false;

  std::vector<double> bins;

  bins.push_back(0.0);
  for( int i = 1; i <= 8; i++ ){
    bins.push_back( bins.back() + tmp.default_width );
  }
  bins.push_back( bins.back() + 2*tmp.default_width );
  bins.push_back( bins.back() + 2*tmp.default_width );
  bins.push_back( bins.back() + 2*tmp.default_width );
  bins.push_back( bins.back() + 3*tmp.default_width );
  bins.push_back( bins.back() + 3*tmp.default_width );

  tmp.bin_edges = bins;
  tmp.nbins = bins.size()-1;
  tmp.min = bins.front();
  tmp.max = bins.back();

  return tmp;
}
//==================================================
axis_binning AnaBinning::GetMuonEnergyUniformBinsGeV( )const{

  axis_binning tmp;

  std::vector<double> bins;
  int nbins = 50;
  double min = 0.0, max = 10.0;  // GeV 
  double binw = (max-min)/(double)nbins;

  for( int i = 0; i <= nbins; i++ ){
    bins.push_back( min + i*binw );
  }
  tmp.bin_edges = bins;
  tmp.nbins = nbins;
  tmp.min = min;
  tmp.max = max;
  tmp.uniform = true;
  tmp.default_width = binw;

  return tmp;
}
//==================================================
axis_binning AnaBinning::GetMuonAngleBinsDegrees( )const{

  axis_binning tmp;

  std::vector<double> bins;
  int nbins = 50;
  double min = 0, max = 25.0;  // degrees
  double binw = (max-min)/(double)nbins;

  for( int i = 0; i <= nbins; i++ ){
    bins.push_back( min + i*binw );
  }
  tmp.bin_edges = bins;
  tmp.nbins = nbins;
  tmp.min = min;
  tmp.max = max;
  tmp.uniform = true;
  tmp.default_width = binw;

  return tmp;
}
//==================================================
axis_binning AnaBinning::GetMuonCosineAngleBins( )const{

  axis_binning tmp;

  std::vector<double> bins;
  int nbins = 50;
  double min = 0.9, max = 1.0;  // degrees
  double binw = (max-min)/(double)nbins;

  for( int i = 0; i <= nbins; i++ ){
    bins.push_back( min + i*binw );
  }
  tmp.bin_edges = bins;
  tmp.nbins = nbins;
  tmp.min = min;
  tmp.max = max;
  tmp.uniform = true;
  tmp.default_width = binw;

  return tmp;
}
//==================================================
axis_binning AnaBinning::GetThetaXBinsDegrees( )const{

  axis_binning tmp;

  std::vector<double> bins;
  int nbins = 50;
  double min = 0.0, max = 25.0;  // degrees
  double binw = (max-min)/(double)nbins;

  for( int i = 0; i <= nbins; i++ ){
    bins.push_back( min + i*binw );
  }
  tmp.bin_edges = bins;
  tmp.nbins = nbins;
  tmp.min = min;
  tmp.max = max;
  tmp.uniform = true;
  tmp.default_width = binw;

  return tmp;
}
//==================================================
axis_binning AnaBinning::GetThetaYBinsDegrees( )const{

  axis_binning tmp;

  std::vector<double> bins;
  int nbins = 50;
  double min = 0.0, max = 25.0;  // degrees
  double binw = (max-min)/(double)nbins;

  for( int i = 0; i <= nbins; i++ ){
    bins.push_back( min + i*binw );
  }
  tmp.bin_edges = bins;
  tmp.nbins = nbins;
  tmp.min = min;
  tmp.max = max;
  tmp.uniform = true;
  tmp.default_width = binw;

  return tmp;
}
//==================================================
axis_binning AnaBinning::GetPhiBinsDegrees( )const{

  axis_binning tmp;

  std::vector<double> bins;
  int nbins = 45;
  double min = -180.0, max = 180.0;  // degrees
  double binw = (max-min)/(double)nbins;

  for( int i = 0; i <= nbins; i++ ){
    bins.push_back( min + i*binw );
  }
  tmp.bin_edges = bins;
  tmp.nbins = nbins;
  tmp.min = min;
  tmp.max = max;
  tmp.uniform = true;
  tmp.default_width = binw;

  return tmp;
}
//==================================================
axis_binning AnaBinning::GetVertexR2Bins( )const{

  axis_binning tmp;

  std::vector<double> bins;
  int nbins = 20;
  double min = 0, max = 10000.0;  // cm^2
  double binw = (max-min)/(double)nbins;

  for( int i = 0; i <= nbins; i++ ){
    bins.push_back( min + i*binw );
  }
  tmp.bin_edges = bins;
  tmp.nbins = nbins;
  tmp.min = min;
  tmp.max = max;
  tmp.uniform = true;
  tmp.default_width = binw;

  return tmp;
}
//==================================================
axis_binning AnaBinning::GetVertexXYBins( )const{

  axis_binning tmp;

  std::vector<double> bins;
  int nbins = 40;
  double min = -100.0, max = 100.0;  // cm 
  double binw = (max-min)/(double)nbins;

  for( int i = 0; i <= nbins; i++ ){
    bins.push_back( min + i*binw );
  }
  tmp.bin_edges = bins;
  tmp.nbins = nbins;
  tmp.min = min;
  tmp.max = max;
  tmp.uniform = true;
  tmp.default_width = binw;

  return tmp;
}
//==================================================
axis_binning AnaBinning::GetVertexZBins( )const{

  axis_binning tmp;

  std::vector<double> bins;
  int nbins = 100;
  double min = 600.0, max = 1100.0;  // cm 
  double binw = (max-min)/(double)nbins;

  for( int i = 0; i <= nbins; i++ ){
    bins.push_back( min + i*binw );
  }
  tmp.bin_edges = bins;
  tmp.nbins = nbins;
  tmp.min = min;
  tmp.max = max;
  tmp.uniform = true;
  tmp.default_width = binw;

  return tmp;
}
//==================================================
axis_binning AnaBinning::GetVertexZBinsModule( )const{

  axis_binning tmp;

  double fullDetectorBins[] = {/*
	400.00,
	405.00,
	410.00,
	415.00,
	420.00,
	424.11224,
	428.22447,
	432.64595,
	437.06743,
	441.48891,
	445.91039,
	450.33188,
	454.75336,
	459.17484,
	463.59632,
	468.01780,
	472.43928,
	476.86077,
	481.28225,
	485.70373,
	490.12521,
	494.54669,
	498.96817, 
	503.38965,
	507.81114,
	512.23262,
	516.44336,
	544.59410,
	549.01558,*/
	553.43706,
	557.85854,
	562.28002,
	566.70151,
	571.12299,
	575.54447,
	579.91524,
	584.43814,
	588.96104,
	593.48394,
	598.00684,
	602.52974,
	607.05264,
	611.57554,
	616.09844,
	620.62134,
	625.14424,
	629.66714,
	634.19004,
	638.71294,
	643.23584,
	647.75874,
	652.28164,
	656.80454,
	661.32744,
	665.85034,
	670.37324,
	674.89614,
	679.41904,
	683.94194,
	688.46484,
	692.98774,
	697.51064,
	702.03354,
	706.55644,
	711.07934,
	715.60224,
	720.12514,
	724.64804,
	729.17094,
	733.69384,
	738.21674,
	742.73964,
	747.26254,
	751.78544,
	756.30834,
	760.83124,
	765.35414,
	769.87704,
	774.39994,
	778.92284,
	783.44574,
	787.96864,
	792.49154,
	797.01444,
	801.53734,
	806.06024,
	810.58314,
	815.10604,
	819.62894,
	824.15184,
	828.67474,
	833.19764,
	837.72054,
	842.24344,
	846.76634,
	851.28924,
	855.81214,
	860.35635,
	864.83663,
	869.31691,
	873.79718,
	878.27746,
	882.75774,
	887.23802,
	891.71830,
	896.19857/*,
	900.67885,
	905.17144,
	909.90470,
	914.63795,
	919.37121,
	924.10446,
	928.83772,
	933.57098,
	938.30423,
	943.03749,
	947.77075,
	952.50400,
	957.23726,
	961.97052,
	966.70377,
	971.43703,
	976.17029,
	980.90354,
	985.63680,
	990.37006,
	995.10331,
	10000.0,
	10050.0,
	10100.0,
	10150.0,
	10200.0,
	10250.0,
	10300.0,
	10350.0,
	10400.0,
	10450.0,
	10500.0,
	10550.0,
	10600.0,
	10650.0,
	10700.0,
	10750.0,
	10800.0,
	10850.0,
	10900.0,
	10950.0,
	11000 */};

  std::vector<double> bins( fullDetectorBins, fullDetectorBins -1 + sizeof(fullDetectorBins)/sizeof(double) );
  int nbins = bins.size()-1;
  double min = bins.front(), max = bins.back();
  double binw = 5.0;

  tmp.bin_edges = bins;
  tmp.nbins = nbins;
  tmp.min = min;
  tmp.max = max;
  tmp.uniform = false;
  tmp.default_width = binw;

  return tmp;
}
//==================================================
axis_binning AnaBinning::GetVertexTBins( )const{

  axis_binning tmp;

  std::vector<double> bins;
  int nbins = 100;
  double min = 0.0, max = 20.0;  // us 
  double binw = (max-min)/(double)nbins;

  for( int i = 0; i <= nbins; i++ ){
    bins.push_back( min + i*binw );
  }
  tmp.bin_edges = bins;
  tmp.nbins = nbins;
  tmp.min = min;
  tmp.max = max;
  tmp.uniform = true;
  tmp.default_width = binw;

  return tmp;
}
//==================================================
axis_binning AnaBinning::GetVertexTBinsZoom( )const{

  axis_binning tmp;

  std::vector<double> bins;
  int nbins = 100;
  double min = 0.0, max = 10.0;  // us 
  double binw = (max-min)/(double)nbins;

  for( int i = 0; i <= nbins; i++ ){
    bins.push_back( min + i*binw );
  }
  tmp.bin_edges = bins;
  tmp.nbins = nbins;
  tmp.min = min;
  tmp.max = max;
  tmp.uniform = true;
  tmp.default_width = binw;

  return tmp;
}
//==================================================
axis_binning AnaBinning::GetMinosVertexXYBins( )const{

  axis_binning tmp;

  std::vector<double> bins;
  int nbins = 40;
  double min = -200.0, max = 200.0;  // cm 
  double binw = (max-min)/(double)nbins;

  for( int i = 0; i <= nbins; i++ ){
    bins.push_back( min + i*binw );
  }
  tmp.bin_edges = bins;
  tmp.nbins = nbins;
  tmp.min = min;
  tmp.max = max;
  tmp.uniform = true;
  tmp.default_width = binw;

  return tmp;
}
//==================================================
axis_binning AnaBinning::GetMinosQPBins( )const{

  axis_binning tmp;

  std::vector<double> bins;
  int nbins = 100;
  double min = -2.0, max = 2.0;
  double binw = (max-min)/(double)nbins;

  for( int i = 0; i <= nbins; i++ ){
    bins.push_back( min + i*binw );
  }
  tmp.bin_edges = bins;
  tmp.nbins = nbins;
  tmp.min = min;
  tmp.max = max;
  tmp.uniform = true;
  tmp.default_width = binw;

  return tmp;
}
//==================================================
void AnaBinning::AddBins( std::vector<double>& bin_edges, const int n_bins, const double bin_width ) const{
  //! If the bins have no defined low entry, assume it should be 0
  double x = 0.;
  if( bin_edges.size() > 0 )
    x = bin_edges[ bin_edges.size() - 1 ];
  else
    bin_edges.push_back( x );

  //! add n_bins spaced by this width
  for( int iBin = 0; iBin < n_bins; ++iBin ) {
    x += bin_width;
    bin_edges.push_back( x );
  }
}
//==================================================
void AnaBinning::AddBins( std::vector<double>& bin_edges, const double bin_width, const double low_bin, const double high_bin ) const{
  double x = low_bin;
  if( bin_edges.size() && low_bin == bin_edges.back() )
    x += bin_width;
  while( x <= high_bin ) {
    bin_edges.push_back( x );
    x += bin_width;
  }
}

#endif


