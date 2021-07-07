#ifndef MNV_HYPERDIMLINEARIZER_h
#define MNV_HYPERDIMLINEARIZER_h 1
#include <vector>
#include <iostream>
#include <utility>
#include "MnvH2D.h"
namespace PlotUtils
{
class HyperDimLinearizer{

 public:
  HyperDimLinearizer(std::vector<std::vector<double> > input, int type);//constructor
  std::pair<int,int> GetBin(std::vector<double> values); // Template get bin for 2,3,4D cases
  std::vector<int> GetValues(int x);//
  std::vector<TH2D*> Get2DHistos(PlotUtils::MnvH2D* result, bool IncludeSys);// This is for type==0
  std::vector<PlotUtils::MnvH2D*> Get2DMnvHistos(PlotUtils::MnvH2D* result, bool IncludeSys);// This is for type==0
  TH2D* Get2DHisto(PlotUtils::MnvH1D* result, bool IncludeSys);// This is for type==1, 2D result only!!
  PlotUtils::MnvH2D* Get2DMnvHisto(PlotUtils::MnvH1D* result, bool IncludeSys);// This is for type==1, 2D result only!!
  void TestFunctionality();// a bunch of prints.
 private:
  int Get1DBin(double value,int el); //Get the bin number for one of the dimensions
  int m_analysis_type; // type 0 is 2D analysis, ignore Y coordinate since that is just used. type1 is full blown 1D.
  std::vector<int> m_el_size;//how many bins are we talking about
  std::vector< std::vector<double> > m_invec; // internal vector of boundaries


};
}
#endif
