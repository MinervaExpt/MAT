#ifndef MNV_HYPERDIMLINEARIZER_cxx
#define MNV_HYPERDIMLINEARIZER_cxx 1
#include "HyperDimLinearizer.h"
#include <string>

namespace PlotUtils
{

  //!Get a set of bins
  HyperDimLinearizer::HyperDimLinearizer(std::vector<std::vector<double> > input, int type)
  {
    int n_dim = input.size();
    m_invec= input;
    m_analysis_type = type;
  
    std::cout << "Contructing class with " << n_dim << " dimensions of type "<< type << std::endl;
    for(unsigned int i=0;i<input.size();i++){
      std::cout << "Bin number " << i << "\t" << input[i].size()+1 << std::endl;
      m_el_size.push_back(input[i].size()+1);//number of bins = vector size-1+2 assuming under/overflow.  
    }
    
  }
  
  //!Giving x,y,z... get the x,y bin//This is templated
  std::pair<int,int> HyperDimLinearizer::GetBin(std::vector<double> values)
  {
    //can I generalize this? YUP
    int y_bin=0;
    int global_x=0;
    std::vector<int> sizes;
    if(m_analysis_type==0){
      for(unsigned int i=0;i<values.size();i++){
        int scale=1;
        int tmp_bin = Get1DBin(values[i],i);
        if(i!=1){
  	for(unsigned int j=0;j<values.size();j++){
  	  if(i==j) break;
  	  if(j!=1) scale *= m_el_size[j];
  	}
  	global_x+=tmp_bin*scale;
        }
        else y_bin = tmp_bin;
      }
    }
    else if(m_analysis_type==1){
      for(unsigned int i=0;i<values.size();i++){
        int tmp_bin = Get1DBin(values[i],i);
        int scale=1;
        for(unsigned int j=0;j<values.size();j++){
  	if(i==j) break;
  	scale *= m_el_size[j];
        }
        global_x+=tmp_bin*scale;
      }
    }
    std::pair<int,int> mypair = std::make_pair(global_x,y_bin);
    return mypair;
    
  }
  
  int HyperDimLinearizer::Get1DBin(double value, int el)
  {
    int b = 0;
    for(unsigned int i=0;i<m_invec[el].size();i++){//loop over bin boundaries
      if(value<m_invec[el][i]){
        break;
      }
      b+=1;//didn't find the bin, add 1. Underflow is 0 and overflow is size()+1
    }
    return b;
  }
  
  //!Giving the x bin get the x,z... bin (or xyz in the schema of a fully linearlized model.
  std::vector<int> HyperDimLinearizer::GetValues(int x)
  {
    //given global x return vector of x,y,z... bin coordinates.
    std::vector<int> ValueCoordinates;
    std::vector<int> ValueReverse;
    if(m_analysis_type==0){
      for(unsigned int i=0;i<m_invec.size();i++){//loop over coordinates
        int scale = 1;
        if(i==1){//skip y
  	ValueCoordinates.push_back(0);
  	continue;
        }
        for(unsigned int j=0;j<m_invec.size();j++){
  	if(j==i) break;//don't scale by bigger super cells
  	if(j==1)continue;//skip y
  	scale*=m_el_size[j];
        }
        ValueCoordinates.push_back((x/scale)%m_el_size[i]);
      }
    }
    else if(m_analysis_type==1){
      for(unsigned int i=0;i<m_invec.size();i++){//loop over coordinates
        int scale = 1;
        for(unsigned int j=0;j<m_invec.size();j++){
  	if(j==i) break;//don't scale by bigger super cells
  	scale*=m_el_size[j];
        }
        ValueCoordinates.push_back((x/scale)%m_el_size[i]);
      }
    }
    return ValueCoordinates;
  }
  
  
  std::vector<TH2D*> HyperDimLinearizer::Get2DHistos(PlotUtils::MnvH2D *result, bool IncludeSys = false){
    //  std::cout <<"Entering Get2DHistos"  << std::endl;
    std::vector<TH2D*> expanded_results;
    TH2D mybigmap;
    if(!IncludeSys) mybigmap = result->GetCVHistoWithStatError();
    else mybigmap = result->GetCVHistoWithError();
    if(m_analysis_type==0){
      //    std::cout << "Starting up get 2D histos with analysis type 0" << std::endl;
      //projected N dims (less Y) come in chunks of X bins (including under/over)
      int num_chunks = 1;
      const int num_x_bins = m_invec[0].size()-1;
      const int num_y_bins = m_invec[1].size()-1;
      //    std::cout << "Master Plot has x " << num_x_bins << "\ty\t" << num_y_bins << std::endl;
      for(unsigned int i=2;i<m_invec.size();i++)num_chunks*=(m_invec[i].size()+1);
      //    std::cout << "I have number of chunks = " << num_chunks << std::endl;
      for(int i=0;i<num_chunks;i++){
        TH2D *tmp_bin = new TH2D(Form("Chunk_%d",i),Form("Chunk_%d",i),num_x_bins,&m_invec[0][0],num_y_bins,&m_invec[1][0]);
        int offset_x = (num_x_bins+2)*i+1;//need to know low bin for chunk
        for(int j=0;j<num_x_bins+2;j++){
  	for(int k=0;k<num_y_bins+2;k++){
  	  double tmpval = mybigmap.GetBinContent(j+offset_x,k);
  	  double tmperr = mybigmap.GetBinError(j+offset_x,k);
  	  tmp_bin->SetBinContent(j,k,tmpval);
  	  tmp_bin->SetBinError(j,k,tmperr);
  	}//end loop over subset y
        }//end loop over subset x
        expanded_results.push_back((TH2D*)tmp_bin->Clone(Form("Clone_%d",i)));//woohoo got a 2D result for one of the chunks!
      }//end loop over chunks
    }
    return expanded_results;
  
  }
  
  
  std::vector<PlotUtils::MnvH2D*> HyperDimLinearizer::Get2DMnvHistos(PlotUtils::MnvH2D *result, bool IncludeSys = false){
    std::cout << "Entering Get2DMnvHistos"  << std::endl;
    std::vector<PlotUtils::MnvH2D*> expanded_results;
    std::vector<TH2D*> CV_vals = Get2DHistos(result,false);//get CV
    std::cout << "I have " << CV_vals.size() << " CV histograms" << std::endl;
    for(uint i=0;i<CV_vals.size();i++) expanded_results.push_back(new PlotUtils::MnvH2D(*CV_vals[i]));
  
    std::vector<std::string> vertnames = result->GetVertErrorBandNames();
    std::vector<std::string> latnames = result->GetLatErrorBandNames();
    
    //Do vert first
    for(uint i=0;i<vertnames.size();i++){
      std::cout << "Working on " << vertnames[i] << std::endl;
      std::vector<std::vector<TH2D*> > unihists;
      PlotUtils::MnvVertErrorBand2D *band = result->GetVertErrorBand(vertnames[i]);
      int bandsize = band->GetNHists();
      for(int uni=0;uni<bandsize;uni++){
        std::vector<TH2D*> tmpbandset = Get2DHistos(new PlotUtils::MnvH2D(*band->GetHist(uni)));//Get the universe hist and spit out the N 2D results.
        unihists.push_back(tmpbandset);
      }
      //Have an N by uni matrix of TH2D. Now time to push back into the primary
      std::cout << "I have created a set of flux hists. This is size of the vector " << unihists.size() << "\t" << unihists[0].size()  << std::endl;
  
      for(int j=0;j<unihists[0].size();j++){//unihists[0].size() is the number of projections needed
        std::vector<TH2D*> tmpband;
        for(int uni=0;uni<bandsize;uni++){
  	tmpband.push_back(unihists[uni][j]);
        }
        expanded_results[j]->AddVertErrorBand(vertnames[i],tmpband);
      }
      
    }
  
    //now lat
     for(uint i=0;i<latnames.size();i++){
      std::cout << "Working on " << latnames[i] << std::endl;
      std::vector<std::vector<TH2D*> > unihists;
      PlotUtils::MnvLatErrorBand2D *band = result->GetLatErrorBand(latnames[i]);
      int bandsize = band->GetNHists();
      for(int uni=0;uni<bandsize;uni++){
        std::vector<TH2D*> tmpbandset = Get2DHistos(new PlotUtils::MnvH2D(*band->GetHist(uni)));//Get the universe hist and spit out the N 2D results.
        unihists.push_back(tmpbandset);
      }
      //Have an N by uni matrix of TH2D. Now time to push back into the primary
  
      std::cout << "I have created a set of flux hists. This is size of the vector " << unihists.size() << "\t" << unihists[0].size()  << std::endl;
  
      for(uint j=0;j<unihists[0].size();j++){//unihists[0].size() is the number of projections needed
        std::vector<TH2D*> tmpband;
        for(int uni=0;uni<bandsize;uni++){
  	tmpband.push_back(unihists[uni][j]);
        }
        expanded_results[j]->AddLatErrorBand(latnames[i],tmpband);
      }
    }
  
    return expanded_results;
  }
  
  
  TH2D* HyperDimLinearizer::Get2DHisto(PlotUtils::MnvH1D *result, bool IncludeSys = false){
    if(m_invec.size()!=2) std::cout << "THIS ONLY WORKS FOR 2D RESULTS.\nIf you are a mapped 3D or more you need to use something different which might not exist." << std::endl;
    std::string myname = Form("Unmapped_%s",result->GetTitle());
    TH1D* mybigmap = NULL;
    if(!IncludeSys) mybigmap = new TH1D(result->GetCVHistoWithStatError());
    else mybigmap = new TH1D(result->GetCVHistoWithError());
    const int num_x_bins = m_invec[0].size()-1;
    const int num_y_bins = m_invec[1].size()-1;
    TH2D* my2D = new TH2D(myname.c_str(),myname.c_str(),num_x_bins,&m_invec[0][0],num_y_bins,&m_invec[1][0]);
    if(m_analysis_type==1){
      std::cout << "Starting up get 2D histo with analysis type 1" << std::endl;
      for(int i=0;i<num_x_bins+2;i++){//includes under/over
        for(int j=0;j<num_y_bins+2;j++){//includes under/over
  	double tmpval = mybigmap->GetBinContent(j*(num_x_bins+2)+(i+1));
  	double tmperr = mybigmap->GetBinError(j*(num_x_bins+2)+(i+1));
  	my2D->SetBinContent(i,j,tmpval);
  	my2D->SetBinError(i,j,tmperr);
        }
      }
    }
  
    return my2D;
  
  }
  
  PlotUtils::MnvH2D* HyperDimLinearizer::Get2DMnvHisto(PlotUtils::MnvH1D *result, bool IncludeSys = false){
    std::cout << "Entering Get2DMnvHisto"  << std::endl;
    if(m_invec.size()!=2) std::cout << "THIS ONLY WORKS FOR 2D RESULTS.\nIf you are a mapped 3D or more you need to use something different which might not exist." << std::endl;
  
    TH2D* CV_vals = Get2DHisto(result,false);
    PlotUtils::MnvH2D *expanded_result = new PlotUtils::MnvH2D(*CV_vals);
    
    std::vector<std::string> vertnames = result->GetVertErrorBandNames();
    std::vector<std::string> latnames = result->GetLatErrorBandNames();
  
      //Do vert first
    for(uint i=0;i<vertnames.size();i++){
      std::cout << "Working on " << vertnames[i] << std::endl;
      std::vector<TH2D*> unihists;
      PlotUtils::MnvVertErrorBand *band = result->GetVertErrorBand(vertnames[i]);
      int bandsize = band->GetNHists();
      for(int uni=0;uni<bandsize;uni++){
        TH2D* tmpbandset = Get2DHisto(new PlotUtils::MnvH1D(*band->GetHist(uni)));//Get the universe hist and spit out the N 2D results.
        unihists.push_back(tmpbandset);
      }
      expanded_result->AddVertErrorBand(vertnames[i],unihists);
    }
  
    for(uint i=0;i<latnames.size();i++){
      std::cout << "Working on " << latnames[i] << std::endl;
      std::vector<TH2D*> unihists;
      PlotUtils::MnvLatErrorBand *band = result->GetLatErrorBand(latnames[i]);
      int bandsize = band->GetNHists();
      for(int uni=0;uni<bandsize;uni++){
        TH2D* tmpbandset = Get2DHisto(new PlotUtils::MnvH1D(*band->GetHist(uni)));//Get the universe hist and spit out the N 2D results.
        unihists.push_back(tmpbandset);
      }
      expanded_result->AddLatErrorBand(latnames[i],unihists);
    }
  
  
    return expanded_result;
  
  }
  
  
  
  void HyperDimLinearizer::TestFunctionality()
  {
    std::cout << "Initializing funcationality test"<<std::endl;
  
    //how many bins
    int n_bins = 1;
    std::cout << "Running with n-dimensions = " << m_invec.size() << std::endl;
    for(unsigned int i=0;i<m_invec.size();i++){
      if(m_analysis_type==0 && i==1) continue;
      std::cout << "Coord " << i << " has " << m_el_size[i] << " bins " << std::endl;
      n_bins*=m_el_size[i];
    }
    std::cout << "This gives us a total of " << n_bins << " bins" << std::endl;
  
  
    for(int i=0;i<n_bins;i++){
      std::vector<int> coordinates = GetValues(i);
      std::cout << i << "\t";
      for(unsigned int j=0;j<coordinates.size();j++){
        std::cout << coordinates[j] << "\t";
      }
      std::cout << std::endl;
    }
  }    

} //namespace PlotUtils
#endif
