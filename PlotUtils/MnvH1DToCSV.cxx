#include <fstream>
#include "PlotUtils/MnvH1DToCSV.h"
#include <string>
#include <iostream>

// function to dump histograms to CSV files



namespace PlotUtils{
  
  
  
  void MnvH1DToCSV(PlotUtils::MnvH1D *hist, std::string name, std::string directory = "./", double scale=1.0, bool fullprecision, bool syserrors){
    
    std::cout << "entering 1DToCSV " << name << std::endl;
    std::ofstream *f_values =new std::ofstream();
    std::ofstream *f_err =new std::ofstream();
    std::ofstream *f_staterr =new std::ofstream();
    std::ofstream *f_syserr =new std::ofstream();
    std::ofstream *f_bins =new std::ofstream();
    std::ofstream *f_corr =new std::ofstream();
    
    
    
    f_values->open((directory+name+"_1d.csv").c_str());
    f_err->open((directory+name+"_errors_1d.csv").c_str());
    f_staterr->open((directory+name+"_staterrors_1d.csv").c_str());
    f_syserr->open((directory+name+"_syserrors_1d.csv").c_str());
    f_bins->open((directory+name+"_bins_1d.csv").c_str());
    f_corr->open((directory+name+"_covariance.csv").c_str());
    
    TH1D stat=hist->GetStatError(); //stat error
    TH1D total=hist->GetCVHistoWithError(); // CV with total error
    TH1D sys=hist->GetTotalError(false); //sys error only
    
    //    *f_bins<<hist->GetXaxis()->GetBinLowEdge(1); //<<std::endl;
    *f_bins<<"Bins";
    
    *f_values << "Values\t";
    *f_err << "err\t";
    *f_staterr << "staterr\t";
    *f_syserr << "syserr\t";
    
    *f_bins << std::endl;
    *f_values << std::endl;
    *f_err << std::endl;
    *f_staterr << std::endl;
    *f_syserr << std::endl;
    *f_bins<<hist->GetXaxis()->GetBinLowEdge(1)<< "\t"; //<<std::endl;
    if(fullprecision){
      
      
      for (int i=1;i<=hist->GetXaxis()->GetNbins();i++)
        {
        
        if (i>1) {
          *f_values << ",\t";
          *f_bins << ",\t";
          *f_err << ",\t";
          *f_staterr << ",\t";
          *f_syserr << ",\t";
        }
        // Bin width normalize if not enu when we want a total x sec
        *f_bins<<hist->GetXaxis()->GetBinUpEdge(i)<< "\t";//<<std::endl;
        *f_values<<Form("%.17e ",total.GetBinContent(i)/(hist->GetXaxis()->GetBinWidth(i))*scale);
        *f_err<<Form("%.17e ",total.GetBinError(i)/(hist->GetXaxis()->GetBinWidth(i))*scale);
        *f_staterr<<Form("%.17e ",stat.GetBinContent(i)/(hist->GetXaxis()->GetBinWidth(i))*scale);
        *f_syserr<<Form("%.17e ",sys.GetBinContent(i)/(hist->GetXaxis()->GetBinWidth(i))*scale);
        
        }
    }
    else{
      //    *f_bins<<hist->GetXaxis()->GetBinLowEdge(1); //<<std::endl;
      /**f_bins<<"Bins\t";
       *f_bins<<hist->GetXaxis()->GetBinLowEdge(1)<< "\t"; //<<std::endl;
       *f_values << "Values\t";
       *f_err << "err\t";
       *f_staterr << "staterr\t";
       *f_syserr << "syserr\t";
       */
      
      for (int i=1;i<=hist->GetXaxis()->GetNbins();i++)
        {
        if (i>1) {
          *f_values << ",\t";
          *f_bins << ",\t";
          *f_err << ",\t";
          *f_staterr << ",\t";
          *f_syserr << ",\t";
        }
        *f_bins<<hist->GetXaxis()->GetBinUpEdge(i)<< "\t";//<<std::endl;
        // Bin width normalize if not enu when we want a total x sec
        *f_values<<Form("%.2f ",total.GetBinContent(i)/(hist->GetXaxis()->GetBinWidth(i))*scale);
        *f_err<<Form("%.2f ",total.GetBinError(i)/(hist->GetXaxis()->GetBinWidth(i))*scale);
        *f_staterr<<Form("%.2f ",stat.GetBinContent(i)/(hist->GetXaxis()->GetBinWidth(i))*scale);
        *f_syserr<<Form("%.2f ",sys.GetBinContent(i)/(hist->GetXaxis()->GetBinWidth(i))*scale);
        
        }
    }
    *f_bins << std::endl;
    *f_values << std::endl;
    *f_err << std::endl;
    *f_staterr << std::endl;
    *f_syserr << std::endl;
    f_values->close();
    f_err->close();
    f_staterr->close();
    f_syserr->close();
    f_bins->close();
    
    //    TMatrixD correlation_matrix= hist->GetTotalCorrelationMatrix();
    TMatrixD correlation_matrix= hist->GetTotalErrorMatrix();
    correlation_matrix *= (scale*scale); // scale by factor of 10^41
    
    int nbins_x=hist->GetNbinsX();
    
    int totalbins=(nbins_x+2);
    
    
    *f_corr<<std::endl;
    for (int x=0;x<totalbins;x++)
      {
      
      double binwidcorri;
      
      binwidcorri = hist->GetXaxis()->GetBinWidth(x);
      
      
      
      if (x==0 ||  x==nbins_x+1 ) continue; // Do not print overflow and underflow
      
      for (int this_x=0;this_x<totalbins;this_x++)
        {
        
        
        if  (this_x==0 || this_x==nbins_x+1 ) continue; // Do not print overflow and underflow
        double binwidcorrj;
        
        binwidcorrj = hist->GetXaxis()->GetBinWidth(this_x);
        if (this_x > 1) *f_corr<< ",\t";
        *f_corr<<Form("%.17e ",correlation_matrix[x][this_x]/binwidcorri/binwidcorrj);
        // need to include bin widths
        }
      *f_corr<<std::endl;
      }
    f_corr->close();
//    if (!syserrors){
//      std::cout << " no systematic errors to consider " << std::endl;
//      return;
//    }
    std::cout << "are we doing systematics " << std::endl;
    std::ofstream * f_errors = new std::ofstream();
    f_errors->open((directory+name+"_sysdump.csv").c_str());
    std::vector<std::string> vert_errBandNames = hist->GetVertErrorBandNames();
    std::vector<std::string> lat_errBandNames  = hist->GetLatErrorBandNames();
    std::vector<std::string> uncorr_errBandNames  = hist->GetUncorrErrorNames();
    std::vector<std::string> cov_errNames = hist->GetCovMatricesNames();
    *f_errors << " vert " << vert_errBandNames.size() <<  std::endl;
    for( std::vector<std::string>::iterator name=vert_errBandNames.begin(); name!=vert_errBandNames.end(); ++name ){
      MnvVertErrorBand* v = hist->GetVertErrorBand(*name);
      unsigned int nunis = v->GetNHists();
      
      for (int i = 0; i< nunis; i++){
        *f_errors << Form("%s_%d",name->c_str(),i) << std::endl;
        TH1* h = v->GetHist(i);
        for (int j=1;j <=h->GetXaxis()->GetNbins();j++)
          {
          if (j>1) {
            *f_errors << ",\t";
          }
          *f_errors<<Form("%.17e ",h->GetBinContent(j)/(h->GetXaxis()->GetBinWidth(j))*scale);
          }
        *f_errors << std::endl;
      }
    }
    *f_errors << " Lat " << lat_errBandNames.size() <<  std::endl;
    for( std::vector<std::string>::iterator name=lat_errBandNames.begin(); name!=lat_errBandNames.end(); ++name ){
      MnvLatErrorBand* v = hist->GetLatErrorBand(*name);
      unsigned int nunis = v->GetNHists();
      
      for (int i = 0; i< nunis; i++){
        *f_errors << Form("%s_%d",name->c_str(),i) << std::endl;
        TH1* h = v->GetHist(i);
        for (int j=1;j <=h->GetXaxis()->GetNbins();j++)
          {
          if (j>1) {
            *f_errors << ",\t";
          }
          *f_errors<<Form("%.17e ",h->GetBinContent(j)/(h->GetXaxis()->GetBinWidth(j))*scale);
          }
        *f_errors << std::endl;
      }
    }
    *f_errors << " covariance " << cov_errNames.size() <<  std::endl;
    for( std::vector<std::string>::iterator name=cov_errNames.begin(); name!=cov_errNames.end(); ++name ){
      TMatrixD  v = hist->GetSysErrorMatrix(*name);
      *f_errors << Form("%s",name->c_str()) << std::endl;
      int nbins_x = hist->GetXaxis()->GetNbins();
      int totalbins=(nbins_x+2);
      
      for (int x=0;x<totalbins;x++)
        {
        double binwidcorri;
        
        binwidcorri = hist->GetXaxis()->GetBinWidth(x);
        if (x==0 ||  x==nbins_x+1 ) continue; // Do not print overflow and underflow
        
        for (int this_x=0;this_x<totalbins;this_x++)
          {
          
          if  (this_x==0 || this_x==nbins_x+1 ) continue; // Do not print overflow and underflow
          double binwidcorrj;
          
          binwidcorrj = hist->GetXaxis()->GetBinWidth(this_x);
          if (this_x > 1) *f_errors<< ",\t";
          *f_errors<<Form("%.17e ",v[x][this_x]/binwidcorri/binwidcorrj*scale*scale);
          // need to include bin widths
          }
        *f_errors<<std::endl;
        }
    }
    f_errors->close();
  }
}
