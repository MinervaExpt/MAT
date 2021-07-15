#include <fstream>
#include "MnvH2DToCSV.h"
#include <string>
#include <iostream>

// function to dump histograms to CSV files



namespace MAT{
 

  
  void MnvH2DToCSV(MAT::MnvH2D *hist, std::string name, std::string directory = "./", double scale = 1.0, bool fullprecision, bool syserrors){
    std::cout << "entering H2DToCsV" << name << std::endl;
    hist->Print();
    std::ofstream *f_values =new std::ofstream();
    std::ofstream *f_err =new std::ofstream();
    std::ofstream *f_staterr =new std::ofstream();
    std::ofstream *f_syserr =new std::ofstream();
    std::ofstream *f_bins =new std::ofstream();
    std::ofstream *f_corr =new std::ofstream();
    
    
    
    f_values->open((directory+name+".csv").c_str());
    f_err->open((directory+name+"_errors.csv").c_str());
    f_staterr->open((directory+name+"_staterrors.csv").c_str());
    f_syserr->open((directory+name+"_syserrors.csv").c_str());
    f_bins->open((directory+name+"_bins.csv").c_str());
    f_corr->open((directory+name+"_covariance.csv").c_str());
    
    
    TH2D stat=hist->GetStatError(); //stat error
    TH2D total=hist->GetCVHistoWithError(); // CV with total error
    TH2D sys=hist->GetTotalError(false); //sys error only
    sys.Print("ALL");
    *f_bins<<"Bins"<<std::endl;
    
    *f_values << "Values\t";
    *f_err << "err\t";
    *f_staterr << "staterr\t";
    *f_syserr << "syserr\t";
    
    *f_values<<std::endl;
    *f_err<<std::endl;
    *f_staterr<<std::endl;
    *f_syserr<<std::endl;
    
    *f_bins<<hist->GetXaxis()->GetBinLowEdge(1)<< "\t"; //<<std::endl;
    
    for (int x=1;x<=hist->GetXaxis()->GetNbins();x++){
      *f_bins<<hist->GetXaxis()->GetBinUpEdge(x)<< "\t";//<<std::endl;
      for (int y=1;y<=hist->GetYaxis()->GetNbins();y++){
        if (y> 1) {
          *f_values<<",\t";
          *f_err<<",\t";
          *f_staterr<<",\t";
          *f_syserr<<",\t";
        }
        if (!fullprecision){
          
          
          *f_values<<Form("%.2f ",total.GetBinContent(x,y)/hist->GetXaxis()->GetBinWidth(x)/hist->GetYaxis()->GetBinWidth(y)*scale);
          *f_err<<Form("%.2f ",total.GetBinError(x,y)/hist->GetXaxis()->GetBinWidth(x)/hist->GetYaxis()->GetBinWidth(y)*scale);
          *f_staterr<<Form("%.2f ",stat.GetBinContent(x,y)/hist->GetXaxis()->GetBinWidth(x)/hist->GetYaxis()->GetBinWidth(y)*scale);
          *f_syserr<<Form("%.2f ",sys.GetBinContent(x,y)/hist->GetXaxis()->GetBinWidth(x)/hist->GetYaxis()->GetBinWidth(y)*scale);
          //std::cout << "syscheck" <<  sys.GetBinContent(x,y) << std::endl;
        }
        else{
          
          *f_values<<Form("%.17e ",total.GetBinContent(x,y)/hist->GetXaxis()->GetBinWidth(x)/hist->GetYaxis()->GetBinWidth(y)*scale);
          *f_err<<Form("%.17e ",total.GetBinError(x,y)/hist->GetXaxis()->GetBinWidth(x)/hist->GetYaxis()->GetBinWidth(y)*scale);
          *f_staterr<<Form("%.17e ",stat.GetBinContent(x,y)/hist->GetXaxis()->GetBinWidth(x)/hist->GetYaxis()->GetBinWidth(y)*scale);
          *f_syserr<<Form("%.17e ",sys.GetBinContent(x,y)/hist->GetXaxis()->GetBinWidth(x)/hist->GetYaxis()->GetBinWidth(y)*scale);
          //std::cout << "syscheck" <<  sys.GetBinContent(x,y) << std::endl;
        }
      }
      *f_values<<std::endl;
      *f_err<<std::endl;
      *f_staterr<<std::endl;
      *f_syserr<<std::endl;
    }
    
    *f_bins<< std::endl;
    for (int y=0;y<=hist->GetYaxis()->GetNbins();y++){
      
      *f_bins<<hist->GetYaxis()->GetBinUpEdge(y) << "\t";//<<std::endl;
    }
    *f_bins << std::endl;
    
    f_values->close();
    f_err->close();
    f_staterr->close();
    f_syserr->close();
    f_bins->close();
    
    //    TMatrixD correlation_matrix= hist->GetTotalCorrelationMatrix();
    TMatrixD correlation_matrix= hist->GetTotalErrorMatrix();
    correlation_matrix *= (scale*scale); // scale by factor of 10^41
    *f_corr << "Covariance" << std::endl;
    
    int nbins_x=hist->GetNbinsX();
    int nbins_y=hist->GetNbinsY();
    int totalbins=(nbins_x+2)*(nbins_y+2);
    
    for (int i=0;i<totalbins;i++)
      {
      int x=i%(nbins_x+2);
      int y=i/(nbins_x+2);
      double binwidcorri;
      
      binwidcorri = hist->GetXaxis()->GetBinWidth(x)*hist->GetYaxis()->GetBinWidth(y);
      
      
      
      if (x==0 || y==0 || x==nbins_x+1 || y== nbins_y+1) continue; // Do not print overflow and underflow
      
      //*f_corr<< "bin_"<<x<<"_"<<y;
      int first = 0;
      for (int j=0;j<totalbins;j++)
        {
        int this_x=j%(nbins_x+2);
        int this_y=j/(nbins_x+2);
        if  (this_x==0 || this_y==0 || this_x==nbins_x+1 || this_y== nbins_y+1) continue; // Do not print overflow and underflow
        double binwidcorrj;
        if ( first != 0) {
          *f_corr<<",\t";
          
        }
        first ++;
       
          binwidcorrj = hist->GetXaxis()->GetBinWidth(this_x)*hist->GetYaxis()->GetBinWidth(this_y);
        
        *f_corr<<Form("%.17e",correlation_matrix[i][j]/binwidcorri/binwidcorrj);   // need to include bin widths
        }
      *f_corr<<std::endl;
      }
    f_corr->close();
    std::ofstream * f_errors = new std::ofstream();
    std::cout << " Lat and Vert bands not implemented yet" << std::endl;
    f_errors->open((directory+name+"_sysdump.csv").c_str());
    std::vector<std::string> cov_errNames = hist->GetCovMatricesNames();
    *f_errors << " covariance " << cov_errNames.size() <<  std::endl;
    for( std::vector<std::string>::iterator name=cov_errNames.begin(); name!=cov_errNames.end(); ++name ){
      TMatrixD  v = hist->GetSysErrorMatrix(*name);
      *f_errors << Form("%s",name->c_str()) << std::endl;
      int nbins_x=hist->GetNbinsX();
      int nbins_y=hist->GetNbinsY();
      int totalbins=(nbins_x+2)*(nbins_y+2);
      
      for (int i=0;i<totalbins;i++)
        {
        int x=i%(nbins_x+2);
        int y=i/(nbins_x+2);
        double binwidcorri;
        
        binwidcorri = hist->GetXaxis()->GetBinWidth(x)*hist->GetYaxis()->GetBinWidth(y);
        
        
        
        if (x==0 || y==0 || x==nbins_x+1 || y== nbins_y+1) continue; // Do not print overflow and underflow
        
        //*f_corr<< "bin_"<<x<<"_"<<y;
        int first = 0;
        for (int j=0;j<totalbins;j++)
          {
          int this_x=j%(nbins_x+2);
          int this_y=j/(nbins_x+2);
          if  (this_x==0 || this_y==0 || this_x==nbins_x+1 || this_y== nbins_y+1) continue; // Do not print overflow and underflow
          double binwidcorrj;
          if ( first != 0) {
            *f_errors<<",\t";
            
          }
          first ++;
          
          binwidcorrj = hist->GetXaxis()->GetBinWidth(this_x)*hist->GetYaxis()->GetBinWidth(this_y);
          
          *f_errors<<Form("%.17e",v[i][j]/binwidcorri/binwidcorrj*scale*scale);   // need to include bin widths
          }
        *f_errors<<std::endl;
        }
    }
    f_errors->close();
  }
}

