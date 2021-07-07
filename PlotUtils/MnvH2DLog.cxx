#ifndef MNV_MnvH2DLOG_cxx
#define MNV_MnvH2DLOG_cxx 1

#include "MnvH2DLog.h"
#include <iostream>

using std::cout;
using std::endl;

using namespace PlotUtils;



//! A helper function which set variables for the deep copy and assignment
void MnvH2DLog::DeepCopy( const MnvH2DLog& h ){
  
  //! Set bin norm width
  //  fNormBinWidthX = h.GetNormBinWidthX();
  //  fNormBinWidthY = h.GetNormBinWidthY();
  
  //! Copy the vert and lat error bands
  std::vector<std::string> vertNames = h.GetVertErrorBandNames();
  for( std::vector<std::string>::iterator name = vertNames.begin(); name != vertNames.end(); ++name )
    fVertErrorBandMap[*name] = new MnvVertErrorBand2D( *h.GetVertErrorBand(*name) );
  
  std::vector<std::string> latNames = h.GetLatErrorBandNames();
  for( std::vector<std::string>::iterator name = latNames.begin(); name != latNames.end(); ++name )
    fLatErrorBandMap[*name] = new MnvLatErrorBand2D( *h.GetLatErrorBand(*name) );
}




void MnvH2DLog::logifyOne(TH2* h, const double lambda)// Logify one histogram in-place
{
  // Loop over bins in the universe histogram
  for(int ibin=0; ibin<h->GetNbinsX()+2; ibin++){
    for(int jbin = 0; jbin < h->GetNbinsY()+2; jbin++){
      double binContent=h->GetBinContent(ibin,jbin);
      
      h->SetBinContent(ibin, jbin, -100.);
      // h->SetBinError(ibin,jbin,0.);
      if(ibin == 0 || ibin == GetNbinsX()+1 || jbin == 0 || jbin ==GetNbinsY()+1) continue;
      if (binContent!=0){
        h->SetBinError(ibin, jbin,h->GetBinError(ibin,jbin)*BoxCoxDeriv(binContent,lambda));
        h->SetBinContent(ibin, jbin, BoxCoxVal(binContent,lambda));
      }
    }
  }
  //cout << "Logify";
  //h->Print("ALL");
}

void MnvH2DLog::logify(const double lambda)
{
  cout << " logify code " << endl;
  //  if (imlogified) return;
  //  imlogified=1;
  // Logify the CV histogram
  // save it first
  TH2D* central = (TH2D*) this->Clone();
  
  logifyOne((TH2D*)this, lambda);
  //this->Print("ALL");
  // Logify all the universes
  
  // Loop over error bands...
  std::vector<std::string> vertErrNames=this->GetVertErrorBandNames();
  for(unsigned int iband=0; iband<vertErrNames.size(); ++iband){
    
    // Loop over universe histograms in the band...
    MnvVertErrorBand2D* band=this->GetVertErrorBand(vertErrNames[iband]);
    // do the band itself - this seems to be a thing for single universes
    logifyOne((TH2*) band, lambda);
    std::vector<TH2D*> univs=band->GetHists();
    for(unsigned int iuniv=0; iuniv<univs.size(); ++iuniv){
      TH2* huniv=univs[iuniv];
      logifyOne(huniv, lambda);
      //huniv->Print();
    }
  }
  // Loop over error bands...
  std::vector<std::string> latErrNames=this->GetLatErrorBandNames();
  for(unsigned int iband=0; iband<latErrNames.size(); ++iband){
    
    // Loop over universe histograms in the band...
    MnvLatErrorBand2D* band=this->GetLatErrorBand(latErrNames[iband]);
    logifyOne((TH2*) band, lambda);
    
    std::vector<TH2D*> univs=band->GetHists();
    for(unsigned int iuniv=0; iuniv<univs.size(); ++iuniv){
      TH2* huniv=univs[iuniv];
      logifyOne(huniv, lambda);
      //huniv->Print();
    }
  }

  cout << "Found it " << HasErrorMatrix("unfoldingCov") << endl;

  // loop over error matrices ...
  std::vector<std::string> errorMatrixNames=this->GetSysErrorMatricesNames();
  for( std::vector<std::string>::const_iterator i = errorMatrixNames.begin(); i != errorMatrixNames.end(); ++i ){
    std::string fname = *i;
    cout << fname << endl;
    if ( HasErrorMatrix( *i ) ){
      std::cout << fSysErrorMatrix.find(fname)->second->GetNrows() << std::endl;
      TMatrixD covmx = *(fSysErrorMatrix.find(fname)->second);
      int nx = covmx.GetNrows();
      for (int ix = 0; ix < nx; ix++){
        double vali = central->GetBinContent(ix);
        for (int jx = 0; jx < nx; jx++){
          double valj = central->GetBinContent(jx);
    
          covmx[ix][jx] *= BoxCoxDeriv(vali,lambda)*BoxCoxDeriv(valj,lambda);
	}
      }
      FillSysErrorMatrix(fname,covmx);
    }
    
  }
}

// may need to implement scaling and operations later...
#endif
