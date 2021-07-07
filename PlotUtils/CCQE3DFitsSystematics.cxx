#ifndef CCQE3DFITS_SYSTEMATICS_CXX
#define CCQE3DFITS_SYSTEMATICS_CXX

#include "CCQE3DFitsSystematics.h"
#include "HyperDimLinearizer.h"
#include "TFile.h"
#include <iostream>

// Helper Functions
namespace PlotUtils{

  //=================================================================================
  // QELike Signal Defintion 
  //=================================================================================
  template <typename T>
  bool IsQELike(const T& universe) {

    int genie_n_muons         = 0;
    int genie_n_mesons        = 0;
    int genie_n_heavy_baryons_plus_pi0s = 0;
    int genie_n_photons       = 0;
  
    int nparticles = universe.GetInt("mc_nFSPart");
    for(int i = 0; i < nparticles; ++i) {
      int pdg = universe.GetVecElemInt("mc_FSPartPDG",i);
      double energy = universe.GetVecElem("mc_FSPartE",i);
  
      if( abs(pdg) == 13 ) genie_n_muons++;
      else if( pdg == 22 && energy > 10 ) genie_n_photons++;
      else if( abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 323 || pdg == 111 || pdg == 130 || pdg == 310 || pdg == 311 || pdg == 313 ) genie_n_mesons++;
      else if( pdg == 3112 || pdg == 3122 || pdg == 3212 || pdg == 3222 || pdg == 4112 || pdg == 4122 || pdg == 4212 || pdg == 4222 || pdg == 411 || pdg == 421 || pdg == 111 ) genie_n_heavy_baryons_plus_pi0s++;
    }


  if(genie_n_muons         == 1 &&
     genie_n_mesons        == 0 &&
     genie_n_heavy_baryons_plus_pi0s == 0 &&
     genie_n_photons       == 0 ) return true;
      
  return false;

  }

}

// Class Definitions
// Constructor
template<typename T>
PlotUtils::CCQE3DFitsUniverse<T>::CCQE3DFitsUniverse( typename T::config_t chw, double nsigma, int variation ) 
  : T(chw, nsigma),m_variation(variation) 
{
  // Set up instance of HyperDimLinearizer
  std::vector<double> pz3Dbins;
  pz3Dbins.push_back(1.5);
  pz3Dbins.push_back(3.5);//added ME
  pz3Dbins.push_back(8.0);
  pz3Dbins.push_back(20.0);

  std::vector<double> pt3Dbins;
  pt3Dbins.push_back(0.0);
  pt3Dbins.push_back(0.075);//added ME
  pt3Dbins.push_back(0.15);
  pt3Dbins.push_back(0.25);//added ME
  pt3Dbins.push_back(0.4);
  pt3Dbins.push_back(0.7);//added ME
  pt3Dbins.push_back(1.0);
  pt3Dbins.push_back(1.5);//added ME
  pt3Dbins.push_back(2.5);

  std::vector<double> recoil3Dbins;
  for(int i=0;i<10;i++)recoil3Dbins.push_back(i*40);
  for(int i=0;i<4;i++)recoil3Dbins.push_back(i*200+400);
  recoil3Dbins.push_back(40000);

  std::vector<std::vector<double> > full3D;
  full3D.push_back(recoil3Dbins);
  full3D.push_back(pt3Dbins);
  full3D.push_back(pz3Dbins);

  PlotUtils::HyperDimLinearizer hyperDimLinearizer = PlotUtils::HyperDimLinearizer(full3D,0);
  
  // Pull hists out of Dan's file and unpack using HyperDimLinearizer 
  TString filename = "root://fndca1.fnal.gov//pnfs/fnal.gov/usr/minerva/persistent/users/finer/highNu/histProcessing/auxiliaryROOTFiles/CrossSection_per_nucleon_3D_pzptreco_iterations_10_CombinedPlaylists.root_big3d.root";
  TFile* CCQE3DFitsWeightFile=TFile::Open(filename,"READ");
  PlotUtils::MnvH2D* dataHist = (PlotUtils::MnvH2D*)CCQE3DFitsWeightFile->Get("h_pzptrec_data_nobck_unfold_effcor_cross_section");
  PlotUtils::MnvH2D* mcHist = (PlotUtils::MnvH2D*)CCQE3DFitsWeightFile->Get("h_pzptrec_mc_nobck_unfold_effcor_cross_section");

  PlotUtils::MnvH2D* ratioHist = dataHist->Clone("3D_Fits_Ratio_hist");
  ratioHist->Divide(dataHist,mcHist);

  m_weightHists = hyperDimLinearizer.Get2DHistos(ratioHist,true);

  CCQE3DFitsWeightFile->Close();
}


template<typename T>
double PlotUtils::CCQE3DFitsUniverse<T>::GetCCQE3DFitsWeight() const {

  if(!PlotUtils::IsQELike(*this)) return T::GetCCQE3DFitsWeight();  // Only reweight QELike events

  // Fetch kinematic variables for this event and map to bin numbers of the unpacked hists
  double pz = T::GetPmuLongitudinalTrue();
  double pt = T::GetPmuTransverseTrue();
  double availableE = T::GetEAvailableTrue();

  int nPzBin;
  if(pz<1.5){nPzBin=1;}
  else if(pz<3.5){nPzBin=2;}
  else if(pz<8.0){nPzBin=3;}
  else{nPzBin=3;}
  int nPtBin = m_weightHists[1]->GetYaxis()->FindBin(pt);
  int nAvailableEBin = m_weightHists[1]->GetXaxis()->FindBin(availableE);

  // Variation type 1 is for the MC excess at High-Pt, Low availableE
  if(m_variation==1&&!((nPtBin==5&&nAvailableEBin==1)||(nPtBin==6&&(nAvailableEBin==1||nAvailableEBin==2)))){
    return 1.0;
  }
  // Variation type 2 is for the MC deficit at low-moderate Pt, Low availableE
  if(m_variation==2&&(nPtBin>4||(nPtBin==1&&nAvailableEBin>5)||(nPtBin==2&&nAvailableEBin>6)||(nPtBin==3&&nAvailableEBin>8)||(nPtBin==4&&nAvailableEBin>10))){
    return 1.0;
  }
  // Variation type 3 is for the MC excess at low-moderate Pt, moderate availableE
  if(m_variation==3&&(nPtBin>4||(nPtBin==1&&nAvailableEBin<=5)||(nPtBin==2&&nAvailableEBin<=6)||(nPtBin==3&&nAvailableEBin<=8)||(nPtBin==4&&nAvailableEBin<=10))){
    return 1.0;
  }

  // Pull weight from unpacked hists 
  double wgt = m_weightHists[nPzBin]->GetBinContent(nAvailableEBin,nPtBin);
  return  1.0 + T::m_nsigma*(wgt-1.);

}


template<typename T>
std::string PlotUtils::CCQE3DFitsUniverse<T>::ShortName() const { return "3D_CCQE_Fits"; }


template<typename T>
std::string PlotUtils::CCQE3DFitsUniverse<T>::LatexName() const { return "3D CCQE Fits"; }


#endif // CCQE3DFITS_SYSTEMATICS_CXX
