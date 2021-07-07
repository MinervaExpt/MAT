#ifndef TRANSWARPEXTRACTION_H
#define TRANSWARPEXTRACTION_H

#include "MinervaUnfold/MnvUnfold.h"
#include "RooUnfold/RooUnfold.h"

#include "PlotUtils/MnvPlotter.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvH3D.h"

#include "TProfile.h"
#include "TRandom3.h"

#include <iostream>

double square_element ( double mean, double element ) { return (element-mean)*(element-mean) ;}

namespace PlotUtils
{
  template <class MnvH> class TransWarpExtraction
  {
    public:
      TransWarpExtraction( std::string data_file,  std::string data_name,  std::string data_truth_file, std::string data_truth_name, std::string reco_file, std::string reco_name, 
                           std::string truth_file, std::string truth_name, std::string migration_file,  std::string migration_name, std::vector<int> iterations,  std::vector<int> exclude_chi2_bins,
                           int random_seed, int stat_universes = 1, bool bIterLogScale = false );
      ~TransWarpExtraction();
      void Setup2DChi2Hists( );
      void ClearInputErrorBands();
      void SetChi2MaxAndStep( double MaxChi2, double Chi2Step ){ m_ChiSquareGuess = MaxChi2; m_ChiSquareStep = Chi2Step; }; 
      void SetStatScale( double stat_scale ){ m_stat_scale = stat_scale; }; 
      void ScalePOT( double pot_scale ); 
      void ScaleDataPOT( double data_pot_scale ); 
      void ScaleMigStatUnc( double stat_scale );
      void UnfoldStatUniverses( );
      void WriteOutput( std::string output_filename, bool b_WriteChi2Verbose = false, bool b_WriteLinearHists = false ); 
      bool UnfoldData( MnvH1D* &h_data_unfolded, MnvH1D* stat_varied, MnvH1D* h_reco, MnvH1D* h_truth, MnvH2D* h_migration, int num_iter ); 
      bool UnfoldData( MnvH2D* &h_data_unfolded, MnvH2D* stat_varied, MnvH2D* h_reco, MnvH2D* h_truth, MnvH2D* h_migration, int num_iter ); 
      bool UnfoldData( MnvH3D* &h_data_unfolded, MnvH3D* stat_varied, MnvH3D* h_reco, MnvH3D* h_truth, MnvH2D* h_migration, int num_iter ); 
      TMatrixD UnfoldDummy( MnvH1D* stat_varied, MnvH1D* h_data_unfolded, MnvH1D* h_reco, MnvH1D* h_truth, MnvH2D* h_migration, int num_iter ); 
      TMatrixD UnfoldDummy( MnvH2D* stat_varied, MnvH2D* h_data_unfolded, MnvH2D* h_reco, MnvH2D* h_truth, MnvH2D* h_migration, int num_iter ); 
      TMatrixD UnfoldDummy( MnvH3D* stat_varied, MnvH3D* h_data_unfolded, MnvH3D* h_reco, MnvH3D* h_truth, MnvH2D* h_migration, int num_iter ); 
      void CalcChi2();
      void MakeBinChi2Dists();
      void MakeMedianChi2Dists( double percentile = 0.5 );
      void MakeTruncatedChi2Dists( double nStdDev = 2.0 );
      void MakeErrorHists( );
      void MakeRatioHists( );
  
    protected:
      MnvH* throwStat(MnvH* hist, TRandom3 *gen, int myuni, std::string prefix);
      MnvH1D* LinearizeHist( MnvH* hist );
      void bookHistos( MnvH1D*& ret_hist, std::string filename, std::string hist_name );
      void bookHistos( MnvH2D*& ret_hist, std::string filename, std::string hist_name );
      void bookHistos( MnvH3D*& ret_hist, std::string filename, std::string hist_name );
      MnvH* ScaleAverageUnfold( MnvH1D* avg_unfold, int num_iter );
      MnvH* ScaleAverageUnfold( MnvH2D* avg_unfold, int num_iter );
      MnvH* ScaleAverageUnfold( MnvH3D* avg_unfold, int num_iter );
      void FillChi2Dists( int stat_uni, MnvH* h_data_unfolded, MnvH* input_truth, MnvH* unfoldtruth, int num_iter );

      double m_stat_scale;
      bool doTailAvgThrows;
      uint m_nStatUniverses;
      double m_ChiSquareGuess;
      double m_ChiSquareStep;
      std::vector<int> m_iterations;
      std::vector<int> m_iterations_bins;

      std::vector<int> m_exclude_chi2_bins; 

      MinervaUnfold::MnvUnfold unfold;
      MnvPlotter *plotter;
      TRandom3 *myrandom;

      MnvH* m_data_truth;
      MnvH* m_data;
      MnvH* m_reco;
      MnvH* m_reco_original;
      MnvH* m_truth;
      MnvH* m_truth_original;
      MnvH2D* m_migration;
      MnvH2D* m_migration_original;
      
      std::map< int, std::vector< MnvH* > > m_data_pull;

      std::map< int, MnvH* > m_avg_data_error;
      std::map< int, MnvH* > m_avg_data_pull;
      std::map< int, MnvH2D* > m_bin_data_errors;
      std::map< int, MnvH2D* > m_bin_data_pulls;
      std::map< int, TProfile* > m_bin_avg_data_errors;
      std::map< int, TProfile* > m_bin_avg_data_pulls;

      std::map< int, std::vector< MnvH* > >   m_ratio_md_td;
      std::map< int, std::vector< MnvH* > >   m_ratio_md_input;
      std::map< int, MnvH2D* >   m_bin_iter_ratio_md_td;
      std::map< int, MnvH2D* >   m_bin_iter_ratio_md_input;
      std::map< int, TProfile* > m_bin_iter_avg_ratio_md_td;
      std::map< int, TProfile* > m_bin_iter_avg_ratio_md_input;

      std::map< int, TMatrixD* > m_avg_unfoldingCovMatrices;
      std::map< int, MnvH* >     m_avg_data_unfolded;

      std::map< int, std::vector<TMatrixD> > m_unfoldingCovMatrices;
      std::map< int, std::vector<MnvH*> > m_outputresults;
      std::map< int, std::vector<MnvH*> > m_inputresults;
      std::map< int, std::vector<TMatrixD*> > m_valmat_vec;
      
      //std::map< int, TMatrixD* > m_avg_valmat_vec;
      std::map< int, MnvH2D* >   m_bin_chi2_iter;
      std::map< int, TProfile* > m_avg_bin_chi2_iter;
    
      TProfile* m_avg_chi2_md_md_iter_chi2;
      TProfile* m_avg_chi2_md_td_iter_chi2;
      TProfile* m_avg_chi2_md_tmc_iter_chi2;
      TProfile* m_avg_chi2_td_tmc_iter_chi2;

      MnvH2D* m_chi2_md_md_iter_chi2;
      MnvH2D* m_chi2_md_td_iter_chi2;
      MnvH2D* m_chi2_md_tmc_iter_chi2;
      MnvH2D* m_chi2_td_tmc_iter_chi2;

      MnvH2D* m_chi2_md_md_iter_stat;
      MnvH2D* m_chi2_md_td_iter_stat;
      MnvH2D* m_chi2_md_tmc_iter_stat;
      MnvH2D* m_chi2_td_tmc_iter_stat;

      TProfile* m_avg_chi2_md_td_iter_chi2_trunc;
      MnvH2D* m_chi2_md_td_iter_chi2_trunc;

      MnvH1D* m_median_chi2_md_td_iter_chi2;
      MnvH2D* m_chi2_md_td_iter_chi2_percentile;
  };
};

#endif
