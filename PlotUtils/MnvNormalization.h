#ifndef MNV_MnvNormalization
#define MNV_MnvNormalization 1

#include <string>
#include <map>

/*
This file defines a namespace MnvNorm and a class MnvNormalizer.

MnvNorm defines normalization corrections applied to observed 
  MC events to account for data/MC discrepancy.


MnvNormalizer is a class that uses the constants in MnvNorm
  and can shift between playlists and analyses.
  Use the class to get your correction factors.

See Tech Note 16: https://neutrino.otterbein.edu/Glaucus/web/view_technote.cgi?technote_id=16

NOTE: reconstruction efficiency affects reconstructed MC only
(the numerator of acceptance); mass model scale affects both
reconstructed MC and truth (the numerator and denominator
of acceptance).

//=====================================================
// Usage
//=====================================================
When loading corrections for a specific playlist and analysis,
  the corrections for the playlist are applied first.

--------------------------------
example - get correction and error
--------------------------------
MnvNormalizer norm;
for( int i = 0; i != nEvents; ++i )
{
  tuple.GetEntry(i);
  const double momentum = tuple.momentum;
  const double normCorrection    = norm.Correction( momentum );
  const double normCorrectionErr = norm.CorrectionErr( momentum );
  hist->Fill( momentum, cv_weight*normCorrection );

  double normErrs = { -normCorrectionErr, normCorrectionErr };
  hist->FillVertErrorBand( "Normalization", momentum, normErrs, cv_weight );
}

--------------------------------
example - specify playlist
--------------------------------
MnvNormalizer norm("minerva5");
--- or ---
MnvNormalizer norm;
norm.LoadPlaylist("minerva5");

--------------------------------
example - specify analysis
--------------------------------
MnvNormalizer norm;
norm.LoadAnalysis("ccqe");

--------------------------------
example - specify playlist and analysis
--------------------------------
MnvNormalizer norm( "minerva5", "ccqe" );
--- or ---
MnvNormalizer norm( "minerva5" );
norm.LoadAnalysis( "ccqe" );
--- or ---
MnvNormalizer norm;
norm.LoadPlaylist("minerva5");
norm.LoadAnalysis( "ccqe" );

*/

namespace PlotUtils
{
  /*! Constants used to account for data/MC discrepancies

    For documentation and explanation of values see docdb7573: 
http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=7573
http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646 - TN016

Abbreviations used are:
<table>
<tr><th> Abbr.</th><th>Meaning</th></tr>
<tr><td>norm</td><td>Normalization</td></tr>
<tr><td>err</td><td>Error</td></tr>
<tr><td>reco</td><td>Reconstruction</td></tr>
<tr><td>eff</td><td>Efficiency</td></tr>
<tr><td>acc</td><td>Acceptance</td></tr>
</table>
   */
  namespace MnvNorm
  {
    //==========================================
    // default corrections and errors
    // note: defaults match playlist minerva1 CCInclusive
    //==========================================
    // Mass model scale + error
    // see: http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=6016
    const double mass_scale     = 1.0;
    const double mass_scale_err = .014;

    // Energy threshold between MINOS tracking efficiency lowP and highP
    const double minosP_minos_trk_eff_threshP = 3000.; //MeV/c

    // Dead time / pile-up in MINERvA + error due to tdead cut
    const double dead_time_tdead     = 1.0;
    const double dead_time_tdead_err = 0.0;

    //==========================================
    // analysis-specific corrections
    //==========================================

    namespace resurrection{ 

      namespace minerva1 {
        
        // MINERvA muon tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double mnv_mu_reco_eff = 1.0;  
        const double mnv_mu_reco_eff_err = 0.002; 
        
        // MINOS muons tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double minos_mu_reco_eff_lowp     = 0.934;
        const double minos_mu_reco_eff_lowp_err  = 0.002;
        const double minos_mu_reco_eff_highp     = 0.982;
        const double minos_mu_reco_eff_highp_err = 0.001;
 
    }
    
    namespace minerva1low {

        // MINERvA muon tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double mnv_mu_reco_eff = 1.0;
//      const double mnv_mu_reco_eff_err = 0.011;
        
        // MINOS muons tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double minos_mu_reco_eff_lowp      = .961;
//      const double minos_mu_reco_eff_lowP_err  = .008;
        const double minos_mu_reco_eff_highp     = .994;
//      const double minos_mu_reco_eff_highP_err = .002;
 
    }

    namespace minerva2 {
        
        // MINERvA muon tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double mnv_mu_reco_eff = 1.0;     
//      const double mnv_mu_reco_eff_err = 0.006; 
        
        // MINOS muons tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double minos_mu_reco_eff_lowp      = 0.930;
//      const double minos_mu_reco_eff_lowP_err  = 0.005;
        const double minos_mu_reco_eff_highp     = 0.979;
//      const double minos_mu_reco_eff_highP_err = 0.002;
 
    }

    namespace minerva3 {

        // MINERvA muon tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double mnv_mu_reco_eff = 1.0;  
//      const double mnv_mu_reco_eff_err = 0.007;
        
        // MINOS muons tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double minos_mu_reco_eff_lowp      = 0.931;
//      const double minos_mu_reco_eff_lowP_err  = 0.008;
        const double minos_mu_reco_eff_highp     = 0.988;
//      const double minos_mu_reco_eff_highP_err = 0.002;
        
    }

    namespace minerva4 {

        // MINERvA muon tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double mnv_mu_reco_eff = 1.0;    
//      const double mnv_mu_reco_eff_err = 0.005; 
        
        // MINOS muons tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double minos_mu_reco_eff_lowp      = 0.912;
//      const double minos_mu_reco_eff_lowP_err  = 0.005;
        const double minos_mu_reco_eff_highp     = 0.961;
//      const double minos_mu_reco_eff_highP_err = 0.002;
 
    }

    namespace minerva5 {

        // MINERvA muon tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double mnv_mu_reco_eff = 1.0;     
//      const double mnv_mu_reco_eff_err = 0.01; 
        
        // MINOS muons tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double minos_mu_reco_eff_lowp      = 0.956;
//      const double minos_mu_reco_eff_lowP_err  = 0.002;
        const double minos_mu_reco_eff_highp     = 0.989;
//      const double minos_mu_reco_eff_highP_err = 0.001;
        
    }

    namespace minerva6 {
        
        // MINERvA muon tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double mnv_mu_reco_eff = 1.0;     
//      const double mnv_mu_reco_eff_err = 0.007; 
        
        // MINOS muons tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double minos_mu_reco_eff_lowp      = 0.987;
//      const double minos_mu_reco_eff_lowP_err  = 0.009;
        const double minos_mu_reco_eff_highp     = 0.991;
//      const double minos_mu_reco_eff_highP_err = 0.002;
        
    }

    namespace minerva7 {

        // MINERvA muon tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double mnv_mu_reco_eff = 1.0;     
//      const double mnv_mu_reco_eff_err = 0.006; 
        
        // MINOS muons tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double minos_mu_reco_eff_lowp      = 0.969;
//      const double minos_mu_reco_eff_lowP_err  = 0.009;
        const double minos_mu_reco_eff_highp     = 0.994;
//      const double minos_mu_reco_eff_highP_err = 0.002;
 
    }

    namespace minerva8 {

        // MINERvA muon tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double mnv_mu_reco_eff = 1.0;     
//      const double mnv_mu_reco_eff_err = 0.005; 
        
        // MINOS muons tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double minos_mu_reco_eff_lowp      = 0.941;
//      const double minos_mu_reco_eff_lowP_err  = 0.004;
        const double minos_mu_reco_eff_highp     = 0.972;
//      const double minos_mu_reco_eff_highP_err = 0.001;
 
    }

    namespace minerva9 {

        // MINERvA muon tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double mnv_mu_reco_eff = 1.0;    
//      const double mnv_mu_reco_eff_err = 0.007; 
        
        // MINOS muons tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double minos_mu_reco_eff_lowp      = 0.951;
//      const double minos_mu_reco_eff_lowP_err  = 0.006;
        const double minos_mu_reco_eff_highp     = 0.994;
//      const double minos_mu_reco_eff_highP_err = 0.002;
 
    }

    namespace minerva10 {

        // MINERvA muon tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double mnv_mu_reco_eff = 1.0;     
//      const double mnv_mu_reco_eff_err = 0.004; 
        
        // MINOS muons tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double minos_mu_reco_eff_lowp      = 0.959;
//      const double minos_mu_reco_eff_lowP_err  = 0.003;
        const double minos_mu_reco_eff_highp     = 0.990;
//      const double minos_mu_reco_eff_highP_err = 0.001;
 
    }

    namespace minerva11 {

        // MINERvA muon tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double mnv_mu_reco_eff = 1.0;     
//      const double mnv_mu_reco_eff_err = 0.005; 
        
        // MINOS muons tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double minos_mu_reco_eff_lowp      = 0.921;
//      const double minos_mu_reco_eff_lowP_err  = 0.005;
        const double minos_mu_reco_eff_highp     = 0.975;
//      const double minos_mu_reco_eff_highP_err = 0.002;
 
    }

    namespace minerva12 {

        // MINERvA muon tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double mnv_mu_reco_eff = 1.0;     
//      const double mnv_mu_reco_eff_err = 0.009;
        
        // MINOS muons tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double minos_mu_reco_eff_lowp      = 0.928;
//      const double minos_mu_reco_eff_lowP_err  = 0.008;
        const double minos_mu_reco_eff_highp     = 0.988;
//      const double minos_mu_reco_eff_highP_err = 0.002;
 
    }
  namespace minerva13 {

        // MINERvA muon tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double mnv_mu_reco_eff = 1.0;     
//      const double mnv_mu_reco_eff_err = 0.003; 
        
        // MINOS muons tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double minos_mu_reco_eff_lowp      = 0.942;
//      const double minos_mu_reco_eff_lowP_err  = 0.002;
        const double minos_mu_reco_eff_highp     = 0.987;
//      const double minos_mu_reco_eff_highP_err = 0.001;
 
    }
    
    namespace minerva13A {

        // MINERvA muon tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double mnv_mu_reco_eff = 1.0;     
//      const double mnv_mu_reco_eff_err = 0.003; 
        
        // MINOS muons tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double minos_mu_reco_eff_lowp      = 0.942;
//      const double minos_mu_reco_eff_lowP_err  = 0.002;
        const double minos_mu_reco_eff_highp     = 0.987;
//      const double minos_mu_reco_eff_highP_err = 0.001;
 
    }
 namespace minerva13B {

        // MINERvA muon tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double mnv_mu_reco_eff = 1.0;     
//      const double mnv_mu_reco_eff_err = 0.003; 
        
        // MINOS muons tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double minos_mu_reco_eff_lowp      = 0.942;
//      const double minos_mu_reco_eff_lowP_err  = 0.002;
        const double minos_mu_reco_eff_highp     = 0.987;
//      const double minos_mu_reco_eff_highP_err = 0.001;
 
    }
    
  namespace minerva13C {

        // MINERvA muon tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double mnv_mu_reco_eff = 1.0;     
//      const double mnv_mu_reco_eff_err = 0.003; 
        
        // MINOS muons tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double minos_mu_reco_eff_lowp      = 0.942;
//      const double minos_mu_reco_eff_lowP_err  = 0.002;
        const double minos_mu_reco_eff_highp     = 0.987;
//      const double minos_mu_reco_eff_highP_err = 0.001;
 
    }

 namespace minerva13D {

        // MINERvA muon tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double mnv_mu_reco_eff = 1.0;     
//      const double mnv_mu_reco_eff_err = 0.003; 
        
        // MINOS muons tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double minos_mu_reco_eff_lowp      = 0.942;
//      const double minos_mu_reco_eff_lowP_err  = 0.002;
        const double minos_mu_reco_eff_highp     = 0.987;
//      const double minos_mu_reco_eff_highP_err = 0.001;
 
    }

 namespace minerva13E {

        // MINERvA muon tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double mnv_mu_reco_eff = 1.0;     
//      const double mnv_mu_reco_eff_err = 0.003; 
        
        // MINOS muons tracking efficiency. http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8646, TN016
        const double minos_mu_reco_eff_lowp      = 0.942;
//      const double minos_mu_reco_eff_lowP_err  = 0.002;
        const double minos_mu_reco_eff_highp     = 0.987;
//      const double minos_mu_reco_eff_highP_err = 0.001;
 
 }

    }//end Resurrection namespace
    namespace eroica{
    //MINOS muon tracking efficiency: http://minerva-docdb.fnal.gov:8080/cgi-bin/RetrieveFile?docid=11367&filename=MINOS%20Tracking%20Efficiencies%20in%20Eroica.pdf&version=8
    //MINERvA muon tracking efficiency: http://minerva-docdb.fnal.gov:8080/cgi-bin/RetrieveFile?docid=11624&filename=MinervaTrackMatchingEff1.pdf&version=1

    namespace minerva1 {
    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp   = 0.963; 
    const double minos_mu_reco_eff_highp     = 0.990;    

    }
    namespace minerva1low {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.963;     
    const double minos_mu_reco_eff_highp     = 0.990;      
    
    }

    namespace minerva2 {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.930;
    const double minos_mu_reco_eff_highp     = 0.979;

    } 
    namespace minerva3 {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.931;
    const double minos_mu_reco_eff_highp     = 0.988;

    }
    
    namespace minerva4 {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.912;
    const double minos_mu_reco_eff_highp     = 0.961;

    }
    namespace minerva5 {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.975;
    const double minos_mu_reco_eff_highp     = 0.995;

    }
   
   namespace minerva6 {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.987;
    const double minos_mu_reco_eff_highp     = 0.991;

    }
 
   namespace minerva7 {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.975;
    const double minos_mu_reco_eff_highp     = 0.995;

    }

   namespace minerva8 {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.941;
    const double minos_mu_reco_eff_highp     = 0.972;

    }
  
   namespace minerva9 {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.972;
    const double minos_mu_reco_eff_highp     = 0.996;

    }

   namespace minerva10 {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.975;
    const double minos_mu_reco_eff_highp     = 0.996;

    }

   namespace minerva11 {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.921;
    const double minos_mu_reco_eff_highp     = 0.975;

    }

   namespace minerva12 {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.928;
    const double minos_mu_reco_eff_highp     = 0.988;

    }

    namespace minerva13A {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.971;
    const double minos_mu_reco_eff_highp     = 0.994;

    }

    namespace minerva13B {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.971;
    const double minos_mu_reco_eff_highp     = 0.994;

    }

    namespace minerva13C {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.971;
    const double minos_mu_reco_eff_highp     = 0.994;

    }

    namespace minerva2p2h {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.971;
    const double minos_mu_reco_eff_highp     = 0.994;

    }

    namespace minervanonreweightables {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.971;
    const double minos_mu_reco_eff_highp     = 0.994;

    }

    namespace minerva13D {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.971;
    const double minos_mu_reco_eff_highp     = 0.994;

    }

    namespace minerva13E {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.971;
    const double minos_mu_reco_eff_highp     = 0.994;

    }

    namespace minerva13 {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.971;
    const double minos_mu_reco_eff_highp     = 0.994;

    }

    namespace minervame1A {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.96;
    const double minos_mu_reco_eff_highp     = 0.99;

    }

    namespace minervame1B {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.96;
    const double minos_mu_reco_eff_highp     = 0.99;

    }

    namespace minervame1C {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.97;
    const double minos_mu_reco_eff_highp     = 0.99;

    }

    namespace minervame1D {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.96;
    const double minos_mu_reco_eff_highp     = 0.99;

    }

    namespace minervame1E {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.96;
    const double minos_mu_reco_eff_highp     = 0.99;

    }

    namespace minervame1F {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.96;
    const double minos_mu_reco_eff_highp     = 0.99;

    }

    namespace minervame1G {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.96;
    const double minos_mu_reco_eff_highp     = 0.99;

    }

    namespace minervame1L {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.96;
    const double minos_mu_reco_eff_highp     = 0.99;

    }

    namespace minervame1M {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.96;
    const double minos_mu_reco_eff_highp     = 0.99;

    }

    namespace minervame1A2p2h {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.971;
    const double minos_mu_reco_eff_highp     = 0.994;

    }

    namespace minervame1B2p2h {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.971;
    const double minos_mu_reco_eff_highp     = 0.994;

    }

    namespace minervame1C2p2h {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.971;
    const double minos_mu_reco_eff_highp     = 0.994;

    }

    namespace minervame1D2p2h {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.971;
    const double minos_mu_reco_eff_highp     = 0.994;

    }

    namespace minervame1E2p2h {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.971;
    const double minos_mu_reco_eff_highp     = 0.994;

    }

    namespace minervame1F2p2h {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.971;
    const double minos_mu_reco_eff_highp     = 0.994;

    }

    namespace minervame1G2p2h {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.971;
    const double minos_mu_reco_eff_highp     = 0.994;

    }

    namespace minervame1L2p2h {

    const double mnv_mu_reco_eff =  0.995;
    const double minos_mu_reco_eff_lowp     = 0.971;
    const double minos_mu_reco_eff_highp     = 0.994;

    }

   }//end Eroica namespace
    /////////////////////////////////////
    namespace inextinguishable{
      //Add documents
    namespace minervame1A {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.95;
    const double minos_mu_reco_eff_highp     = 0.99;

    }

    namespace minervame1B {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.95;
    const double minos_mu_reco_eff_highp     = 0.991;

    }

    namespace minervame1C {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.95;
    const double minos_mu_reco_eff_highp     = 0.991;

    }

    namespace minervame1D {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.95;
    const double minos_mu_reco_eff_highp     = 0.991;

    }

    namespace minervame1E {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.94;
    const double minos_mu_reco_eff_highp     = 0.98;

    }

    namespace minervame1F {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.93;
    const double minos_mu_reco_eff_highp     = 0.99;

    }

    namespace minervame1G {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.94;
    const double minos_mu_reco_eff_highp     = 0.99;

    }

    namespace minervame1L {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.93;
    const double minos_mu_reco_eff_highp     = 0.99;

    }

    namespace minervame1M {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.93;
    const double minos_mu_reco_eff_highp     = 0.99;

    }

    namespace minervame1N {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.92;
    const double minos_mu_reco_eff_highp     = 0.98;

    }

    namespace minervame1O {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.92;
    const double minos_mu_reco_eff_highp     = 0.98;

    }

    namespace minervame1P {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.92;
    const double minos_mu_reco_eff_highp     = 0.98;

    }
    namespace minervame5A {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.92;
    const double minos_mu_reco_eff_highp     = 0.98;

    }
    namespace minervame6A {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.92;
    const double minos_mu_reco_eff_highp     = 0.98;

    }
    namespace minervame6B {

    const double mnv_mu_reco_eff =  1;
    const double minos_mu_reco_eff_lowp     = 0.92;
    const double minos_mu_reco_eff_highp     = 0.98;

    }


   }//end Inextinguishable namespace
  }//end of MnvNormalization namespace


  //==========================================
  // MnvNormalizer class
  //==========================================
  //! Class to apply normalizations to MC due to data/MC discrepancies
  class MnvNormalizer
  {

     
    struct correction_struct {

        double mnv_mu_reco_eff;
        double mnv_mu_reco_eff_err;

        double minos_mu_reco_eff_lowp;
        double minos_mu_reco_eff_lowp_err;
        double minos_mu_reco_eff_highp;
        double minos_mu_reco_eff_highp_err;

    };

    public:
      //! Default constructor
     // MnvNormalizer();

      //! Specify playlist and possibly analysis
      explicit MnvNormalizer( const std::string& processing_name, const std::string& playlist );

      //! Default destructor
      virtual ~MnvNormalizer();


      //=========================================================
      // publicly accessible normalization corrections and errors
      // todo: make them read-only to outsiders
      //=========================================================
      // Mass model scale + error
      double mass_scale;
      double mass_scale_err;

      //Dead time / pile-up in MINERvA + error due to tdead cut.
      double dead_time_tdead;
      double dead_time_tdead_err;

      /*! Load default corrections/errors for this playlist.
        @param[in] playlist name of the playlist to load
        @return nothing.  Throws if playlist name is not found.
       */
      void LoadPlaylist( const std::string& playlist_name );

      /*! Load default corrections/errors for this analysis.
        @param[in] analysis name of the analysis to load
        @return nothing.  Throws if analysis name is not found.
       */
      void LoadAnalysis( const std::string& analysis );

      /*! Load corrections specific for the this playlist-analysis pair.
        @param[in] playlist name of the playlist to load
        @param[in] analysis name of the analysis to load
        @return nothing.  Throws if analysis name is not found.
        */
      void Load( const std::string& playlist_name, const std::string& analysis );

      /*! Get the total correction factor
      @param[in] minosP momentum of muon at face of MINOS in MeV/c
      @return total correction factor
      */
      double Correction( double minosP ) const;
      double GetCorrection(double minosP) const;
      
      /*! Get the error on the total correction factor
      @param[in] minosP momentum of muon at face of MINOS in MeV/c
      @return absolute error on the correction factor
      */
      double CorrectionErr( double minosP ) const;
      double GetCorrectionErr(double minosP) const;

      /*! Get the mass model scale
      @return mass model scale
      */
      double MassModelScale() const;
      double GetMassModelScale() const;

      /*! Get the error on the mass model scale
      @return error on the mass model scale
      */
      double MassModelErr() const;
      double GetMassModelErr() const;

      void ReadCorrectionTable();
    private:

      std::string m_processing;

      std::string playlist_;
      std::string analysis_;

      correction_struct* __the_correction_struct;

      std::map<std::string, std::map<std::string, correction_struct*> > __correction_table;
 

 }; // end class MnvNormalizer

}//end namespace PlotUtils

#endif //MNV_MnvNormization
