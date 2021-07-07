## import common python modules
import ROOT,PlotUtils,os,sys
from argparse import ArgumentParser
## modify $PYTHONPATH to know about custom python modules
PLOTUTILSROOT = os.path.expandvars("$PLOTUTILSROOT")
sys.path.append("{0}/NSF_ValidationSuite/py_classes".format(PLOTUTILSROOT))
## import custom python modules
from plottingClasses import * 
from errorMaps import error_bands

## Configure MnvPlotter`  
plotter = PlotUtils.MnvPlotter()
# Manually override default error summary groups
plotter.error_summary_group_map.clear()
for group in error_bands:
  for error in error_bands[group]:
    plotter.error_summary_group_map[group].push_back(error)

## Miscellaneous Plotting prep
plotter.axis_maximum_group = 0.01

##set ROOT to batch mode
ROOT.gROOT.SetBatch()
#Load and implement Phil's plot style header file
ROOT.gROOT.ProcessLine(".L ../style/myPlotStyle.h")
ROOT.myPlotStyle()

# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)

############################################################################################## Preamble above

def produceComparisonPlots(hist_user,hist_reference,histString,lowerBoundX,upperBoundX,horizontalAxisLabel,plotDir,loopWeights):

  # For ratio plots
  lineAt1 = ROOT.TLine(lowerBoundX,1,upperBoundX,1)
  lineAt1.SetLineColor(ROOT.kRed+1)
  lineAt1.SetLineWidth(2)
  
  UNIQUE_USER_VERT_ERROR_BANDS = [
    'Michel_Efficiency',
    'Proton_TrackEff',
    'Reweight_Neutron',
    'Reweight_Pion',
    'Reweight_Proton',
    'Target_Mass',
  ]
  
  UNIQUE_USER_LAT_ERROR_BANDS = [
    'Bethe_Bloch',
    'Birks_Response_Proton',
    'Crosstalk',
    'MEU_Proton',
    'Mass_Model_Proton'
  ]
  
  for eb in UNIQUE_USER_VERT_ERROR_BANDS:
    hist_user.PopVertErrorBand(eb)
  for eb in UNIQUE_USER_LAT_ERROR_BANDS:
    hist_user.PopLatErrorBand(eb)
  
  #############################################################################################################
  ### Plot CV #################################################################################################
  #############################################################################################################
  
  # Pull out CV
  cvHist_user = hist_user.GetCVHistoWithStatError()
  cvHist_reference = hist_reference.GetCVHistoWithStatError()
      
  with makeEnv_TCanvas('{0}/CV_inclusive_Tmu.png'.format(plotDir)):
    cvHist_user.SetMarkerColor(ROOT.kRed)
    cvHist_user.SetLineColor(ROOT.kRed)
    cvHist_user.GetXaxis().SetTitle(horizontalAxisLabel)
    #cvHist_user.Draw()
    cvHist_reference.SetMarkerColor(ROOT.kBlue)
    cvHist_reference.SetLineColor(ROOT.kBlue)
    #cvHist_reference.Draw("same")
    cvHist_reference.Draw()
    cvHist_user.Draw("same") 
   
  hist_ratioCV = cvHist_user.Clone("ratioCV")
  hist_ratioCV.Divide(hist_ratioCV,cvHist_reference)
  hist_ratioCV.GetYaxis().SetRangeUser(0.99,1.01)
  hist_ratioCV.GetXaxis().SetRangeUser(lowerBoundX,upperBoundX)
  
  hist_ratioCV.GetXaxis().SetTitle(horizontalAxisLabel)
  hist_ratioCV.GetYaxis().SetTitle("CV Ratio (user/reference)")
  
  with makeEnv_TCanvas('{0}/CV-Comparison_ratio_inclusive_Tmu_user_to_reference.png'.format(plotDir)):
    hist_ratioCV.Draw()
    lineAt1.Draw()
  
  hist_differenceCV = cvHist_user.Clone("differenceCV")
  hist_differenceCV.Add(cvHist_reference,-1)
  hist_differenceCV.GetYaxis().SetRangeUser(-5,5)
  #hist_differenceCV.GetXaxis().SetRangeUser(lowerBoundX,upperBoundX)
  
  hist_differenceCV.GetXaxis().SetTitle(horizontalAxisLabel)
  hist_differenceCV.GetYaxis().SetTitle("CV Difference (user - reference)")
  
  with makeEnv_TCanvas('{0}/CV-Comparison_difference_inclusive_Tmu_user_minus_reference.png'.format(plotDir)):
    hist_differenceCV.Draw()
  
  with makeEnv_TCanvas('{0}/CV_errorSummary_user.png'.format(plotDir)) as canvas:
    hist_user_local = hist_user.Clone("hist_user_local")
    hist_user_local.GetXaxis().SetRangeUser(lowerBoundX,upperBoundX)
    hist_user_local.GetXaxis().SetTitle("{0} [user]".format(horizontalAxisLabel))
    localDrawErrorSummary(plotter,hist_user_local) 
   
  with makeEnv_TCanvas('{0}/CV_errorSummary_reference.png'.format(plotDir)):
    hist_reference_local = hist_reference.Clone("hist_reference_local")
    hist_reference_local.GetXaxis().SetRangeUser(lowerBoundX,upperBoundX)
    hist_reference_local.GetXaxis().SetTitle("{0} [reference]".format(horizontalAxisLabel))
    localDrawErrorSummary(plotter,hist_reference_local)
  
  #############################################################################################################
  ### Compare sets of error bands #############################################################################
  #############################################################################################################
  
  ## Read set of error bands for each
  print "reference's Error Bands: "
  errorBands_reference = hist_reference.GetVertErrorBandNames()
  for eb in errorBands_reference:
    print "eb: " , eb
  
  print "user's Vertical Error Bands: "
  errorBands_user = hist_user.GetVertErrorBandNames()
  for eb in errorBands_user:
    print "eb: " , eb
  print "user's Lateral Error Bands: "
  errorBands_user = hist_user.GetLatErrorBandNames()
  for eb in errorBands_user:
    print "eb: " , eb
  
  #############################################################################################################
  ### Loop over vertical error bands ##########################################################################
  #############################################################################################################
  
  #for fakeFold in [0]:
  if True:
    for errorBand in ['Flux',
                      'GENIE_AGKYxF1pi',
                      'GENIE_AhtBY',
                      'GENIE_BhtBY',
                      'GENIE_CCQEPauliSupViaKF',
                      'GENIE_CV1uBY',
                      'GENIE_CV2uBY',
                      'GENIE_EtaNCEL',
                      'GENIE_FrAbs_N',
                      'GENIE_FrAbs_pi',
                      'GENIE_FrCEx_N',
                      'GENIE_FrCEx_pi',
                      'GENIE_FrElas_N',
                      'GENIE_FrElas_pi',
                      'GENIE_FrInel_N',
                      'GENIE_FrPiProd_N',
                      'GENIE_FrPiProd_pi',
                      'GENIE_MFP_N',
                      'GENIE_MFP_pi',
                      'GENIE_MaCCQEshape',
                      'GENIE_MaNCEL',
                      'GENIE_MaRES',
                      'GENIE_MvRES',
                      'GENIE_NormCCQE',
                      'GENIE_NormDISCC',
                      'GENIE_NormNCRES',
                      'GENIE_RDecBR1gamma',
                      'GENIE_Rvn1pi',
                      'GENIE_Rvn2pi',
                      'GENIE_Rvp1pi',
                      'GENIE_Rvp2pi',
                      'GENIE_Theta_Delta2Npi',
                      'GENIE_VecFFCCQEshape',
                      'MINOS_Reconstruction_Efficiency',
                      'Low_Recoil_2p2h_Tune',
                      'RPA_HighQ2',
                      'RPA_LowQ2'
                      ]:
    
      exec("nUniverses = hist_user.GetVertErrorBand(\"{0}\").GetNHists()".format(errorBand))
    
      #for hist,nameString in [[hist_user,"user"],[hist_reference,"reference"]]:
      for hist,nameString in [[hist_reference,"reference"],[hist_user,"user"]]:
        
        # Set hist specs
        exec("hist.GetVertErrorBand(\"{0}\").GetXaxis().SetTitle(horizontalAxisLabel)".format(errorBand))
        exec("hist.GetVertErrorBand(\"{0}\").GetXaxis().SetRangeUser(lowerBoundX,upperBoundX)".format(errorBand))
      
        # Pull out error bands and save to use in this loop as well as outside of this loop
        exec("eb_{0}  = hist.GetVertErrorBand(\"{1}\").GetErrorBand(True)".format(nameString,errorBand))
    
        for i in range(nUniverses):
          exec("eb{0}_{1} = hist.GetVertErrorBand(\"{2}\").GetHist({0})".format(i,nameString,errorBand))
    
        exec("hist_local = hist.Clone(\"hist_local_{0}\")".format(nameString))
        
        with makeEnv_TCanvas('{0}/errorBand_inclusive_{1}_{2}.png'.format(plotDir,nameString,errorBand)):
          exec("hist_local.GetVertErrorBand(\"{0}\").DrawAll(\"\",True)".format(errorBand))
          exec("cvHist_local = cvHist_{0}.Clone(\"cvHist_local_{0}\")".format(nameString))
          cvHist_local.SetMarkerColor(ROOT.kGreen)
          cvHist_local.SetLineColor(ROOT.kGreen)
          cvHist_local.Draw("same")
    
        with makeEnv_TCanvas('{0}/errorBand_inclusive_{1}_{2}_fractional.png'.format(plotDir,nameString,errorBand)):
          exec("fracEB = hist_local.GetVertErrorBand(\"{0}\").GetErrorBand(True)".format(errorBand))
          fracEB.Draw()
  
        exec("hist_local.DivideSingle(hist_local,cvHist_{0})".format(nameString))
        hist_local.GetXaxis().SetTitle(horizontalAxisLabel)
        hist_local.GetXaxis().SetRangeUser(lowerBoundX,upperBoundX)
        hist_local.GetYaxis().SetRangeUser(0.99,1.01)
    
        with makeEnv_TCanvas('{0}/errorBand_inclusive_ratio_{1}_{2}.png'.format(plotDir,nameString,errorBand)):
          hist_local_cv = hist_local.GetCVHistoWithStatError()
          hist_local_cv.Draw()
          exec("hist_local.GetVertErrorBand(\"{0}\").DrawAll(\"same\",True)".format(errorBand))
    
      hist_ratio = eb_user.Clone("ratio_{0}".format(errorBand))
      hist_ratio.Divide(eb_user,eb_reference)
      hist_ratio.GetYaxis().SetRangeUser(0.99,1.01)
      
      with makeEnv_TCanvas('{0}/ratio_inclusive_user_to_reference_{1}.png'.format(plotDir,errorBand)):
        hist_ratio.Draw()
    
      for i in range(nUniverses):
        exec("hist_ratio_universe{0} = eb{0}_user.Clone(\"ratio_{1}_universe0\")".format(i,errorBand))
        exec("hist_ratio_universe{0}.Divide(eb{0}_user,eb{0}_reference)".format(i))
        exec("hist_ratio_universe{0}.GetYaxis().SetRangeUser(0.99,1.01)".format(i))
        exec("hist_ratio_universe{0}.GetXaxis().SetRangeUser(lowerBoundX,upperBoundX)".format(i))
      
        with makeEnv_TCanvas('{0}/ratio_inclusive_user_to_reference_{1}_universe{2}.png'.format(plotDir,errorBand,i)):
          exec("hist_ratio_universe{0}.Draw()".format(i))
          lineAt1.Draw()
  
  #############################################################################################################
  ### Loop over lateral error bands ###########################################################################
  #############################################################################################################
  
  #for fakeFold in [0]:
  if True:
    for errorBand in ['Muon_Energy_MINERvA',
                      'Muon_Energy_MINOS',
                      'Muon_Energy_Resolution',
                      'BeamAngleX',
                      'BeamAngleY',
                      'Proton_Response',
                      'Pion_Response',
                      'EM_Response',
                      'Other_Response' 
                     ]:
    
      exec("nUniverses = hist_user.GetLatErrorBand(\"{0}\").GetNHists()".format(errorBand))
      
      for hist,nameString in [[hist_user,"user"],[hist_reference,"reference"]]:
    
        errorBandType = "Lat" if nameString == "user" else "Vert"
    
        # Set hist specs
        exec("hist.Get{0}ErrorBand(\"{1}\").GetXaxis().SetTitle(horizontalAxisLabel)".format(errorBandType,errorBand))
        exec("hist.Get{0}ErrorBand(\"{1}\").GetXaxis().SetRangeUser(lowerBoundX,upperBoundX)".format(errorBandType,errorBand))
      
        # Pull out error bands and save to use in this loop as well as outside of this loop
        exec("eb_{0}  = hist.Get{1}ErrorBand(\"{2}\").GetErrorBand(True)".format(nameString,errorBandType,errorBand))
    
        for i in range(nUniverses):
          exec("eb{0}_{1} = hist.Get{2}ErrorBand(\"{3}\").GetHist({0})".format(i,nameString,errorBandType,errorBand))
     
        exec("hist_local = hist.Clone(\"hist_{0}_local\")".format(nameString))
        
        with makeEnv_TCanvas('{0}/errorBand_inclusive_{1}_{2}.png'.format(plotDir,nameString,errorBand)):
          if nameString == 'user':
            exec("hist_local.GetLatErrorBand(\"{1}\").DrawAll(\"\",True)".format(nameString,errorBand))
          else:
            exec("hist_local.GetVertErrorBand(\"{1}\").DrawAll(\"\",True)".format(nameString,errorBand))
          exec("cvHist_local = cvHist_{0}.Clone(\"cvHist_local\")".format(nameString))
          cvHist_local.SetMarkerColor(ROOT.kGreen)
          cvHist_local.SetLineColor(ROOT.kGreen)
          cvHist_local.Draw("same")
    
        exec("hist_local.DivideSingle(hist_local,cvHist_{0})".format(nameString))
        hist_local.GetXaxis().SetTitle(horizontalAxisLabel)
        hist_local.GetXaxis().SetRangeUser(lowerBoundX,upperBoundX)
        hist_local.GetYaxis().SetRangeUser(0.99,1.01)
    
        with makeEnv_TCanvas('{0}/errorBand_inclusive_ratio_{1}_{2}.png'.format(plotDir,nameString,errorBand)):
          hist_local_cv = hist_local.GetCVHistoWithStatError()
          hist_local_cv.Draw()
          if nameString == 'user':
            exec("hist_local.GetLatErrorBand(\"{0}\").DrawAll(\"same\",True)".format(errorBand))
          else:
            exec("hist_local.GetVertErrorBand(\"{0}\").DrawAll(\"same\",True)".format(errorBand))
    
      hist_ratio = eb_user.Clone("ratio_{0}".format(errorBand))
      hist_ratio.Divide(eb_user,eb_reference)
      hist_ratio.GetYaxis().SetRangeUser(0.99,1.01)
      
      with makeEnv_TCanvas('{0}/ratio_inclusive_user_to_reference_{1}.png'.format(plotDir,errorBand)):
        hist_ratio.Draw()
    
      for i in range(nUniverses):
        exec("hist_ratio_universe{0} = eb{0}_user.Clone(\"ratio_{1}_universe0\")".format(i,errorBand))
        exec("hist_ratio_universe{0}.Divide(eb{0}_user,eb{0}_reference)".format(i))
        exec("hist_ratio_universe{0}.GetYaxis().SetRangeUser(0.99,1.01)".format(i))
        exec("hist_ratio_universe{0}.GetXaxis().SetRangeUser(lowerBoundX,upperBoundX)".format(i))
      
        with makeEnv_TCanvas('{0}/ratio_inclusive_user_to_reference_{1}_universe{2}.png'.format(plotDir,errorBand,i)):
          exec("hist_ratio_universe{0}.Draw()".format(i))
          lineAt1.Draw()
  
  #############################################################################################################
  ### Loop over weight histograms #############################################################################
  #############################################################################################################
  
  if loopWeights:
  
    for wgtDef,danHistName in [['WgtVsPmuMinos_MINOS','h_mwgt_minosp_CC'],
                               ['WgtVsTmu_MINOS','h_mwgt_mut_CC'],
                               ['WgtVsTmu_GENIE',''],
                               ['WgtVsTmu_RPA','h_rpa_mut_CC'],
                               ['WgtVsTmu_2p2h','h_lowrecoil2p2h_mut_CC'],
                               ['WgtVsTmu_Flux','h_frw_mut_CC']
                              ]:
  
      exec("hist_reference_{0} = histFile_reference.Get(\"h_inclusive_{0}\")".format(wgtDef))
      with makeEnv_TCanvas('{0}/{1}_reference.png'.format(plotDir,wgtDef)):
        exec("hist_reference_{0}.Draw(\"colz\")".format(wgtDef))
  
      if not danHistName == '':
        exec("hist_user_{0} = histFile_user.Get(\"{1}\")".format(wgtDef,danHistName))
        with makeEnv_TCanvas('{0}/{1}_user.png'.format(plotDir,wgtDef)):
          exec("hist_user_{0}.Draw(\"colz\")".format(wgtDef))
  
        exec("hist_user_{0}.Add(hist_reference_{0},-1)".format(wgtDef))
        with makeEnv_TCanvas('{0}/WgtComparison_difference_{1}_userMinusreference.png'.format(plotDir,wgtDef)):
          exec("hist_user_{0}.GetZaxis().SetRangeUser(-3,3)".format(wgtDef))
          exec("hist_user_{0}.Draw(\"colz\")".format(wgtDef))
  
    # MINOS wgt
    hist_reference_Wgt_MINOS = histFile_reference.Get("h_inclusive_WgtVsPmuMinos_MINOS").ProjectionX()
    hist_user_Wgt_MINOS = histFile_user.Get("h_mwgt_minosp_CC").ProjectionX()
    
    with makeEnv_TCanvas('{0}/Wgt_MINOS_reference.png'.format(plotDir)):
      hist_reference_Wgt_MINOS.Draw("colz")
    with makeEnv_TCanvas('{0}/Wgt_MINOS_user.png'.format(plotDir)):
      hist_user_Wgt_MINOS.Draw("colz")
    
    hist_user_Wgt_MINOS.Add(hist_reference_Wgt_MINOS,-1)
    with makeEnv_TCanvas('{0}/WgtComparison_difference_Wgt_MINOS_userMinusreference.png'.format(plotDir)):
      hist_user_Wgt_MINOS.GetYaxis().SetRangeUser(-10,10)
      hist_user_Wgt_MINOS.Draw("colz")

def main():

  print "I'm inside main!"

  #############################################################################################################
  ### User customizations #####################################################################################
  #############################################################################################################
  
  user_Tmu_histogram_name = "h_mut_CC"
  user_Ptmu_histogram_name = "h_ptmu_CC"
  user_recoil_histogram_name = "h_recoil_CC"
  
  #############################################################################################################
  ### Parse user arguments and set filepaths ##################################################################
  #############################################################################################################
  
  ## Parse user args
  parser = ArgumentParser(description='Process optional inputs')
  parser.add_argument('--input', dest='inputFilePath', action='store')
  parser.add_argument('--outdir', dest='outDir', action='store')
  parser.add_argument('--compare', dest='compareTo', action='store')
  OPTS_VEC = parser.parse_args()
  
  ## Output directory
  plotDir = OPTS_VEC.outDir
  if not os.path.isdir(plotDir):
    print "Making plot directory {0}".format(plotDir)
    os.system( "mkdir %s" % plotDir )
  
  ## Reference histograms
  refDir = '/minerva/data/NSF_Validation/referenceHists/'
  refFile = 'NSF_compareToDan_CCQENu_mc_minervame1L_2020-04-09.root ' ## Specific to the original test against CCQENu
  filePath_reference = "{0}/{1}".format(refDir,refFile)
 
  ## User input histogram file
  filePath_user = OPTS_VEC.inputFilePath
  
  #############################################################################################################
  ### Pull histograms out of input files ######################################################################
  #############################################################################################################
  
  histFile_user = ROOT.TFile(filePath_user)
  histFile_reference = ROOT.TFile(filePath_reference)
  
  hist_Tmu_user = histFile_user.Get(user_Tmu_histogram_name)
  hist_Ptmu_user = histFile_user.Get(user_Ptmu_histogram_name)
  hist_recoil_user = histFile_user.Get(user_recoil_histogram_name)
  hist_recoil_user.Rebin(10)
  
  hist_Tmu_reference = histFile_reference.Get("h_inclusive_Emu")
  hist_Ptmu_reference = histFile_reference.Get("h_inclusive_Pt")
  hist_recoil_reference = histFile_reference.Get("h_inclusive_Nu")
  hist_recoil_reference.Rebin(10)
  
  for [histString,lowerBoundX,upperBoundX,horizontalAxisLabel] in [['Tmu',0,120,"Reconstructed E^{kin}_{#mu} (GeV)"],
                                                                   ['Ptmu',0,2.5,"Reconstructed P^{T}_{#mu} (GeV)"],
                                                                   ['recoil',0,5000,"Recoil Energy (MeV)"]
                                                                  ]:
    exec("hist_user = hist_{0}_user".format(histString))
    exec("hist_reference = hist_{0}_reference".format(histString))
  
    ## Output subdirectory
    plotSubdir = "{0}/{1}".format(plotDir,histString)
    if not os.path.isdir(plotSubdir):
      print "Making plot directory {0}".format(plotSubdir)
      os.system( "mkdir %s" % plotSubdir )
  
    produceComparisonPlots(hist_user,hist_reference,histString,lowerBoundX,upperBoundX,horizontalAxisLabel,plotSubdir,False)

if __name__ == "__main__":
  main()

