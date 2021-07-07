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
ROOT.gROOT.ProcessLine(".L {0}/NSF_ValidationSuite/style/myPlotStyle.h".format(PLOTUTILSROOT))
ROOT.myPlotStyle()

# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)

############################################################################################## Preamble above

def produceComparisonPlots(hist_user,hist_reference,histString,lowerBoundX,upperBoundX,horizontalAxisLabel,plotDir):

  ## Set vertical axis bounds
  lowerBoundY = 0.9
  upperBoundY = 1.1

  # For ratio plots
  lineAt1 = ROOT.TLine(lowerBoundX,1,upperBoundX,1)
  lineAt1.SetLineColor(ROOT.kRed+1)
  lineAt1.SetLineWidth(2)
  
  #############################################################################################################
  ### Plot CV #################################################################################################
  #############################################################################################################
  
  # Pull out CV
  cvHist_user = hist_user.GetCVHistoWithStatError()
  cvHist_reference = hist_reference.GetCVHistoWithStatError()
      
  with makeEnv_TCanvas('{0}/CV_inclusive_{1}.png'.format(plotDir,histString)):
    cvHist_user.SetMarkerColor(ROOT.kRed)
    cvHist_user.SetLineColor(ROOT.kRed)
    cvHist_user.GetXaxis().SetTitle(horizontalAxisLabel)
    cvHist_reference.SetMarkerColor(ROOT.kBlue)
    cvHist_reference.SetLineColor(ROOT.kBlue)
    cvHist_reference.Draw()
    cvHist_user.Draw("same") 
   
  hist_ratioCV = cvHist_user.Clone("ratioCV")
  hist_ratioCV.Divide(hist_ratioCV,cvHist_reference)
  hist_ratioCV.GetYaxis().SetRangeUser(lowerBoundY,upperBoundY)
  hist_ratioCV.GetXaxis().SetRangeUser(lowerBoundX,upperBoundX)
  
  hist_ratioCV.GetXaxis().SetTitle(horizontalAxisLabel)
  hist_ratioCV.GetYaxis().SetTitle("CV Ratio (user/reference)")
  
  with makeEnv_TCanvas('{0}/CV-Comparison_ratio_inclusive_{1}_user_to_reference.png'.format(plotDir,histString)):
    hist_ratioCV.Draw()
    lineAt1.Draw()
  
  hist_differenceCV = cvHist_user.Clone("differenceCV")
  hist_differenceCV.Add(cvHist_reference,-1)
  hist_differenceCV.GetYaxis().SetRangeUser(-5,5)
  #hist_differenceCV.GetXaxis().SetRangeUser(lowerBoundX,upperBoundX)
  
  hist_differenceCV.GetXaxis().SetTitle(horizontalAxisLabel)
  hist_differenceCV.GetYaxis().SetTitle("CV Difference (user - reference)")
  
  with makeEnv_TCanvas('{0}/CV-Comparison_difference_inclusive_{1}_user_minus_reference.png'.format(plotDir,histString)):
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
  ### Loop over vertical error bands ##########################################################################
  #############################################################################################################
  
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
                    'RPA_LowQ2',
                    'Muon_Energy_MINERvA',
                    'Muon_Energy_MINOS',
                    'Muon_Energy_Resolution',
                    'BeamAngleX',
                    'BeamAngleY',
                    'Proton_Response',
                    'Pion_Response',
                    'EM_Response',
                    'Other_Response' 
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
      hist_local.GetYaxis().SetRangeUser(lowerBoundY,upperBoundY)
  
      with makeEnv_TCanvas('{0}/errorBand_inclusive_ratio_{1}_{2}.png'.format(plotDir,nameString,errorBand)):
        hist_local_cv = hist_local.GetCVHistoWithStatError()
        hist_local_cv.Draw()
        exec("hist_local.GetVertErrorBand(\"{0}\").DrawAll(\"same\",True)".format(errorBand))
  
    hist_ratio = eb_user.Clone("ratio_{0}".format(errorBand))
    hist_ratio.Divide(eb_user,eb_reference)
    hist_ratio.GetYaxis().SetRangeUser(lowerBoundY,upperBoundY)
    
    with makeEnv_TCanvas('{0}/ratio_inclusive_user_to_reference_{1}.png'.format(plotDir,errorBand)):
      hist_ratio.Draw()
  
    for i in range(nUniverses):
      exec("hist_ratio_universe{0} = eb{0}_user.Clone(\"ratio_{1}_universe0\")".format(i,errorBand))
      exec("hist_ratio_universe{0}.Divide(eb{0}_user,eb{0}_reference)".format(i))
      exec("hist_ratio_universe{0}.GetYaxis().SetRangeUser(lowerBoundY,upperBoundY)".format(i))
      exec("hist_ratio_universe{0}.GetXaxis().SetRangeUser(lowerBoundX,upperBoundX)".format(i))
    
      with makeEnv_TCanvas('{0}/ratio_inclusive_user_to_reference_{1}_universe{2}.png'.format(plotDir,errorBand,i)):
        exec("hist_ratio_universe{0}.Draw()".format(i))
        lineAt1.Draw()

def main():

  print "I'm inside main!"

  #############################################################################################################
  ### User customizations #####################################################################################
  #############################################################################################################
  
  user_Emu_histogram_name = "h_inclusive_Emu"
  user_Ptmu_histogram_name = "h_inclusive_Pt"
  user_recoil_histogram_name = "h_inclusive_Nu"
  
  #############################################################################################################
  ### Parse user arguments and set filepaths ##################################################################
  #############################################################################################################
  
  ## Parse user args
  parser = ArgumentParser(description='Process optional inputs')
  parser.add_argument('--input', dest='inputFilePath', action='store')
  parser.add_argument('--outdir', dest='outDir', action='store')
  parser.add_argument('--compare', dest='compareTo', action='store')
  parser.add_argument('--refHists', dest='userRefHists', action='store',default='')
  OPTS_VEC = parser.parse_args()
  
  ## Output directory
  plotDir = OPTS_VEC.outDir
  if not os.path.isdir(plotDir):
    print "Making plot directory {0}".format(plotDir)
    os.system( "mkdir %s" % plotDir )
  
  ## Reference histograms
  if not OPTS_VEC.userRefHists == '': filePath_reference = OPTS_VEC.userRefHists
  else:
    refDir = '/minerva/data/NSF_Validation/referenceHists/'
    if   OPTS_VEC.compareTo == 'CCQENu':  refFile = 'NSF_MnvGENIEv1_CCQENu_mc_minervame1L_2020-04-17.root'
    elif OPTS_VEC.compareTo == 'NukeCC':  refFile = 'NSF_MnvGENIEv1_NukeCC_mc_minervame1L_2020-04-17.root'
    else:
      print "You asked to compare your histograms to reference histograms that don't exist. I'm exiting gracefully."
      sys.exit()
    
    filePath_reference = "{0}/{1}".format(refDir,refFile)
  
  ## User input histogram file
  filePath_user = OPTS_VEC.inputFilePath
  
  #############################################################################################################
  ### Pull histograms out of input files ######################################################################
  #############################################################################################################
  
  histFile_user = ROOT.TFile(filePath_user)
  histFile_reference = ROOT.TFile(filePath_reference)
  
  hist_Emu_user = histFile_user.Get(user_Emu_histogram_name)
  hist_Ptmu_user = histFile_user.Get(user_Ptmu_histogram_name)
  hist_recoil_user = histFile_user.Get(user_recoil_histogram_name)
  hist_recoil_user.Rebin(10)
  
  hist_Emu_reference = histFile_reference.Get("h_inclusive_Emu")
  hist_Ptmu_reference = histFile_reference.Get("h_inclusive_Pt")
  hist_recoil_reference = histFile_reference.Get("h_inclusive_Nu")
  hist_recoil_reference.Rebin(10)
  
  for [histString,lowerBoundX,upperBoundX,horizontalAxisLabel] in [['Emu',0,120,"Reconstructed E_{#mu} (GeV)"],
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
  
    produceComparisonPlots(hist_user,hist_reference,histString,lowerBoundX,upperBoundX,horizontalAxisLabel,plotSubdir)

if __name__ == "__main__":
  main()

