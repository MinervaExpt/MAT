############# This is a script to rebin flux MnvH1Ds and to reorder the "Flux" universes therein by the weight
############# prescribed by the nue constraint on a per-universe basis (in 'weightsFilePath').
############# Authored by RDF in Oct 2019 -- see DocDB #24738 for context and details
#############
#############################################################################################################
#############################################################################################################
#############################################################################################################

import ROOT,PlotUtils,os
from array import *

##set ROOT to batch mode
ROOT.gROOT.SetBatch()

# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)

plotter = PlotUtils.MnvPlotter()

# Manually define the "Flux" error group within the MnvPlotter config so that DrawErrorSummary does something pretty
plotter.error_summary_group_map.clear()
plotter.error_summary_group_map["Flux"].push_back("Flux")

def writeHist(hist,outFile):
  outFile.cd()
  print "Writing {0} to output file".format(hist.GetName())
  hist.Write()

class makeEnv_TCanvas(object):

  def __init__(self,plotName,logy=False):
    self.plotName = plotName
    self.logy = logy

  def __enter__(self):
    self.canvas = ROOT.TCanvas( "canvas" , "canvas" , 10 , 10 , 1000, 750 )
    if self.logy: self.canvas.SetLogy(0)
    return self

  def __exit__(self,*exc):
    # If directory for plotName doesn't exist yet, make it
    plotNameComponents = self.plotName.split('/')
    plotDir = '/'.join(plotNameComponents[:-1]) 
    if not os.path.isdir(plotDir):
      print "Making plot directory {0}".format(plotDir)
      os.system( "mkdir %s" % plotDir )
    self.canvas.Print(self.plotName)
    del self.canvas

############################################################################################## Preamble above
#############################################################################################################

#############################################################################################################
### Specify output pths #####################################################################################
#############################################################################################################

USER = os.path.expandvars("$USER")
MPARAMFILESROOT = os.path.expandvars("$MPARAMFILESROOT")
# path for new weights file
newWeightsFilePath = "{0}/data/FluxConstraints/nu+e_ME_spectrum_temp.txt".format(MPARAMFILESROOT)
# paths to output directories for the new flux files and supporting plots
individualHistDir = "/minerva/data/users/{0}/FluxRearranger_hists".format(USER)
plotDir = "/minerva/data/users/{0}/FluxRearranger_plots".format(USER)
# create the output dirs if they don't exist
if not os.path.isdir(plotDir):
  print "Making plot directory {0}".format(plotDir)
  os.system( "mkdir {0}".format(plotDir) )
if not os.path.isdir(individualHistDir):
  print "Making plot directory {0}".format(individualHistDir)
  os.system( "mkdir {0}".format(individualHistDir) )

#############################################################################################################
### Rebinning ###############################################################################################
#############################################################################################################

# Provide binning breakout and assemble list of bin boundaries
##############################################################

binMap = [#[min,max,size] (all in units of GeV)
          [0,1.5,0.5],
          [1.5,10,0.1],
          [10,20,0.5],
          [20,30,1],
          [30,60,5],
          [60,100,10] # The flux generation only goes up to 100 GeV
]
bins_nuE = []
for minBin,maxBin,binSize in binMap:
  nBins = (maxBin-minBin)/binSize
  for iBin in range(int(nBins)):
    binBoundary = minBin + iBin*binSize
    bins_nuE.append(binBoundary)
bins_nuE.append(100.) # The above prescription doesn't add the last bin
nBins_nuE = len(bins_nuE)-1

print "This is the binning I'm going to use!: {0}".format(bins_nuE)
print "That's {0} bins total: ".format(nBins_nuE)

# Make reference hist that has the binning we want to adopt
###############################################################

referenceHist = PlotUtils.MnvH1D( 'h_flux_reference' , 'h_flux_reference' , nBins_nuE , array('d',bins_nuE))

#############################################################################################################
### Rearrange universes to sort by nu-e constraint weight ###################################################
#############################################################################################################

# Fetch current mapping of universe to wgt and store in dict
############################################################

weightsFilePath = "{0}/data/FluxConstraints/nu+e_ME_spectrum_originalOrdering.txt".format(MPARAMFILESROOT)
weightsFile = open(weightsFilePath,"r")
wgtList = weightsFile.readlines()
wgtDict = {}

for i,entry in enumerate(wgtList):
  elements = entry.split(' ')
  universe = int(elements[1])
  wgt = float(elements[2].rstrip('\n'))
  wgtDict[universe] = wgt

print 'old weight mapping: ' , wgtDict

## # Re-order universe mapping to sort by wgt
## #############################################################################
## 
## sortedWgtDict = {}
## 
## for i,universe in enumerate(sorted(wgtDict,key=wgtDict.get,reverse=True)):
##   #print i,universe,wgtDict[universe]
##   sortedWgtDict[i]=[wgtDict[universe],universe]
## 
## print 'new weight mapping: ' , sortedWgtDict
## 
## # Write new mapping to text file
## ############################################################
## 
## newWeightsFile = open(newWeightsFilePath,"w+")
## 
## for i,universe in enumerate(sorted(wgtDict,key=wgtDict.get,reverse=True)):
##   newWeightsFile.write('Flux {0} {1}\n'.format(i,wgtDict[universe]))
## 
## print 'writing new weights file to {0}'.format(newWeightsFilePath)
## newWeightsFile.close()

# Re-order universe mapping to preserve current randomness, but through out bottom 500 universes by weight
###############################################################################

sortedWgtDict = {}
index = 0
index2 = 500

for universe in wgtDict:
  wgt = wgtDict[universe]
  if wgt < 0.0684478529886:
    print index2,universe,wgtDict[universe]
    sortedWgtDict[index2]=[wgtDict[universe],universe]
    index2+=1
  else:
    print index,universe,wgtDict[universe]
    sortedWgtDict[index]=[wgtDict[universe],universe]
    index+=1

print 'new weight mapping: ' , sortedWgtDict

# Write new mapping to text file
############################################################

newWeightsFile = open(newWeightsFilePath,"w+")

for i,universe in enumerate(sortedWgtDict):
  newWeightsFile.write('Flux {0} {1}\n'.format(i,sortedWgtDict[universe][0]))

print 'writing new weights file to {0}'.format(newWeightsFilePath)
newWeightsFile.close()

#############################################################################################################
### Implement the rebinning and universe shuffle ############################################################
#############################################################################################################

# Fetch original and rebinned flux MnvH1Ds using FluxReweighter
###############################################################

## User should specify playlists and fluxType
############################################

# Loop over all permutations of playlist, reweighted/generate, and PDG code 
#for playlist in ['minervame1D','minervame1M','minervame1N','minervame1D1M1NWeightedAve']:
for playlist in ['minervame5A','minervame6A']:
  for fluxType in ['Reweighted','Generated']:
    for flavor in ['14','-14','12','-12']:
#Debug
#for playlist in ['minervame1D']:
#  for fluxType in ['Reweighted']:
#    for flavor in ['14']:

      exec("fluxReweighter = PlotUtils.FluxReweighter({0},False,PlotUtils.FluxReweighter.{1},PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv6,1000)".format(flavor,playlist))
      exec("fluxHist_originalBinning = fluxReweighter.GetFlux{0}({1})".format(fluxType,flavor))
      fluxHist_originalBinning.SetName('flux_{0}_{1}_{2}_originalBinning'.format(playlist,fluxType,flavor))
      exec("fluxHist_rebinned = fluxReweighter.GetRebinnedFlux{0}({1},referenceHist)".format(fluxType,flavor))
      if fluxType == 'Reweighted':
        fluxHist_rebinned.SetName('flux_E_cvweighted')
      else:
        fluxHist_rebinned.SetName('flux_E_unweighted')
 
      #############################################################################################################
      ### Write rebinned flux MnvH1D to output file ###############################################################
      #############################################################################################################
  
      fluxTag = 'g4numiv6' if fluxType == 'Generated' else 'gen2thin' 
      histOutputFilePath = '{0}/flux-{1}-pdg{2}-{3}.root'.format(individualHistDir,fluxTag,flavor,playlist)
      histOutputFile = ROOT.TFile(histOutputFilePath,"recreate")
      histOutputFile.cd()
      
      writeHist(fluxHist_rebinned,histOutputFile)

      #############################################################################################################
      ### Reorder universes #######################################################################################
      #############################################################################################################
  
      ## Pop these error bands because nobody uses them and everyone has to get rid of them...
      ## Actually, we decided not to do this
      #fluxHist_rebinned.PopVertErrorBand("Flux_BeamFocus")
      #fluxHist_rebinned.PopVertErrorBand("ppfx1_Total")
      
      # This is the error band we wnat to rearrange
      fluxErrorBand_oldOrder = fluxHist_rebinned.GetVertErrorBand('Flux')
      ## Debug
      ebs = fluxHist_rebinned.GetErrorBandNames()     
      for eb in ebs:
        print "eb: ", eb
      nUniverses = fluxErrorBand_oldOrder.GetNHists()
      print "nUniverses: " , nUniverses
 
      # Pop the old Flux error band (with old universe ordering)
      fluxHist_rebinned.PopVertErrorBand("Flux")
      
      # This is the duplicate error band that will hold the rearranged universes
      fluxHist_rebinned.AddVertErrorBandAndFillWithCV("Flux",1000)
      fluxErrorBand_newOrder = fluxHist_rebinned.GetVertErrorBand("Flux")
     
      for newUniverseNumber,wgtList in sortedWgtDict.items(): 
        for i in range(nBins_nuE+2):
          fluxErrorBand_newOrder.GetHist(newUniverseNumber).SetBinContent(i,fluxErrorBand_oldOrder.GetHist(wgtList[1]).GetBinContent(i))
        # While we're here, remove the stat. error from each universe, because it doesn't get used and is just extra memory footprint
        fluxErrorBand_newOrder.GetHist(newUniverseNumber).Sumw2(False)  
   
      #############################################################################################################
      ### Write new flux MnvH1D to output file ####################################################################
      #############################################################################################################
  
      histOutputFilePath_constrained = '{0}/flux-{1}-pdg{2}-{3}_constrained.root'.format(individualHistDir,fluxTag,flavor,playlist)
      histOutputFile_constrained = ROOT.TFile(histOutputFilePath_constrained,"recreate")
      histOutputFile_constrained.cd()
      
      writeHist(fluxHist_rebinned,histOutputFile_constrained)
      
      #############################################################################################################
      ### Make plots to compare before and after ##################################################################
      #############################################################################################################
  
      plotter.axis_maximum = 0.5
   
      with makeEnv_TCanvas('{0}/{1}_{2}_{3}_CVHist_beforeRearrange.png'.format(plotDir,playlist,fluxType,flavor)):
        cvHist_originalBinning = fluxHist_originalBinning.GetCVHistoWithError()
        cvHist_originalBinning.Scale()
        cvHist_originalBinning.GetYaxis().SetRangeUser(-0.5e-5,1.6e-4)
        cvHist_originalBinning.Draw()
 
      with makeEnv_TCanvas('{0}/{1}_{2}_{3}_errorSummary_beforeRearrange.png'.format(plotDir,playlist,fluxType,flavor)):
        plotter.DrawErrorSummary(fluxHist_originalBinning,"TL",False,True,0.00001,False,"Flux")#,True,"",True)
      
      with makeEnv_TCanvas('{0}/{1}_{2}_{3}_CVHist_afterRearrange.png'.format(plotDir,playlist,fluxType,flavor)):
        cvHist_rebinned = fluxHist_rebinned.GetCVHistoWithError()
        cvHist_rebinned.GetYaxis().SetRangeUser(-0.5e-5,1.6e-4)
        cvHist_rebinned.Draw()
      
      with makeEnv_TCanvas('{0}/{1}_{2}_{3}_errorSummary_afterRearrange.png'.format(plotDir,playlist,fluxType,flavor)):
        plotter.DrawErrorSummary(fluxHist_rebinned,"TL",False,True,0.00001,False,"Flux")#,True,"",True)
  
      #############################################################################################################
      ### Close output file; other business #######################################################################
      #############################################################################################################
      histOutputFile.Close()
      histOutputFile_constrained.Close()

#############################################################################################################
### Produce official flux files (root and txt) ##############################################################
#############################################################################################################

if 0:

  officialFluxDir = "/minerva/data/users/{0}/FluxRearranger_officialStyleFluxPlots".format(USER)
  # create the output dir if it doesn't exist
  if not os.path.isdir(officialFluxDir):
    print "Making local official-style flux plot directory {0}".format(officialFluxDir)
    os.system( "mkdir {0}".format(officialFluxDir) )

  officialFluxRootFilePath = "{0}/officialMinervaFlux_rebinned_nuEConstrainted_2019-10-08.root".format(officialFluxDir)
  officialFluxTextFilePath = "{0}/officialMinervaFlux_rebinned_nuEConstrainted_2019-10-08.dat".format(officialFluxDir)

  fluxReweighter_withNuEConstraint = PlotUtils.FluxReweighter(14,True,PlotUtils.FluxReweighter.minervame1D1M1NWeightedAve,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv6,False)
  fluxHist_rebinned_withNuEConstraint = fluxReweighter_withNuEConstraint.GetRebinnedFluxReweighted(14,referenceHist)
  
  officialFlux = fluxHist_rebinned_withNuEConstraint.GetCVHistoWithStatError()
  officialFlux.SetName("reweightedflux_rebinned_nueconstrained_CV_WithStatErr")
  
  officialFlux.GetXaxis().SetTitle("Energy (GeV)")
  officialFlux.GetYaxis().SetTitle("#nu_{mus}/m^{2}/P.O.T./GeV")
  
  officialFluxRootFile = ROOT.TFile(officialFluxRootFilePath,"recreate")
  officialFluxRootFile.cd()
  
  writeHist(officialFlux,officialFluxRootFile)
  officialFluxRootFile.Close()
  
  officialFluxTextFile = open(officialFluxTextFilePath,"w+")
  officialFluxTextFile.write("Bin_low[GeV],Bin_high[GeV],FluxVal[nu_mu/m^2/POT/GeV]\n")
  
  for binNo,binLowEdge in enumerate(bins_nuE):
    if binLowEdge == 100.: continue
    binVal = officialFlux.GetBinContent(binNo+1)
    binHighEdge = bins_nuE[binNo+1]
    officialFluxTextFile.write("{0},{1},{2}\n".format(binLowEdge,binHighEdge,binVal))
  
  officialFluxTextFile.close()

