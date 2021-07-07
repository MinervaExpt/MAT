## import common python modules
import ROOT,os,sys,array
from argparse import ArgumentParser
## modify $PYTHONPATH to know about custom python modules
PLOTUTILSROOT = os.path.expandvars("$PLOTUTILSROOT")
sys.path.append("{0}/NSF_ValidationSuite/py_classes".format(PLOTUTILSROOT))
## import custom python modules
from PlotUtils.HistWrapper import HistWrapper
from validationConfigCard import *
from validationFunctions import *

################################################################ Preamble above

def makeRecoHistsForMCTuple( OPTS_VEC , tuplePath , outputFile , toolName ):

  print "I'm inside makeRecoHistsForMCTuple"

  # Declare chain
  chain = PlotUtils.ChainWrapper( toolName )
  chain.Add( os.path.expandvars(tuplePath) )

  # Get list of systematic universes to loop through 
  systematicUniverses = defineSystematicUniverses(chain,"mc") 

  # Declare hists 
  
  hist_inclusiveEmu = HistWrapper('inclusive_Emu',nBins_Emu,array.array('d',bins_Emu),systematicUniverses)
  hist_inclusivePt= HistWrapper('inclusive_Pt',nBins_pt,array.array('d',bins_pt),systematicUniverses)
  hist_inclusiveNu = HistWrapper('inclusive_Nu',nBins_nu,array.array('d',bins_nu),systematicUniverses)
  
  hist_inclusive_WgtVsPmuMinos_MINOS = HistWrapper('inclusive_WgtVsPmuMinos_MINOS',40,0,20,100,0.9,1.1,systematicUniverses)
  hist_inclusive_WgtVsEmu_MINOS = HistWrapper('inclusive_WgtVsEmu_MINOS',40,0,20,100,0.9,1.1,systematicUniverses)
  hist_inclusive_WgtVsEmu_GENIE = HistWrapper('inclusive_WgtVsEmu_GENIE',40,0,20,200,0.0,2.0,systematicUniverses)
  hist_inclusive_WgtVsEmu_RPA = HistWrapper('inclusive_WgtVsEmu_RPA',40,0,20,200,0.0,2.0,systematicUniverses)
  hist_inclusive_WgtVsEmu_2p2h = HistWrapper('inclusive_WgtVsEmu_2p2h',40,0,20,200,0.0,2.0,systematicUniverses)
  hist_inclusive_WgtVsEmu_Flux = HistWrapper('inclusive_WgtVsEmu_Flux',40,0,20,200,0.0,2.0,systematicUniverses)
 
  # Loop over entries in the chain
  for entry in range(chain.GetEntries()):
  
    if entry % 1000 == 0:
      print "entry #: " , entry

    # Implement pre-cuts
    if not passCuts(chain,toolName,entry): continue 
    
    # Loop over systematics
    for systematicUniverseClass in systematicUniverses:
      for i,systematicUniverse in enumerate(systematicUniverses[systematicUniverseClass]):

        systematicUniverse.SetEntry(entry)

        TMu = systematicUniverse.GetRecoMuonEKinetic()/1000.
        EMu = systematicUniverse.GetEmu()/1000.
        nu = systematicUniverse.GetRecoilEnergyLocal(toolName)
        PmuT = systematicUniverse.GetPmuTransverse()/1000.
       
        wgt = systematicUniverse.getWeight()

        # Wgt Plots
        PmuMinos = systematicUniverse.GetPmuMinos()       
        wgt_MINOS = systematicUniverse.GetMinosEfficiencyWeight()
        wgt_GENIE = systematicUniverse.GetGenieWeight()
        wgt_RPA = systematicUniverse.GetRPAWeight()
        wgt_2p2h = systematicUniverse.GetLowRecoil2p2hWeight()
        wgt_Flux = systematicUniverse.GetFluxAndCVWeight()

        # Individual hists for particular studies
        hist_inclusiveEmu.FillUniverse(systematicUniverse,EMu,wgt)
        hist_inclusivePt.FillUniverse(systematicUniverse,PmuT,wgt)
        hist_inclusiveNu.FillUniverse(systematicUniverse,nu,wgt)

        hist_inclusive_WgtVsPmuMinos_MINOS.FillUniverse(systematicUniverse,PmuMinos/1000,wgt_MINOS,1.0)
        hist_inclusive_WgtVsEmu_MINOS.FillUniverse(systematicUniverse,EMu,wgt_MINOS,1.0)
        hist_inclusive_WgtVsEmu_GENIE.FillUniverse(systematicUniverse,EMu,wgt_GENIE,1.0) 
        hist_inclusive_WgtVsEmu_RPA.FillUniverse(systematicUniverse,EMu,wgt_RPA,1.0) 
        hist_inclusive_WgtVsEmu_2p2h.FillUniverse(systematicUniverse,EMu,wgt_2p2h,1.0)
        hist_inclusive_WgtVsEmu_Flux.FillUniverse(systematicUniverse,EMu,wgt_Flux,1.0)

  print 'tuplePath: ', tuplePath
  
  outputFile.cd()

  # Sync CV
  hist_inclusiveEmu.SyncCVHistos()
  hist_inclusivePt.SyncCVHistos()
  hist_inclusiveNu.SyncCVHistos()
  hist_inclusive_WgtVsPmuMinos_MINOS.SyncCVHistos()
  hist_inclusive_WgtVsEmu_MINOS.SyncCVHistos()
  hist_inclusive_WgtVsEmu_GENIE.SyncCVHistos()
  hist_inclusive_WgtVsEmu_RPA.SyncCVHistos()
  hist_inclusive_WgtVsEmu_2p2h.SyncCVHistos()
  hist_inclusive_WgtVsEmu_Flux.SyncCVHistos()
 
  # Write out MnvHXDs 
  hist_inclusiveEmu.hist.Write()
  hist_inclusivePt.hist.Write()
  hist_inclusiveNu.hist.Write()
  hist_inclusive_WgtVsPmuMinos_MINOS.hist.Write()
  hist_inclusive_WgtVsEmu_MINOS.hist.Write()
  hist_inclusive_WgtVsEmu_GENIE.hist.Write()
  hist_inclusive_WgtVsEmu_RPA.hist.Write()
  hist_inclusive_WgtVsEmu_2p2h.hist.Write()
  hist_inclusive_WgtVsEmu_Flux.hist.Write()

def main():
  
  print "I'm inside main"

  ## Parse user args
  parser = ArgumentParser(description='Process optional inputs')
  parser.add_argument('--outdir', dest='outDir', action='store', required=True)
  OPTS_VEC = parser.parse_args()

  toolName = validationSampleName 

  # Loop over tuples (defined in validationConfigCard.py)
  for tuplePath in tuplePaths:
  
    print "I'm processing this tuple: {0}".format(tuplePath)

    outPath = constructOutputFilePath(tuplePath,OPTS_VEC.outDir,toolName)
    outputFile = ROOT.TFile( outPath , "recreate" )
    print "I'm going to write histograms here: {0}".format(outPath)

    # These variables are all set in validationConfigCard.py
    PlotUtils.DefaultCVUniverse.SetNuEConstraint(useNuEScatteringConstraint)
    PlotUtils.DefaultCVUniverse.SetNonResPiReweight(useNonResPiReweight)
    PlotUtils.DefaultCVUniverse.SetPlaylist(playlist)
    PlotUtils.DefaultCVUniverse.SetAnalysisNuPDG(analysisNuPDG)
    PlotUtils.DefaultCVUniverse.SetNFluxUniverses(nFluxUniverses)

    makeRecoHistsForMCTuple(OPTS_VEC,tuplePath,outputFile,toolName)
    copyMetaTreeToOutput(OPTS_VEC,tuplePath,outputFile) 

    outputFile.Close()

if __name__ == "__main__":
  main()

