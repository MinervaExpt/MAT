import os
import sys
import ROOT
import glob
from optparse import OptionParser
import logging

# Configure option parser
parser = OptionParser()
parser.add_option("-m", "--mergedFile", dest="mergedFileName",
                  help="full path to merged File")
parser.add_option("-u", "--unmergedFile", dest="unmergedFileName",
                  help="full path to top-level unmerged file directory")
parser.add_option("-p", "--playlist", dest="inputPlaylist",
                  help="full path to input playlist")
parser.add_option("-t", "--tool", dest="toolName",
                  help="name of ana tool")
parser.add_option("-g", "--grid", action="store_true",dest="onGrid",default=False,
                  help="are you auditing on the grid? Default is not")
(options, args) = parser.parse_args()

# Store options to local variables
onGrid = options.onGrid
unmergeddirectory = options.unmergedFileName
if onGrid: 
  mergedfile = ROOT.TFile.Open(options.mergedFileName)
else:
  mergedfile = ROOT.TFile(options.mergedFileName)
base = os.path.basename(options.mergedFileName)
myrecotree = options.toolName
playlistFilePath = options.inputPlaylist
if playlistFilePath: playlistFile = open(playlistFilePath)

# Assume the user is running over MC. If it's actually data, we'll determine that by the absence of a truth tree and adjust this variable accordingly
isMC = True 

# Check local directory to see if this audit has already been done
localfiles = os.popen("find ./ -type f").readlines()
found = False

for e in localfiles:
    if(e.find(base)!=-1):
        found = True

if(found):
    print "This file (%s) has already been checked"%(base)
    sys.exit(0)

# Configure logfile writer to write to stdout and log file
logger = logging.getLogger()
logger.setLevel(logging.INFO)

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)

if onGrid:
  fh = logging.FileHandler('%s.txt'%(base))
  fh.setLevel(logging.INFO)
  logger.addHandler(fh)

# If an input playlist is specified, use that 
if playlistFilePath:
  unmergedfiles = []
  # Populate subPlaylist list
  for fileName in playlistFile.readlines():
    unmergedfiles.append(fileName.rstrip('\n'))
# Otherwise, use the '-u' argument
else:
  unmergedfiles = os.popen("find %s -type f "%unmergeddirectory).readlines()

#Check POT, n events in truth and reco in merged
merg_POT_Used = 0
merg_POT_Total = 0
merg_reco_events = 0
merg_truth_events = 0

merg_meta = mergedfile.Get("Meta")
if not merg_meta:
    logger.info("Your merged file doesn't have a meta tree!")
    sys.exit(1)

if not merg_meta.GetListOfBranches().FindObject("POT_Used"):
    logger.info("This set of files has no POT_Used branches.  Cannot test POT. Does your anatuples have Meta?")

    metachain = ROOT.TChain("Meta")
    for f in unmergedfiles:
        filename = f.rstrip("\n")
        metachain.Add(filename)
    if not metachain or not metachain.GetListOfBranches().FindObject("POT_Used"): 
        logger.info("Looks like there's a problem with your anatuples Meta trees")
        logger.info("%s CHECK ANATUPLES"%(base))
    else:
        logger.info("Your anatuples have Meta trees") 
        logger.info("%s BAD"%(base))
    sys.exit(0)
    
if not merg_meta.GetListOfBranches().FindObject("POT_Total"):
    logger.info("This set of files has no POT_Total branches.  Cannot test POT")
    metachain = ROOT.TChain("Meta")
    for f in unmergedfiles:
        filename = f.rstrip("\n")
        metachain.Add(filename)
    if not metachain or not metachain.GetListOfBranches().FindObject("POT_Used"): 
        logger.info("Looks like there's a problem with your anatuples Meta trees")
        logger.info("%s CHECK ANATUPLES"%(base))
    else:
        logger.info("Your anatuples have Meta trees") 
        logger.info("%s BAD"%(base))
    sys.exit(0)

for e in merg_meta:
    merg_POT_Used+=e.POT_Used
    merg_POT_Total+=e.POT_Total

merg_reco = mergedfile.Get(myrecotree)
if not merg_reco:
    logger.info("Your merged file doesn't have a %s tree!!"%(myrecotree))
    sys.exit(1)

merg_reco_events = merg_reco.GetEntries()
merg_truth = mergedfile.Get("Truth")

if not merg_truth:
    logger.info("Your merged file doesn't have a truth tree. Assuming this is data")
    isMC = False
else:
    merg_truth_events=merg_truth.GetEntries()

logger.info("I have checked your merged file and found:")
logger.info("POT_Used,POT_Total = {:.3e},{:.3e}".format(merg_POT_Used,merg_POT_Total))
logger.info("Number of %s tree entries = %s"%(myrecotree,merg_reco_events))
logger.info("Number of Truth tree entries = %s"%merg_truth_events)

#Check POT, n events in truth and reco in unmerged
POT_Used = 0
POT_Total = 0
reco_events = 0
truth_events = 0

metachain = ROOT.TChain("Meta")
recochain = ROOT.TChain(myrecotree)
truthchain = ROOT.TChain("Truth")

for f in unmergedfiles:
    filename = f.rstrip("\n")
    metachain.Add(filename)
    recochain.Add(filename)
    if isMC:
        truthchain.Add(filename)

if not metachain.GetListOfBranches().FindObject("POT_Used"):
    logger.info("This set of files has no POT_Used branches.  Cannot test POT")
    logger.info("%s CHECK ANATUPLES"%(base))
    sys.exit(0)
    
if not metachain.GetListOfBranches().FindObject("POT_Total"):
    logger.info("This set of files has no POT_Total branches.  Cannot test POT")
    logger.info("%s CHECK ANATUPLES"%(base))
    sys.exit(0)

for e in metachain:
    POT_Used+=e.POT_Used
    POT_Total+=e.POT_Total
reco_events = recochain.GetEntries()
truth_events = truthchain.GetEntries()

logger.info("-----------------------------------------")
logger.info("I have checked your unmerged files and found:")
logger.info("POT_Used,POT_Total = {:.3e},{:.3e}".format(POT_Used,POT_Total))
logger.info("Number of %s tree entries = %s"%(myrecotree,reco_events))
logger.info("Number of Truth tree entries = %s"%truth_events)

good = True
good_tracker = False
frac_used_pot = 1-merg_POT_Used/POT_Used if POT_Used > 0 else 1
frac_total_pot = 1-merg_POT_Total/POT_Total if POT_Total > 0 else 1

if POT_Used > 0:
  frac_used_pot = 1-merg_POT_Used/POT_Used
else:
  if merg_POT_Used > 0:
    frac_used_pot = 1
  else:
    frac_used_pot = 0

if POT_Total > 0:
  frac_total_pot = 1-merg_POT_Total/POT_Total if POT_Total > 0 else 1
else:
  if merg_POT_Total > 0:
    frac_total_pot = 1
  else:
    frac_total_pot = 0

if(frac_used_pot < -0.0001 or frac_used_pot>0.0001): good = False
if(frac_total_pot < -0.0001 or frac_total_pot>0.0001): good = False
if(merg_reco_events!=reco_events): good = False
if(merg_truth_events!=truth_events):
    #It is possible that this merged is done by --tracker_only option, hence number of truth tree entries is different.
    tracker_truth_events = truthchain.Draw("1","mc_vtx[2]>5891 && mc_vtx[2]<8439")
    logger.info("Number of Truth tree entries in tracker = %s"%tracker_truth_events)
    if (merg_truth_events==tracker_truth_events): good_tracker = good 
    good = False

if good:
    logger.info("This set of files merged successfully")
    logger.info("%s GOOD"%(base))
    sys.exit(0)

elif good_tracker:
    logger.info("This set of files merged successfully if you merged with --tracker_only")
    logger.info("%s GOOD if tracker only"%(base))
    sys.exit(0)

else:
    logger.info("This set of files fails audit checks!!!")
    logger.info("%s BAD"%(base))
    sys.exit(1)
