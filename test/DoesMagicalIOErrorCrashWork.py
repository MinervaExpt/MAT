"""File: DoesMagicalIOErrorCrashWork.py
   Brief: Attempts to read a MINERvA AnaTuple from a .root file and
          loop over the values in a specific branch.  This program should
          crash with a nonzero return code to the operating system if
          the file can't be opened or the branch doesn't exist.

          This should be a straightforward translation of
          DoesMagicalIOErrorCrashWork.cpp.
   Author: Andrew Olivier aolivier@ur.rochester.edu"""

import sys
import ROOT
import PlotUtils

#TODO: Write just a function that registers python warnings/exceptions instead
#forceErrorHandling = ROOT.ROOT.detail.beforeMain()
ROOT.PlotUtils.HandleErrorsInPython()

#Read arguments from the command line.
#No pipes allowed because I'm not going to try to read from STDIN
#in python.
treeName = sys.argv[1]
branchName = sys.argv[2]
filesToRead = sys.argv[3:]

#Process each file.  Any failed file opens or branch lookups should
#stop the script with a nonzero exit code.
for fileName in filesToRead:
  inFile = ROOT.TFile.Open(fileName)

  tree = inFile.Get(treeName)

  #Strangely, ROOT doesn't warn users about this case.  So, the error
  #handling machinery can't help.  I'll put in a warning myself.
  if tree is None:
    print "Failed to read a TTree named " + tree.GetName() + " from " + inFile.GetName() + ".  The error handling code in PlotUtils isn't able to catch this failure mode, so you'll have to handle it yourself."
    sys.exit(1)

  tree.SetBranchStatus("*", 0)
  tree.SetBranchStatus(branchName, 1)
 
  whichEntry = 0 
  for entry in tree:
    if(whichEntry % 1000 == 0):
      print "Processing entry " + str(whichEntry) + "..."
    whichEntry += 1

#If everything worked, make certain that this script returns 0 to the
#operating system.
sys.exit(0)
