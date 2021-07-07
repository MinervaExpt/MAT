"""
  NIMStyle.py:
   Loads the NIM paper plot style from PlotUtils.
  
   Original author: J. Wolcott (jwolcott@fnal.gov)
                    March 2013
"""

import os, os.path
import sys

import ROOT

# the ROOT interpreter gobbles up any arguments
# that might be in sys.argv, whether you like it or not.
# we need to tuck them away for later.
args = sys.argv[:]
sys.argv = sys.argv[:0]

def SetStyle(path_to_plotutils):
	ROOT.gROOT.SetMacroPath( "%s:%s" % (ROOT.gROOT.GetMacroPath(), path_to_plotutils) )
	ROOT.gROOT.Macro("NIMStyle.C")


if "PLOTUTILSROOT" in os.environ:
	SetStyle("%s/PlotUtils" % os.path.expandvars("$PLOTUTILSROOT"))
else:
	print >> sys.stderr, "WARNING: $PLOTUTILSROOT is not defined in the current environment.  NIMStyle cannot be loaded automatically."

# now that ROOT has done its thing,
# we can restore the arguments...
sys.argv = args
