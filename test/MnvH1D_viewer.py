import sys,os,string
import math

# plot an MNvH1D CV and the ratios of its error bands to their CV's.

from ROOT import *

from PlotUtils import *

from array import array
full = False
binwidth = True
norm = True
if len(sys.argv)< 3:
  print ("viewer args are filename histname")
  sys.exit(1)
file = sys.argv[1]
hist = sys.argv[2]
localdir = os.path.dirname(file)
dir= os.path.join(localdir,os.path.basename(file).replace(".root","_csvdump"))
if not os.path.exists(dir):
  os.makedirs(dir)


f = TFile.Open(file,"READONLY")
#f.ls()

h = MnvH1D()
h = f.Get(hist)
# hack to fix escape characters
h.GetXaxis().SetTitle(h.GetXaxis().GetTitle().replace("\nu","\\nu"))
h.GetYaxis().SetTitle(h.GetXaxis().GetTitle().replace("\nu","\\nu"))
h.SetTitle(h.GetTitle().replace("\nu","\\nu"))

fname = file[0:-5]
o = TFile.Open(fname+"_"+hist+"_cv_bands_TH.root","RECREATE")

o.cd()
cv = h.GetCVHistoWithStatError()
xtra = "_binwidth"

if not binwidth:
 xtra = ""

if not full:
  xtra += "_short"

h.MnvH1DToCSV(h.GetName()+xtra,dir,1.E39,full,True,False,binwidth)
cv.Write()


n = h.GetNVertErrorBands()

names = h.GetErrorBandNames()

for name in names:
  #print ("error band:", name)
  band = h.GetVertErrorBand(name)
  hists = band.GetHists()
  bcv = MnvH1D()
  bcv = band.GetErrorBand()
  bcv.SetName(bcv.GetName()+ " error band CV")
  if norm:
    bcv.Divide(bcv,cv)
  bcv.Write()
  
  for i in range (0,band.GetNHists()):
    if(norm):
      hists[i].Divide(cv)
    
    hists[i].Write()

  

names = h.GetLatErrorBandNames()

for name in names:
  band = h.GetLatErrorBand(name)
  hists = band.GetHists()
  for i in range (0,band.GetNHists()):
    if norm:
      hists[i].Divide(cv)
    
    hists[i].Write()
o.Close()




