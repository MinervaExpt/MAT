import sys,os,string
import math

# plot an MNvH1D CV and the ratios of its error bands to their CV's.

from ROOT import *

from PlotUtils import *

from array import array

norm = True
if len(sys.argv)< 3:
  print ("viewer args are filename histname")
  sys.exit(1)
file = sys.argv[1]
hist = sys.argv[2]
 
dir = os.path.dirname(file)

f = TFile.Open(file,"READONLY")
f.ls()

h = MnvH1D()
h = f.Get(hist)
h.Scale(1.E42)


fname = file[0:-5]
o = TFile.Open(fname+"_"+hist+"_cv_bands_TH.root","RECREATE")

o.cd()
cv = h.GetCVHistoWithStatError()
xtra = "_binwidth"
binwidth = True
if not binwidth:
 xtra = ""
full = False

h.MnvH1DToCSV(h.GetName()+xtra,dir,1.,full,True,False,binwidth)
cv.Write()

#o = TFile.Open(fname+"_"+hist+"_cv_bands_TH.root","RECREATE")

n = h.GetNVertErrorBands()

names = h.GetErrorBandNames()

for name in names:
  print ("error band:", name)
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




