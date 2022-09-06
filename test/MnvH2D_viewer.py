import sys,os,string
import math

from ROOT import *

from PlotUtils import *

from array import array

def SyncBands2(hist):
 #print (hist.GetName())
 
 theCVhisto = MnvH2D()
 theCVHisto = hist.Clone()
 theCVHisto.SetDirectory(0);
 bandnames = hist.GetErrorBandNames();
 # "Synching Error Band CV's with MnvH1D's CV" << std::endl;
 # loop this MnvH1D's bands
 for bandname in bandnames:
   #std::cout << *bandname << std::endl;
   # band is a reference to the error band's CV histo
   band = MnvVertErrorBand()
   band = hist.GetVertErrorBand(bandname)
   # set the BAND's CV histo to the MnvH1D's CV histo
   for i in range(0,band.GetNbinsX()+2):
      for j in range(0,band.GetNbinsY()+2):
     #if(i < 4):
     #  print "band ", hist.GetName(), bandname, i,theCVHisto.GetBinContent(i), band.GetBinContent(i)
        band.SetBinContent(i,j,theCVHisto.GetBinContent(i,j))
        band.SetBinError(i,j,theCVHisto.GetBinError(i,j))
 return theCVHisto


binwidth = True  # correct for bindwidth

full = False # full precision
Fractional = True # show individual error bands as fractions
xtra = ""

if binwidth:
  xtra = "_binwidth"
if not full:
  xtra += "_short"
def isTrivial(c,h):
  eps = 0.01
  trivial = True
  nx = h.GetNbinsX()
  ny = h.GetNbinsY()
  
  # ignore the last bin
  for i in range(1,nx):
    for j in range(1,ny):
      z0 = c.GetBinContent(i,j)
      z = h.GetBinContent(i,j)
      if z0 == 0:
        continue
      #print h.GetName(), i,j,z,z0,z/z0-1, eps
      if (z/z0 -1. > eps or z/z0 -1  < -eps):
      #  print h.GetName(), " is ok"
        trivial = False
        return trivial
  #print h.GetName(), " is trivial"
  return True


if len(sys.argv)< 3:
  print ("viewer args are filename histname")
  sys.exit(1)
file = sys.argv[1]
localdir = os.path.dirname(file)
dir= os.path.join(localdir,os.path.basename(file).replace(".root","_csvdump"))
if not os.path.exists(dir):
  os.makedirs(dir)
hist = sys.argv[2]
histname = hist
newname = sys.argv[2]
oldname = newname
if len(sys.argv)>3:
  oldname = sys.argv[3]
if len(sys.argv)>4:
  newname = sys.argv[4]
#print ("----",hist,"-----")
xtra = ""
if binwidth:
  xtra = "_binwidth"
if not full:
  xtra += "_short"
f = TFile.Open(file,"READONLY")

h = f.Get(hist)
if h == None:
  print ("no ",sys.argv[2]," in ", sys.argv[1])
  sys.exit(1)

# hack to fix escape characters
h.GetXaxis().SetTitle(h.GetXaxis().GetTitle().replace("\nu","\\nu"))
h.GetYaxis().SetTitle(h.GetYaxis().GetTitle().replace("\nu","\\nu"))
h.SetTitle(h.GetTitle().replace("\nu","\\nu"))
print ("name is ",h.GetName(),h.GetTitle())


central = h.GetCVHistoWithStatError()
central.SetDirectory(0)
cx = central.ProjectionX()
cy = central.ProjectionY()

h.MnvH2DToCSV(h.GetName()+xtra,dir,1.E44,full,True,Fractional,binwidth)

b0 = h.Clone()
b0.SetDirectory(0)
#b = SyncBands2(b0)
#b0.Print("ALL")
bx = b0.ProjectionX()
bx.GetXaxis().SetTitle(h.GetXaxis().GetTitle())
by = b0.ProjectionY()
by.GetXaxis().SetTitle(h.GetYaxis().GetTitle())
bx.MnvH1DToCSV(bx.GetName()+xtra,dir,1.e44,full,True,Fractional,binwidth)
by.MnvH1DToCSV(by.GetName()+xtra,dir,1.e44,full,True,Fractional,binwidth)
fname = file[0:-5]

if not binwidth:
  o = TFile.Open(fname+"_"+hist+"_cv_bands_TH.root","RECREATE")
if binwidth:
  o = TFile.Open(fname+"_"+hist+"_cv_bands_binwidth_TH.root","RECREATE")

o.cd()


cv = TH2D(h.GetCVHistoWithStatError())
cv.SetDirectory(0)
if binwidth:
  cv.Scale(1.0,"width")
  cx.Scale(1.0,"width")
  cy.Scale(1.0,"width")
cv.Write()
cx.SetName(cx.GetName().replace(oldname,newname))
cx.Write()
cy.SetName(cy.GetName().replace(oldname,newname))
cy.Write()

#mx = TMatrixD()
#mx = b0.GetTotalCorrelationMatrix()
##mx = TMatrixD(cx.GetTotalCorrelationMatrix)
#mx.ResizeTo(30,30)
##mx.Print()




names = h.GetVertErrorBandNames()
#print "Vert Error band names ", names
for name in names:
  band = h.GetVertErrorBand(name)
  vcv = MnvH2D(band)
  
  cx = band.ProjectionX()
  cy = band.ProjectionY()
  
  hists = band.GetHists()
  for i in range (0,min(10,band.GetNHists())):
    if isTrivial(vcv, hists[i]):
      continue
      #hists[i].Write()
    hists[i].SetName(hists[i].GetName().replace(oldname,newname))
    hx = (hists[i].ProjectionX())
    hx.Divide(cx)
    nx = hx.GetName()
    hx.SetName(nx.replace(oldname,newname))
    hx.Write()
    hy = (hists[i].ProjectionY())
    hy.Divide(cy)
    ny = hy.GetName()
    hy.SetName(ny.replace(histname,newname))
    hy.Write()
    hists[i].Divide(central)
    hists[i].Write()



names = h.GetLatErrorBandNames()
for name in names:
  band = h.GetLatErrorBand(name)
  vcv = MnvH2D(band)
  hists[i].SetName(hist[i].GetName().replace(oldname,newname))
  cx = band.ProjectionX()
  cy = band.ProjectionY()
  hists = band.GetHists()
  for i in range (0,min(10,band.GetNHists())):
    if isTrivial(vcv, hists[i]):
      continue

    hx = (hists[i].ProjectionX())
    hx.Divide(cx)
    hx.Write()
    hy = (hists[i].ProjectionY())
    hy.Divide(cy)
    hy.Write()
    hists[i].Divide(central)
    hists[i].Write()
o.Close()




