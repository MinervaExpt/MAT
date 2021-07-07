import sys,os,string
import math

from ROOT import *

from PlotUtils import *

def doit(keys,file,output,dim,canvas,pdfname):
 
  #print keys
  for key in keys:
    object = key.ReadObj()
    #print "the object", object.IsA(),type(object.IsA())
    #print object.GetName(),object.IsA().GetName()
    if "TDirectory" in object.IsA().GetName():
      print (" a directory ", object.GetName())
      newkeys = object.GetListOfKeys()
      #for newkey in newkeys:
      #  print "in the directory there is ", newkey.ReadObj().GetName()
      doit(newkeys,file,output,dim,canvas,pdfname)
    
    if dim == "1D" and ("MnvH1D" in object.IsA().GetName() or "TH1D" in object.IsA().GetName()):
      TH1 = TH1D(object)
      #TH1.Scale(1.,"width")
      #TH1.Print()
      TH1.Write()
      TH1.SetTitle(TH1.GetName())
      TH1.SetMinimum(0.8)
      TH1.SetMaximum(1.2)
      TH1.Draw("hist")
      canvas.Print(pdfname,"Title:"+TH1.GetName())
      continue
    if dim == "2D" and ("MnvH2D" in object.IsA().GetName() or "TH2D" in object.IsA().GetName()):
        TH2temp = TH2D(object)

        if dim == "px":
          TH2X = TH2temp.ProjectionX()
        else:
          TH2X = TH2temp.ProjectionY()

        #TH2X.Scale(1.,"width")

        TH2X.Write()
        TH2X.Draw()
        canvas.Print(pdfname,"Title:"+TH2X.GetName())
        TH2temp.SetStats(0)
        canvas.SetLogy()
        TH2temp.SetTitle(TH2temp.GetName())
        TH2temp.SetMinimum(0.8)
        TH2temp.SetMaximum(1.2)
        TH2temp.Draw("COLZ")
        
        canvas.Print(pdfname,"Title:"+TH2temp.GetName())
        continue
    
    else:
      continue

#-------------------------

filename = sys.argv[1]

thefile = TFile.Open(filename,"READONLY")

thedir = thefile.GetListOfKeys()
output = TFile.Open(filename.replace(".root","_sum.root"),"RECREATE")

#print thedir, thefile,output

pdfname = filename.replace(".root","_1D.pdf")
pdfstart = pdfname.replace(".pdf",".pdf(")
pdfend = pdfname.replace(".pdf",".pdf)")
canvas = TCanvas("1D")
canvas.SetLeftMargin(0.2);
canvas.SetBottomMargin(0.2);
canvas.Print(pdfstart,"pdf")
doit(thedir,thefile,output,"1D",canvas,pdfname)
canvas.Print(pdfend,"pdf")

pdfname = filename.replace(".root","_2D.pdf")
pdfstart = pdfname.replace(".pdf",".pdf(")
pdfend = pdfname.replace(".pdf",".pdf)")
canvas = TCanvas("2D")
canvas.SetLeftMargin(0.2);
canvas.SetBottomMargin(0.2);
canvas.Print(pdfstart,"pdf")
doit(thedir,thefile,output,"2D",canvas,pdfname)
canvas.Print(pdfend,"pdf")





