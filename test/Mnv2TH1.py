import sys,os,string
import math

from ROOT import *

from PlotUtils import *

def doit(keys,file,output,dim):

  #print keys  
  for key in keys:
    object = key.ReadObj()
    #print "the object", object.IsA(),type(object.IsA())
    if "TDirectory" in object.IsA().GetName():
      print (" a directory ", object.GetName())
      newkeys = object.GetListOfKeys()
      #for newkey in newkeys:
      #  print "in the directory there is ", newkey.ReadObj().GetName()
      doit(newkeys,file,output,dim)
    
    if "MnvH1D" in object.IsA().GetName():
      TH1 = TH1D(object)
      #TH1.Scale(1.,"width")
      TH1.Print()
      TH1.Write()
    if "MnvH2D" in object.IsA().GetName():
        TH2temp = TH2D(object)
        if dim == "px":
          TH2X = TH2temp.ProjectionX()
        else:
          TH2X = TH2temp.ProjectionY()
       
        #TH2X.Scale(1.,"width")
      
        TH2X.Write()
        
    else:
      continue

#-------------------------

filename = sys.argv[1]

thefile = TFile.Open(filename,"READONLY")

thedir = thefile.GetListOfKeys()
output = TFile.Open(filename.replace(".root","_py_TH.root"),"RECREATE")
print thedir, thefile,output
doit(thedir,thefile,output,"py")
output = TFile.Open(filename.replace(".root","_px_TH.root"),"RECREATE")
print thedir, thefile,output
doit(thedir,thefile,output,"px")






