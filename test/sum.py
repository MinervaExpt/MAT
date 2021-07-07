import sys,os,string
import math

# plot an MNvH1D CV and the ratios of its error bands to their CV's.

from ROOT import *

from PlotUtils import *

from array import array

def dot(a,b):
  n = len(a)
  if len(b)!= n:
    return None
  d = 0
  for i in range(1,n):
    d[i] += a[i] * b[i]
  return d
  
def getvals(h):
  n = h.GetNbinsX()
  t = 0
  s = 0
  y = 0
  
  hs = TH1D()
  hs =h.GetStatError()
  
  hy = TH1D()
  hy = h.GetTotalError(False)
 

  for i in range(1,n+1):
    t += h.GetBinContent(i)
    s += hs.GetBinContent(i)*hs.GetBinContent(i)
    y += hy.GetBinContent(i)*hy.GetBinContent(i)
 
  return [t,math.sqrt(s),math.sqrt(y)]

  
def getvals2(h):
  n = h.GetNbinsX()
  m = h.GetNbinsY()
  t = 0
  s = 0
  y = 0
  hs = TH1D()
  hs = h.GetStatError()
  hy = TH1D()
  hy = h.GetTotalError(False)
  
  #hy.Print("ALL")
  for i in range(1,n+1):
    for j in range(1,m+1):
      t += h.GetBinContent(i,j)
      s += hs.GetBinContent(i,j)*hs.GetBinContent(i,j)
      y += hy.GetBinContent(i,j)*hy.GetBinContent(i,j)
 
  return [t,math.sqrt(s),math.sqrt(y)]
 
  
def getErrorMatrix(h):
  m = TMatrixD()
  m = h.GetTotalErrorMatrix(False)
  n = h.GetNbinsX()
  v = TVectorD(n)
  l = []
  for i in range(1,n+1):
    v[i-1] = m[i][i]
    l.append(m[i][i])
  print (l)
  #m.Print()
  
  
def pretty(n,v):
  s = "%30s : %g +- %g +- %g"%(n, v[0],v[1],v[2])
  return s
def doit(keys):

 #print keys
 for key in keys:
   object = key.ReadObj()
   
   if "TDirectory" in object.IsA().GetName():
     print (" a directory ", object.GetName())
     newkeys = object.GetListOfKeys()
     
  
     doit(newkeys,file,output,dim,canvas,pdfname)
   
   if "MnvH1D" in object.IsA().GetName():
     h = MnvH1D()
     h = object.Clone()
     if "migration" in h.GetName():
        continue
     print (pretty(h.GetName(), getvals(h)))
     continue
   if "MnvH2D" in object.IsA().GetName():
      h = MnvH2D()
      h = object.Clone()
      if "migration" in h.GetName():
        continue
      print (pretty(h.GetName(), getvals2(h)))
      continue
   
   else:
     continue
     
filename = sys.argv[1]

thefile = TFile.Open(filename,"READONLY")

thedir = thefile.GetListOfKeys()

doit(thedir)


