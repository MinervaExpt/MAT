import sys,os,string
import math

from ROOT import *

from PlotUtils import *

from array import array

trials = 100000

def func(x):
  return 1.0



def Test1D():
  TH1.AddDirectory(False)
  
  yscale = array('f', [0.98,1.02])
  xoffset = array('f', [-0.5,0.5])
  xscale = array('f', [0.95,1.05])
  xbins = array('d',[0.0,0.5,1.0,2.0,3.0,4.0])
  print (xbins,len(xbins))
  h_1 = MnvH1D("h_1","Central Value;x",len(xbins)-1,xbins)
  h_1.SetStrict(1)
  f = TFile.Open("test_output.root",'RECREATE')
  r = TRandom()

  
  h_1.AddVertErrorBand("yscale",2)
  h_1.AddLatErrorBand("xoffset",2)
  h_1.AddLatErrorBand("xscale",2)


  for i in range(0,trials):
    x = r.Gaus(5,2)
    wt = func(x)
    h_1.Fill(x,wt)
    h_1.FillVertErrorBand("yscale",x,yscale[0],yscale[1],wt)
    h_1.FillLatErrorBand("xscale",x,x*xscale[0]-x,x*xscale[1]-x,wt)
    h_1.FillLatErrorBand("xoffset",x,xoffset[0],xoffset[1],wt)
  nx = h_1.GetNbinsX()+2
  cov = TMatrixD(nx,nx)
  
  for i in range(0,cov.GetNrows()):
    for j in range(0,i+1):
      cov[i][j] = r.Gaus(0,.1);
      cov[j][i] = cov[i][j]
      if( i == j):
        cov[i][j] = 1;

  h_1.FillSysErrorMatrix("covariance",cov)
  
  h_4 = h_1.Clone("double")
  h_4.Add(h_1,2.)
  
  m2 = TMatrixD
  m2 = h_1.GetSysErrorMatrix("covariance")
  m4 = TMatrixD
  m4 = h_4.GetSysErrorMatrix("covariance")
  
  m2.Print()
  m4.Print()

  h_1.Print("ALL")
  h_1.Scale(1,"width")
  f.cd()
  h_1.Write()
  h_1.MnvH1DToCSV(h_1.GetName(),"./",1.0,True,True)
  h_2 = h_1.Clone()
  h_2.SetStrict(1)
  print ("h_1 errors", h_1.GetSysErrorMatricesNames())
  print ("h_2 errors", h_2.GetSysErrorMatricesNames())

  for e in h_1.GetSysErrorMatricesNames():
    print ("h_1 e",e)

  for e in h_2.GetSysErrorMatricesNames():
    print ("h_2 e",e)


def Test2D():
  TH1.AddDirectory(False)
  xbins = array('d',[0.0,1.0,2.0,3.0,4.0,6.0])
  ybins = array('d',[-5.,-2.,-1.,1, 2.,5.])
  
  h_2 = MnvH2D("h_2","Central Value 2D;x;y;z",5,xbins,5,ybins)
  h_2.SetStrict(1)
  yscale = array('f', [0.98,1.02])
  xoffset = array('f', [-0.5,0.5])
  xscale = array('f', [0.95,1.05])
  
  nx = h_2.GetNbinsX()+2
  ny = h_2.GetNbinsY()+2
  
  cov = TMatrixD(nx*ny,nx*ny)
  
  for i in range(0,cov.GetNrows()):
    for j in range(0,cov.GetNrows()):
      cov[i][j] = i*j

  
  
  f = TFile.Open("test2_output.root",'RECREATE')
  r = TRandom()
  
  h_2.AddVertErrorBand("yscale",2)
  #h_2.AddLatErrorBand("xoffset",2)
  #h_2.AddLatErrorBand("xscale",2)
  
  
  for i in range(0,trials):
    x = r.Gaus(0,3)
    y = r.Gaus(0,3)
    wt = func(x)*func(y)
    
    h_2.Fill(x,y,wt)
    h_2.FillVertErrorBand("yscale",x,yscale[0],yscale[1],wt)



  names = h_2.GetSysErrorMatricesNames()
  print (names)
  h_2.PushCovMatrix("testcov",cov)

  h_4 = h_2.Clone("double")
  h_4.Add(h_2,2.)
  
  m2 = TMatrixD
  m2 = h_2.GetSysErrorMatrix("testcov")
  m4 = TMatrixD
  m4 = h_4.GetSysErrorMatrix("testcov")
  m2.ResizeTo(4,4)
  m4.ResizeTo(4,4)
  m2.Print()
  m4.Print()
  
  h_3 = h_2.Clone("h_3_scale")
  
  h_3.Scale(1.,"width")

  h_3.Print()

  h_4 = h_2.Clone("ratio")

  print ("h_3 errors", h_3.GetVertErrorBandNames())
  h_3.Print()
  h_3.SetName("scaled")
  print ("h_2 error names ", h_2.GetSysErrorMatricesNames())
  print ("h_3 error names ", h_3.GetSysErrorMatricesNames())
  print ("h_4 error names ", h_3.GetSysErrorMatricesNames())
  h_4.Divide(h_3,h_2)
  
  newcov = h_2.GetSysErrorMatrix("testcov")
#print "newcov"
#newcov.Print()
  ratcov = h_3.GetSysErrorMatrix("testcov")
  print ("ratcov")
  ratcov.Print()

#h_2.Print("ALL")
  print (" make a projection")
  h_2.Print("ALL")
  h_2.MnvH2DToCSV(h_2.GetName(),"./",1.0,True)
  h_x = h_2.ProjectionX()
  h_x.MnvH1DToCSV(h_x.GetName(),"./",1.0,True)
  h_x.Print("ALL")
  h_x2 = h_x.Clone("h_x2")
  h_x2.Scale(1.,"width")
  h_y = h_2.ProjectionY()
  h_y.Print("ALL")
  h_y2 = h_y.Clone("h_y2")
  h_y2.Scale(1.0,"width")
  h_y2.Print("ALL")
  f.cd()
  print (" try to write out h_2")
  h_2.Write()
  h_x2.Write()
  h_x.Write()
  h_y.Write()
  h_y2.Write()
  h_3.Write()
  h_4.Write()
  print ("end of Test2D()")
  return



Test1D()
Test2D()
sys.exit(0)

