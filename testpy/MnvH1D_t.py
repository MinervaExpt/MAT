import ROOT, array, sys
import PlotUtils
from PlotUtils import *
from math import sqrt   

h_test = MnvH1D("test","title",10,0.0,10.0)
h_test.Sumw2()

for i in range(0,h_test.GetNbinsX()+2):
    h_test.SetBinContent(i,i*2+1.)
    h_test.SetBinError(i,sqrt(h_test.GetBinContent(i)))

h_test.Print("ALL")
# Draw the correlation matrix for kicks
c = ROOT.TCanvas()

