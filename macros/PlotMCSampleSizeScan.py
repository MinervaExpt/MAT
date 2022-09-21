import os,sys
import ROOT
from ROOT import PlotUtils

ROOT.TH1.AddDirectory(False)


def sortInputFiles(filelist):
    newlist = []
    factors = []
    for f in filelist:
        fn = f.rstrip("\n")
        mcfactor = float(fn.split(".txt")[-2].split("_")[-1])
        factors.append(mcfactor)
        

    factors = sorted(factors)
    for f in factors:
        for fi in filelist:
            mcfactor = float(fi.split(".txt")[-2].split("_")[-1])
            if(f==mcfactor): newlist.append(fi)

    return newlist

colors = [1,20,30,40,50,60,70,80,90,100,51,1,20,30,40,50,60,70,80,90,100,51]#22 allowed inputs
styles = [1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2]#22 allowed inputs


inputfiles = sys.argv[1:]
newlist = sortInputFiles(inputfiles)
for el in newlist:
    print el.rstrip("\n")
legend = ROOT.TLegend(0.2,0.7,0.82,0.9)
legend.SetFillStyle(0)
legend.SetNColumns(3)

ndf=6;
var=newlist[0].split("_")[1]
prefix=""
if var=="MnvGenie":
    var==newlist[0].split("_")[2]
    prefix=newlist[0].split("_")[0]
if var=="Enu":
    ndf=8;



chi2_by_iter_by_file = []
for f in newlist:
    fn = f.rstrip("\n")
    temp_file = open(fn,"r").readlines()
    n_uni = 0
    current_iter = 0
    chi2_by_iter = []
    for l in temp_file:
        temp_line = l.split()
        chi2 = float(temp_line[0])
        iteration = int(temp_line[1])
        universe = int(temp_line[2])
        if(universe>n_uni):
            n_uni=universe+1
        if(iteration!=current_iter):
            if(current_iter!=0):
                chi2_by_iter[-1][1]/=n_uni
            chi2_by_iter.append([iteration,chi2])
            current_iter=iteration
        else:
            chi2_by_iter[-1][1]+=chi2
        if l == temp_file[-1]:
            chi2_by_iter[-1][1]/=n_uni
#            chi2_by_iter.append([iteration,chi2])

    print chi2_by_iter
    chi2_by_iter_by_file.append(chi2_by_iter)

mygraphs = []
for i in range(0,len(chi2_by_iter_by_file)):
    tmpgraph = ROOT.TGraph()
    for j,el in enumerate(chi2_by_iter_by_file[i]):
        tmpgraph.SetPoint(j,el[0],el[1])
        
    tmpgraph.SetLineColor(colors[i])
    tmpgraph.SetLineStyle(styles[i])
    tmpgraph.SetLineWidth(2)
    tmpgraph.SetMarkerColor(colors[i])

    mygraphs.append(tmpgraph)
        
c1 = ROOT.TCanvas("c1","c1", 1200,800)
c1.cd()
c1.SetLogx()
F=1e10
for i,g in enumerate(mygraphs):
    F=float(newlist[i].split("_")[-2])
    legend.AddEntry(g,"%s"%(float(newlist[i].split(".txt")[-2].split("_")[-1])),"l")
    if(i==0):
        g.SetTitle(var +", uncfactor=" + str(F))
        g.GetYaxis().SetTitle("#chi^{2} (ndf="+str(ndf)+")")
        g.GetXaxis().SetTitle("(Unfolded Data: True Data) # of Iterations")
        g.GetXaxis().SetRangeUser(1,50)
        g.GetYaxis().SetRangeUser(0,15)
        g.Draw("ALP")
        g.SetMinimum(0)
    else:
        g.Draw("SAMELP")

legend.Draw("SAME")
target=newlist[0].split("MC_")[1].split("Nu")[0]
c1.Print(prefix+"statunc_" +target+var+ "_F"+str(F)+".png")
raw_input("DONE")
