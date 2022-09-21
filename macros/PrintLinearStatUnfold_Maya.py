from ROOT import *
from ROOT import PlotUtils
from ROOT.PlotUtils import *
from ROOT.PlotUtils import MnvH1D
from collections import namedtuple
from OrderedDict import OrderedDict
from array import array
import sys


itcollect  =namedtuple("itcollect","iteration legtitle color")
multihist  =namedtuple("multihist","collect legend")
histcollect=namedtuple("histcollect","filename histname title color mstyle msize")

def set_palette():
  red   = array('d',[0   ,1.0   ,1.0  ])
  green = array('d',[0   ,1.0   ,0   ])
  blue  = array('d',[1.0 ,1.0   ,0   ])
  stops = array('d',[0.0 ,0.5 ,1.0])

  TColor.CreateGradientColorTable(3,stops,red,green,blue,50)  

def bookHistos( filename, histnames, titles, newnames =[]):
  if type(histnames) == list:
    ret = {}
    for i,name in enumerate(histnames):
      getfile = TFile(filename[i])
      tmphist=getfile.Get(name)
      tmphist.Print("all")
      #if tmphist is None:
      print "Hist doesn't exist: " + name+ "in file: "+filename[i]
      #if titles[i]!="":
      #  tmphist.SetTitle(titles[i])
      tmphist.SetDirectory(0)
      if len(newnames)==0:
        ret[name]=tmphist
      else:
        ret[newnames[i]]=tmphist
      getfile.Close()
    return ret
  else:
    getfile = TFile(filename)
    tmphist=getfile.Get(histnames)
    if titles!="":
      tmphist.SetTitle(titles)
    tmphist.SetDirectory(0)
    getfile.Close()
    return tmphist

def LinearizeHist(hist):
  nBins = hist.fN
  rethist=TH1D(hist.GetName()+"_Linear",hist.GetTitle(),nBins,0,nBins)
  for b in range(nBins):
    rethist.SetBinContent(b+1,hist.GetBinContent(b))
    rethist.SetBinError(b+1,hist.GetBinError(b))
  return rethist

#filenames = ["MayaUnfold_Closure_BjorkenX_IT1000_SU1000","MayaUnfold_Closure_Enu_IT1000_SU1000","MayaUnfold_BjorkenX_IT1000_SU1000","MayaUnfold_Enu_IT1000_SU1000"]
#/minerva/data/users/afilkins/NukeHists/TransWarpStudies_210117/AllNuME
filenames = ["MayaUnfold_New_Closure_BjorkenX_IT1000_SU1000","MayaUnfold_New_Closure_Enu_IT1000_SU1000","MayaUnfold_New_BjorkenX_IT1000_SU1000","MayaUnfold_New_Enu_IT1000_SU1000"]
#filenames = ["StupidBinErrorMayaUnfold_New_Closure_BjorkenX_IT1000_SU1000","StupidBinErrorMayaUnfold_New_Closure_Enu_IT1000_SU1000","StupidBinErrorMayaUnfold_New_BjorkenX_IT1000_SU1000","StupidBinErrorMayaUnfold_New_Enu_IT1000_SU1000"]
#filenames = ["StupidBinErrorMayaUnfold_New_Enu_IT1000_SU1000"]
#filenames = ["MayaUnfold_Closure_Enu_IT1000_SU1000"]
itercollect = [itcollect(1   ,"1 Iteration(s)",    kRed+2),
               itcollect(2   ,"2 Iteration(s)",    kOrange+7),
               itcollect(5   ,"5 Iteration(s)",    kYellow+2),
               itcollect(10  ,"10 Iteration(s)",   kSpring-7),
               itcollect(20  ,"20 Iteration(s)",   kGreen+2),
               itcollect(50  ,"50 Iteration(s)",   kTeal-7),
               itcollect(100 ,"100 Iteration(s)",  kCyan+3),
               itcollect(200 ,"200 Iteration(s)",  kAzure+7),
               itcollect(500 ,"500 Iteration(s)",  kBlue+1 ), 
               itcollect(1000,"1000 Iteration(s)", kViolet+3)] 

iterbins = []
for it in itercollect:
  iterbins.append(it.iteration) 
iterbins.append(iterbins[-1]*10)

gStyle.SetLegendFillColor(0)
gStyle.SetLegendBorderSize(0)

set_palette()
###NOTE!  NEED TO WORRY ABOUT CHI2 avg (derived quanity

#migfile="/minerva/data/users/$USER/NukeHists/v1_TargetCode/minervame1D/Hists_Migration_t4_z82_Nu_v1_TargetCode.root"
migfile=sys.argv[1]
#migfile="/minerva/data/users/${USER}/NukeHists/v1_Unfolding_FullSample/AllNuME/Hists_TrueEnergy_MC_t4_z82_Nu_v1_Unfolding_FullSample.root"
migcollect=OrderedDict()
migcollect["Enu"]=histcollect(migfile,"reco_dis_v_true_dis_MC_Enu_t4_z82",";Reco Migration;True",kBlack,1,1)
migcollect["x"  ]=histcollect(migfile,"reco_dis_v_true_dis_MC_x_t4_z82",";Reco Migration;True",kBlack, 1,1)
tmpmatrix=bookHistos([migcollect[_].filename for _ in migcollect],[migcollect[_].histname for _ in migcollect],[migcollect[_].title for _ in migcollect])
migmatrix={}
for mig in migcollect:
  migmatrix[mig]=tmpmatrix[migcollect[mig].histname]

for filename in filenames:
  if filename.find("Enu")>-1:
    stackTitle="E_{#nu}" 
    titleKey = "Enu"
  elif filename.find("BjorkenX")>-1:  
    stackTitle="Bjorken X"
    titleKey = "x"

  dt=filename
  if filename.find("Closure")>-1:
    dt="Closure"+filename
  #InputFile="/minerva/data/users/$USER/NukeHists/v1_TargetCode/minervame1D/Hists_TrueEnergy_MC_t4_z82_Nu_v1_TargetCode.root"
  InputFile=migfile
  InputCollect=OrderedDict()
  InputCollect["DataTruthClosure"+filename]= histcollect(InputFile,"sample_true_mat_MC_{0}_t4_z82".format(titleKey),"",kGreen+3,26,1)
  InputCollect["DataTruth"+filename]= histcollect(InputFile,"true_dis_rw_MC_{0}_t4_z82".format(titleKey),"",kGreen+3,26,1)
  InputCollect["MCReco"+filename   ]= histcollect(InputFile,"sample_dis_reco_MC_{0}_t4_z82".format(titleKey),"",kViolet+10,27,.75)
  InputCollect["MCTruth"+filename  ]= histcollect(InputFile,"sample_true_mat_MC_{0}_t4_z82".format(titleKey),"",kOrange+1,28,.75)
  
  tmpInput=bookHistos([InputCollect[_].filename for _ in InputCollect],[InputCollect[_].histname for _ in InputCollect],[InputCollect[_].title for _ in InputCollect],[_ for _ in InputCollect])
  InputHists = {}
  for inp in InputCollect:
    InputHists[inp]=LinearizeHist(tmpInput[inp])
    InputHists[inp].GetXaxis().SetTitle(InputCollect[inp].title)
    InputHists[inp].SetMarkerColor(InputCollect[inp].color)
    InputHists[inp].SetLineColor(InputCollect[inp].color)
    InputHists[inp].SetMarkerSize(InputCollect[inp].msize)
    InputHists[inp].SetMarkerStyle(InputCollect[inp].mstyle)

  outhists  =bookHistos([filename+".root" for _ in range(2)],["m_avg_chi2_md_td_iter_chi2","m_avg_chi2_md_tmc_iter_chi2"],["" for _ in range(2)]).values()
  outhists2D=bookHistos([filename+".root" for _ in range(2)],["h_chi2_modelData_trueData_iter_chi2","h_chi2_modelData_trueMC_iter_chi2"],["" for _ in range(2)]).values()

  avgcomp=[]
  inputcomp=[]
  effcomp=[]
  effRcomp=[]
  matcomp=[]
  
  #(reco-true)/true vs iteration
  dist_bins = []
  tmp_hbin=LinearizeHist(bookHistos(filename+".root","Unfolded_Data/Stat_0_Iter_{0}".format(itercollect[0].iteration),""))
  tmp_hbin.SetName("tmp_hbin")
  for b in range(1,tmp_hbin.GetNbinsX()+2):
     dist_bins.append(tmp_hbin.GetBinLowEdge(b))
  h_pull_truth = TH2D("h_pull_truth",";Iterations (Scale #frac{{reco-true}}{{true}});Average Stat Uni {0} Bin Number".format(stackTitle),len(iterbins)-1,array('d',iterbins),len(dist_bins)-1,array('d',dist_bins))
  h_pull_truth.SetDirectory(0)

  #(reco-true)/unc vs distribution
  tmppull=[]
  tmperror =[]
  legend=TLegend(0.55,0.7,0.9,0.9)
  legend.SetFillStyle(0)
  legerr=TLegend(0.55,0.7,0.9,0.9)
  legerr.SetFillStyle(0)
  legend.SetHeader("Pull ")
  legerr.SetHeader("Error ")
  
  minmax=[0,0]  
  for i,it in enumerate(itercollect):
    tmphist=LinearizeHist(bookHistos(filename+".root","Pull_Histograms/avg_Iter_{0}_pull".format(it.iteration),""))
    tmperr =LinearizeHist(bookHistos(filename+".root","Pull_Histograms/avg_Iter_{0}_error".format(it.iteration),""))
    tmphist.SetDirectory(0); tmphist.SetLineColor(it.color); tmphist.SetMarkerColor(it.color)
    tmperr.SetDirectory(0); tmperr.SetLineColor(it.color); tmperr.SetMarkerColor(it.color)
    tmphist.SetMinimum(-1.5); tmphist.SetMaximum(1.5)
    legend.AddEntry(tmphist,it.legtitle,"l")
    legerr.AddEntry(tmphist,it.legtitle,"l")
    tmppull.append(tmphist)
    tmperr.GetYaxis().SetTitle("error")
    tmperror.append(tmperr)

  pullhists=multihist(tmppull,legend)
  errhists=multihist(tmperror,legerr)

  for it in itercollect:
    WarpCollect=OrderedDict()
    WarpCollect["Unfold_Stat0"]=histcollect(filename,"Unfolded_Data/Stat_0_Iter_{0}".format(it.iteration),";{0};".format(stackTitle),kBlack,1,1.25)
    WarpCollect["Unfold_Stat1"]=histcollect(filename,"Unfolded_Data/Stat_1_Iter_{0}".format(it.iteration),";{0}".format(stackTitle),kBlue,1,1)
    WarpCollect["Unfold_Avg"]  =histcollect(filename,"Unfolded_Data/Avg_Iter{0}_Scaled".format(it.iteration)    ,";{0}".format(stackTitle),kRed,1,.5)
    WarpCollect["Input_Stat0"]=histcollect(filename,"Stat_Varied_Input_Data/Input_Stat_0_Iter_{0}".format(it.iteration),";{0};".format(stackTitle),kBlack,1,1)
    WarpCollect["Input_Stat1"]=histcollect(filename,"Stat_Varied_Input_Data/Input_Stat_1_Iter_{0}".format(it.iteration),";{0};".format(stackTitle),kBlue,1,0.75)

    tmpWarp=bookHistos([WarpCollect[_].filename+".root" for _ in WarpCollect],[WarpCollect[_].histname for _ in WarpCollect],[WarpCollect[_].title for _ in WarpCollect])
    WarpHists=OrderedDict()
    MnvWarpHists=OrderedDict()
    for warp in WarpCollect:
      WarpHists[warp]=LinearizeHist(tmpWarp[WarpCollect[warp].histname])
      MnvWarpHists[warp]=tmpWarp[WarpCollect[warp].histname]
    #XXX get stats
    StatBinContent=0
    nUniverses=1001
    #for i in range(nUniverses):
    #  tmpStat = bookHistos(filename+".root","Unfolded_Data/Stat_{0}_Iter_{1}".format(i,it.iteration),"")
    #  tmpCov  = TH2D(tmpStat.GetTotalErrorMatrix())
    #  tmpCov.SetName("tmp_{0}_{1}_{2}".format(filename,i,it.iteration))
    #  StatBinContent+=tmpCov.GetBinContent(3,3)
    #statAvg = float(StatBinContent)/nUniverses
    statAvg=0

    for i in range(1,InputHists["DataTruth"+dt].GetNbinsX()+2):
      if InputHists["DataTruth"+dt].GetBinContent(i)==0:
        h_pull_truth.Fill(it.iteration,InputHists["DataTruth"+dt].GetBinLowEdge(i),0)
      else:
        h_pull_truth.Fill(it.iteration,InputHists["DataTruth"+dt].GetBinLowEdge(i),(WarpHists["Unfold_Avg"].GetBinContent(i)-InputHists["DataTruth"+dt].GetBinContent(i))/InputHists["DataTruth"+dt].GetBinContent(i))
        print

    for warp in WarpHists:
        WarpHists[warp].SetTitle(WarpCollect[warp].title)
        WarpHists[warp].SetMarkerColor(WarpCollect[warp].color)
        WarpHists[warp].SetLineColor(WarpCollect[warp].color)
        WarpHists[warp].SetMarkerSize(WarpCollect[warp].msize)
        #WarpHists[warp].SetMarkerStyle(WarpCollect[warp].mstyle)

    tmp_unfold_legend=TLegend(0.6,0.7,0.8,0.9)
    tmp_unfold_legend.SetFillStyle(0)
    tmp_unfold_legend.SetHeader("Unfolded Iteration {0}".format(it.iteration))    
    tmp_unfold_legend.AddEntry(WarpHists["Unfold_Stat0"],"Stat Uni 0","p")    
    tmp_unfold_legend.AddEntry(WarpHists["Unfold_Stat1"],"Stat Uni 1","p")    
    tmp_unfold_legend.AddEntry(WarpHists["Unfold_Avg"]  , "Average","p")    
    avgcomp.append(multihist([WarpHists["Unfold_Stat0"],WarpHists["Unfold_Stat1"],WarpHists["Unfold_Avg"]],tmp_unfold_legend)) 

    tmp_input_legend=TLegend(0.5,0.7,0.8,0.9)
    tmp_input_legend.SetFillStyle(0)
    tmp_input_legend.SetHeader("Input Iteration {0}".format(it.iteration))    
    tmp_input_legend.AddEntry(WarpHists["Input_Stat0"],"Stat Uni 0","p")    
    tmp_input_legend.AddEntry(WarpHists["Input_Stat1"],"Stat Uni 1","p")    
    inputcomp.append(multihist([WarpHists["Input_Stat0"],WarpHists["Input_Stat1"]],tmp_input_legend)) 
                         
    tmp_input_stat2 =WarpHists["Input_Stat1"].Clone( WarpHists["Input_Stat1"].GetName()+"_2")
    tmp_unfold_stat2=WarpHists["Unfold_Stat1"].Clone(WarpHists["Unfold_Stat1"].GetName()+"_2")
    tmp_unfold_avg2 =WarpHists["Unfold_Avg"].Clone(  WarpHists["Unfold_Avg"].GetName()+"_2")

    tmp_input_stat2.SetMarkerColor(kBlack)
    tmp_unfold_stat2.SetMarkerColor(kPink+10)
    tmp_input_stat2.SetMarkerStyle(4)
    tmp_unfold_stat2.SetMarkerStyle(30)
    tmp_unfold_avg2.SetMarkerStyle(25)

    tmp_uneff_legend=TLegend(0.6,0.7,0.8,0.9)
    tmp_uneff_legend.SetFillStyle(0)
    tmp_uneff_legend.SetHeader("Unfolding Effects Iteration {0}".format(it.iteration))
    tmp_uneff_legend.AddEntry(tmp_input_stat2,"Input Stat Uni 1","p")    
    tmp_uneff_legend.AddEntry(tmp_unfold_stat2,"Unfolded Stat Uni 1","p")    
    tmp_uneff_legend.AddEntry(tmp_unfold_avg2, "Average Unfolding","p")    
    tmp_uneff_legend.AddEntry(InputHists["DataTruth"+dt], "Data Truth","p")    
    tmp_uneff_legend.AddEntry(InputHists["MCTruth"+filename], "MC Truth","p")    
    tmp_uneff_legend.AddEntry(InputHists["MCReco"+filename], "MC Reco","p")    
    effcomp.append(multihist([tmp_input_stat2,tmp_unfold_stat2,tmp_unfold_avg2,InputHists["DataTruth"+dt],InputHists["MCTruth"+filename],InputHists["MCReco"+filename]],tmp_uneff_legend))
  
    tmp_uneff_legend2=TLegend(0.5,0.7,0.8,0.9)
    tmp_uneff_legend2.SetFillStyle(0)
    
    tmp_unfold_stat3DT=WarpHists["Unfold_Stat1"].Clone(WarpHists["Unfold_Stat1"].GetName()+"_3DT")
    tmp_unfold_stat3MT=WarpHists["Unfold_Stat1"].Clone(WarpHists["Unfold_Stat1"].GetName()+"_3MT")
    tmp_unfold_avg3DT =WarpHists["Unfold_Avg"].Clone(  WarpHists["Unfold_Avg"].GetName()+"_3DT")
    tmp_unfold_avg3MT =WarpHists["Unfold_Avg"].Clone(  WarpHists["Unfold_Avg"].GetName()+"_3MT")
    tmp_unfold_stat3DT.SetLineColor(kGreen+3)
    tmp_unfold_stat3MT.SetLineColor(kGreen-5)
    tmp_unfold_avg3DT.SetLineColor(kOrange+1) 
    tmp_unfold_avg3MT.SetLineColor(kOrange+3) 

    tmp_unfold_stat3DT.Add(InputHists["DataTruth"+dt],-1)
    tmp_unfold_stat3MT.Add(InputHists["MCTruth"+filename],-1)
    tmp_unfold_avg3DT.Add(InputHists["DataTruth"+dt],-1)
    tmp_unfold_avg3MT.Add(InputHists["MCTruth"+filename],-1)

    tmp_unfold_stat3DT.Divide(tmp_unfold_stat3DT,InputHists["DataTruth"+dt])
    tmp_unfold_stat3MT.Divide(tmp_unfold_stat3MT,InputHists["MCTruth"+filename])
    tmp_unfold_avg3DT.Divide(tmp_unfold_avg3DT,InputHists["DataTruth"+dt])
    tmp_unfold_avg3MT.Divide(tmp_unfold_avg3MT,InputHists["MCTruth"+filename])
  
    tmp_unfold_stat3DT.SetTitle(";Linear Bin Number;#frac{reco-true}{true}")
    tmp_unfold_stat3MT.SetTitle(";Linear Bin Number;#frac{reco-true}{true}")
    tmp_unfold_avg3DT.SetTitle(";Linear Bin Number;#frac{reco-true}{true}")
    tmp_unfold_avg3MT.SetTitle(";Linear Bin Number;#frac{reco-true}{true}")
  
    tmp_uneff_legend2.SetHeader("Unfolding Ratio Iteration {0}".format(it.iteration))
    tmp_uneff_legend2.AddEntry(tmp_unfold_stat3DT,"Unfolded Stat Uni 1/Data Truth","l")
    tmp_uneff_legend2.AddEntry(tmp_unfold_stat3MT,"Unfolded Stat Uni 1/MC Truth","l")
    tmp_uneff_legend2.AddEntry(tmp_unfold_avg3DT,"Unfolded Average/Data Truth","l")
    tmp_uneff_legend2.AddEntry(tmp_unfold_avg3MT,"Unfolded Average/MC Truth","l")
    effRcomp.append(multihist([tmp_unfold_stat3DT,tmp_unfold_stat3MT,tmp_unfold_avg3DT,tmp_unfold_avg3MT],tmp_uneff_legend2))

    #Matrices
    stat1_matrix = TH2D(MnvWarpHists["Unfold_Stat1"].GetTotalErrorMatrix())
    stat1_matrix.SetName("stat1_matrix_Iter{0}".format(it.iteration))
    stat1_matrix.GetXaxis().SetTitle("Iteration {0} Pushed Stat Universe 1 Cov".format(it.iteration))

    avg_matrix  = TH2D(MnvWarpHists["Unfold_Avg"].GetTotalErrorMatrix())
    avg_matrix.SetName("avg_matrix_Iter{0}".format(it.iteration))
    avg_matrix.GetXaxis().SetTitle("Iteration {0} Pushed Average Universe 0 Cov".format(it.iteration))
    print"Avg", statAvg,avg_matrix.GetBinContent(3,3)

    minz=min(stat1_matrix.GetMinimum(),avg_matrix.GetMinimum())
    maxz=max(stat1_matrix.GetMaximum(),avg_matrix.GetMaximum())
    stat1_matrix.SetMinimum(minz); stat1_matrix.SetMaximum(maxz)
    avg_matrix.SetMinimum(minz);  avg_matrix.SetMaximum(maxz)

    ##stat_matrix.SetMinimum(min(stat_matrix.GetMinimum(),stat_cov.GetMinimum())); 
    ##stat_matrix.SetMaximum(max(stat_matrix.GetMaximum(),stat_cov.GetMaximum()))
    ##stat_cov.SetMinimum(min(stat_matrix.GetMinimum(),stat_cov.GetMinimum()));    
    ##stat_cov.SetMaximum(max(stat_matrix.GetMaximum(),stat_cov.GetMaximum()))

    ##stat1_matrix.SetMinimum(min(stat1_matrix.GetMinimum(),stat1_cov.GetMinimum())); 
    ##stat1_matrix.SetMaximum(max(stat1_matrix.GetMaximum(),stat1_cov.GetMaximum()))
    ##stat1_cov.SetMinimum(min(stat1_matrix.GetMinimum(),stat1_cov.GetMinimum()));    
    ##stat1_cov.SetMaximum(max(stat1_matrix.GetMaximum(),stat1_cov.GetMaximum()))

    ##avg_matrix.SetMinimum(min(avg_matrix.GetMinimum(),avg_cov.GetMinimum())); 
    ##avg_matrix.SetMaximum(max(avg_matrix.GetMaximum(),avg_cov.GetMaximum()))
    ##avg_cov.SetMinimum(min(avg_matrix.GetMinimum(),avg_cov.GetMinimum()));    
    ##avg_cov.SetMaximum(max(avg_matrix.GetMaximum(),avg_cov.GetMaximum()))
    ###stat1_matrix.SetMinimum(minz); stat1_matrix.SetMaximum(maxz)
    ###avg_matrix.SetMinimum(minz);  avg_matrix.SetMaximum(maxz)

    #matcomp.append([stat_matrix,stat_cov,stat1_matrix,stat1_cov,avg_matrix,avg_cov])
    if filename.find("Enu")>-1:
      matcomp.append([migmatrix["Enu"],stat1_matrix,avg_matrix])
    else:
      matcomp.append([migmatrix["x"],stat1_matrix,avg_matrix])

  outname="{0}NoSquare_Linear.pdf".format(filename)
  c=TCanvas()
  c.Print(outname+"[")
  c.SetLogx()
  for outhist in outhists2D:
    outhist.Draw("colz")
    c.Print(outname)
  for outhist in outhists:
    outhist.Draw("hist")
    c.Print(outname)
    if filename.find("Bjorken")>-1:
      outhist.SetMaximum(5000)
      outhist.Draw("hist")
      c.Print(outname)

  h_pull_truth.Draw("colz")
  c.Print(outname)
  c.SetLogx(False)

  pullmax=[_.GetMaximum() for _ in pullhists.collect]
  pullmin=[_.GetMinimum() for _ in pullhists.collect]
  for i,hist in enumerate(pullhists.collect):
    hist.SetMaximum(1.5*max(pullmax))
    hist.SetMinimum(1.5*max(pullmin))
    if i==0:
      hist.Draw("hist")
    else:
      hist.Draw("hist same")
  pullhists.legend.Draw()
  c.Print(outname)

  errmax=[_.GetMaximum() for _ in errhists.collect]
  for i,hist in enumerate(errhists.collect):
    hist.SetMaximum(1.5*max(errmax))
    if i==0:
      hist.Draw("hist")
    else:
      hist.Draw("hist same")
  errhists.legend.Draw()
  c.Print(outname)

  for comp in effRcomp:
    for i,hist in enumerate(comp.collect):
      hist.SetMaximum(1.5)
      hist.SetMinimum(-1.5)
      if i==0:
        hist.Draw("hist")
      else:
        hist.Draw("hist same")
    comp.legend.Draw()
    c.Print(outname)
  for comp in effcomp:
    for i,hist in enumerate(comp.collect):
      hist.SetMaximum(1.5*hist.GetMaximum())
      if i==0:
        hist.Draw("hist p")
      else:
        hist.Draw("hist same p")
    comp.legend.Draw()
    c.Print(outname)
  for comp in avgcomp:
    for i,hist in enumerate(comp.collect):
      hist.SetMaximum(1.5*hist.GetMaximum())
      if i==0:
        hist.Draw("E0")
      else:
        hist.Draw("E0 same")
    comp.legend.Draw()
    c.Print(outname)
  for j,comp in enumerate(inputcomp):
    if j>0:
      break
    for i,hist in enumerate(comp.collect):
      hist.SetMaximum(1.5*hist.GetMaximum())
      if i==0:
        hist.Draw("E0")
      else:
        hist.Draw("E0 same")
    comp.legend.Draw()
    c.Print(outname)
  zmin = min([h.GetMinimum() for m in matcomp for h in m ])
  zmax = max([h.GetMaximum() for m in matcomp for h in m ])
  zmin = -abs(min(zmin,zmax))
  zmax = abs(min(zmin,zmax))
  for i,comp in enumerate(matcomp):
    for j,mat in enumerate(comp):
      if j==0 and i!=0:
        continue       
      if j!=0:
        mat.SetMinimum(zmin)
        mat.SetMaximum(zmax)
      mat.Draw("colz")
      c.Print(outname)
  c.Print(outname+"]")
