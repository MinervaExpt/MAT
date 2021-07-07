#! /bin/env python
# Filename: PrintWarpingStudy.py
# Created: 08/28/2018
# Author: Aaron Bercellie
#
# Print out nice plots taking TransWarpExtraction input

from ROOT import *
from PlotUtils import *
from collections import namedtuple
from OrderedDict import OrderedDict
from array import array
import colorsys, re, sys, getopt

drawhist =namedtuple("drawhist","hist legend opt")
iterColorCode = lambda iCol, nColors: TColor.GetColor(*colorsys.hls_to_rgb((float(iCol)/nColors)*(2.0/3),0.35,1))
binColorCode = lambda iCol, nColors: TColor.GetColor(*colorsys.hls_to_rgb((1.618*iCol)%1,0.55-(0.4*float(iCol)/nColors),1))#Using some golden ratio things

def set_palette():
  red   = array('d',[0.1098 ,0   ,0  ])
  green = array('d',[0.7098 ,0   ,0  ])
  blue  = array('d',[0.87843,0.26,0  ])
  stops = array('d',[0.0    ,0.8 ,1.0])

  TColor.CreateGradientColorTable(3,stops,red,green,blue,50)

#def set_palette_covariance():
#  red   = array('d',[0   ,1.0   ,1.0  ])
#  green = array('d',[0   ,1.0   ,0   ])
#  blue  = array('d',[1.0 ,1.0   ,0   ])
#  stops = array('d',[0.0 ,0.5 ,1.0])
#
#  TColor.CreateGradientColorTable(3,stops,red,green,blue,50)  

def myPlotStyle():
  gStyle.SetPadTopMargin(0.05);
  gStyle.SetPadRightMargin(0.1);
  gStyle.SetPadBottomMargin(0.14);
  gStyle.SetPadLeftMargin(0.15);

  gStyle.SetLabelFont(42, "XYZ");
  gStyle.SetTitleFont(42, "XYZ");

  gStyle.SetTextSize(0.08);

  gStyle.SetLabelSize(0.05,"x");
  gStyle.SetTitleSize(0.045,"x");
  gStyle.SetLabelSize(0.05,"y");
  gStyle.SetTitleSize(0.045,"y");
  gStyle.SetLabelSize(0.05,"z");
  gStyle.SetTitleSize(0.04,"z");
  gStyle.SetTitleFillColor(0);
  gStyle.SetTitleX(0.25);
  gStyle.SetTitleFontSize(0.08);
  gStyle.SetTitleOffset(1.2, "x");
  gStyle.SetTitleOffset(1.2, "Y");

  gStyle.SetMarkerStyle(20);
  gStyle.SetHistLineWidth(2);
  gStyle.SetLineStyleString(2,"[12 12]");

  #gStyle.SetErrorX(0.001);

  #gStyle.SetOptTitle(0);
  gStyle.SetOptStat(0);
  gStyle.SetOptFit(0);

  gStyle.SetTickLength(0.01, "Y");
  gStyle.SetTickLength(0.02, "X");

  gStyle.SetNdivisions(505, "XYZ");
  gStyle.SetStripDecimals(false);

def GetHist(inputfile,dirname,histname):
  if dirname!="":
    rethist = inputfile.Get("{0}/{1}".format(dirname,histname))
  else:
    rethist = inputfile.Get("{0}".format(histname))

  if not rethist:
    print("Cannot find {0}/{1}".format(dirname,histname))
  rethist.SetDirectory(0)
  return rethist

def GetHistograms(inputfile,dirname="",regex=""):
  ret=OrderedDict()
  if dirname=="":
    hists=iter(inputfile.GetListOfKeys())
  else:
    inputfile.cd(dirname)
    hists=iter(gDirectory.GetListOfKeys())

  while True:
    try:
      hist=hists.next()
      name=hist.GetName()
      rename=re.search(regex,name)

      if regex!="": #If the field isn't there, always get it, put it into the directory
        if not rename:
          continue
        else:
          ngroups=re.compile(regex).groups
          if ngroups==0:
            name=rename.group(0)
          else:
            name=""
            for _ in range(ngroups-1):
              name+="{0}_".format(rename.group(_+1))
            name+=rename.group(ngroups)

      h=hist.ReadObj()

      if not h.IsA().InheritsFrom("TH1") and not h.IsA().InheritsFrom("TMatrixT<double>"):
        del hist; del h;
        continue
      else:
        if h.IsA().InheritsFrom("TH1"):
          h.SetDirectory(0)
        ret[name]=h
        #if name in ret:
        #  ret[name].append(h)
        #else:
        #  ret[name]=[]
        #  ret[name].append(h)

        del h; del hist
 
    except StopIteration:
      break

  if hists:
    del hists
  return ret

def GetDebugStatHists(inputfile,iter_stat,iteration_bins,addLinear,maxRatio,minmaxPull):
  stats=[]
  for it in iter_stat:
    for stat in iter_stat[it]:
      if stat[0]>1e-2:
        stats.append(stat[1])
  uniq_stat=set(sorted(stats))

  persistent_hists=[]
  iterations=sorted(iter_stat.keys())
  debug_page=OrderedDict()
  for su in uniq_stat:
    smeared_hist = GetHist(inputfile,"Stat_Varied_Smeared_Data","Input_Stat_{0}_Iter_{1}".format(su,iterations[0])+addLinear)

    input_hist  = GetHist(inputfile,"Input_Hists","h_data")
    input_ratio = smeared_hist.Clone("Input_Stat_{0}_Iter_{1}_Ratio".format(su,iterations[0])+addLinear)
    input_ratio.Divide(input_ratio,input_hist)
    input_ratio.GetYaxis().SetTitle("Poisson Varied Smeared Data/Input Smeared Data")
    input_ratio.SetMinimum(0.0);input_ratio.SetMaximum(maxRatio)
    for b in range(input_ratio.GetNbinsX()+2):
      input_ratio.SetBinError(0,0)

    unfolded_hists       = GetHistograms(inputfile,"Unfolded_Data","Stat_{0}_Iter_([0-9]+)".format(su)+addLinear+"$")
    ratio_unfolded_true  = GetHistograms(inputfile,"Ratio_Unfolded_True_Histograms","Ratio_modelData_trueData_Stat_{0}_Iter_([0-9]+)".format(su)+addLinear+"$")
    ratio_unfolded_input = GetHistograms(inputfile,"Ratio_Unfolded_Input_Histograms","Ratio_modelData_input_Stat_{0}_Iter_([0-9]+)".format(su)+addLinear+"$")

    stack_xtitle=unfolded_hists.values()[0].GetXaxis().GetTitle()
    iterations = sorted(unfolded_hists.keys(),key=lambda x: int(x))

    palsize = len(unfolded_hists)
    Unfolded_legend = TLegend(0.7,0.7,0.9,0.94)
    Unfolded_legend.SetTextSize(17); Unfolded_legend.SetFillStyle(0)
    Unfolded_Stack  = THStack("hs_Unfolded_Stat_{0}".format(su),";#bf{{Stat Universe {0}}} {1};{2}".format(su,stack_xtitle,unfolded_hists.values()[0].GetYaxis().GetTitle()))
    for color,it in enumerate(iterations):
      hist = unfolded_hists[it]
      if color==0:
        stack_xtitle=hist.GetXaxis().GetTitle()
      hist.SetLineColor(iterColorCode(color,palsize))
      Unfolded_legend.AddEntry(hist,"Iteration {0}".format(it),"l")
      Unfolded_Stack.Add(hist,"hist")
    Unfolded_Stack.Add(smeared_hist,"E0"); Unfolded_legend.AddEntry(smeared_hist,"Poisson Varied Smeared Data","p")
      
    palsize = len(ratio_unfolded_true)
    UT_legend = TLegend(0.7,0.7,0.9,0.94)
    UT_legend.SetTextSize(17); UT_legend.SetFillStyle(0)
    UT_Stack  = THStack("hs_Ratio_Unfolded_True_Stat_{0}".format(su),";#bf{{Stat Universe {0}}} {1};{2}".format(su,stack_xtitle,"Unfolded Data/True Data"))
    UT_Stack.SetMinimum(0.0);UT_Stack.SetMaximum(maxRatio)
    for color,it in enumerate(iterations):
      hist = ratio_unfolded_true[it]
      if color==0:
        stack_xtitle=hist.GetXaxis().GetTitle()
      hist.SetLineColor(iterColorCode(color,palsize))
      UT_legend.AddEntry(hist,"Iteration {0}".format(it),"l")
      UT_Stack.Add(hist)
      
    UI_legend = TLegend(0.7,0.7,0.9,0.94)
    UI_legend.SetTextSize(17); UI_legend.SetFillStyle(0)
    UI_Stack  = THStack("hs_Ratio_Unfolded_Input_Stat_{0}".format(su),";#bf{{Stat Universe {0}}} {1};{2}".format(su,stack_xtitle,"Unfolded Data/Poisson Varied Input Data"))
    UI_Stack.SetMinimum(0.0);UI_Stack.SetMaximum(maxRatio)
    for color,it in enumerate(iterations):
      hist = ratio_unfolded_input[it]
      if color==0:
        stack_xtitle=hist.GetXaxis().GetTitle()
      hist.SetLineColor(iterColorCode(color,palsize))
      UI_legend.AddEntry(hist,"Iteration {0}".format(it),"l")
      UI_Stack.Add(hist)

    debug_page["Stat Uni {0} P1".format(su)]=OrderedDict() 
    debug_page["Stat Uni {0} P1".format(su)]["Unfold Input"] = drawhist( Unfolded_Stack,Unfolded_legend,"nostack")
    debug_page["Stat Uni {0} P1".format(su)]["Input Ratio"]  = drawhist( input_ratio,None,"hist") 
    debug_page["Stat Uni {0} P1".format(su)]["Unfold/True"]  = drawhist( UT_Stack,UT_legend,"nostack")
    debug_page["Stat Uni {0} P1".format(su)]["Unfold/Input"] = drawhist( UI_Stack,UI_legend,"nostack")
    
    pull_hists = GetHistograms(inputfile,"Pull_Histograms","Stat_{0}_Iter_([0-9]+)_pull".format(su)+addLinear+"$")

    pull_legend = TLegend(0.7,0.7,0.9,0.94)
    pull_legend.SetTextSize(17); pull_legend.SetFillStyle(0)
    pull_Stack  = THStack("hs_pull_Stat_{0}".format(su),";#bf{{Stat Universe {0}}} {1};{2}".format(su,stack_xtitle,"#frac{reco-true}{unc}"))
    pull_Stack.SetMinimum(minmaxPull[0]);pull_Stack.SetMaximum(minmaxPull[1])
    for color,it in enumerate(iterations):
      hist = pull_hists[it]
      if color==0:
        stack_xtitle=hist.GetXaxis().GetTitle()
      hist.SetLineColor(iterColorCode(color,palsize))
      pull_legend.AddEntry(hist,"Iteration {0}".format(it),"l")
      pull_Stack.Add(hist)

    chi2_matrices = GetHistograms(inputfile,"Chi2Maps","Chi2Map_Stat_{0}_Iter_([0-9]+)".format(su))
    bin_chi2_iter_legend = TLegend(0.7,0.7,0.9,0.94)
    bin_chi2_iter_legend.SetTextSize(17); bin_chi2_iter_legend.SetFillStyle(0)
    bin_chi2_iter_Stack  = THStack("hs_bin_chi2_iter_Stat_{0}".format(su),";#bf{{Stat Universe {0}}} Iterations;{1}".format(su,"#chi^{2}"))
    for b in range(1,input_hist.GetNbinsX()+1):
      bin_chi2_iterations = TH1D("Bin_{0}_Stat_{1}_Iter_Chi2".format(b,su),";#bf{{ Bin {0} }} Iterations;#chi^{{2}}".format(b),len(iteration_bins)-1,array('d',iteration_bins))
      bin_chi2_iterations.SetDirectory(0)

      #PyRoot does cleanup after we return.  Returning this stops cleanup on these histograms
      persistent_hists.append(bin_chi2_iterations)
      
      bin_chi2 = dict()
      chi2_nCols = chi2_matrices.values()[0].GetNcols()
      bin_chi2_iterations.SetLineColor(binColorCode(b,input_hist.GetNbinsX()))

      bin_chi2 = dict()
      chi2_nCols = chi2_matrices.values()[0].GetNcols()
      if int(b)-1 < chi2_nCols:
        for it in chi2_matrices:
          tmp_bin_chi2=0
          for col in range(chi2_nCols):
            tmp_bin_chi2 = tmp_bin_chi2 + chi2_matrices[it](int(b)-1,col)
          bin_chi2_iterations.SetBinContent(int(it),tmp_bin_chi2)
      bin_chi2_iter_legend.AddEntry(bin_chi2_iterations,"Bin {0}".format(b),"l")
      bin_chi2_iter_Stack.Add(bin_chi2_iterations,"hist")

    bin_chi2_iter_Stack.SetMaximum(1.5*bin_chi2_iter_Stack.GetMaximum())
 
    debug_page["Stat Uni {0} P2".format(su)]=OrderedDict() 
    debug_page["Stat Uni {0} P2".format(su)]["Pull"] = drawhist( pull_Stack, pull_legend,"nostack")
    debug_page["Stat Uni {0} P2".format(su)]["Chi2"] = drawhist( bin_chi2_iter_Stack, bin_chi2_iter_legend,"nostack") 

  return debug_page,persistent_hists

def DebugStatUniverses(inputfile,nTopUniverses):
  h_chi2_stat = GetHist(inputfile,"Chi2_Iteration_Dists","h_chi2_modelData_trueData_iter_stat")
  nStats    = h_chi2_stat.GetNbinsY()
  nIterBins = h_chi2_stat.GetNbinsX()

  iter_stat_debug = OrderedDict()
  for it in range(1,nIterBins+1):
    if int(h_chi2_stat.GetXaxis().GetBinLowEdge(it))==0:
      continue
    chi2_stat = [h_chi2_stat.GetBinContent(it,_) for _ in range(1,nStats+1)]
    top_chi2 = sorted(zip(chi2_stat,range(nStats)),reverse=True)[:nTopUniverses]
    iter_stat_debug[int(h_chi2_stat.GetXaxis().GetBinLowEdge(it))] = top_chi2

  return iter_stat_debug 

def DrawTruncate(canvas,outname,page,bLogIteration):
  canvas.Divide(2,2)
  for i,histname in enumerate(page):
    dhist = page[histname]
    canvas.cd(i+1)
    if bLogIteration and dhist.opt.find("LOG")>-1:
      if dhist.opt.find("LOGX")>-1:
        canvas.SetLogx(1)
      if dhist.opt.find("LOGY")>-1:
        canvas.SetLogy(1)
    dhist.hist.Draw(dhist.opt)
    if dhist.legend:
      dhist.legend.Draw()
  canvas.Print("{0}.pdf".format(outname))
  for i in range(1,5):
    canvas.cd(i)
    canvas.SetLogx(0); canvas.SetLogy(0);
  canvas.Clear()

def DrawRowNormalized( h_migration ):
  print h_migration.GetName()
  nbinsX = h_migration.GetNbinsX()+2
  nbinsY = h_migration.GetNbinsY()+2
  m_migration = TMatrixD( nbinsX, nbinsY )
  for yBin in range(0,h_migration.GetNbinsY()+2):
    norm = 0
    for xBin in range(0,h_migration.GetNbinsX()+2):
      norm += h_migration.GetBinContent(xBin,yBin);
    if fabs(norm) > 1E-8:
      for xBin in range(0,h_migration.GetNbinsX()+2):
        percentage =  100 * h_migration.GetBinContent(xBin,yBin) / norm;
        m_migration[yBin][xBin] = percentage; 
  tmp = TH2D( m_migration)
  tmp.SetDirectory(0)
  tmp.GetXaxis().SetTitle("Reco "+h_migration.GetXaxis().GetTitle() + " Bins")
  tmp.GetYaxis().SetTitle("True "+h_migration.GetYaxis().GetTitle() + " Bins")

  #gStyle.SetHistMinimumZero(False)
  tmp.SetMinimum(1)
  gStyle.SetPaintTextFormat("2.0f")
  tmp.SetMarkerSize(2)
  tmp.DrawCopy("colz text")
  #gStyle.SetHistMinimumZero(True)

def DrawOutput(outname,dhists,bTruncate,bLogIteration):
  canvas = TCanvas("canvas"+outname,"canvas",1600,1200)
  canvas.Print("{0}.pdf[".format(outname))
  for pagename in dhists:
    page = dhists[pagename]
    if bTruncate and pagename!="Migration":
      DrawTruncate(canvas,outname,page,bLogIteration)
    else:
      for histname in page:
        dhist = page[histname]
        canvas.cd(0)
        if bLogIteration and dhist.opt.find("LOG")>-1:
          if dhist.opt.find("LOGX")>-1:
            canvas.SetLogx(1)
          if dhist.opt.find("LOGY")>-1:
            canvas.SetLogy(1)
        if dhist.opt.find("matrix")>-1:#Draw the migration matrix with row norm and numbers
          DrawRowNormalized( dhist.hist )
        else: #Regular draw
          dhist.hist.Draw(dhist.opt)
        if dhist.legend:
          dhist.legend.Draw()
        canvas.Print("{0}.pdf".format(outname))
        canvas.SetLogx(0); canvas.SetLogy(0);
        canvas.Clear()
         
  canvas.Print("{0}.pdf]".format(outname))
  del canvas

def main(filename,outname,bTruncate,bPrintLinear,bLogIteration,nDebug,maxRatio,minmaxPull):
  myPlotStyle()
  set_palette()
  gStyle.SetLegendFont(43)
  gStyle.SetLegendFillColor(0)
  gStyle.SetLegendBorderSize(0)

  #Print out linearized hists
  if bPrintLinear:
    addLinear="_Linear"
  else:
    addLinear=""

  infile = TFile(filename)

  pagehist=OrderedDict()

  ###########################################################################################
  # Overview histograms
  ###########################################################################################

  #Input histograms
  pagehist["Input"] = OrderedDict() 
  pagehist["Input"]["Input Smeared Data"]   = drawhist( GetHist(infile,"Input_Hists","h_data"+addLinear      ),None,"E0") 
  pagehist["Input"]["True Data"]      = drawhist( GetHist(infile,"Input_Hists","h_data_truth"+addLinear),None,"E0") 
  pagehist["Input"]["Reco MC"]        = drawhist( GetHist(infile,"Input_Hists","h_mc_reco"+addLinear   ),None,"E0") 
  pagehist["Input"]["Truth MC"]       = drawhist( GetHist(infile,"Input_Hists","h_mc_truth"+addLinear  ),None,"E0") 
  for histname in pagehist["Input"]:
    pagehist["Input"][histname].hist.GetXaxis().SetTitle("#bf{{ {0} }} {1}".format(histname,pagehist["Input"][histname].hist.GetXaxis().GetTitle()))

  #Migration Matrix (want to draw separately)
  pagehist["Migration"] = OrderedDict() 
  pagehist["Migration"]["Migration"] = drawhist( GetHist(infile,"Input_Hists","h_migration_matrix"),None,"colz") 
  pagehist["Migration"]["MigrationMatrix"] = drawhist( GetHist(infile,"Input_Hists","h_migration_matrix"),None,"colz matrix") 

  #Chi2 vs iterations   
  chi2_avg =    GetHist(infile,"Chi2_Iteration_Dists","m_avg_chi2_modelData_trueData_iter_chi2")
  chi2_median = GetHist(infile,"Chi2_Iteration_Dists","h_median_chi2_modelData_trueData_iter_chi2")
  chi2_avg_trunc = GetHist(infile,"Chi2_Iteration_Dists","m_avg_chi2_modelData_trueData_iter_chi2_truncated")
  chi2_avg_trunc.SetMarkerColor(kRed); chi2_avg_trunc.SetLineColor(kRed);  
  chi2_avg.SetMarkerSize(1.5); chi2_avg_trunc.SetMarkerSize(1.5);  
  chi2_avg_trunc.SetTitle(";Iterations;Truncated Avg #chi^{2}")
  #Try to get a truncated view with the median as a way of setting the maximum
  #trunc_avg_bincontents = [ chi2_avg_trunc.GetBinContent(_) for _ in range(chi2_avg_trunc.GetNbinsX()) ]
  #trunc_avg_bincontents.sort()
  #chi2_avg_trunc.SetMaximum(1.2*trunc_avg_bincontents[len(trunc_avg_bincontents)//2])

  avg_chi2_legend = TLegend(0.7,0.7,0.9,0.94)
  avg_chi2_legend.SetTextSize(17); avg_chi2_legend.SetFillStyle(0)
  avg_chi2_Stack  = THStack("hs_avg_chi2",";Iterations;{1}".format(chi2_avg.GetXaxis().SetTitle(),"Average #chi^{2}"))

  avg_chi2_Stack.Add(chi2_avg,"E0")
  avg_chi2_Stack.Add(chi2_avg_trunc,"E0")
  avg_chi2_legend.AddEntry(chi2_avg,"Average","p")
  avg_chi2_legend.AddEntry(chi2_avg_trunc,"Truncated Average","p")
  avg_chi2_Stack.SetMinimum(0.0); avg_chi2_Stack.SetMaximum(1.5*chi2_avg.GetMaximum())

  pagehist["Chi2 Ratio"] = OrderedDict()
  pagehist["Chi2 Ratio"]["2D"]  = drawhist( GetHist(infile,"Chi2_Iteration_Dists","h_chi2_modelData_trueData_iter_chi2")    ,None,"colz LOGX")
  pagehist["Chi2 Ratio"]["Avg"] = drawhist( avg_chi2_Stack,avg_chi2_legend,"nostack LOGX")
  pagehist["Chi2 Ratio"]["TruncAvg"] = drawhist( chi2_avg_trunc,None,"LOGX LOGY")
  pagehist["Chi2 Ratio"]["Median"] = drawhist( chi2_median,None,"LOGX LOGY")

  #ratio by iterations
  ratio_unfolded_true  = GetHistograms(infile,"Ratio_Unfolded_True_Histograms","Ratio_modelData_trueData_Stat_0_Iter_([0-9]+)"+addLinear+"$")
  ratio_unfolded_input = GetHistograms(infile,"Ratio_Unfolded_Input_Histograms","Ratio_modelData_input_Stat_0_Iter_([0-9]+)"+addLinear+"$")

  stack_xtitle=ratio_unfolded_true.values()[0].GetXaxis().GetTitle()

  palsize = len(ratio_unfolded_true)
  UT_legend = TLegend(0.7,0.7,0.9,0.94)
  UT_legend.SetTextSize(17); UT_legend.SetFillStyle(0)
  UT_Stack  = THStack("hs_Ratio_Unfolded_True_Stat_0",";#bf{{Stat Universe 0}} {0};{1}".format(stack_xtitle,"Unfolded Data/True Data"))
  UT_Stack.SetMinimum(0.0);UT_Stack.SetMaximum(maxRatio)
  for color,it in enumerate(sorted(ratio_unfolded_true, key=lambda x: int(x))):
    hist = ratio_unfolded_true[it]
    if color==0:
      stack_xtitle=hist.GetXaxis().GetTitle()
    hist.SetLineColor(iterColorCode(color,palsize))
    UT_legend.AddEntry(hist,"Iteration {0}".format(it),"l")
    UT_Stack.Add(hist)
    
  UI_legend = TLegend(0.7,0.7,0.9,0.94)
  UI_legend.SetTextSize(17); UI_legend.SetFillStyle(0)
  UI_Stack  = THStack("hs_Ratio_Unfolded_TInput_Stat_0",";#bf{{Stat Universe 0}} {0};{1}".format(stack_xtitle,"Unfolded Data/Poisson Varied Input Data"))
  UI_Stack.SetMinimum(0.0);UI_Stack.SetMaximum(maxRatio)
  for color,it in enumerate(sorted(ratio_unfolded_input, key=lambda x: int(x))):
    hist = ratio_unfolded_input[it]
    if color==0:
      stack_xtitle=hist.GetXaxis().GetTitle()
    hist.SetLineColor(iterColorCode(color,palsize))
    UI_legend.AddEntry(hist,"Iteration {0}".format(it),"l")
    UI_Stack.Add(hist)
    
  pagehist["Chi2 Ratio"]["Unfold/True"] = drawhist( UT_Stack,UT_legend,"nostack")
  pagehist["Chi2 Ratio"]["Unfold/Input"]= drawhist( UI_Stack,UI_legend,"nostack")

  #####################################V######################################################
  # Pull histograms
  ###########################################################################################
  avg_pull_hists = GetHistograms(infile,"Average_Pull_Histograms","avg_Iter_([0-9]+)_pull"+addLinear+"$")

  for i,it in enumerate(sorted(avg_pull_hists,key=lambda x: int(x))):
    if i%4 == 0:
      pullpage = i/4
      pagehist["Pull Hists {0}".format(pullpage)]    = OrderedDict()
    avg_pull_hists[it].SetMinimum(minmaxPull[0]); avg_pull_hists[it].SetMaximum(minmaxPull[1])
    avg_pull_hists[it].GetXaxis().SetTitle("#bf{{ Iteration {0} }} {1}".format(it,avg_pull_hists[it].GetXaxis().GetTitle()))
    avg_pull_hists[it].SetMarkerSize(1.5)
    pagehist["Pull Hists {0}".format(pullpage)][i] = drawhist( avg_pull_hists[it], None, "E0" )

  #####################################V######################################################
  # Bin histograms
  ###########################################################################################
  bin_pull_hists                  = GetHistograms(infile,"Bin_Pull_Histograms","^Bin_([0-9]+)_Iter_pull")
  avg_bin_pull_hists              = GetHistograms(infile,"Bin_Pull_Histograms","Avg_Bin_([0-9]+)_Iter_pull")
  bin_ratio_unfold_true_hists     = GetHistograms(infile,"Bin_Ratio_Unfolded_True_Histograms","^Bin_([0-9]+)_Iter_Ratio_modelData_trueData")
  avg_bin_ratio_unfold_true_hists = GetHistograms(infile,"Bin_Ratio_Unfolded_True_Histograms","Avg_Bin_([0-9]+)_Iter_Ratio_modelData_trueData")
  bin_chi2_iterations             = GetHistograms(infile,"Bin_Chi2_Iteration_Dists","^Bin_([0-9]+)_Iter_Chi2")
  avg_bin_chi2_iterations         = GetHistograms(infile,"Bin_Chi2_Iteration_Dists","Avg_Bin_([0-9]+)_Iter_Chi2")

  maxbin = max(*[int(_) for _ in bin_pull_hists])
  #print bin_pull_hists
  #print avg_bin_chi2_iterations
  #Get some binning
  ndf = pagehist["Input"]["True Data"].hist.GetNbinsX()
  iteration_bins = [pagehist["Chi2 Ratio"]["2D"].hist.GetXaxis().GetBinLowEdge(_) for _ in range(1,pagehist["Chi2 Ratio"]["2D"].hist.GetNbinsX()+2)] 
  tmp_chi2_bins = [pagehist["Chi2 Ratio"]["2D"].hist.GetYaxis().GetBinLowEdge(_) for _ in range(1,pagehist["Chi2 Ratio"]["2D"].hist.GetNbinsY()+2)] 
  chi2_bins =[float(_)/10 for _ in tmp_chi2_bins] 

  for i, b in enumerate(sorted(bin_pull_hists,key=lambda x:int(x))):
  #for i, b in enumerate(sorted(avg_bin_chi2_iterations,key=lambda x:int(x))):
    avg_bin_pull_hists[b].SetMarkerSize(1.5)              
    avg_bin_ratio_unfold_true_hists[b].SetMarkerSize(1.5)
    avg_bin_pull_hists[b].SetMinimum(minmaxPull[0]); avg_bin_pull_hists[b].SetMaximum(minmaxPull[1]);
    avg_bin_ratio_unfold_true_hists[b].SetMinimum(0.0); avg_bin_ratio_unfold_true_hists[b].SetMaximum(maxRatio);
    avg_bin_pull_hists[b].GetYaxis().SetTitle("Average {0}".format(avg_bin_pull_hists[b].GetYaxis().GetTitle()))
    avg_bin_ratio_unfold_true_hists[b].GetYaxis().SetTitle("Average {0}".format(avg_bin_ratio_unfold_true_hists[b].GetYaxis().GetTitle()))

    pagehist["Bin {0} Hists".format(b)]=OrderedDict()
    pagehist["Bin {0} Hists".format(b)]["Bin Pull"]                  = drawhist( bin_pull_hists[b]                 , None, "colz LOGX" ) 
    pagehist["Bin {0} Hists".format(b)]["Avg Bin Pull"]              = drawhist( avg_bin_pull_hists[b]             , None, "E0 LOGX"   ) 
    pagehist["Bin {0} Hists".format(b)]["Bin Ratio Unfold/True"]     = drawhist( bin_ratio_unfold_true_hists[b]    , None, "colz LOGX" ) 
    pagehist["Bin {0} Hists".format(b)]["Avg Bin Ratio Unfold/True"] = drawhist( avg_bin_ratio_unfold_true_hists[b], None, "E0 LOGX"   ) 

    if int(b) == 0 or int(b) == maxbin:#Chi2 info does not exist for this
      continue
    avg_bin_chi2_iterations[b].SetMarkerSize(1.5)
    pagehist["Bin {0} Chi2 Hists".format(b)] = OrderedDict() 
    pagehist["Bin {0} Chi2 Hists".format(b)]["Bin Chi2"]             = drawhist( bin_chi2_iterations[b]     , None, "colz LOGX" ) 
    pagehist["Bin {0} Chi2 Hists".format(b)]["Avg Bin Chi2"]         = drawhist( avg_bin_chi2_iterations[b] , None, "E0 LOGX"   ) 

  if nDebug>0:
    iter_stat_debug = DebugStatUniverses(infile,nDebug)
    debug_hists,dummy_pyroot = GetDebugStatHists(infile,iter_stat_debug,iteration_bins,addLinear,maxRatio,minmaxPull)

  infile.Close()

  DrawOutput(outname,pagehist,bTruncate,bLogIteration)

  if nDebug>0:
    DrawOutput(outname+"_Debug",debug_hists,bTruncate,bLogIteration)
    print
    print "Largest Chi2 Per Iteration" 
    print "|{0}|".format("-"*41)
    print "| Iteration | Stat Universe | Chi2        |"
    for it in iter_stat_debug:
      print "|{0}|".format("-"*41)
      if len(iter_stat_debug[it])>0:
        for su in range(nDebug):
          if su==0:
            itprint = it
          else:
            itprint = "         "
          print "| {0:>9} | {1:>13} | {2:11.2f} |".format(itprint,iter_stat_debug[it][su][1],iter_stat_debug[it][su][0])
    print "|{0}|".format("-"*41)

############################################################################################
#Input/Output/Options
############################################################################################
if __name__ == "__main__":
  helptext='PrintWarpingStudy.py [-i <input file name> -o <output file name>]\n'\
           'Prints out some nicely formatted plots for warping studies       \n'\
           'Required Arguments\n'\
           '                 -i --input      Input File Name\n'\
           '                 -o --output     Output file name, no file extension\n'\
           'Optional Arguments\n'\
           '                 -d --debug      (Default: 0) The number of stat universes per iteration looked at when looking for\n'\
           '                                 high chi^2. If 0, do not search for universes with high chi^2.  Outputs extra file\n'\
           '                                 <output_file>_Debug.pdf\n'\
           '                 -t --truncated  Puts more plots on each page\n'\
           '                 -l --linearbins Uses linearized histograms outputted by TransWarpExtraction\n'\
           '                 -L --logscale   Displays iterations on a log scale\n'\
           '                 -h --help       Displays this message\n'

  try:
    opts,args=getopt.getopt(sys.argv[1:],"i:o:d:tlLh",["input=","output=","debug=","truncated","linearbins","logscale","help"])
  except getopt.GetoptError:
    print helptext
    sys.exit(1)

  gROOT.SetBatch(True)

  input_name  = ""
  output_name = ""
  bTruncate   = False
  bLinearBins = False
  bLogScale   = False
  nDebug      = 0

  for opt, arg in opts:
    if opt in ("-h","--help"):
      print helptext
      sys.exit(0)
    if opt in ("-i","--input"):
      input_name = arg
    if opt in ("-o","--output"):
      output_name = arg
    if opt in ("-t","--truncated"):
      bTruncate = True
    if opt in ("-l","--linearbins"):
      bLinearBins = True
    if opt in ("-L","--logscale"):
      bLogScale = True
    if opt in ("-d","--debug"):
      try:
        nDebug = int(arg)
      except ValueError:
        o#rint "Invalid Debug Numter: ",arg
        sys.exit(1)

  #Possible user options in the future, putting them here
  maxRatio = 2.0
  minmaxPull = [-3.0,3.0]

  main(input_name,output_name,bTruncate,bLinearBins,bLogScale,nDebug,maxRatio,minmaxPull)
