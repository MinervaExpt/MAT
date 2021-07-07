import ROOT
import PlotUtils
import sys

isFHC = bool(int(sys.argv[1]))
doLog = bool(int(sys.argv[2]))
doFillStyle = bool(int(sys.argv[3]))
fluxstring = "Joint_LE_ME_Flux_isFHC_%d_doLog_%d_fillStyle_%d"%(isFHC,doLog,doFillStyle)
if(isFHC):
    numu_flux = PlotUtils.FluxReweighter(14,True,PlotUtils.FluxReweighter.minervame1D1M1NWeightedAve,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv6,100)
    antinumu_flux = PlotUtils.FluxReweighter(-14,True,PlotUtils.FluxReweighter.minervame1D1M1NWeightedAve,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv6,100)
    numu_flux_LE = PlotUtils.FluxReweighter(14,True,PlotUtils.FluxReweighter.minerva13,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv5,100)
    antinumu_flux_LE = PlotUtils.FluxReweighter(-14,True,PlotUtils.FluxReweighter.minerva13,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv5,100)
elif(not isFHC):
    numu_flux = PlotUtils.FluxReweighter(14,False,PlotUtils.FluxReweighter.minervame6A,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv6,100)
    antinumu_flux = PlotUtils.FluxReweighter(-14,False,PlotUtils.FluxReweighter.minervame6A,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv6,100)
    numu_flux_LE = PlotUtils.FluxReweighter(14,True,PlotUtils.FluxReweighter.minerva5,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv5,100)
    antinumu_flux_LE = PlotUtils.FluxReweighter(-14,True,PlotUtils.FluxReweighter.minerva5,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv5,100)
else:
    sys.exit("Something weird")
    


y_max = 1e-3
y_min = 1e-6

y_max_scaled = 1.6e17
y_min_scaled = 0
if(doLog): y_min_scaled = 9e13

pot = 12.1e20
pot_LE = 4e20
if(not isFHC):
    pot = 12.4e20
    pot_LE = 1.7e20






numu = numu_flux.GetFluxReweighted(14)
antinumu = antinumu_flux.GetFluxReweighted(-14)

numu_LE = numu_flux_LE.GetFluxReweighted(14)
antinumu_LE = antinumu_flux_LE.GetFluxReweighted(-14)


#fix titl
numu.GetYaxis().CenterTitle(True)
numu.GetYaxis().SetTitleOffset(1.3)
numu.GetYaxis().SetTitleSize(0.05)
numu.GetYaxis().SetLabelSize(0.05)

numu.GetXaxis().CenterTitle(True)
numu.GetXaxis().SetTitle("Neutrino Energy (GeV)")
numu.GetXaxis().SetTitleSize(0.05)
numu.GetXaxis().SetLabelSize(0.05)

antinumu.GetYaxis().CenterTitle(True)
antinumu.GetYaxis().SetTitleOffset(1.3)
antinumu.GetYaxis().SetTitleSize(0.05)
antinumu.GetYaxis().SetLabelSize(0.05)

antinumu.GetXaxis().CenterTitle(True)
antinumu.GetXaxis().SetTitle("Neutrino Energy (GeV)")
antinumu.GetXaxis().SetTitleSize(0.05)
antinumu.GetXaxis().SetLabelSize(0.05)

numu.GetYaxis().SetTitle("#nu / m^{2} / P.O.T/ GeV")
antinumu.GetYaxis().SetTitle("#nu / m^{2} / P.O.T/ GeV")


me = ROOT.TColor.GetColor("#005b8e")
le = ROOT.TColor.GetColor("#DD7E33")

numu.SetLineColor(me)
antinumu.SetLineStyle(2)
numu_LE.SetLineColor(le)
antinumu_LE.SetLineStyle(2)
antinumu_LE.SetLineColor(le)


if(doFillStyle):
    if(isFHC):
        numu.SetFillColor(me)
        numu_LE.SetFillColor(le)
        numu.SetFillStyle(3001)
        numu_LE.SetFillStyle(3001)
    else:
        antinumu.SetFillColor(me)
        antinumu_LE.SetFillColor(le)
        antinumu.SetFillStyle(3001)
        antinumu_LE.SetFillStyle(3001)
    


leg = ROOT.TLegend(0.55,0.6,0.82,0.9)
leg.SetFillColor(0)
leg.SetLineColor(0)
if(doFillStyle):
    if(isFHC):
        leg.AddEntry(numu,"Medium Energy #scale[1.7]{#nu_{#mu}}","f")
        leg.AddEntry(antinumu,"Medium Energy #scale[1.7]{#bar{#nu}_{#mu}}","l")
        leg.AddEntry(numu_LE,"Low Energy #scale[1.7]{#nu_{#mu}}","f")
        leg.AddEntry(antinumu_LE,"Low Energy #scale[1.7]{#bar{#nu}_{#mu}}","l")
    else:
        leg.AddEntry(numu,"Medium Energy #scale[1.7]{#nu_{#mu}}","l")
        leg.AddEntry(antinumu,"Medium Energy #scale[1.7]{#bar{#nu}_{#mu}}","f")
        leg.AddEntry(numu_LE,"Low Energy #scale[1.7]{#nu_{#mu}}","l")
        leg.AddEntry(antinumu_LE,"Low Energy #scale[1.7]{#bar{#nu}_{#mu}}","f")
  
else:
    if(isFHC):
        leg.AddEntry(numu,"Medium Energy #scale[1.7]{#nu_{#mu}}","l")
        leg.AddEntry(antinumu,"Medium Energy #scale[1.7]{#bar{#nu}_{#mu}}","l")
        leg.AddEntry(numu_LE,"Low Energy #scale[1.7]{#nu_{#mu}}","l")
        leg.AddEntry(antinumu_LE,"Low Energy #scale[1.7]{#bar{#nu}_{#mu}}","l")
    else:
        leg.AddEntry(numu,"Medium Energy #scale[1.7]{#nu_{#mu}}","l")
        leg.AddEntry(antinumu,"Medium Energy #scale[1.7]{#bar{#nu}_{#mu}}","l")
        leg.AddEntry(numu_LE,"Low Energy #scale[1.7]{#nu_{#mu}}","l")
        leg.AddEntry(antinumu_LE,"Low Energy #scale[1.7]{#bar{#nu}_{#mu}}","l")
    

can = ROOT.TCanvas()
if(isFHC):
    numu.Draw("HIST")
    antinumu.Draw("SAMEHIST")
    numu_LE.Draw("SAMEHIST")
    antinumu_LE.Draw("SAMEHIST")
else:
    antinumu.Draw("HIST")
    numu.Draw("SAMEHIST")
    numu_LE.Draw("SAMEHIST")
    antinumu_LE.Draw("SAMEHIST")
    
leg.Draw("SAME")

can.SetLogy(doLog)
numu.GetXaxis().SetRangeUser(0,20)
antinumu.GetXaxis().SetRangeUser(0,20)
if(doLog): numu.GetYaxis().SetRangeUser(y_min,y_max)

can.Print("%s.eps"%(fluxstring))
can.Print("%s.png"%(fluxstring))
can.Print("%s.C"%(fluxstring))

if(isFHC): print "ME",numu.Integral(0,-1,"width"),"LE",numu_LE.Integral(0,-1,"width"),numu.Integral(0,-1,"width")/numu_LE.Integral(0,-1,"width")
else:
    print "ME",antinumu.Integral(0,-1,"width"),"LE",antinumu_LE.Integral(0,-1,"width"),antinumu.Integral(0,-1,"width")/antinumu_LE.Integral(0,-1,"width")
print "---------------------------_"


can2 = ROOT.TCanvas()
numu.Scale(pot)
antinumu.Scale(pot)
numu_LE.Scale(pot_LE)
antinumu_LE.Scale(pot_LE)

numu.GetYaxis().SetTitle("#nu / m^{2}/ GeV")
antinumu.GetYaxis().SetTitle("#nu / m^{2}/ GeV")
numu.GetXaxis().SetRangeUser(0,20)
antinumu.GetXaxis().SetRangeUser(0,20)
numu.GetYaxis().SetRangeUser(y_min_scaled,y_max_scaled)
antinumu.GetYaxis().SetRangeUser(y_min_scaled,y_max_scaled)


fluxstring = fluxstring+"_POT_%e"%(pot)

can2.SetLogy(doLog)
if(isFHC):
    numu.Draw("HIST")
    antinumu.Draw("SAMEHIST")
    numu_LE.Draw("SAMEHIST")
    antinumu_LE.Draw("SAMEHIST")
else:
    antinumu.Draw("HIST")
    numu.Draw("SAMEHIST")
    numu_LE.Draw("SAMEHIST")
    antinumu_LE.Draw("SAMEHIST")

leg.Draw("SAME")
can2.Print("%s.eps"%(fluxstring))
can2.Print("%s.png"%(fluxstring))
can2.Print("%s.C"%(fluxstring))

if(isFHC): print "ME",numu.Integral(0,-1,"width"),"LE",numu_LE.Integral(0,-1,"width"),numu.Integral(0,-1,"width")/numu_LE.Integral(0,-1,"width")
else:
    print "ME",antinumu.Integral(0,-1,"width"),"LE",antinumu_LE.Integral(0,-1,"width"),antinumu.Integral(0,-1,"width")/antinumu_LE.Integral(0,-1,"width")
raw_input("Done")
