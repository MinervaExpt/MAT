import ROOT
import PlotUtils
import sys

isFHC = bool(int(sys.argv[1]))
applyNuE = bool(int(sys.argv[2]))
isME = bool(int(sys.argv[3]))

fluxstring = "Flux_isFHC_%d_applyNuE_%d_isME_%d"%(isFHC,applyNuE,isME)

print "You will be running this plotting with isFHC=%d and useNue=%d and isME=%d"%(isFHC,applyNuE,isME)
raw_input("Is this correct? If not exit now")

if(isFHC and isME):
    numu_flux = PlotUtils.FluxReweighter(14,applyNuE,PlotUtils.FluxReweighter.minervame1D1M1NWeightedAve,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv6,100)
#    nue_flux = PlotUtils.FluxReweighter(12,applyNuE,PlotUtils.FluxReweighter.minervame1D1M1NWeightedAve,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv6,100)
    antinumu_flux = PlotUtils.FluxReweighter(-14,applyNuE,PlotUtils.FluxReweighter.minervame1D1M1NWeightedAve,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv6,100)
#    antinue_flux = PlotUtils.FluxReweighter(-12,applyNuE,PlotUtils.FluxReweighter.minervame1D1M1NWeightedAve,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv6,100)
elif(isFHC and not isME):
    numu_flux = PlotUtils.FluxReweighter(14,applyNuE,PlotUtils.FluxReweighter.minerva13,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv5,100)
#    nue_flux = PlotUtils.FluxReweighter(12,applyNuE,PlotUtils.FluxReweighter.minerva13,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv5,100)
    antinumu_flux = PlotUtils.FluxReweighter(-14,applyNuE,PlotUtils.FluxReweighter.minerva13,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv5,100)
#    antinue_flux = PlotUtils.FluxReweighter(-12,applyNuE,PlotUtils.FluxReweighter.minervame13,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv5,100)
elif(not isFHC and isME):
    if(applyNuE):
        sys.exit("You chose to run with NuE for anti-nu which is NOT support for ME. Exiting!")
    numu_flux = PlotUtils.FluxReweighter(14,applyNuE,PlotUtils.FluxReweighter.minervame6A,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv6,100)
#    nue_flux = PlotUtils.FluxReweighter(12,applyNuE,PlotUtils.FluxReweighter.minervame6A,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv6,100)
    antinumu_flux = PlotUtils.FluxReweighter(-14,applyNuE,PlotUtils.FluxReweighter.minervame6A,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv6,100)
#    antinue_flux = PlotUtils.FluxReweighter(-12,applyNuE,PlotUtils.FluxReweighter.minervame6A,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv6,100)
elif(not isFHC and not isME):
    numu_flux = PlotUtils.FluxReweighter(14,applyNuE,PlotUtils.FluxReweighter.minerva5,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv5,100)
#    nue_flux = PlotUtils.FluxReweighter(12,applyNuE,PlotUtils.FluxReweighter.minerva5,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv5,100)
    antinumu_flux = PlotUtils.FluxReweighter(-14,applyNuE,PlotUtils.FluxReweighter.minerva5,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv5,100)
#    antinue_flux = PlotUtils.FluxReweighter(-12,applyNuE,PlotUtils.FluxReweighter.minerva5,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv5,100)
else:
    sys.exit("Something weird")
    

y_max = 1e-3
y_min = 1e-6

y_max_scaled = 2e17
y_min_scaled = 1e14

pot = 12.1e20

if(not isFHC and isME):
    pot = 12.4e20
elif(isFHC and not isME):
    pot = 4e20
elif(not isFHC and not isME):
    pot = 1.7e20







numu = numu_flux.GetFluxReweighted(14)
#nue = nue_flux.GetFluxReweighted(12)
antinumu = antinumu_flux.GetFluxReweighted(-14)
#antinue = antinue_flux.GetFluxReweighted(-12)


#fix title
numu.GetYaxis().CenterTitle(True)
numu.GetYaxis().SetTitleOffset(1.3)
numu.GetYaxis().SetTitleSize(0.05)
numu.GetYaxis().SetLabelSize(0.05)

numu.GetXaxis().CenterTitle(True)
numu.GetXaxis().SetTitle("Neutrino Energy (GeV)")
numu.GetXaxis().SetTitleSize(0.05)
numu.GetXaxis().SetLabelSize(0.05)

numu.GetYaxis().SetTitle("#nu / m^{2} / P.O.T/ GeV")


antinumu.SetLineStyle(2)
#nue.SetLineStyle(3)
#antinue.SetLineStyle(4)


leg = ROOT.TLegend(0.65,0.75,0.82,0.9)
leg.SetFillColor(0)
leg.SetLineColor(0)
leg.AddEntry(numu,"#nu_{#mu}","l")
leg.AddEntry(antinumu,"#bar{#nu}_{#mu}","l")
#leg.AddEntry(nue,"#nu_{e}","l")
#leg.AddEntry(antinue,"#bar{#nu}_{e}","l")


can = ROOT.TCanvas()
numu.Draw("HIST")
#nue.Draw("SAMEHIST")
antinumu.Draw("SAMEHIST")
#antinue.Draw("SAMEHIST")
leg.Draw("SAME")

can.SetLogy(True)
numu.GetXaxis().SetRangeUser(0,20)
numu.GetYaxis().SetRangeUser(y_min,y_max)

can.Print("%s.eps"%(fluxstring))
can.Print("%s.png"%(fluxstring))
can.Print("%s.C"%(fluxstring))



can2 = ROOT.TCanvas()
numu.Scale(pot)
#nue.Scale(pot)
antinumu.Scale(pot)
#antinue.Scale(pot)
numu.GetYaxis().SetTitle("#nu / m^{2}/ GeV")
numu.GetXaxis().SetRangeUser(0,20)
numu.GetYaxis().SetRangeUser(y_min_scaled,y_max_scaled)
mylatex = ROOT.TLatex()
mylatex.SetTextSize(0.035)
mylatex.SetTextFont(42)
mystring = "Exposure : %1.2e P.O.T"%pot


fluxstring = fluxstring+"_POT_%e"%(pot)
if(applyNuE):
    mystring = mystring + "                 #nu+e constraint applied"
else:
    mystring = mystring + "                 #nu+e constraint not applied"

can2.SetLogy(True)

numu.Draw("HIST")
#nue.Draw("SAMEHIST")
antinumu.Draw("SAMEHIST")
#antinue.Draw("SAMEHIST")
leg.Draw("SAME")
mylatex.DrawLatexNDC(0.15,0.92,mystring)
can2.Print("%s.eps"%(fluxstring))
can2.Print("%s.png"%(fluxstring))
can2.Print("%s.C"%(fluxstring))


raw_input("Done")
