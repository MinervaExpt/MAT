import array
import os.path

import ROOT
ROOT.gROOT.SetBatch(True)

OUT_DIR = "/scratch/minerva/jwolcott"

import PlotUtils

pl = PlotUtils.MnvPlotter()

h = PlotUtils.MnvH1D("test", "test", 1, 0, 1)
h.AddVertErrorBand("Flux", 4)
h.AddVertErrorBand("dummy", 4)

wgts = array.array('d', [0.5, 0.75, 1.25, 1.5])
h.Fill(0.5, 100)
h.FillVertErrorBand("Flux", 0.5, wgts, 100)
h.GetVertErrorBand("Flux").SetUseSpreadError(False)

wgts = array.array('d', [i*0.2 for i in range(1,5)])
h.FillVertErrorBand("dummy", 0.5, wgts, 100)
h.GetVertErrorBand("dummy").SetUseSpreadError(False)

c = ROOT.TCanvas()
pl.DrawAllUniverses(h)
c.Print( os.path.join(OUT_DIR, "test.png") )
c.Clear()

pl.DrawErrorSummary(h, "TR", False, True, 0.00001, False, "", False) # ... lots of arguments just to get the last one (non-fractional errs)
c.Print( os.path.join(OUT_DIR, "test_errs.png") )
c.Clear()
pl.DrawErrorSummary(h, "TR", False, False) # now do fractional errs
c.Print( os.path.join(OUT_DIR, "test_errs_frac.png") )

constrainer = PlotUtils.MnvHistoConstrainer()
constrainer.LoadConstraint("test", "constraint_test_wgts.txt")

constrainer.SetSpectatorCorrectionStrategy("dummy", PlotUtils.MnvHistoConstrainer.PRESERVE_ABSOLUTE_ERR)
h_constr = constrainer.ConstrainHisto(h)

pl.DrawAllUniverses(h_constr)
c.Print(  os.path.join(OUT_DIR, "test_constrained.png") )
c.Clear()

pl.DrawErrorSummary(h_constr, "TR", False, True, 0.00001, False, "", False)
c.Print( os.path.join(OUT_DIR, "test_errs_constrained.png") )
c.Clear()
pl.DrawErrorSummary(h_constr, "TR", False, False)
c.Print( os.path.join(OUT_DIR, "test_errs_frac_constrained.png") )



