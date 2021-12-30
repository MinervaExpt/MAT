import os
import math

import ROOT
ROOT.gROOT.SetBatch(True)
import PlotUtils

# getntargets(material, isMC, NNucleons = True)
# material = tracker, carbon, water, iron, iron_t1, lead, lead_t5

# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)

DO_UNIVERSES = False


plotter = PlotUtils.MnvPlotter()
plotter.SetRedHeatPalette()


# Calculate flux using the technique in GENIEXSecExtract/src/XSecLooper.cxx's getFluxHist()
# Flux = rate / xsec
# xsec is integral of spline from 
# $MPARAMFILESROOT/data/GENIE/spline_files/gxspl-nuclear-MINERVA_Full_v2126.root

supported_materials = ("carbon", "iron", "lead", "oxygen")



def getspline(material):
    assert material in supported_materials, "Can't find material %s in supported_materials" % material
    
    tf = ROOT.TFile("/minerva/app/users/kleykamp/cmtuser/Minerva_v22r1p1/MParamFiles/data/GENIE/spline_files/gxspl-nuclear-MINERVA_Full_v2126.root")
    splinename = None
    if material == "carbon":
        splinename = "nu_mu_C12/tot_cc"
    if material == "lead":
        splinename = "nu_mu_Pb207/tot_cc"
    if material == "iron":
        splinename = "nu_mu_Fe56/tot_cc"
    if material == "oxygen":
        splinename = "nu_mu_O16/tot_cc"
    assert splinename != None, "Wasn't able to get a spline name from material name %s" % material
    
    out = tf.Get(splinename)
    assert out != None, "Spline is none %s for spline name %s" % (out, splinename)
    ## Convert to actual function so that integrals are smoother
    funcname = "func_" + material
    func = ROOT.TF1(funcname, "pol6", 0, 100)
    out.Fit(funcname)
    return func
    
def unit_test_getspline():
    carbon = getspline("carbon")
    iron = getspline("iron")
    lead = getspline("lead")
    oxygen = getspline("oxygen")
    c = ROOT.TCanvas('c', 'c', 800, 800)
    carbon.SetLineColor(ROOT.kRed)
    iron.SetLineColor(ROOT.kBlue)
    oxygen.SetLineColor(ROOT.kGreen)
    lead.SetLineColor(ROOT.kBlack)
    lead.Draw()
    carbon.Draw("same")
    iron.Draw("same")
    oxygen.Draw("same")
    c.Print("unit_tests/getspline.png")
    return True

'''
  KEY: PlotUtils::MnvH1D        truehydrogen;1  Hydrogen
  KEY: PlotUtils::MnvH1D        truecarbon;1    Carbon
  KEY: PlotUtils::MnvH1D        trueoxygen;1    Oxygen
  KEY: PlotUtils::MnvH1D        truesteel;1     Iron
  KEY: PlotUtils::MnvH1D        truelead;1      Lead
'''    
def gethistname(material):
    if material == "tracker" or material == "carbon":
        return "truecarbon"
    if "iron" in material:
        return "truesteel"
    if "lead" in material:
        return "truelead"
    if material == "water":
        return "trueoxygen"

def getsplinename(material):
    if material == "tracker" or material == "carbon":
        return "carbon"
    if "iron" in material:
        return "iron"
    if "lead" in material:
        return "lead"
    if "water" in material:
        return "oxygen"
        
def make_clean_copy(hist, add_error_bands = True):
    # Try to copy the histogram and make a blank one
    new = type(hist)(hist)
    # Removes seg fault by taking direct control of this histogram.
    # See https://root.cern.ch/root/htmldoc/guides/users-guide/ObjectOwnership.html
    new.SetDirectory(0)
    #new = copy(hist)
    new.Reset()
    if add_error_bands:
        new.AddMissingErrorBandsAndFillWithCV(hist)
    return new
    

def getntargets(material, isMC, NNucleons = True):
    apothem = 850.0
    
    # TODO do we want this at 0 or the default 25mm? Right now the truth cuts don't cut in x or y (ie leave a buffer between targets) so 0 is correct.
    PlotUtils.TargetUtils.Get().SetDistToDivCut(0)
    
    usual_uncertainty = 0.005 # most targets have 0.005% uncertainty except water
    #water_uncertainty = 0.043 # 4.3% uncertainty. See studies/water_target_mass.txt
    water_uncertainty = 0.01 # 1% uncertainty as of 2/3. Based on Howard's study
    
    if NNucleons:
        trackerfunc = PlotUtils.TargetUtils.Get().GetTrackerNNucleons
        targetfunc = PlotUtils.TargetUtils.Get().GetPassiveTargetNNucleons
    else:
        trackerfunc = PlotUtils.TargetUtils.Get().GetTrackerNNeutrons
        targetfunc = PlotUtils.TargetUtils.Get().GetPassiveTargetNNeutrons
    
    if material == "tracker":
        return trackerfunc(5890.0, 8467.0, isMC, apothem), usual_uncertainty
        #return PlotUtils.TargetUtils.Get().GetTrackerNAtoms(5890.0, 8467.0, isMC, 850.0)
    if "iron_" in material and material != "iron_t15":
        i = int(material[-1:])
        return targetfunc(i, 26, isMC, apothem), usual_uncertainty
    if material == "iron":
        out = 0
        for i in [2, 3, 5]:  # Don't include T1
            out += targetfunc(i, 26, isMC, apothem)
        return out, usual_uncertainty
    if material == "iron_t15":
        out = 0
        for i in [1, 2, 3, 5]:  # include T1
            out += targetfunc(i, 26, isMC, apothem)
        return out, usual_uncertainty
    if "lead_" in material and material != "lead_t15":
        i = int(material[-1:])
        return targetfunc(i, 82, isMC, apothem), usual_uncertainty
    if material == "lead":
        out = 0
        for i in [2, 3, 4, 5]:  # Don't include T1
            out += targetfunc(i, 82, isMC, apothem)
        return out, usual_uncertainty
    if material == "lead_t15":
        out = 0
        for i in [1, 2, 3, 4, 5]:  # include T1
            out += targetfunc(i, 82, isMC, apothem)
        return out, usual_uncertainty
    if material == "carbon":
        return targetfunc(3, 6, isMC, apothem), usual_uncertainty
    water_cut = None
    if "water" in material:
        if "apothem" in material: water_cut = "apothem"
        if "radius" in material: water_cut = "radius"
        material = "water"
    assert water_cut in ("radius", "apothem", None)
    if material == "water":
        if isMC:
            if water_cut == None: WATER_MASS = 638.94 # kg, calculated from Kevin's parameterization.
            if water_cut == "radius": WATER_MASS = 459.12 # kg, calculated from Kevin's parameterization with a radius of 900mm
            if water_cut == "apothem": WATER_MASS = 452.90 # kg, calculated from 900mm radius but corrected for 850mm apothem
        else:
            if water_cut == None: WATER_MASS = 643.400 # kg, from ECL 50679
            if water_cut == "radius": WATER_MASS = 424.2 # kg, from Howard's calculation with a radius of 900mm
            if water_cut == "apothem": WATER_MASS = 428.9 # kg, from Howard's calculation using an apothem of 850mm
        WATER_MASS *= 1000 # convert to grams 
        if NNucleons:
            n_targets = 16 + 2 # protons/neutrons in oxygen + 2 protons in hydrogen
        else:
            n_targets = 8 # neutrons
        # n molecules = Water mass / (grams/mol) * (n molecules / mol [N_A])
        n_water_molecules = (WATER_MASS / 18.015) * 6.0221409e+23
        return n_targets * n_water_molecules, water_uncertainty
    raise NotImplemented("getntargets: %s is not implemented yet" % material)
        
def safe_getntargets(ntargets_material):
    if ntargets_material == "tracker":
        PlotUtils.TargetUtils.Get().SetDistToDivCut(0)
        out = PlotUtils.TargetUtils.Get().GetTrackerNAtoms(5890.0, 8467.0, True, 850.0)
        temp, temp_err = getntargets(ntargets_material, True, True) 
        temp /= 13.0
        print out, temp
        out = temp
        out = 2.5162271154514E+029
        # That's the function we should use.
        out = PlotUtils.TargetUtils.Get().GetTrackerNCarbonAtoms(5890.0, 8467.0, True, 850.0)
        return out, 0.005
    
    ntargets, ntargets_err = getntargets(ntargets_material, True, True)
    # N nucleons but I want n atoms
    if "lead" in ntargets_material:
        ntargets /= 207.0
    if "iron" in ntargets_material:
        ntargets /= 56.0
    if "water" in ntargets_material:
        ntargets /= 18.0
    if "carbon" in ntargets_material or "tracker" in ntargets_material:
        ntargets /= 12.0
        
    return ntargets, ntargets_err

def main(output_dir, input_filename, histname, pot, material):
    is2d = "_vs" in output_dir
    
    c = ROOT.TCanvas('c', 'c', 800, 800)
    
    if True:
        print "Working on", material
        filename = input_filename
        fullfilename = input_filename
        tf = ROOT.TFile(fullfilename)
        assert not tf.IsZombie(), "Wasn't able to open the file %s" % fullfilename

        hist = tf.Get(histname) 
        assert hist != None, "Wasn't able to load the histogram from %s" % fullfilename
        
        hist_2d = None
        if is2d:
            print "Getting projection X"
            hist_2d = hist
            hist = hist_2d.ProjectionX("hist_projection_" + material)
            print "Done getting projection X"
        
        fluxhist = make_clean_copy(hist)
        N_flux_universes = 500
        flux_universe_name = "Flux"
        fluxhist.AddVertErrorBand(flux_universe_name, N_flux_universes)
        xsechist = make_clean_copy(hist)
        
        splinename = getsplinename(material)
        spline = getspline(splinename)
        
        # Now we have the CC rate hist and the spline
        units = 1e-38 # converts to 1e-38cm
        units *= 1 / float(100 * 100)
        
        potnorm = pot
            
        bins = range(1, hist.GetNbinsX() + 1)
        for b in bins: 
            binlowedge = hist.GetBinLowEdge(b)
            binwidth = hist.GetBinWidth(b)
            
            # Splines go from 0 to 100 and there are 1000 bins. So 10x the binlowedge
            #spline_low = int(binlowedge * 10)
            #spline_high = int((binlowedge + binwidth) * 10)
            spline_low = binlowedge
            spline_high = binlowedge + binwidth
            xsecInt = spline.Integral(spline_low, spline_high)
            xsec = xsecInt / binwidth * units
            #print "For {:0.1f} to {:0.1f}, xsec was {:0.1f} which is {:0.1f}. Spline int from {} to {}".format(binlowedge, binlowedge + binwidth, xsecInt * 1000, xsec, spline_low, spline_high)
            
            ntargets_material = material
            ntargets, ntargets_err = safe_getntargets(ntargets_material)
            
            rate = hist.GetBinContent(b)
            rate_err = hist.GetBinError(b)
            
            flux = rate/(xsec * ntargets * potnorm)
            flux_err = math.sqrt(rate_err**2/(xsec * ntargets * potnorm)**2 + ntargets_err**2 / ntargets**2 * rate**2/(xsec * ntargets * potnorm)**2)
            
            fluxhist.SetBinContent(b, flux)
            fluxhist.SetBinError(b, flux_err)
            
            if DO_UNIVERSES:
                vert = fluxhist.GetVertErrorBand(flux_universe_name)
                assert vert != None, "Wasn't able to get flux universe from flux hist"
                vert_hist = hist.GetVertErrorBand(flux_universe_name)
                assert vert_hist != None, "Wasn't able to get flux universe from rate hist"
                for i in range(N_flux_universes):
                    rate_hist = vert_hist.GetHist(i)
                    assert rate_hist != None, "Wasn't able to get the %sth rate hist" % i
                    rate_i = rate_hist.GetBinContent(b)
                    rate_err_i = rate_hist.GetBinError(b)
                    flux_i = rate_i/(xsec * ntargets * potnorm)
                    flux_err_i = math.sqrt(rate_err_i**2/(xsec * ntargets * potnorm)**2 + ntargets_err**2 / ntargets**2 * rate_i**2/(xsec * ntargets * potnorm)**2)
                    vert.GetHist(i).SetBinContent(b, flux_i)
                    vert.GetHist(i).SetBinError(b, flux_err_i)
                
            
            print b, binlowedge, binwidth, xsecInt, xsec, rate, flux, ntargets
            
            xsechist.SetBinContent(b, xsec)
            xsechist.SetBinError(b, 0)
            
        if is2d:
            flux_2d = make_clean_copy(hist_2d)
            ybins = range(1, flux_2d.GetNbinsY() + 1)
            # Loop over each flux bin and scale by flux/hist
            for b in bins:
                hist_content = hist.GetBinContent(b)
                hist_error = hist.GetBinError(b)
                flux_content = fluxhist.GetBinContent(b)
                flux_error = fluxhist.GetBinError(b)
                if hist_content > 0:
                    scale_factor = flux_content / hist_content
                else: scale_factor = 1
                #print "bin {}, hist content = {} +/- {}, flux content = {} +/- {}, scale factor {}".format(b, hist_content, hist_error, flux_content, flux_error, scale_factor)
                for by in ybins:
                    global_bin = flux_2d.GetBin(b, by)
                    current_content = hist_2d.GetBinContent(global_bin)
                    current_error = hist_2d.GetBinError(global_bin)
                    new_content = current_content * scale_factor
                    # TODO this scale factor is only good for content, not error
                    new_error = current_error * scale_factor
                    #print "For {} {}, new content = {} +/- {}, scale factor {}".format(b, by, new_content, new_error, scale_factor)
                    flux_2d.SetBinContent(global_bin, new_content)
                    flux_2d.SetBinError(global_bin, new_error)
                    # TODO add vert error universes
                    
                
            
        hist.Scale(1, "width")
        hist.GetXaxis().SetRangeUser(0, 20)
        hist.Draw("e")
        c.Print(os.path.join(output_dir, "%s_rate.png" % material))
        fluxhist.Scale(1, "width")
        fluxhist.GetXaxis().SetRangeUser(0, 20)
        fluxhist.Draw("e")
        c.Print(os.path.join(output_dir, "%s_flux.png" % material))
        #xsechist.Scale(1, "width")
        xsechist.GetXaxis().SetRangeUser(0, 20)
        xsechist.Draw("hist")
        c.Print(os.path.join(output_dir, "%s_xsec.png" % material))
        
        
        if is2d:
            flux_2d.GetXaxis().SetRangeUser(0, 20)
            flux_2d.Draw("colz")
            c.Print(os.path.join(output_dir, "%s_flux_2d.png" % material))
        
        tf2 = ROOT.TFile(os.path.join(output_dir, "flux_%s.root" % material), "recreate")
        tf2.cd()
        hist.Write("rate")
        fluxhist.Write("flux")
        xsechist.Write("xsec")
        if is2d:
            flux_2d.Write("flux_2d")
        tf2.Close()
        tf3 = ROOT.TFile(os.path.join(output_dir, "simple_flux_%s.root" % material), "recreate")
        tf3.cd()
        hist.GetCVHistoWithStatError().Write("rate")
        fluxhist.GetCVHistoWithStatError().Write("cv_flux")
        fluxhist.GetCVHistoWithError().Write("cv_flux_with_sys_error")
        xsechist.GetCVHistoWithStatError().Write("xsec")
        tf3.Close()
        print "Done working on", material
        print "\n\n\n\n"
   
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1: 
        print "python flux_calculation.py output_directory input_filename histname pot material"
        print "histname: name of the histogram inside input_filename, which comprises of only events on specified material"
        print """materials:
tracker (assumes specific Z constraints, which you can edit in the script)\niron/lead (assumes T2-5)\niron_t15/lead_t15 (assumes T1-5)\ncarbon
water (assumes no x/y cuts), water_apothem (assumes 850mm apothem cut), water_radius (assumes 900mm radius cut)"""
        exit()
    output_directory = sys.argv[1]
    input_filename = sys.argv[2]
    histname = sys.argv[3]
    pot = float(sys.argv[4])
    material = sys.argv[5]
    main(output_directory, input_filename, histname, pot, material)  























