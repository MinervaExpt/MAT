import os
import sys
import math
import glob
import ROOT
ROOT.gROOT.SetBatch(True)
import PlotUtils

# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)

c = ROOT.TCanvas("c", "c", 800, 800)



N_TARGET_UNCERTAINTY = 0.005 # most targets have 0.005% uncertainty except water
N_TARGET_STARTING_Z = 5890.0
N_TARGET_ENDING_Z = 8467.0
N_TARGET_APOTHEM = 850.0
N_TARGET_N_NUCLEONS = True

def copy(hist, isROOT = False):
    #print "copy: hist n vert errors bars", hist.GetNVertErrorBands()
    #print "copy hist n lat errors bars", hist.GetNLatErrorBands()
    #new = type(hist)(hist)
    #new.TransferErrorBands(hist, False)
    new = make_clean_copy(hist, not isROOT)
    new.Add(hist)
    #print "copy: new n vert errors bars", new.GetNVertErrorBands()
    #print "copy: new n lat errors bars", new.GetNLatErrorBands()
    return new
    
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



def unit_test_simple_average():
    filename = "xsec/XSec_2021-05-25_test_daisy/tki_muon_p_flux_daisy/default_flux/tracker_daisy_*/out.root"
    
    files = glob.glob(filename)
    assert len(files) > 0
    
    hist = None
    hists = []
    hist_name = "xsec_no_width_scale_mc"
    hist_name = "xsec_data"
    scale = 1 / float(len(files))
    for f in files:
        tf = ROOT.TFile(f)
        assert not tf.IsZombie(), "Couldn't load %s" % f
        temp = tf.Get(hist_name)
        assert temp != None, "Couldn't get hist %s" % hist_name
        # Temp fix because N targets is too large by a factor 12
        # TODO remove this
        temp.Scale(12.0)
        if hist == None:
            hist = make_clean_copy(temp)
        hist.Add(temp, scale)
        hists.append(temp)
    hist.Draw("e1")
    for h in hists:
        h.SetLineColor(ROOT.kRed)
        h.Draw("hist same")
    filename2 = "xsec/XSec_2021-05-25_test_daisy/tki_muon_p_flux_daisy/default_flux/tracker/out.root"
    tf2 = ROOT.TFile(filename2)
    assert not tf2.IsZombie()
    hist_original = tf2.Get(hist_name)
    assert hist_original != None
    #hist_original.Scale(1.1)
    hist_original.SetLineColor(ROOT.kBlue)
    hist_original.Draw("e1 same")
    head = hist.GetMaximum() * 1.2
    head = max(head, hist_original.GetMaximum() * 1.2)
    hist.GetYaxis().SetRangeUser(0, head)
    hist.Draw("e1 same")
    hist_original.Draw("e same")
    c.Print("test/unit_test_simple_average.png")
    
    
import re
def unit_test_carbon():
    filename = "xsec/XSec_2021-05-25_test_daisy/tki_muon_p_flux_daisy/default_flux/tracker_daisy_*/out.root"
    
    files = glob.glob(filename)
    assert len(files) > 0
    
    files2 = []
    match = r"tracker_daisy_(\d\d)"
    for f in files:
        number = re.search(match, f).group(1)
        number = int(number)
        if number <= 5 or number == 11:
            files2.append(f)
    files = files2
    
    hist = None
    hists = []
    hist_name = "xsec_no_width_scale_mc"
    hist_name = "xsec_data"
    scale = 1 / float(len(files))
    for f in files:
        tf = ROOT.TFile(f)
        assert not tf.IsZombie(), "Couldn't load %s" % f
        temp = tf.Get(hist_name)
        assert temp != None, "Couldn't get hist %s" % hist_name
        # Temp fix because N targets is too large by a factor 12
        # TODO remove this
        temp.Scale(12.0)
        if hist == None:
            hist = make_clean_copy(temp)
        hist.Add(temp, scale)
        hists.append(temp)
    hist.Draw("e1")
    for h in hists:
        h.SetLineColor(ROOT.kRed)
        h.Draw("hist same")
    filename2 = "xsec/XSec_2021-05-25_test_daisy/tki_muon_p_flux_daisy/default_flux/carbon/out.root"
    filename2 = "xsec/XSec_Final_2021-04-25/tki_muon_p/default_flux/carbon/out.root"
    tf2 = ROOT.TFile(filename2)
    assert not tf2.IsZombie()
    hist_original = tf2.Get(hist_name)
    assert hist_original != None
    #hist_original.Scale(1.1)
    hist_original.SetLineColor(ROOT.kBlue)
    hist_original.Draw("e1 same")
    
    
    filename3 = "xsec/XSec_2021-05-25_test_daisy/tki_muon_p_flux_daisy/default_flux/tracker/out.root"
    tf3 = ROOT.TFile(filename3)
    assert not tf3.IsZombie()
    tracker_original = tf3.Get(hist_name)
    assert hist_original != None
    #hist_original.Scale(1.1)
    tracker_original.SetLineColor(ROOT.kMagenta)
    tracker_original.Draw("e1 same")
    
    
    head = hist.GetMaximum() * 1.2
    head = max(head, hist_original.GetMaximum() * 1.2)
    hist.GetYaxis().SetRangeUser(0, head)
    hist.Draw("e1 same")
    hist_original.Draw("e1 same")
    hist_original_errors = hist_original.GetCVHistoWithError(True, False)
    hist_original_errors.Draw("e1 same")
    c.Print("test/unit_test_carbon.png")
    
def add_hists(params, hists, out = None):
    if out == None:
        out = make_clean_copy(hists[0])
    nbins = params.GetNbinsX()
    for b in range(1, nbins + 1):
        par = params.GetBinContent(b)
        i = b - 1
        hist = hists[i]
        out.Add(hist, par)
    return out
    
def add_hists_with_param_hist(param_hist, hists):
    print "Adding CV hist"
    out = add_hists(param_hist, hists)
    for vert_name in param_hist.GetVertErrorBandNames():
        print "Working on", vert_name
        error_band = param_hist.GetVertErrorBand(vert_name)
        n_hists = error_band.GetNHists()
        out.AddVertErrorBand(vert_name, n_hists)
        out_error_band = out.GetVertErrorBand(vert_name)
        for i in range(n_hists):
            error_band_hist = error_band.GetHist(i)
            out_error_band_hist = out_error_band.GetHist(i)
            add_hists(error_band_hist, hists, out_error_band_hist)
    return out
    
def get_and_add_hists(param_hist, filename, hist_name):
    
    files = glob.glob(filename)
    assert len(files) > 0
    
    hists = []
    scale = 1 / float(len(files))
    for f in files:
        tf = ROOT.TFile(f)
        assert not tf.IsZombie(), "Couldn't load %s" % f
        temp = tf.Get(hist_name)
        assert temp != None, "Couldn't get hist %s" % hist_name
        # Temp fix because N targets is too large by a factor 12
        # TODO remove this
        #temp.Scale(12.0)
        hists.append(temp)
        
    added = add_hists_with_param_hist(param_hist, hists)
    return added
    
def unit_test_add_hists_with_param_hist():
    param_filename = "flux_prism_2d_output/lambda_100/unit_test_simple_carbon_parameters.root"
    project_name = "lambda_100"
    #param_filename = "flux_prism_2d_output/no_weight/unit_test_simple_carbon_parameters.root"
    #project_name = "no_weight"
    #param_filename = "flux_prism_2d_output/simple_carbon/unit_test_simple_carbon_parameters.root"
    #project_name = "simple_carbon"

    tf2 = ROOT.TFile(param_filename)
    assert not tf2.IsZombie()
    param_hist = tf2.Get("param_hist")
    
    filename = "xsec/XSec_2021-06-02_test_daisy/muon_p_flux_daisy/default_flux/tracker_daisy_*/out.root"
    output_directory = "xsec/XSec_2021-06-02_test_daisy/muon_p_flux_daisy/default_flux/tracker_combined_{}".format(project_name)
    filename = "xsec/XSec_2021-07-14_daisy/tki_alpha_flux_daisy/default_flux/tracker_daisy_*/out.root"
    output_directory = "xsec/XSec_2021-07-14_daisy/tki_alpha_flux_daisy/default_flux/tracker_combined_{}".format(project_name)
    
    hist_names = ["xsec_mc", "xsec_data", "xsec_mc_ccqelike_and_ccqe", "xsec_mc_ccqelike_and_type2p2h", "xsec_mc_ccqelike_and_res", "xsec_mc_ccqelike_and_dis", "xsec_mc_ccqelike"]
    out = dict()
    for hist_name in hist_names:
        hist = get_and_add_hists(param_hist, filename, hist_name)
        out[hist_name] = hist
    
    
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    tf = ROOT.TFile(os.path.join(output_directory, "out.root"), "recreate")
    assert not tf.IsZombie()
    tf.cd()
    for hist_name in hist_names:
        hist = out[hist_name]
        hist.Write(hist_name)
    tf.Close()
    

#unit_test_add_hists_with_param_hist()
    
    
#unit_test_carbon()

def getflux(hist, material, project_dir):
    filename_base = "../data/flux_daisy/%s/flux_with_errors/flux_%s.root"
    filename_target = filename_base % (project_dir, material)
    filename_tracker = filename_base % (project_dir, "tracker")
    
    tf = ROOT.TFile(filename_target)
    if tf.IsZombie():
        raise ValueError("Cannot resolve flux for material '%s'" % material)
    flux_target = tf.Get("flux")
    assert flux_target != None, "Wasn't able to load histogram 'flux' from %s" % filename_target
    
    
    tf = ROOT.TFile(filename_tracker)
    if tf.IsZombie():
        raise ValueError("Cannot resolve flux for material '%s'" % material)
    flux_tracker = tf.Get("flux")
    assert flux_target != None, "Wasn't able to load histogram 'flux' from %s" % filename_tracker
    
    #  (int nuPDG, bool applyNuEConstraint, enum EPlaylist playlist, enum EFluxVersion fluxVersion, enum EG4NumiVersion g4NumiVersion, int nUniverses=200)
    print "Warning: Getting flux with the following settings: PlotUtils.FluxReweighter(14, True, 19, 2, 1, 500)"
    print "Please adjust the code to get the specific flux you want."
    frw = PlotUtils.FluxReweighter(14, True, 19, 2, 1, 500)
    flux = frw.GetFluxReweighted(14);
    
    # Make a copy of the flux histogram, but fill with ratio of target / tracker
    ratio_hist = make_clean_copy(flux, False) # Want same bins
    nbins = ratio_hist.GetXaxis().GetNbins()
    for b in range(0, nbins + 2):
        bin_min_x = ratio_hist.GetBinLowEdge(b)
        bin_center_x = ratio_hist.GetBinCenter(b)
        bin_max_x = ratio_hist.GetBinLowEdge(b+1)
        
        bin_tracker = flux_tracker.FindBin(bin_center_x)
        content_tracker = flux_tracker.GetBinContent(bin_tracker)
        bin_target = flux_target.FindBin(bin_center_x)
        content_target = flux_tracker.GetBinContent(bin_target)
        
        ratio = 1.0
        if content_tracker > 0:
            ratio = content_target / content_tracker
        ratio_hist.SetBinContent(b, ratio)
        ratio_hist.SetBinError(b, 0)
        
    # Now fill all universes with the same ratio of target / tracker
    ratio_hist.AddMissingErrorBandsAndFillWithCV(flux)
    # Now make the correct flux by multipling the ratio hist
    flux.Multiply(flux, ratio_hist)
    
    out, foo = getconstflux_fromhist(hist, flux)
    
    return out, foo
        
    
    

def getflux_old(hist, material, project_dir):
    filename = "../data/flux_daisy/%s/flux_with_errors/flux_%s.root" % (project_dir, material)
    tf = ROOT.TFile(filename)
    if tf.IsZombie():
        raise ValueError("Cannot resolve flux for material '%s'" % material)
    f = tf.Get("flux")
    vert = f.GetVertErrorBand("Flux")
    # Remove sys universes
    for name in f.GetVertErrorBandNames():
        if "flux" in name.lower():
            print "Found a flux universe:", name
        else:
            print "Removing:", name
            f.PopVertErrorBand(name)
    for name in f.GetLatErrorBandNames():
        if "flux" in name.lower():
            print "Found a flux universe:", name
        else:
            print "Removing:", name
            f.PopLatErrorBand(name)
    f.AddMissingErrorBandsAndFillWithCV(hist)
    
    vert = f.GetVertErrorBand("Flux")
    # Takes care of flux universes
    # Also converts from m^-2 to cm^-2
    # h_flux_integrated *= 1.0e-4 
    out, foo = getconstflux_fromhist(hist, f)
    vert = out.GetVertErrorBand("Flux")
    return out, foo
    
def getconstflux_fromhist(hist, h_flux):
    
    # Get clean histogram without systematic universes
    h_flux_integrated = make_clean_copy(hist, False)
    Emin = 0
    Emax = 1000.0
    
    b_min = h_flux.FindBin( Emin );
    b_max = h_flux.FindBin( Emax );
    flux_cv = h_flux.Integral( b_min, b_max , "width" );
    nbins = hist.GetNbinsX() + 2
    if hist.GetNbinsY() > 1: nbins *= hist.GetNbinsY() + 2
    for i in range(nbins):
        h_flux_integrated.SetBinContent(i,flux_cv)
        
    # Add the remaining error bands
    print list(hist.GetVertErrorBandNames())
    h_flux_integrated.AddMissingErrorBandsAndFillWithCV(hist)
    
    # Now change the error bands that are in h_flux (ie flux universes)
    names = h_flux.GetVertErrorBandNames()
    for name in names:
        if name not in h_flux_integrated.GetVertErrorBandNames(): 
            print "Skipping", name, "because it's not in output"
            print list(h_flux_integrated.GetVertErrorBandNames())
            continue
        vert = h_flux.GetVertErrorBand(name)
        print "Working on", name, "which has", vert.GetNHists(), "histograms"
        assert vert.GetNHists() >= h_flux_integrated.GetVertErrorBand(name).GetNHists(), "Flux histogram is missing n hists: " + str(h_flux_integrated.GetVertErrorBand(name).GetNHists()) + " vs " + str(vert.GetNHists())
        nhists = h_flux_integrated.GetVertErrorBand(name).GetNHists()
        for j in range(nhists):
            #print "On histogram", j
            flux_universe = vert.GetHist(j).Integral( b_min, b_max , "width" )
            vert_new = h_flux_integrated.GetVertErrorBand(name)
            for i in range(nbins):
                vert_new.GetHist(j).SetBinContent(i,flux_universe)
                
    # TODO add flux wiggle universes here
        
    # Convert from m^-2 to cm^-2
    h_flux_integrated *= 1.0e-4 
    print "CV flux is", flux_cv * 1.0e-4, "cm^-2"
    
    return h_flux_integrated, None

def get_n_targets(isMC):
    print "WARNING: N targets function assumes the following. If not true, then adjust script."
    print "Tracker Z %s-%s" % (N_TARGET_STARTING_Z, N_TARGET_ENDING_Z)
    print "Apothem:", N_TARGET_APOTHEM
    if N_TARGET_N_NUCLEONS: print "Counting N nucleons"
    else: print "Counting N neutrons"
    if N_TARGET_N_NUCLEONS:
        trackerfunc = PlotUtils.TargetUtils.Get().GetTrackerNNucleons
    else:
        trackerfunc = PlotUtils.TargetUtils.Get().GetTrackerNNeutrons
    return trackerfunc(N_TARGET_STARTING_Z, N_TARGET_ENDING_Z, isMC, N_TARGET_APOTHEM), N_TARGET_UNCERTAINTY
    
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
    
def copy(hist, isROOT = False):
    new = make_clean_copy(hist, not isROOT)
    new.Add(hist)
    return new
    
def add_hists_with_flux(hists, material, project_name, normalization):
    if material == "tracker":
        param_hist = PlotUtils.MnvH1D("tracker_param_hist", "tracker_param_hist;Flux Bin;Weight", 12, 0, 12)
        for i in range(12):
            b = i + 1
            param_hist.SetBinContent(b, 1)
            param_hist.SetBinError(b, 0)
    else:
        param_hist_filename = None
        if project_name == "targets_2345_temp":
            param_hist_filename = "../data/flux_daisy/targets_2345_temp/out_%s_000100.root" % (material)
          
        assert param_hist_filename != None, "Don't have a param hist for these targets: %s" % project_name
        tf = ROOT.TFile(param_hist_filename)
        assert not tf.IsZombie(), "Can't load param hist"
        param_hist = tf.Get("param_hist")
        assert param_hist != None, "Can't get param_hist from %s" % param_hist_filename
    
    added_together = add_hists_with_param_hist(param_hist, hists)
    added_together.Scale(normalization)
    flux, foo = getflux(added_together, material, project_name)
    added_together.Divide(added_together, flux)
    
    return added_together

if __name__ == "__main__":
    output_filename = sys.argv[1]
    material = sys.argv[2]
    other_filenames = sys.argv[3:]
    assert len(other_filenames) == 12, "I expected 12 tracker segments"
    
    n_targets_data, n_targets_data_uncertainty = get_n_targets(False)
    n_targets_mc, n_targets_mc_uncertainty = get_n_targets(True)
    
    hist_names = ["eff_corrected_data", "eff_corrected_mc", "eff_corrected_mc_ccqelike_and_ccqe", "eff_corrected_mc_ccqelike_and_type2p2h", "eff_corrected_mc_ccqelike_and_res", "eff_corrected_mc_ccqelike_and_dis", "eff_corrected_mc_ccqelike"]
    
    pot = 1.05139333022e+21 
    N_TARGET_STARTING_Z = 5890.0
    N_TARGET_ENDING_Z = 8467.0
    N_TARGET_APOTHEM = 850.0
    ntargets_mc = PlotUtils.TargetUtils.Get().GetTrackerNNucleons(N_TARGET_STARTING_Z, N_TARGET_ENDING_Z, True, N_TARGET_APOTHEM)
    ntargets_data = PlotUtils.TargetUtils.Get().GetTrackerNNucleons(N_TARGET_STARTING_Z, N_TARGET_ENDING_Z, True, N_TARGET_APOTHEM)
    normalization_mc = 1 / (pot * ntargets_mc)
    normalization_data = 1 / (pot * ntargets_data)
    project_name = "targets_2345_temp"
    
    final_hists = dict()
    for hist_name in hist_names:
        hists = []
        for f in other_filenames:
            tf = ROOT.TFile(f)
            assert not tf.IsZombie(), "Couldn't load %s" % f
            temp = tf.Get(hist_name)
            assert temp != None, "Couldn't get hist %s" % hist_name
            hists.append(temp)
        normalization = normalization_mc
        if "data" in hist_name:
            normalization = normalization_data
        xsec = add_hists_with_flux(hists, material, project_name, normalization)
        
        width_scaled = copy(xsec)
        width_scaled.Scale(1, "width")
        
        new_name = hist_name.replace("eff_corrected", "xsec") + "_no_width_scale"
        final_hists[new_name] = xsec
        new_name2 = hist_name.replace("eff_corrected", "xsec")
        final_hists[new_name2] = width_scaled
        
    # Now save the files
    base_directory = os.path.dirname(output_filename)
    if not os.path.exists(base_directory):
        os.makedirs(base_directory)
    tf = ROOT.TFile(output_filename, "recreate")
    tf.cd()
    for name in final_hists:
        hist = final_hists[name]
        hist.Write(name)
        print name, hist.Integral(), hist.Integral("width")
    tf.Close()
        
    
    
























