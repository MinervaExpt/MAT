import ROOT
ROOT.gROOT.SetBatch(True)

import flux_prism_add
import PlotUtils


# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)



def aaron():
    filename = "/minerva/data/users/abercell/hists/xsec_inputs/Merge_BkgdSub_Unfold_POTNorm_daisyAll_MinervaME1ABCDEFGLMNOP.root"
    
    tf = ROOT.TFile(filename)
    assert not tf.IsZombie(), "Cannot open %s" % filename
    
    tf.ls()
    
    
    ## Calculate the normalization
    # The xsec is the weighted sum of the eff corrected hists
    # scaled by the normalization
    isMC = False
    pot = 1.05139333022e+21 
    N_TARGET_STARTING_Z = 5890.0
    N_TARGET_ENDING_Z = 8467.0
    N_TARGET_APOTHEM = 850.0
    ntargets = PlotUtils.TargetUtils.Get().GetTrackerNNucleons(N_TARGET_STARTING_Z, N_TARGET_ENDING_Z, isMC, N_TARGET_APOTHEM)
    normalization = 1 / (pot * ntargets)
    
    # Figure out the material
    # tracker weights everything with 1. This should give you the actual tracker xsec
    # other options are carbon, water, iron and lead
    material = "carbon"
    
    # Decide on a project directory
    # This directory was tuned with only targets 2345
    project_directory = "targets_2345_temp"
    
    ## You could add additional loops here to loop over all possible variables and data/mc
    
    # Load your twelve hists however you want
    variable_name = "muon_p"
    histname_template = "h_{}_daisy{}_data_eff_corrected"
    twelve_hists = []
    for daisy_number in range(12):
        histname = histname_template.format(variable_name, daisy_number)
        hist = tf.Get(histname)
        assert hist != None, "Could not load hist '%s' from file %s" % (histname, filename)
        twelve_hists.append(hist)
        
    # This is the function that gets the xsec for you weighted by the prism weights
    # It's not width normalized.
    xsec = flux_prism_add.add_hists_with_flux(twelve_hists, material, project_directory, normalization)
    
    # Save xsec here
    
    
aaron()
    
