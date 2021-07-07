from systematicsFunctions import *
import glob

def defineSystematicUniverses(chain,dataSwitch):

  truthSwitch = True if dataSwitch == 'truth' else False

  systematicUniverses = {} 
  systematicUniverses["CV"] = []
  systematicUniverses["CV"].append( CVUniverse(chain) )

  # In the case of data, we don't assess any of the systematics on the data directly. We'll fill in the CV for all SUs and they will become meaningful as we extract the cross section/flux
  if dataSwitch == 'data': return systematicUniverses

  ## Flux Universes
  #nFluxUniverses = 1000
  nFluxUniverses = 100
  #nFluxUniverses = 2
  systematicUniverses["Flux"] = []
  for i in range(nFluxUniverses):
    systematicUniverses["Flux"].append( FluxUniverse(chain,1,i) )


  ## ##  GENIE Universes
  ## ## I switch between the below two sets of GENIE knobs
  ## ## depending I whether I'm testing or want the full suite

  ## ## GENIE_UNIVERSES = [
  ## ##   "AGKYxF1pi", 
  ## ##   "AhtBY",
  ## ## ]

  GENIE_UNIVERSES = [
    "AGKYxF1pi", 
    "AhtBY",
    "BhtBY",
    "CCQEPauliSupViaKF",
    "CV1uBY",
    "CV2uBY",
    "EtaNCEL",
    "FrAbs_N",
    "FrAbs_pi",
    "FrCEx_N",
    "FrCEx_pi",
    "FrElas_N",
    "FrElas_pi",
    "FrInel_N",
    "FrPiProd_N",
    "FrPiProd_pi",
    "MFP_N",
    "MFP_pi",
    "MaCCQEshape",
    "MaNCEL",
    "MaRES",
    "MvRES",
    "NormCCQE",
    "NormDISCC",
    "NormNCRES",
    "RDecBR1gamma",
    "Rvn1pi",
    "Rvn2pi",
    "Rvp1pi",
    "Rvp2pi",
    "Theta_Delta2Npi",
    "VecFFCCQEshape" 
  ]

  for univ in GENIE_UNIVERSES:
    systematicUniverses["GENIE_{0}".format(univ)] = []
    systematicUniverses["GENIE_{0}".format(univ)].append( GenieUniverse(chain,-1,univ) )
    systematicUniverses["GENIE_{0}".format(univ)].append( GenieUniverse(chain,+1,univ) )

  ## 2p2h, RPA Universes
  systematicUniverses["Low_Recoil_2p2h_Tune"] = []
  systematicUniverses["Low_Recoil_2p2h_Tune"].append( Universe2p2h(chain,1,1) )
  systematicUniverses["Low_Recoil_2p2h_Tune"].append( Universe2p2h(chain,1,2) )
  systematicUniverses["Low_Recoil_2p2h_Tune"].append( Universe2p2h(chain,1,3) )

  systematicUniverses["RPA_HighQ2"] = []
  systematicUniverses["RPA_HighQ2"].append( RPAUniverse(chain,1,1,"HighQ2") )
  systematicUniverses["RPA_HighQ2"].append( RPAUniverse(chain,1,2,"HighQ2") )
  systematicUniverses["RPA_LowQ2"] = []
  systematicUniverses["RPA_LowQ2"].append( RPAUniverse(chain,1,3,"LowQ2") )
  systematicUniverses["RPA_LowQ2"].append( RPAUniverse(chain,1,4,"LowQ2") )
 
  # Reco-only systematics ------------------------------------------------------------------------------------------------------
  if not truthSwitch:

    ## Muon Energy Scale Universes
    systematicUniverses["Muon_Energy_MINERvA"] = []
    systematicUniverses["Muon_Energy_MINERvA"].append( MuonUniverseMinerva(chain,-1) )
    systematicUniverses["Muon_Energy_MINERvA"].append( MuonUniverseMinerva(chain,+1) )
    systematicUniverses["Muon_Energy_MINOS"] = []
    systematicUniverses["Muon_Energy_MINOS"].append( MuonUniverseMinos(chain,-1) )
    systematicUniverses["Muon_Energy_MINOS"].append( MuonUniverseMinos(chain,+1) )

    systematicUniverses["Muon_Energy_Resolution"] = []
    systematicUniverses["Muon_Energy_Resolution"].append( MuonResolutionUniverse(chain,-1) )
    systematicUniverses["Muon_Energy_Resolution"].append( MuonResolutionUniverse(chain,+1) )

    ## Minos Efficiency Universes
    systematicUniverses["MINOS_Reconstruction_Efficiency"] = []
    systematicUniverses["MINOS_Reconstruction_Efficiency"].append( MinosEfficiencyUniverse(chain,-1) )
    systematicUniverses["MINOS_Reconstruction_Efficiency"].append( MinosEfficiencyUniverse(chain,+1) )

    ## Muon Angle Universes
    systematicUniverses["BeamAngleX"] = []
    systematicUniverses["BeamAngleX"].append( BeamAngleXUniverse(chain,-1) )
    systematicUniverses["BeamAngleX"].append( BeamAngleXUniverse(chain,+1) )
    systematicUniverses["BeamAngleY"] = []
    systematicUniverses["BeamAngleY"].append( BeamAngleYUniverse(chain,-1) )
    systematicUniverses["BeamAngleY"].append( BeamAngleYUniverse(chain,+1) )

  RESPONSE_UNIVERSES = {
    "proton": "Proton",
    "meson": "Pion",
    "em": "EM",
    "other": "Other"
  }

  for univ,Univ in RESPONSE_UNIVERSES.items():
    systematicUniverses["{0}_Response".format(Univ)] = []
    systematicUniverses["{0}_Response".format(Univ)].append( ResponseUniverse(chain,-1,univ) )
    systematicUniverses["{0}_Response".format(Univ)].append( ResponseUniverse(chain,+1,univ) )

  return systematicUniverses
 
bins_Emu = [0,1,2,3,4,5,7,9,12,15,18,22,36,50,75,100,120]
nBins_Emu = len(bins_Emu)-1

bins_nu = range(0,5050,50)
nBins_nu = len(bins_nu)-1

bins_pt = [0.0,0.075,0.15,0.25,0.325,0.4,0.475,0.55,0.7,0.85,1.0,1.25,1.5,2.5,4.5]
nBins_pt = len(bins_pt)-1
 
analysisNuPDG = 14 # muon-neutrino
nFluxUniverses = 100
useNonResPiReweight = True
useNuEScatteringConstraint = True
playlist = "minervame1L"

## By default the tuples are taken from /minerva/data/NSF_Validation/referenceTuples/
#validationSampleName = "CCQENu"
validationSampleName = "NukeCC"
refDir = "/minerva/data/NSF_Validation/referenceTuples/{0}".format(validationSampleName)
tuplePaths = glob.glob("{0}/*.root".format(refDir))

