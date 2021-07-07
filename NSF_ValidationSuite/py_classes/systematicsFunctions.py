import PlotUtils,sys
import math

class CVUniverse(PlotUtils.DefaultCVUniverse,object):

  def __init__(self,chain):

    PlotUtils.DefaultCVUniverse.__init__(self,chain)

  def getWeight(self):

    ## # MINOS
    wgt_MINOS = self.GetMinosEfficiencyWeight()

    # GENIE
    wgt_GENIE = self.GetGenieWeight()
 
    # Flux
    wgt_flux = self.GetFluxAndCVWeight()

    # 2p2h, RPA
    wgt_2p2h = self.GetLowRecoil2p2hWeight()
    #wgt_2p2h = 1.
    wgt_RPA = self.GetRPAWeight()

    return wgt_MINOS * wgt_GENIE * wgt_flux * wgt_2p2h * wgt_RPA 
 
  def GetPmuTransverse(self):
    Pmu = self.GetPmu()
    Thetamu = self.GetThetamu()
    return Pmu * math.sin(Thetamu)

  def GetRecoMuonEKinetic(self):

    return self.GetEmu()-105.6583

  def GetRecoilEnergyCCQENu(self):
    Ehad = self.GetVecElem("recoil_summed_energy",0)

    return Ehad

  def GetRecoilEnergyNukeCC(self):
    Ehad = self.GetDouble("NukeCC_nu_energy_recoil")

    return Ehad

  def GetRecoilEnergyLocal(self,toolName):
    if toolName == "CCQENu": return self.GetRecoilEnergyCCQENu()
    elif toolName == "NukeCC": return self.GetRecoilEnergyNukeCC()
    else:
      print "I don't know about a recoil energy definition for your chosen toolName, so I'm going to exit gracefully"
      sys.exit()

# Explicitly tell Flux universe that it inherits from CVUniverse
class FluxUniverse(PlotUtils.FluxUniverse(PlotUtils.DefaultCVUniverse), CVUniverse, object):

  def __init__(self,chain,nsigma,universe_number):
    PlotUtils.FluxUniverse(PlotUtils.DefaultCVUniverse).__init__(self,chain,nsigma,universe_number)

  pass

# Explicitly tell GENIE universe that it inherits from CVUniverse
class GenieUniverse(PlotUtils.GenieUniverse(PlotUtils.DefaultCVUniverse), CVUniverse, object):

  def __init__(self,chain,nsigma,universe_name):
    PlotUtils.GenieUniverse(PlotUtils.DefaultCVUniverse).__init__(self,chain,nsigma,universe_name)

  pass

# Explicitly tell GENIE Rvx1pi universe that it inherits from CVUniverse
class GenieRvx1piUniverse(PlotUtils.GenieRvx1piUniverse(PlotUtils.DefaultCVUniverse), CVUniverse, object):

  def __init__(self,chain,nsigma,universe_name):
    PlotUtils.GenieRvx1piUniverse(PlotUtils.DefaultCVUniverse).__init__(self,chain,nsigma,universe_name)

  pass

# Explicitly tell GENIE FaCCQE universe that it inherits from CVUniverse
class GenieFaCCQEUniverse(PlotUtils.GenieFaCCQEUniverse(PlotUtils.DefaultCVUniverse), CVUniverse, object):

  def __init__(self,chain,nsigma,universe_number):
    PlotUtils.GenieFaCCQEUniverse(PlotUtils.DefaultCVUniverse).__init__(self,chain,nsigma,universe_number)

  pass

# Explicitly tell 2p2h universe that it inherits from CVUniverse
class Universe2p2h(PlotUtils.Universe2p2h(PlotUtils.DefaultCVUniverse), CVUniverse, object):

  def __init__(self,chain,nsigma,universe_number):
    PlotUtils.Universe2p2h(PlotUtils.DefaultCVUniverse).__init__(self,chain,nsigma,universe_number)

  pass

# Explicitly tell RPA universe that it inherits from CVUniverse
class RPAUniverse(PlotUtils.RPAUniverse(PlotUtils.DefaultCVUniverse), CVUniverse, object):

  def __init__(self,chain,nsigma,universe_number,q2_region):
    PlotUtils.RPAUniverse(PlotUtils.DefaultCVUniverse).__init__(self,chain,nsigma,universe_number,q2_region)

  pass

# Explicitly tell MINOS Efficiency universe that it inherits from CVUniverse
class MinosEfficiencyUniverse(PlotUtils.MinosEfficiencyUniverse(PlotUtils.DefaultCVUniverse), CVUniverse, object):

  def __init__(self,chain,nsigma):
    PlotUtils.MinosEfficiencyUniverse(PlotUtils.DefaultCVUniverse).__init__(self,chain,nsigma)

  pass

# Explicitly tell Minerva (Energy Scale) universe that it inherits from CVUniverse
class MuonUniverseMinerva(PlotUtils.MuonUniverseMinerva(PlotUtils.DefaultCVUniverse), CVUniverse, object):

  def __init__(self,chain,nsigma):
    PlotUtils.MuonUniverseMinerva(PlotUtils.DefaultCVUniverse).__init__(self,chain,nsigma)

  pass

# Explicitly tell MINOS (Energy Scale) universe that it inherits from CVUniverse
class MuonUniverseMinos(PlotUtils.MuonUniverseMinos(PlotUtils.DefaultCVUniverse), CVUniverse, object):

  def __init__(self,chain,nsigma):
    PlotUtils.MuonUniverseMinos(PlotUtils.DefaultCVUniverse).__init__(self,chain,nsigma)

  pass

# Explicitly tell muon resolution universe that it inherits from CVUniverse
class MuonResolutionUniverse(PlotUtils.MuonResolutionUniverse(PlotUtils.DefaultCVUniverse), CVUniverse, object):

  def __init__(self,chain,nsigma):
    PlotUtils.MuonResolutionUniverse(PlotUtils.DefaultCVUniverse).__init__(self,chain,nsigma)

  pass

# Explicitly tell Beam Angle universes that they inherit from CVUniverse
class BeamAngleXUniverse(PlotUtils.BeamAngleXUniverse(PlotUtils.DefaultCVUniverse), CVUniverse, object):

  def __init__(self,chain,nsigma):
    PlotUtils.BeamAngleXUniverse(PlotUtils.DefaultCVUniverse).__init__(self,chain,nsigma)

  pass

class BeamAngleYUniverse(PlotUtils.BeamAngleYUniverse(PlotUtils.DefaultCVUniverse), CVUniverse, object):

  def __init__(self,chain,nsigma):
    PlotUtils.BeamAngleYUniverse(PlotUtils.DefaultCVUniverse).__init__(self,chain,nsigma)

  pass

# Explicitly tell Response universes that they inherit from CVUniverse
class ResponseUniverse(PlotUtils.ResponseUniverse(PlotUtils.DefaultCVUniverse), CVUniverse, object):

  def __init__(self,chain,nsigma,universe_name):
    PlotUtils.ResponseUniverse(PlotUtils.DefaultCVUniverse).__init__(self,chain,nsigma,universe_name) 

  # Overwrite GetRecoilEnergyCCQENu() method of the base class to instead get the recoil energy corresponding to this systematic universe variation
  def GetRecoilEnergyCCQENu(self):
    shift_val = super(ResponseUniverse,self).GetRecoilShift()
    Ehad = self.GetVecElem("recoil_summed_energy",0)

    return shift_val+Ehad
  
