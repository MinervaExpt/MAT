import PlotUtils
from array import *

class HistWrapper(object):
 
  def __init__(self,*args):

    if len(args) == 3: self.initFromMnvHXD(*args)
    if len(args) == 4: self.init1DFromVec(*args)
    if len(args) == 5: self.init1D(*args)
    if len(args) == 6: self.init2DFromVecs(*args)
    if len(args) == 8: self.init2D(*args)
    if len(args) == 11: self.init3D(*args)
  
    # CV is a special case, because we'll store the MnvHXD itself here, as opposed to the various error bands of the MnvHXD  
    self.AddUniverses(args[-1])

  def initFromMnvHXD(self, mnvHXD, clear_errorband, systEvents):
    self.hist = mnvHXD.Clone()
    if clear_errorband and self.hist.GetNErrorSources()>0:
      print("HistWrapper: Clearing error bands from MnvH1D {} with existing error bands\n".format(self.hist.GetName()))
      self.hist.ClearAllErrorBands();

    self.title = self.hist.GetTitle()
    self.eventToUnivMap = {}
    self.hist.SetDirectory(0)

  def init1D(self,title,nBins,xMin,xMax,systEvents):

    self.title = title
    self.dimensionality = 1
    self.nBins = nBins
    self.xMin = xMin
    self.xMax = xMax
    #self.systematicUniverses = systEvents
    self.eventToUnivMap = {}

    self.hist = PlotUtils.MnvH1D('h_%s'%self.title,self.title,self.nBins,self.xMin,self.xMax)

  def init1DFromVec(self,title,nXBins,xBins,systEvents):
    self.title = title
    self.dimensionality = 1
    self.nXBins = nXBins
    self.xBins = xBins
    #self.systematicUniverses = systEvents
    self.eventToUnivMap = {}

    self.hist = PlotUtils.MnvH1D('h_%s'%self.title,self.title,self.nXBins,self.xBins)
    self.hist.SetDirectory(0)

  def init2D(self,title,nXBins,xMin,xMax,nYBins,yMin,yMax,systEvents):
  
    self.title = title
    self.dimensionality = 2
    self.nXBins = nXBins
    self.xMin = xMin
    self.xMax = xMax
    self.nYBins = nYBins
    self.yMin = yMin
    self.yMax = yMax
    #self.systematicUniverses = systEvents
    self.eventToUnivMap = {}

    self.hist = PlotUtils.MnvH2D('h_%s'%self.title,self.title,self.nXBins,self.xMin,self.xMax,self.nYBins,self.yMin,self.yMax)
    self.hist.SetDirectory(0)
   
  def init2DFromVecs(self,title,nXBins,xBins,nYBins,yBins,systEvents):
  
    self.title = title
    self.dimensionality = 2
    self.nXBins = nXBins
    self.xBins = xBins
    self.nYBins = nYBins
    self.yBins = yBins
    #self.systematicUniverses = systEvents
    self.eventToUnivMap = {}

    self.hist = PlotUtils.MnvH2D('h_%s'%self.title,self.title,self.nXBins,self.xBins,self.nYBins,self.yBins)
    self.hist.SetDirectory(0)

  def init3D(self,title,nXBins,xMin,xMax,nYBins,yMin,yMax,nZBins,zMin,zMax,systEvents):
    self.title = title
    self.dimensionality = 3
    self.nXBins = nXBins
    self.xMin = xMin
    self.xMax = xMax
    self.nYBins = nYBins
    self.yMin = yMin
    self.yMax = yMax
    self.nZBins = nZBins
    self.zMin = zMin
    self.zMax = zMax
    #self.systematicUniverses = systEvents
    self.eventToUnivMap = {}

    self.hist = PlotUtils.MnvH3D('h_%s'%self.title,self.title,self.nXBins,self.xMin,self.xMax,self.nYBins,self.yMin,self.yMax,self.nZBins,self.zMin,self.zMax)
    self.hist.SetDirectory(0)
 
  def SyncCVHistos(self):
    for error_band_name in self.hist.GetErrorBandNames():
      if self.dimensionality == 1:
        super(PlotUtils.MnvVertErrorBand,self.hist.GetVertErrorBand(error_band_name)).__assign__(self.hist) 
      if self.dimensionality == 2:
        super(PlotUtils.MnvVertErrorBand2D,self.hist.GetVertErrorBand(error_band_name)).__assign__(self.hist) 


  def FillUniverse(self,universe,*arg):
    self.eventToUnivMap[universe].Fill(*arg)

  #def univHist(self,systematicUniverseClass,i):
  #  return self.eventToUnivMap[systematicUniverseClass][i]   

  def AddUniverses(self,systEvents):
    for systematicUniverseClass,univVec in systEvents.iteritems():
      # CV is special
      if systematicUniverseClass.upper() == 'CV':
        self.eventToUnivMap[univVec[0]] = self.hist
        continue

      # Create error band if it doesn't exist (it shouldn't), and add corresponding entry to univMap dict
      if not self.hist.HasVertErrorBand(systematicUniverseClass):
        # The Flux gets special treatment in the case that the nueConstraint is being used
        if systematicUniverseClass == 'Flux' and univVec[0].UseNuEConstraint():
          PlotUtils.flux_reweighter(univVec[0].GetPlaylist(),univVec[0].GetAnalysisNuPDG(),True,univVec[0].GetNFluxUniverses()).AddFluxErrorBand(self.hist)
        else:
          self.hist.AddVertErrorBand(systematicUniverseClass,len(systEvents[systematicUniverseClass]))
      else:
        assert (self.hist.GetVertErrorband(systematicUniverseClass).GetNHists()>= len(systEvents[systematicUniverseClass])),"There is more {} universes than the histogram read".format(systematicUniverseClass)

      for i,systematicUniverse in enumerate(univVec):
        self.eventToUnivMap[systematicUniverse] = self.hist.GetVertErrorBand(systematicUniverseClass).GetHist(i)
