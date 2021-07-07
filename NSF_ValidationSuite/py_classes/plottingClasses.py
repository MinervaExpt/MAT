import ROOT,os

class makeEnv_TCanvas(object):

  def __init__(self,plotName,logy=False):
    self.plotName = plotName
    self.logy = logy

  def __enter__(self):
    self.canvas = ROOT.TCanvas( "canvas" , "canvas" , 10 , 10 , 1000, 750 )
    if self.logy: self.canvas.SetLogy(0)
    return self

  def __exit__(self,*exc):
    # If directory for plotName doesn't exist yet, make it
    plotNameComponents = self.plotName.split('/')
    plotDir = '/'.join(plotNameComponents[:-1]) 
    if not os.path.isdir(plotDir):
      print "Making plot directory {0}".format(plotDir)
      os.system( "mkdir %s" % plotDir )
    self.canvas.Print(self.plotName)
    del self.canvas

def localDrawErrorSummary( plotter , hist , errorBand="" ):

  plotter.axis_maximum = 0.20
  plotter.DrawErrorSummary(hist,"TR",True,True,0.00001,False,errorBand,True,"",True)

