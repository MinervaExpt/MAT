from ROOT import *
from PlotUtils import *
import sys
from array import array
import getopt

def getChi2StatUni(infile, iteration, n_universes = 1):
  #Get the top universe
  chi2_stat_limit = infile.Get("Chi2_Iteration_Dists/h_chi2_modelData_trueData_iter_stat")
  max_stat_uni = []
  max_stat_chi2 = 0

  iterBin = chi2_stat_limit.GetXaxis().FindBin(iteration)

  for iBinY in range(1,chi2_stat_limit.GetNbinsY()+1):
    if chi2_stat_limit.GetBinContent(iterBin,iBinY) > max_stat_chi2:
      max_stat_chi2 = chi2_stat_limit.GetBinContent(iterBin,iBinY)
      max_stat_uni.insert(0,iBinY-1)
      if len(max_stat_uni)>n_universes:
        max_stat_uni.pop(n_universes) 
  return max_stat_uni

def main(input_file, lower_limit):
  f = TFile(input_file)

  truth_hist = f.Get("Input_Hists/h_data_truth")
  #Want to get the overflow/underflow bins
  #Everything inherits from TH2
  nOverflowX = 1
  nOverflowY = int( truth_hist.IsA().InheritsFrom("TH2") ) 
  nOverflowZ = int( truth_hist.IsA().InheritsFrom("TH3") )
  
  exclude_bins = {}
  #Different iterations?
  iterations = [1]
  for iteration in iterations:
    #Possibly expand this for many different universes, for now stick with 1
    max_stat_uni = getChi2StatUni(f, iteration)
    for stat_uni in max_stat_uni:
      unfold_hist = f.Get("Unfolded_Data/Stat_{0}_Iter_{1}".format(stat_uni,iteration))

      #Get the lower limit bins
      nBinsZ = unfold_hist.GetNbinsZ()
      nBinsY = unfold_hist.GetNbinsY()
      nBinsX = unfold_hist.GetNbinsX()
      for iBinZ in range(nBinsZ+2*nOverflowZ):
        for iBinY in range(nBinsY+2*nOverflowY):
          for iBinX in range(nBinsX+2*nOverflowX):
            if iBinX == nOverflowX-1 or iBinY == nOverflowY-1 or iBinZ == nOverflowZ-1:
              continue
            if iBinX == nBinsX+1 or iBinY == nBinsY+1 or iBinZ == nBinsZ+1:
              continue

            iGlobalBin = unfold_hist.GetBin(iBinX,iBinY,iBinZ)

            nEntryStat  = unfold_hist.GetBinContent(iGlobalBin)
            nEntryTruth = truth_hist.GetBinContent(iGlobalBin)

            #Bins that are very both small are supressed
            if nEntryStat <= 1e-6 and nEntryTruth <= 1e-6:
              continue
            if nEntryStat <= 300 or nEntryTruth <= 300:
              #Binning without under/overflow
              iBin       = ((iBinZ-nOverflowZ)*nBinsY+(iBinY-nOverflowY))*nBinsX + (iBinX - nOverflowX)
              exclude_bins[iBin] = iGlobalBin
              #print " {0:>3} {1:>4} {2:>4} {3:>4} {4:>4} {5:6.2f} {6:6.2f}".format(iBin, iGlobalBin, iBinX, iBinY, iBinZ, nEntryStat, nEntryTruth)

      #Now, let's get the chi2 map 
      #chi2maps are using no under/overflow. First seen bin is 0
      chi2map = TH2D(f.Get("Chi2Maps/Chi2Map_Stat_{0}_Iter_{1}".format(stat_uni,iteration)))
      chi2map.SetName("Chi2Map{0}".format(iteration))
     
      print "Bins to exclude (Stat {0} Iter {1}): {2}".format(stat_uni,iteration, ",".join([str(_) for _ in sorted(exclude_bins.values())]) )      


      tmp_exclude_bins = [] 
      #Now, let's look at each row
      #Modifier: look at high chi2 bins?
      #Let's calculate the absolute chi2 of each row then, look at the components that make up the chi2
      # if a bin provides ~0.25 of the absolute chi2, AND it's an excluded bin, let's remove it from exclusion
      # If the component comes from the bin itself, skip over
      for iBin in range(1,chi2map.GetNbinsX()+1):
        #Calculate absolute chi2
        abs_chi = 0
        for jBin in range(1,chi2map.GetNbinsY()+1):
          abs_chi += abs(chi2map.GetBinContent(iBin,jBin))

        #check the percentage
        for jBin in range(1,chi2map.GetNbinsY()+1):
          if iBin == jBin:
            continue
          if abs_chi<1e-6:
            continue

          if chi2map.GetBinContent(iBin,jBin)/abs_chi >= 0.25:
            if jBin-1 in exclude_bins.keys() and jBin-1 not in tmp_exclude_bins:
              tmp_exclude_bins.append(jBin-1)          

  for b in tmp_exclude_bins:
    exclude_bins.pop(b)

  #need to return the number of excluded bins
  #These need to be translated to over/underflow bins
  print "Excluded bins (no overflow) {0}".format( ",".join([str(_) for _ in sorted(exclude_bins.keys())]) )
  print " "
  print "(Use this)         --exclude_bins {0}".format( ",".join([str(_) for _ in sorted(exclude_bins.values())]) )

if __name__ == "__main__":
  helptext='FindExcludeChi2Bins.py [-i <input file name> -L <lower limit>]\n'\
           'Find some excluded bins \n'\
           'Required Arguments\n'\
           '                 -i --input         Input File Name\n'\
           '                 -L --lower_limit   Lower limit\n'\
           'Optional Arguments\n'\
           '                 -h --help       Displays this message\n'

  try:
    opts,args=getopt.getopt(sys.argv[1:],"i:L:h",["input=","lower_limit=","help"])
  except getopt.GetoptError:
    print helptext
    sys.exit(1)

  gROOT.SetBatch(True)

  input_name  = ""
  lower_limit = 0

  for opt, arg in opts:
    if opt in ("-h","--help"):
      print helptext
      sys.exit(0)
    if opt in ("-i","--input"):
      input_name = arg
    if opt in ("-L","--lower_limit"):
      lower_limit = arg

  main(input_name,lower_limit)
