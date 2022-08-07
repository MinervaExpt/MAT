import os,sys
import PlotUtils
import ROOT
import optparse

def applyCorrection(covmatrix,histobj,factor):
    #This will add a correction of (1+1/factor) to the covariance matrix
    #This will add a correction of sqrt(1+1/factor) to the diagonal uncertainties stored in your mnvh object
    
    correction = (1+1/factor)
    sqrtcorrection = ROOT.TMath.Sqrt(correction)
    histobj.PopSysErrorMatrix("unfoldingCov")
    m_cov = covmatrix
    m_cov*=correction
    if(histobj.GetNbinsY()==1): #1D
        for i in range(0,histobj.GetNbinsX()+2):
            histobj.SetBinError(i,histobj.GetBinError(i)*sqrtcorrection)
    else: #2D
        for i in range(0,histobj.GetNbinsX()+2):
            for j in range(0,histobj.GetNbinsY()+2):
                histobj.SetBinError(i,j,histobj.GetBinError(i,j)*sqrtcorrection)
    histobj.PushCovMatrix("unfoldingCov",m_cov)
    return histobj.Clone()

def getIterations(mystring):
    myret = []
    mylist_string =mystring.split(",")
    for el in mylist_string: myret.append(int(el))
    
    return myret

    
parser = optparse.OptionParser()
parser.add_option("--n_uni",dest="n_uni",help="Number of universes")
parser.add_option("--iters",dest="iters",type="str",help="comma separate iteration list")
parser.add_option("--input",dest="input_file",help="Input file name")
parser.add_option("--uncfactor",dest="uncfactor",help="Extra unc. factor")
parser.add_option("--f_option_used_transwarp",dest="mcfrac",help="What -f option was run for the input file?")


(options,args) = parser.parse_args()


ROOT.TH1.AddDirectory(False)
plotter = PlotUtils.MnvPlotter()
f_input = ROOT.TFile(options.input_file)
n_uni = int(options.n_uni) + 1
myiters = getIterations(options.iters)
factor = float(options.uncfactor)
mcfactor = float(options.mcfrac)

outputfile = open("%s_%f_%f.txt"%(os.path.basename(options.input_file).rstrip("\n"),factor,mcfactor),"w")

for i in myiters:
    for u in range(0,n_uni):
        
        h_truth = f_input.Get("Input_Hists/h_data_truth")#Fake data truth
        h_unfolded = f_input.Get("Unfolded_Data/Stat_%d_Iter_%d"%(u,i))#unfolded 
        cov = f_input.Get("Cov_Matrix/CovMatrix_Stat_%d_Iter_%d"%(u,i))#Covariance
        print h_truth,h_unfolded,cov
        m_unfolded = applyCorrection(cov,h_unfolded,factor)
#        sys.exit()
        chi2 = plotter.Chi2DataMC(m_unfolded,h_truth,1.0,True,False,False,None)
        outputfile.write("%e\t%d\t%d\n"%(chi2,i,u))
