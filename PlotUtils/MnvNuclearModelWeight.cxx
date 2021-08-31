#include "PlotUtils/MnvNuclearModelWeight.h"


#include <iostream>
#include <cstdlib>
#include <cmath>
#include <TString.h>
#include <TError.h>
#include <TSystem.h>

namespace PlotUtils
{

  //================================
  // Constructors
  //================================
  MnvNuclearModelWeight::MnvNuclearModelWeight( ):
    model_("NOMODEL"), 
    GenieSFLoc("/pnfs/minerva/persistent/users/wospakrk/NukeFiles/GENIESF.root"),
    minQ2(0.0),  
    maxQ2(100.0),
    minX(0.0),
    maxX(1.1)
  
  {
    //This is a work around until I can find out how to use TH2F copy constructors
    
    //Error("MnvNuclearModelWeight::()",GenieSFLoc.c_str());
    //TFile *file = new TFile(GenieSFLoc.c_str(), "READ");
     
      
    //Take you're stinking paws off my histograms you damn dirty ROOT!
    TH1::AddDirectory(kFALSE); // do this instead setting each directory manually
      TSystem* tSystem_ = new TSystem();
      TString* condorDirInput_ = new TString( std::getenv("CONDOR_DIR_INPUT") );
  
    /*
    GF1FeP->SetDirectory(0);
    GF2FeP->SetDirectory(0);
    GxF3FeP->SetDirectory(0);
    */
    TFile *file;
    TString *directory_str;
    if(std::getenv("CONDOR_DIR_INPUT")){
      TString basename = tSystem_->BaseName( GenieSFLoc.c_str() );
      TString filename = *condorDirInput_+ "/" + basename;
      *directory_str = *condorDirInput_;
      file = new TFile(filename,"READ");
    }  
    else{
      file = new TFile(GenieSFLoc.c_str(), "READ");
      *directory_str = std::getenv("FILES");
    }  
      if(file){
      GF1FeP  = (TH2F*)file->Get("F1p0");  
      GF2FeP  = (TH2F*)file->Get("F2p0");
      GxF3FeP  = (TH2F*)file->Get("xF3p0");
      
      file->Close();
      delete file; 
    
    }// end of if file
    
    else
      Error("MnvNuclearModelWeight::()", "Input file for GENIE SFs is NULL!");
   
  }
  
  MnvNuclearModelWeight::MnvNuclearModelWeight( const std::string& model) :
    model_(model), 
    GenieSFLoc("/pnfs/minerva/persistent/users/wospakrk/NukeFiles/GENIESF.root"),
    ModelSFLoc(GenieSFLoc),
    minQ2(0.01),
    maxQ2(29.9),
    minX(0.001),
    maxX(0.75)
  
  {
  
    //Take you're stinking paws off my histograms you damn dirty ROOT!
    TH1::AddDirectory(kFALSE); //do this instead of manually setting each directory
      TSystem * tSystem_ = new TSystem();
      TString* condorDirInput_ = new TString( std::getenv("CONDOR_DIR_INPUT") );
      TString* directory_str;
      if(std::getenv("CONDOR_DIR_INPUT")){
        directory_str = new TString(std::getenv("CONDOR_DIR_INPUT") );
      }
      else{
        directory_str = new TString(std::getenv("FILES"));
      }
      TString filename = *directory_str + "/GENIESF.root";
      TFile* fileG = new TFile(filename,"READ");
    //First Load in all the base GENIE structure Functions
    if(fileG){
      GF1FeP  = (TH2F*)fileG->Get("F1p0");  
      GF2FeP  = (TH2F*)fileG->Get("F2p0");
      GxF3FeP  = (TH2F*)fileG->Get("xF3p0");
      
      GF1FeN  = (TH2F*)fileG->Get("F1n0");
      GF2FeN  = (TH2F*)fileG->Get("F2n0");
      GxF3FeN = (TH2F*)fileG->Get("xF3n0");
      
      GF1P = (TH2F*)fileG->Get("F1p1"); 
      GF2P = (TH2F*)fileG->Get("F2p1");
      GxF3P = (TH2F*)fileG->Get("xF3p1");
      
      GF1N = (TH2F*)fileG->Get("F1n2"); 
      GF2N = (TH2F*)fileG->Get("F2n2");
      GxF3N = (TH2F*)fileG->Get("xF3n2");
  
    /*    
      //Take you're stinking paws off my histograms you damn dirty ROOT!
      //GF1FeP->SetDirectory(0);
      //GF2FeP->SetDirectory(0);
      //GxF3FeP->SetDirectory(0);
      
      //GF1FeN->SetDirectory(0);
      //GF2FeN->SetDirectory(0);
      //GxF3FeN->SetDirectory(0);
      
      //GF1P->SetDirectory(0);
      //GF2P->SetDirectory(0);
      //GxF3P->SetDirectory(0);
      
      //GF1N->SetDirectory(0);
      //GF2N->SetDirectory(0);
      //GxF3N->SetDirectory(0);    
    */  
      fileG->Close();
      delete fileG; 
    
    }// end of if file
    
    else
      Error("MnvNuclearModelWeight::()", "Input file for GENIE SFs is NULL!");
    
    
    //Now select what model you want by defining isOverallScaling, and the overall scale functions or the structure function ratios
    
    
    if(model_ == "BY13") {
       isOverallScaling = true;
       ModelSFLoc = *directory_str + "/FreeNucleonCloseout.root"; //Not Used in this case. It's a dummy
       minQ2 = 0.0;
       maxQ2 = 200.0;
       minX  = 0.0;
       maxX  = 1.5;
       
       f_Fe = TF1(Form("f_Fe_%s", model.c_str()), "1.096 - 0.38*x - 0.3*TMath::Exp(-23.0*x) + 8*TMath::Power(x,15)" ,0, 1);
       f_Pb = TF1(Form("f_Pb_%s", model.c_str()), "(0.932 + 2.461*x - 24.23*x**2 + 101.03*x**3 - 203.47*x**4 + 193.85*x**5 - 69.82*x**6)*(1.096 - 0.38*x - 0.3*TMath::Exp(-23.0*x) + 8*TMath::Power(x,15) )" ,0, 1);
       f_C = TF1(Form("f_C_%s", model.c_str()), "(1.096 - 0.38*x - 0.3*TMath::Exp(-23.0*x) + 8*TMath::Power(x,15) )/(0.919 + 1.844*x - 12.73*x**2 + 36.89*x**3 - 46.77*x**4 + 21.22*x**5)" ,0, 1);
       f_D = TF1(Form("f_D_%s", model.c_str()), "0.985*(1 + 0.422*x - 2.745*x**2 + 7.570*x**3 - 10.335*x**4 + 5.422*x**5)", 0, 1);                
    }
    
    else if(model_ == "Closure"){
      isOverallScaling = false;
      minQ2 = 0.01;
      maxQ2 = 29.9;
      minX  = 0.01;
      maxX  = 0.99;
      ModelSFLoc =*directory_str +"/Closeout.root";
  
      f_Fe = TF1(Form("f_Fe_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_Pb = TF1(Form("f_Pb_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_C =  TF1(Form("f_C_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_D =  TF1(Form("f_D_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
    }
  
    else if(model_ == "Kulagin"){
      isOverallScaling = false;
      minQ2 = 0.50;
      maxQ2 = 0.50;
      minX  = 0.001;
      maxX  = 0.958;
      ModelSFLoc =*directory_str +"/KPSF.root";
      f_Fe = TF1(Form("f_Fe_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_Pb = TF1(Form("f_Pb_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_C =  TF1(Form("f_C_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_D =  TF1(Form("f_D_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
    
    }
  
    else if(model_ == "KulaginNoIso"){
      isOverallScaling = false;
      minQ2 = 0.50;
      maxQ2 = 0.50;
      minX  = 0.001;
      maxX  = 0.958;
      ModelSFLoc =*directory_str +"/KPSFNoIso.root";
      f_Fe = TF1(Form("f_Fe_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_Pb = TF1(Form("f_Pb_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_C =  TF1(Form("f_C_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_D =  TF1(Form("f_D_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
    
    }
  
    else if(model_ == "KulaginTwistZero"){
      isOverallScaling = false;
      minQ2 = 0.1;
      maxQ2 = 29.9;
      minX  = 0.01;
      maxX  = 0.99;
      ModelSFLoc =*directory_str +"/HT0KPSF.root";
      f_Fe = TF1(Form("f_Fe_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_Pb = TF1(Form("f_Pb_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_C =  TF1(Form("f_C_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_D =  TF1(Form("f_D_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
    
    }
  
    else if(model_ == "KulaginTwistOne"){
      isOverallScaling = false;
      minQ2 = 0.50;
      maxQ2 = 9.9;
      minX  = 0.005;
      maxX  = 0.958;
      ModelSFLoc =*directory_str +"/HT1KPSF.root";
      f_Fe = TF1(Form("f_Fe_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_Pb = TF1(Form("f_Pb_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_C =  TF1(Form("f_C_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_D =  TF1(Form("f_D_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
    
    }
    
    else if (model_ == "CloetQ21"){
      isOverallScaling = false;
      minQ2 = 1.0;
      maxQ2 = 1.0;
      minX  = 0.05;
      maxX  = 1.0;
      ModelSFLoc =*directory_str +"/Cloet.root";
      f_Fe = TF1(Form("f_Fe_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_Pb = TF1(Form("f_Pb_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_C =  TF1(Form("f_C_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_D =  TF1(Form("f_D_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
    
    
    }
  
    else if (model_ == "CloetQ25"){
      isOverallScaling = false;
      minQ2 = 0.01;
      maxQ2 = 29.9;
      minX  = 0.01;
      maxX  = 0.99;
      ModelSFLoc =*directory_str +"/Cloet.root";
      f_Fe = TF1(Form("f_Fe_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_Pb = TF1(Form("f_Pb_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_C =  TF1(Form("f_C_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_D =  TF1(Form("f_D_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
    
    
    }
  
    else if (model_ == "CloetNoFreeP"){
      isOverallScaling = false;
      minQ2 = 0.01;
      maxQ2 = 29.9;
      minX  = 0.01;
      maxX  = 0.99;
      ModelSFLoc =*directory_str +"/CloetNoFreeP.root";
      f_Fe = TF1(Form("f_Fe_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_Pb = TF1(Form("f_Pb_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_C =  TF1(Form("f_C_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_D =  TF1(Form("f_D_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
    
    
    }
    
    //No model
    else if(model_ == "NoModel"){
      isOverallScaling = true;
      ModelSFLoc =*directory_str +"/FreeNucleonCloseout.root";
      f_Fe = TF1(Form("f_Fe_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_Pb = TF1(Form("f_Pb_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_C =  TF1(Form("f_C_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_D =  TF1(Form("f_D_%s", model.c_str()), "1.0+0.0*x" ,0, 1);               
    }
    
    //If nothing specified, go with the CV
    else{
      isOverallScaling = false;
      //minQ2 = 0.8;
      ModelSFLoc =*directory_str +"/FreeNucleonCloseout.root";
      f_Fe = TF1(Form("f_Fe_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_Pb = TF1(Form("f_Pb_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_C =  TF1(Form("f_C_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
      f_D =  TF1(Form("f_D_%s", model.c_str()), "1.0+0.0*x" ,0, 1);
    
    
    
    }
    
    //These are selected by input file
    TFile *fileM = new TFile(ModelSFLoc.c_str(), "READ");
    
    if(fileM && !isOverallScaling){
      //F1_P   = (TH2F*)fileM->Get("F1P");
      //F1_N   = (TH2F*)fileM->Get("F1N");
      F1_CH   = (TH2F*)fileM->Get("F1CH");
      F1_C   = (TH2F*)fileM->Get("F1C");  
      F1_Fe  = (TH2F*)fileM->Get("F1Fe");
      F1_Pb  = (TH2F*)fileM->Get("F1Pb");
  
      /*
      F1_P->SetDirectory(0);
      F1_N->SetDirectory(0);
      F1_CH->SetDirectory(0);
      F1_C->SetDirectory(0);
      F1_Fe->SetDirectory(0);
      F1_Pb->SetDirectory(0);
      */   
   
     // F2_P   = (TH2F*)fileM->Get("F2P");
     // F2_N   = (TH2F*)fileM->Get("F2N");
      F2_CH   = (TH2F*)fileM->Get("F2CH");
      F2_C   = (TH2F*)fileM->Get("F2C");  
      F2_Fe  = (TH2F*)fileM->Get("F2Fe");
      F2_Pb  = (TH2F*)fileM->Get("F2Pb");
  
      /*
      F2_P->SetDirectory(0);
      F2_N->SetDirectory(0);
      F2_CH->SetDirectory(0);
      F2_C->SetDirectory(0);
      F2_Fe->SetDirectory(0);
      F2_Pb->SetDirectory(0);
      */
  
     // xF3_P   = (TH2F*)fileM->Get("xF3P");
     // xF3_N   = (TH2F*)fileM->Get("xF3N");
      xF3_CH   = (TH2F*)fileM->Get("xF3CH");
      xF3_C   = (TH2F*)fileM->Get("xF3C");  
      xF3_Fe  = (TH2F*)fileM->Get("xF3Fe");
      xF3_Pb  = (TH2F*)fileM->Get("xF3Pb");
    
      /*  
      xF3_P->SetDirectory(0);
      xF3_N->SetDirectory(0);
      xF3_CH->SetDirectory(0);
      xF3_C->SetDirectory(0);
      xF3_Fe->SetDirectory(0);
      xF3_Pb->SetDirectory(0);
      */
  
      fileM->Close();
      delete fileM;
    }
    
    else
      Error("MnvNuclearModelWeight::()", "Input file for Model SFs is NULL!");
                      
  
  }
  //================================
  // Destructor
  //================================
  MnvNuclearModelWeight::~MnvNuclearModelWeight() {
  
    /*if(F1_P)
     delete F1_P;
    if(F1_N)
     delete F1_N;*/ 
    if(F1_CH)
     delete F1_CH; 
    if(F1_C)
     delete F1_C;
    if(F1_Fe)
     delete F1_Fe;
    if(F1_Pb)
     delete F1_Pb;
    
    /*if(F2_P)
     delete F2_P;
    if(F2_N)
     delete F2_N; */
    if(F2_CH)
     delete F2_CH; 
    if(F2_C)
     delete F2_C;
    if(F2_Fe)
     delete F2_Fe;
    if(F2_Pb)
     delete F2_Pb;      
    
    /*if(xF3_P)
     delete xF3_P;
    if(xF3_N)
     delete xF3_N;*/
    if(xF3_CH)
     delete xF3_CH;
    if(xF3_C)
     delete xF3_C;
    if(xF3_Fe)
     delete xF3_Fe;
    if(xF3_Pb)
     delete xF3_Pb;
  
    if(GF1FeP)
      delete GF1FeP;
    if(GF2FeP)
      delete GF2FeP;
    if(GxF3FeP)
      delete GxF3FeP;
      
    if(GF1FeN)
      delete GF1FeN;
    if(GF2FeN)
      delete GF2FeN;
    if(GxF3FeN)
      delete GxF3FeN;
    
    if(GF1P)
      delete GF1P;
    if(GF2P)
      delete GF2P;
    if(GxF3P)
      delete GxF3P; 
    
    if(GF1N)
      delete GF1N;
    if(GF2N)
      delete GF2N;
    if(GxF3P)
      delete GxF3N;                 
      
  
  };
  
  double MnvNuclearModelWeight::CalculateWeight(double Q2, double x, double y, double Enu, int tgtNucleon, int Target_A){
    bool Debug = false;
   
    // Compiler warning for unused variable tgtNucleon
    tgtNucleon = 0;
   
    //double Weight = 1.0;
    const double TargetMass =(.9396 + .9383 ) / 2.;
    double ScaledX = x;
    
    //BY13 uses a different scaling variable.
    if(model_ == "BY13"){
       double Ehad = Q2 / (2*TargetMass*x);
       ScaledX = Q2 / (TargetMass*Ehad*(1+sqrt(1+Q2/pow(Ehad, 2))));
    
    }
    
    //Check to see if we are in a kinematic region allowed by the model.
    //If not, freeze the calculation at the minimum / maximum point
    if(x < minX ){
      // Warning("MnvNuclearModelWeight::CalculateWeight", Form("x = %.2f less than x min: %.2f  using x min", x, minX) );
       x = minX;
    
    }
    
    if(x > maxX ){
      // Warning("MnvNuclearModelWeight::CalculateWeight", Form("x = %.2f greater than x max: %.2f  using x max", x, maxX) );
       x = maxX;
    
    }
    
    if(Q2 < minQ2 ){
      // Warning("MnvNuclearModelWeight::CalculateWeight", Form("Q2 = %.2f less than Q2 min: %.2f  using Q2 min", Q2, minQ2) );
       Q2 = minQ2;
    
    }
    
    if(Q2 > maxQ2 ){
      // Warning("MnvNuclearModelWeight::CalculateWeight", Form("Q2 = %.2f greater than Q2 max: %.2f  using Q2 max", Q2, maxQ2) );
       Q2 = maxQ2;
    
    }  
      
    //Is this a flat weight?
    if(isOverallScaling){
       
       //Hyrdrogen is 1
       if(Target_A == 1)
          return 1.0;
       else if(Target_A == 12)
          return f_C.Eval(ScaledX)*f_D.Eval(x)/NuclModUtils::GENIE_BY_bug(x, Q2, Target_A);
       else if(Target_A == 56)
          return f_Fe.Eval(ScaledX)*f_D.Eval(x)/NuclModUtils::GENIE_BY_bug(x, Q2, Target_A);    
       else if(Target_A == 207)
          return f_Pb.Eval(ScaledX)*f_D.Eval(x)/NuclModUtils::GENIE_BY_bug(x, Q2, Target_A);
       else
          return f_Fe.Eval(ScaledX)*f_D.Eval(x)/NuclModUtils::GENIE_BY_bug(x, Q2, Target_A);
    }//end of isOverallScaling
    
    else{
       
       int xBinM = -1;
       int q2BinM = -1;
       
       //xBinM = x*100; // KP
       xBinM = x*1000; //Cloet
       //q2BinM = Q2*100;
       //xBinM  = F1_Fe->GetXaxis()->FindBin(x);     
       q2BinM = F1_Fe->GetYaxis()->FindBin(Q2);
       
      // if(Debug)
      // Warning("MnvNuclearModelWeight", Form("x bin: %d, y bin: %d", xBin, q2Bin) );
      
       //const double Num = CalculateDifferentialCrossSection(Q2, x, y, Enu, Target_A );
       const double Num = CalculateDifferentialCrossSection(q2BinM, xBinM, x, y, Enu, Target_A );
       
       int xBin = -1;
       int q2Bin = -1;
       //xBin  = GF1FeP->GetXaxis()->FindBin(x);     
       q2Bin = GF1FeP->GetYaxis()->FindBin(Q2);
       //q2Bin = Q2*100;
       xBin = x*100;
       //xBin = xBin - 1;
       
      
       //const double Denom = CalculateGENIEDifferentialCrossSection(Q2, x, y, Enu, Target_A);
       const double Denom = CalculateGENIEDifferentialCrossSection(q2Bin, xBin, x, y, Enu, Target_A);
       
       Debug = Debug && (Q2 >= 1.0 && Q2 < 1.5 && x > 0.1 && x < 0.15);
       
       if(Debug){
          
  	Warning("MnvNuclearModelWeight::CalculateWeight", Form("x: %.3f y: %.2f Enu: %.2f Q2: %.2f F1: %.2f Rat: %.2f", x, y, Enu, Q2, F1_C->GetBinContent(xBinM, q2BinM), F1_Pb->GetBinContent(xBinM, q2BinM) )); 
          Warning("MnvNuclearModelWeight::CalculateWeight", Form("x: %.3f y: %.2f Enu: %.2f Q2: %.2f F2: %.2f Rat: %.2f", x, y, Enu, Q2, F2_C->GetBinContent(xBinM, q2BinM), F2_Pb->GetBinContent(xBinM, q2BinM) ) );
         Warning("MnvNuclearModelWeight::CalculateWeight", Form("x: %.3f y: %.2f Enu: %.2f Q2: %.2f xF3: %.2f Rat: %.2f", x, y, Enu, Q2, xF3_C->GetBinContent(xBinM, q2BinM), xF3_Pb->GetBinContent(xBinM, q2BinM) ) );
         if(Num == 0){
            Warning("MnvNuclearModelWeight::CalculateWeight", Form("x: %.7f bin: %d Q2: %.2f bin: %d xF3 %.2f", x, xBin, Q2,  q2Bin, xF3_Fe->GetBinContent(xBinM, q2BinM) ) );
         
         }
         
     }
       
       if(Debug){
         Warning("MnvNuclearModelWeight::CalculateWeight", Form("x: %.3f y: %.2f Enu: %.2f Q2: %.2f GF1P Fe: %.2f GF1N Fe: %.2f ", x, y, Enu, Q2, GF1FeP->GetBinContent(xBin, q2Bin), GF1FeN->GetBinContent(xBin, q2Bin) ) ); 
         Warning("MnvNuclearModelWeight::CalculateWeight", Form("x: %.3f y: %.2f Enu: %.2f Q2: %.2f GF2P Fe: %.2f GF2N Fe: %.2f", x, y, Enu, Q2, GF2FeP->GetBinContent(xBin, q2Bin), GF2FeN->GetBinContent(xBin, q2Bin) ) );
         Warning("MnvNuclearModelWeight::CalculateWeight", Form("x: %.3f y: %.2f Enu: %.2f Q2: %.2f GxF3P Fe: %.2f GxF3N Fe: %.2f", x, y, Enu, Q2, GxF3FeP->GetBinContent(xBin, q2Bin), GxF3FeN->GetBinContent(xBin, q2Bin) ) );
       }
       
       
       if(Denom == 0.0)
         return -1.0;
       
       
       if(Debug){              	
          
  	Warning("MnvNuclearModelWeight::CalculateWeight", Form("Weight: %.2f", (Num/Denom) ) );
          Warning("MnvNuclearModelWeight::CalculateWeight", Form("Num: %.2f", (Num) ) );
          Warning("MnvNuclearModelWeight::CalculateWeight", Form("Denom: %.2f", (Denom) ) );
  		                  
        }
           
         
         
         
         return Num/Denom;					  
    }//end of else
  
  
  }//end of CalculateWeight
  
  
  double  MnvNuclearModelWeight::CalculateDifferentialCrossSection(double Q2, double x, double y, double Enu, int Target_A){
     double dSigma = 0.0;   
     double M = (.9396 + .9383 ) / 2.;
  
     TH2F *F1;
     TH2F *F2;
     TH2F *xF3;
     if(Target_A == 0 ){
       F1  = F1_CH;
       F2  = F2_CH;
       xF3 = xF3_CH;
     
     }
     
     else if(Target_A == 1 ){
       F1  = F1_P;
       F2  = F2_P;
       xF3 = xF3_P;
     
     }
     
     else if(Target_A == 12){
       F1  = F1_C;
       F2  = F2_C;
       xF3 = xF3_C;
     
     }
  
     else if(Target_A == 56){
       F1  = F1_Fe;
       F2  = F2_Fe;
       xF3 = xF3_Fe;
     
     }
     
     else if(Target_A == 207){
       F1  = F1_Pb;
       F2  = F2_Pb;
       xF3 = xF3_Pb;
     
     }
     
     else{
       F1  =  F1_Fe;
       F2  =  F2_Fe;
       xF3 = xF3_Fe;
    
    }
    
       dSigma = 0.5*pow(y,2)*F1->Interpolate(x, Q2) + (1 - y - ((M*x*y)/(2*Enu)) )*F2->Interpolate(x, Q2) + y*(1-y/2)*xF3->Interpolate(x, Q2); 
                   
     return dSigma; //fudge it
  }//end of CalculateCrossSection
  
  double  MnvNuclearModelWeight::CalculateGENIEDifferentialCrossSection(double Q2, double x, double y, double Enu, int Target_A){
     double dSigma = 0.0;
     double dSigmaP = 0.0;
     double dSigmaN = 0.0;
     double M = (.9396 + .9383 ) / 2.;
     int Z = 0;
     if(Target_A == 0 ){
       Z = 7;
       Target_A = 13;
     
     }
     
     else if(Target_A == 1){
       Z = 1;
     
     }
     
     else if(Target_A == 12){
       Z = 6;
     
     }
  
     else if(Target_A == 56){
       Z = 26;
     
     }
     
     else if(Target_A == 207){
       Z = 82;
     
     }
     
     else{
       Z = Target_A / 2;
    
    }
  
     
    if(Target_A != 1){
      //Proton Piece       
      M = 0.9383;
      dSigmaP = 0.5*pow(y,2)*GF1FeP->Interpolate(x, Q2) + (1 - y - ((M*x*y)/(2*Enu)) )*GF2FeP->Interpolate(x, Q2) + y*(1-y/2)*GxF3FeP->Interpolate(x, Q2);
    dSigmaP *= Z;
    
      //Neutron Piece
      M = 0.9396;
      dSigmaN = 0.5*pow(y,2)*GF1FeN->Interpolate(x, Q2) + (1 - y - ((M*x*y)/(2*Enu)) )*GF2FeN->Interpolate(x, Q2) + y*(1-y/2)*GxF3FeN->Interpolate(x, Q2);
    dSigmaN *= (Target_A - Z);
    }
    
    else {
      //Proton Piece       
      M = 0.9383;
      dSigmaP = 0.5*pow(y,2)*GF1P->Interpolate(x, Q2) + (1 - y - ((M*x*y)/(2*Enu)) )*GF2P->Interpolate(x, Q2) + y*(1-y/2)*GxF3P->Interpolate(x, Q2);
    dSigmaP *= Z;
    
      //Neutron Piece
      M = 0.9396;
      dSigmaN = 0.5*pow(y,2)*GF1N->Interpolate(x, Q2) + (1 - y - ((M*x*y)/(2*Enu)) )*GF2N->Interpolate(x, Q2) + y*(1-y/2)*GxF3N->Interpolate(x, Q2);
    dSigmaN *= (Target_A - Z);
    }  
    
    //sum and divide by the number of nucleons
    if(Target_A != 0)
      dSigma = (dSigmaP + dSigmaN)/(Target_A);
    else   
      dSigma = (dSigmaP + dSigmaN);    
     
       
     return dSigma;
  }//end of CalculateCrossSection
  
  double  MnvNuclearModelWeight::CalculateGENIEDifferentialCrossSection(int q2Bin, int xBin, double x, double y, double Enu, int Target_A){
     double dSigma = 0.0;
     double dSigmaP = 0.0;
     double dSigmaN = 0.0;
     double M = (.9396 + .9383 ) / 2.;
     int Z = 0;
     if(Target_A == 0 ){
       Z = 7;
       Target_A = 13;
     
     }
     
     else if(Target_A == 1){
       Z = 1;
     
     }
     
     else if(Target_A == 12){
       Z = 6;
     
     }
  
     else if(Target_A == 56){
       Z = 26;
     
     }
     
     else if(Target_A == 207){
       Z = 82;
     
     }
     
     else{
       Z = Target_A / 2;
    
    }
  
     
    if(Target_A != 1){
      //Proton Piece       
      //M = 0.9383;
      dSigmaP = 0.5*pow(y,2)*GF1FeP->GetBinContent(xBin, q2Bin) + (1 - y - ((M*x*y)/(2*Enu)) )*GF2FeP->GetBinContent(xBin, q2Bin) + y*(1-y/2)*GxF3FeP->GetBinContent(xBin, q2Bin);
    dSigmaP *= Z;
    
      //Neutron Piece
      //M = 0.9396;
      dSigmaN = 0.5*pow(y,2)*GF1FeN->GetBinContent(xBin, q2Bin) + (1 - y - ((M*x*y)/(2*Enu)) )*GF2FeN->GetBinContent(xBin, q2Bin) + y*(1-y/2)*GxF3FeN->GetBinContent(xBin, q2Bin);
    dSigmaN *= (Target_A - Z);
    }
    
    else{
      //Proton Piece       
      //M = 0.9383;
      dSigmaP = 0.5*pow(y,2)*GF1P->GetBinContent(xBin, q2Bin) + (1 - y - ((M*x*y)/(2*Enu)) )*GF2P->GetBinContent(xBin, q2Bin) + y*(1-y/2)*GxF3P->GetBinContent(xBin, q2Bin);
    dSigmaP *= Z;
    
      //Neutron Piece
      //M = 0.9396;
      dSigmaN = 0.5*pow(y,2)*GF1N->GetBinContent(xBin, q2Bin) + (1 - y - ((M*x*y)/(2*Enu)) )*GF2N->GetBinContent(xBin, q2Bin) + y*(1-y/2)*GxF3N->GetBinContent(xBin, q2Bin);
    dSigmaN *= (Target_A - Z);
    }  
    
    //sum and divide by the number of nucleons
    if(Target_A != 0)
      dSigma = (dSigmaP + dSigmaN)/(Target_A);
    else   
      dSigma = (dSigmaP + dSigmaN);    
     
       
     return dSigma;
  }//end of CalculateCrossSection
  
  double  MnvNuclearModelWeight::CalculateDifferentialCrossSection(int q2Bin, int xBin, double x, double y, double Enu, int Target_A){
     double dSigma = 0.0;   
     double M = (.9396 + .9383 ) / 2.;
  
     TH2F *F1;
     TH2F *F2;
     TH2F *xF3;
     if(Target_A == 0 ){
       F1  = F1_CH;
       F2  = F2_CH;
       xF3 = xF3_CH;
     
     }
     
     else if(Target_A == 1 ){
       F1  = F1_P;
       F2  = F2_P;
       xF3 = xF3_P;
     
     }
     
     else if(Target_A == 12){
       F1  = F1_C;
       F2  = F2_C;
       xF3 = xF3_C;
     
     }
  
     else if(Target_A == 56){
       F1  = F1_Fe;
       F2  = F2_Fe;
       xF3 = xF3_Fe;
     
     }
     
     else if(Target_A == 207){
       F1  = F1_Pb;
       F2  = F2_Pb;
       xF3 = xF3_Pb;
     
     }
     
     else{
       F1  =  F1_Fe;
       F2  =  F2_Fe;
       xF3 = xF3_Fe;
    
    }
    
       dSigma = 0.5*pow(y,2)*F1->GetBinContent(xBin, q2Bin) + (1 - y - ((M*x*y)/(2*Enu)) )*F2->GetBinContent(xBin, q2Bin) + y*(1-y/2)*xF3->GetBinContent(xBin, q2Bin); 
                   
     return dSigma;
  }//end of CalculateCrossSection
  
  //Placeholder for Whitlow
  double  MnvNuclearModelWeight::CalcR(double Q2, double x) const{
    const double C2 = TMath::Power(0.125, 2);
    const double B1 =  0.0635;
    const double B2 =  0.5747;
    const double B3 = -0.3534;
  
    double Q2R   = TMath::Max(Q2, 0.35);
    double ss    = TMath::Log(Q2R/.04);
  
    double x2    = TMath::Power(x,   2.);
    double Q4R   = TMath::Power(Q2R, 2.);
    double Q4    = TMath::Power(Q2,  2.);
  
    double theta = 1. + (12.*Q2R/(Q2R+1.)) * (C2/(C2+x2));
    double R     = (B1/ss)*theta + B2/Q2R + B3/(Q4R+.09);
  
    R = TMath::Max(R,0.);
  
    if(Q2 < 0.35)  {
        R *= ( 3.207*(Q2/(Q4+1.)) );
    }
    return R;
  }//end of CalcR
  
  double  MnvNuclearModelWeight::CalcR(double Q2, double x, int fit) const{
  //R is taken from Abe et al, Physics Letters B 452 1999 194 - 200
  //It's similar to Whitlow, but extends further down in x
  const double FitA[] = {0.0485, 0.547, 2.062, -0.3804, 0.509, -0.0285};
  const double FitB[] = {0.0481, 0.6114, -0.3509, -0.4611, 0.7172, -0.037};
  const double FitC[] = {0.0577, 0.4644, 1.8288, 12.3708, -43.1043, 41.7415};
  std::vector<double> Coeffs;
  
  if(fit == 1)
    Coeffs.assign(FitA, FitA+6);
  else if(fit == 2)
    Coeffs.assign(FitB, FitB+6);
  else if(fit == 3)
    Coeffs.assign(FitC, FitC+6);
  else{
    Error( "MnvNuclearModelWeight::LoadModel", "Need to specify a fit number (1, 2 or 3)  !");
    return -1.0; 
    }
    
  const double  ThetaTerm = Coeffs[0]/log10(Q2/0.04)*Theta(Q2, x);
  
  double SecondTerm = 0.0;
  
  if(fit == 1)
    SecondTerm = Coeffs[1]*pow(x, Coeffs[5])*(1 + Coeffs[3]*x + Coeffs[4]*pow(x, 2))/(pow(Q2, 4)+pow(Coeffs[2], 4));
  else if(fit == 2)
    SecondTerm = pow(x, Coeffs[5])*(1 + Coeffs[3]*x + Coeffs[4]*pow(x, 2))*(Coeffs[1]/Q2 + Coeffs[2]/(pow(Q2, 2) + pow(0.3, 2)));
  else if(fit == 3)
   SecondTerm = Coeffs[1]/sqrt(pow( (Q2 - Coeffs[3] + Coeffs[4]*pow(x,2) + Coeffs[5]*pow(x, 3) ), 2) + pow(Coeffs[2], 2));
  else{
    Error("MnvNuclearModelWeight::CalcR", "Need to specify a fit number (1, 2 or 3)  !");
    return -1.0; 
    }
    
  return ThetaTerm + SecondTerm;      
  
  }//end of CalcR
  
  double  MnvNuclearModelWeight::Theta(double Q2, double x) const{
  
     return 1 + 12*(Q2/(Q2 + 1))*pow(0.125, 2)/(pow(0.125, 2) + pow(x, 2));
  
  
  }//end of Theta
} //namespace PlotUtils
