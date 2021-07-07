#ifndef MNVHADRONREWEIGHT_CXX 
#define MNVHADRONREWEIGHT_CXX 

#include <TF1.h>
#include <TLeaf.h>
#include "MnvHadronReweight.h"
#include "NSFDefaults.h"
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TGraph.h>

#include <stdlib.h>
//#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdio.h>

using std::vector;


using namespace PlotUtils;

namespace XSecFunctions 
{
  // Carbon
  TGraph* carbon_pion_tgraph;
  TGraph* carbon_proton_tgraph;
  TGraph* carbon_neutron_tgraph;
  TGraph* carbon_pion_elastic_tgraph;
  TGraph* carbon_proton_elastic_tgraph;
  TGraph* carbon_neutron_elastic_tgraph;
  TGraph* carbon_pion_total_tgraph;
  TGraph* carbon_proton_total_tgraph;
  TGraph* carbon_neutron_total_tgraph;
  // Lead
  TGraph* lead_pion_tgraph;
  TGraph* lead_proton_tgraph;
  TGraph* lead_neutron_tgraph;
  TGraph* lead_pion_elastic_tgraph;
  TGraph* lead_proton_elastic_tgraph;
  TGraph* lead_neutron_elastic_tgraph;
  TGraph* lead_pion_total_tgraph;
  TGraph* lead_proton_total_tgraph;
  TGraph* lead_neutron_total_tgraph;
  // Iron
  TGraph* iron_pion_tgraph;
  TGraph* iron_proton_tgraph;
  TGraph* iron_neutron_tgraph;
  TGraph* iron_pion_elastic_tgraph;
  TGraph* iron_proton_elastic_tgraph;
  TGraph* iron_neutron_elastic_tgraph;
  TGraph* iron_pion_total_tgraph;
  TGraph* iron_proton_total_tgraph;
  TGraph* iron_neutron_total_tgraph;
  // Hydrogen
  TGraph* hydrogen_neutron_total_tgraph;
  TGraph* hydrogen_neutron_elastic_tgraph;
  TGraph* hydrogen_neutron_inelastic_tgraph;
  // New Carbon neutron
  TGraph* new_carbon_neutron_elastic_tgraph;
  TGraph* new_carbon_neutron_inelastic_tgraph;
  TGraph* new_carbon_neutron_total_tgraph;
  double changeLowEElasticXSec = -1;
  double useHDXSec;
  //const double scintillator_hydrogen_fraction = 0.077418;
  TGraph* loadTGraph(const char* filename, const char* graphname)
  {
    std::cout<<"Trying to load filename: "<<filename<<", graphname: "<<graphname<<std::endl; // TODO remove
    std::string plotutils = std::string(getenv("PLOTUTILSROOT"));
    if (plotutils == std::string("")) 
    {
      std::cout<<"PlotUtils.MnvHadronReweight.XSecFunctions.loadTGraph: ERROR: $PLOTUTILSROOT is not defined"<<graphname<<std::endl;
      throw 66;
    }
    std::string fullfilename_temp = plotutils + std::string(filename);
    const char* fullfilename = fullfilename_temp.c_str();
    TFile f(fullfilename);
    if (f.IsZombie()) 
    {
      std::cout<<"PlotUtils.MnvHadronReweight.XSecFunctions.loadTGraph: ERROR: File is missing or is zombie: "<<fullfilename<<std::endl;
      throw 66;
    }
    TGraph* out = (TGraph*)f.Get(graphname);
    if (!out)
    {
      std::cout<<"PlotUtils.MnvHadronReweight.XSecFunctions.loadTGraph: ERROR: TGraph is missing in file: "<<graphname<<std::endl;
      throw 66;
    }
    out->Sort(); // Sort points so it hopefully does a binary search, but this version might not have that feature
    return out;
  }
  double pip_s(double* x, double* par)
  {
    double ke = x[0];
    
    double gamS = 350.0;
    double wS   = 100.0;
    double mPi  = pion_mass;
 
    double f = sqrt( (pow(ke+mPi,2.0) - mPi*mPi)/(pow(180.0+mPi,2.0) - mPi*mPi) );
    f       *= 120.0*gamS*gamS/( pow(ke-wS,2.0) + gamS*gamS/4.0 );
    f       += 200.0*(1.0 - exp(-1.0*ke/900.0) );
    
    par = NULL;
    
    return f*(1.0e-27);  //units of cm^2
  }
  
  double pip_d(double* x, double* par)
  {
    double ke = x[0]; 
    
    double gamD = 200.0;
    double wD   = 140.0;
    double mPi  = pion_mass;
 
    double f = sqrt( (pow(ke+mPi,2.0) - mPi*mPi)/(pow(180.0+mPi,2.0) - mPi*mPi) );
    f       *= 65.0*gamD*gamD/( pow(ke-wD,2.0) + gamD*gamD/4.0 );
    f       += 100.0*(1.0 - exp(-1.0*ke/900.0) );
    
    par = NULL;
    
    return f*(1.0e-27);  //units of cm^2
  }
  
  double pip_t(double* x, double* par)
  {
    return pip_s(x,par) + pip_d(x,par);
  }

  double pro_t(double* x, double* par)
  {  
    double ke = x[0];
    double f = 250.0*exp( -pow(ke-25.0,2.0)/500.0) + 80.0*exp( -pow(ke-80.0,2.0)/2500.0) + 200.0*(1.0 + 0.0003*ke);
    
    par = NULL;
    return f*(1.0e-27);
  }

  double kaon_t(double* x, double* par)
  { 
    double ke = x[0];  
    
    par = NULL;
    return (100.0 + (ke/8.0))*(1.0e-27);
  }
  
  double tgrapheval(TGraph* graph, double x, double min, double max, bool lerpzero = true)
  {
    // It's safer to round into range than to extrapolate past the bounds
    if (x > max) x = max;
    if (x < min) // HD XSec goes all the way down to 1 MeV, so this isn't needed
    {
      if (lerpzero)
      {
        double foo = graph->Eval(min) * 1.0e-27;
        if (x < 0.5) x = 0.5;
        return foo * x / min; // do a lerp toward zero
      }
      else
      {
        double out = graph->Eval(min) * 1.0e-27;
        if (changeLowEElasticXSec > 0) 
        {
          // Want 1 at 5 MeV and 0 at 100 MeV
          double lerp = (100.0 - x) / 95.0;
          // At 
          out *= 1 + changeLowEElasticXSec * lerp;
        }
        return out;
      }
    }
    // Assumes tgraph is in mb
    return graph->Eval(x) * 1.0e-27;
  }
  
  double tgrapheval_withHD(TGraph* graph, double x, double min, double max, bool lerpzero = true)
  {
    // It's safer to round into range than to extrapolate past the bounds
    if (x > max) x = max;
    if (!useHDXSec && x < min) // HD XSec goes all the way down to 1 MeV, so this isn't needed
    {
      if (lerpzero)
      {
        double foo = graph->Eval(min) * 1.0e-27;
        if (x < 0.5) x = 0.5;
        return foo * x / min; // do a lerp toward zero
      }
      else
      {
        double out = graph->Eval(min) * 1.0e-27;
        if (changeLowEElasticXSec > 0) 
        {
          // Want 1 at 5 MeV and 0 at 100 MeV
          double lerp = (100.0 - x) / 95.0;
          // At 
          out *= 1 + changeLowEElasticXSec * lerp;
        }
        return out;
      }
    }
    // Assumes tgraph is in mb
    return graph->Eval(x) * 1.0e-27;
  }
  
  // The maxes for the tgraphs
  const double DEFAULT_MIN = 100, DEFAULT_MAX = 5000;
  const double DEFAULT_MIN_ALT = 40, DEFAULT_MAX_ALT = 2000; // For lead I did a smaller min
  // Carbon inelastic
  double neutron_inelastic_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval_withHD(carbon_neutron_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX); 
  }
  double pion_inelastic_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval(carbon_pion_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX); 
  }
  double proton_inelastic_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval(carbon_proton_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX); 
  }
  // Carbon elastic
  double neutron_elastic_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval_withHD(carbon_neutron_elastic_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
  }
  double pion_elastic_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval(carbon_pion_elastic_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
  }
  double proton_elastic_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval(carbon_proton_elastic_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
  }
  // Carbon total
  double neutron_total_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval_withHD(carbon_neutron_total_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
  }
  double pion_total_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval(carbon_pion_total_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
  }
  double proton_total_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval(carbon_proton_total_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
  }
  // Lead inelastic
  double lead_pion_inelastic_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval(lead_pion_tgraph, x[0], DEFAULT_MIN_ALT, DEFAULT_MAX_ALT); 
  }
  double lead_proton_inelastic_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval(lead_proton_tgraph, x[0], DEFAULT_MIN_ALT, DEFAULT_MAX_ALT); 
  }
  double lead_neutron_inelastic_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval_withHD(lead_neutron_tgraph, x[0], DEFAULT_MIN_ALT, DEFAULT_MAX_ALT); 
  }
  // Lead elastic
  double lead_pion_elastic_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval(lead_pion_elastic_tgraph, x[0], DEFAULT_MIN_ALT, DEFAULT_MAX_ALT, false); 
  }
  double lead_proton_elastic_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval(lead_proton_elastic_tgraph, x[0], DEFAULT_MIN_ALT, DEFAULT_MAX_ALT, false); 
  }
  double lead_neutron_elastic_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval_withHD(lead_neutron_elastic_tgraph, x[0], DEFAULT_MIN_ALT, DEFAULT_MAX_ALT, false); 
  }
  // Lead total
  double lead_pion_total_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval(lead_pion_total_tgraph, x[0], DEFAULT_MIN_ALT, DEFAULT_MAX_ALT, false); 
  }
  double lead_proton_total_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval(lead_proton_total_tgraph, x[0], DEFAULT_MIN_ALT, DEFAULT_MAX_ALT, false); 
  }
  double lead_neutron_total_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval_withHD(lead_neutron_total_tgraph, x[0], DEFAULT_MIN_ALT, DEFAULT_MAX_ALT, false); 
  }
  // Iron inelastic
  double iron_pion_inelastic_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval(iron_pion_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX); 
  }
  double iron_proton_inelastic_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval(lead_proton_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX); 
  }
  double iron_neutron_inelastic_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval_withHD(lead_neutron_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX); 
  }
  // Iron inelastic
  double iron_pion_elastic_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval(iron_pion_elastic_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
  }
  double iron_proton_elastic_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval(lead_proton_elastic_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
  }
  double iron_neutron_elastic_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval_withHD(lead_neutron_elastic_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
  }
  // Iron total
  double iron_pion_total_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval(iron_pion_total_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
  }
  double iron_proton_total_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval(lead_proton_total_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
  }
  double iron_neutron_total_func(double* x, double* par)
  { 
    par = NULL;
    return tgrapheval_withHD(lead_neutron_total_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
  }
  
  // Scintillator
  double scint_neutron_elastic_func(double* x, double* par)
  { 
    par = NULL;
    if (!useHDXSec) return tgrapheval_withHD(carbon_neutron_elastic_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
    else
    {
      double carbon = tgrapheval_withHD(carbon_neutron_elastic_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
      double hydrogen = tgrapheval_withHD(hydrogen_neutron_elastic_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
      return hydrogen + carbon;
    }
  }
  double scint_neutron_inelastic_func(double* x, double* par)
  { 
    par = NULL;
    if (!useHDXSec) return tgrapheval_withHD(carbon_neutron_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
    else
    {
      double carbon = tgrapheval_withHD(carbon_neutron_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
      double hydrogen = tgrapheval_withHD(hydrogen_neutron_inelastic_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
      return hydrogen + carbon;
    }
  }
  double scint_neutron_total_func(double* x, double* par)
  { 
    par = NULL;
    if (!useHDXSec) return tgrapheval_withHD(carbon_neutron_total_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
    else
    {
      double carbon = tgrapheval_withHD(carbon_neutron_total_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
      double hydrogen = tgrapheval_withHD(hydrogen_neutron_total_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
      return hydrogen + carbon;
    }
  }
  
  
  
  // Scintillator
  double new_scint_neutron_elastic_func(double* x, double* par)
  { 
    par = NULL;
    {
      double carbon = tgrapheval_withHD(new_carbon_neutron_elastic_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
      double hydrogen = tgrapheval_withHD(hydrogen_neutron_elastic_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
      return hydrogen + carbon;
    }
  }
  double new_scint_neutron_inelastic_func(double* x, double* par)
  { 
    par = NULL;
    {
      double carbon = tgrapheval_withHD(new_carbon_neutron_inelastic_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
      double hydrogen = tgrapheval_withHD(hydrogen_neutron_inelastic_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
      return hydrogen + carbon;
    }
  }
  double new_scint_neutron_total_func(double* x, double* par)
  { 
    par = NULL;
    {
      double carbon = tgrapheval_withHD(new_carbon_neutron_total_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
      double hydrogen = tgrapheval_withHD(hydrogen_neutron_total_tgraph, x[0], DEFAULT_MIN, DEFAULT_MAX, false); 
      return hydrogen + carbon;
    }
  }
  
  // Carbon
  TF1 carbon_pion_inelastic("carbon_pion_inelastic", pion_inelastic_func, 0.0, 5000.0, 0);
  TF1 carbon_proton_inelastic("carbon_proton_inelastic", proton_inelastic_func, 0.0, 5000.0, 0); 
  TF1 carbon_neutron_inelastic("carbon_neutron_inelastic", neutron_inelastic_func, 0.0, 5000.0, 0);
  TF1 carbon_pion_elastic("carbon_pion_elastic", pion_elastic_func, 0.0, 5000.0, 0);
  TF1 carbon_proton_elastic("carbon_proton_elastic", proton_elastic_func, 0.0, 5000.0, 0); 
  TF1 carbon_neutron_elastic("carbon_neutron_elastic", neutron_elastic_func, 0.0, 5000.0, 0);
  TF1 carbon_pion_total("carbon_pion_total", pion_total_func, 0.0, 5000.0, 0);
  TF1 carbon_proton_total("carbon_proton_total", proton_total_func, 0.0, 5000.0, 0); 
  TF1 carbon_neutron_total("carbon_neutron_total", neutron_total_func, 0.0, 5000.0, 0);
  // Scintillator
  TF1 scint_neutron_inelastic("scint_neutron_inelastic", scint_neutron_inelastic_func, 0.0, 5000.0, 0);
  TF1 scint_neutron_elastic("scint_neutron_elastic", scint_neutron_elastic_func, 0.0, 5000.0, 0); 
  TF1 scint_neutron_total("scint_neutron_total", scint_neutron_total_func, 0.0, 5000.0, 0);
  //TF1 carbon_kaon("carbon_kaon", kaon_t, 0.0, 5000.0, 0);
  // Lead
  TF1 lead_pion_inelastic("lead_pion_inelastic", lead_pion_inelastic_func, 0.0, 5000.0, 0);
  TF1 lead_proton_inelastic("lead_proton_inelastic", lead_proton_inelastic_func, 0.0, 5000.0, 0); 
  TF1 lead_neutron_inelastic("lead_neutron_inelastic", lead_neutron_inelastic_func, 0.0, 5000.0, 0);
  TF1 lead_pion_elastic("lead_pion_elastic", lead_pion_elastic_func, 0.0, 5000.0, 0);
  TF1 lead_proton_elastic("lead_proton_elastic", lead_proton_elastic_func, 0.0, 5000.0, 0); 
  TF1 lead_neutron_elastic("lead_neutron_elastic", lead_neutron_elastic_func, 0.0, 5000.0, 0);
  TF1 lead_pion_total("lead_pion_total", lead_pion_total_func, 0.0, 5000.0, 0);
  TF1 lead_proton_total("lead_proton_total", lead_proton_total_func, 0.0, 5000.0, 0); 
  TF1 lead_neutron_total("lead_neutron_total", lead_neutron_total_func, 0.0, 5000.0, 0);
  // Iron
  TF1 iron_pion_inelastic("iron_pion_inelastic", iron_pion_inelastic_func, 0.0, 5000.0, 0);
  TF1 iron_proton_inelastic("iron_proton_inelastic", iron_proton_inelastic_func, 0.0, 5000.0, 0); 
  TF1 iron_neutron_inelastic("iron_neutron_inelastic", iron_neutron_inelastic_func, 0.0, 5000.0, 0);
  TF1 iron_pion_elastic("iron_pion_elastic", iron_pion_elastic_func, 0.0, 5000.0, 0);
  TF1 iron_proton_elastic("iron_proton_elastic", iron_proton_elastic_func, 0.0, 5000.0, 0); 
  TF1 iron_neutron_elastic("iron_neutron_elastic", iron_neutron_elastic_func, 0.0, 5000.0, 0);
  TF1 iron_pion_total("iron_pion_total", iron_pion_total_func, 0.0, 5000.0, 0);
  TF1 iron_proton_total("iron_proton_total", iron_proton_total_func, 0.0, 5000.0, 0); 
  TF1 iron_neutron_total("iron_neutron_total", iron_neutron_total_func, 0.0, 5000.0, 0);
  // New carbon neutron
  TF1 new_carbon_neutron_elastic("new_carbon_neutron_elastic", new_scint_neutron_elastic_func, 0.0, 5000.0, 0);
  TF1 new_carbon_neutron_inelastic("new_carbon_neutron_inelastic", new_scint_neutron_inelastic_func, 0.0, 5000.0, 0);
  TF1 new_carbon_neutron_total("new_carbon_neutron_total", new_scint_neutron_total_func, 0.0, 5000.0, 0);
  
  
  bool loadTGraphs(bool hdXSec)
  {
    useHDXSec = hdXSec;
    if (!hdXSec)
    {
        const char* tgraphname1 = "inelXS";
        // Neutron
        carbon_neutron_tgraph = loadTGraph("/data/hadronReweight/cross_section_neutron_carbon.root", tgraphname1); 
        iron_neutron_tgraph = loadTGraph("/data/hadronReweight/cross_section_neutron_iron.root", tgraphname1);  
        lead_neutron_tgraph = loadTGraph("/data/hadronReweight/cross_section_neutron_lead.root", tgraphname1); 
        // Proton
        carbon_proton_tgraph = loadTGraph("/data/hadronReweight/cross_section_proton_carbon.root", tgraphname1); 
        iron_proton_tgraph = loadTGraph("/data/hadronReweight/cross_section_proton_iron.root", tgraphname1); 
        lead_proton_tgraph = loadTGraph("/data/hadronReweight/cross_section_proton_lead.root", tgraphname1);
        // Pion
        carbon_pion_tgraph = loadTGraph("/data/hadronReweight/cross_section_pion_carbon.root", tgraphname1); 
        iron_pion_tgraph = loadTGraph("/data/hadronReweight/cross_section_pion_iron.root", tgraphname1); 
        lead_pion_tgraph = loadTGraph("/data/hadronReweight/cross_section_pion_lead.root", tgraphname1); 
        const char* tgraphname2 = "elXS";
        // Neutron
        carbon_neutron_elastic_tgraph = loadTGraph("/data/hadronReweight/cross_section_neutron_carbon.root", tgraphname2); 
        iron_neutron_elastic_tgraph = loadTGraph("/data/hadronReweight/cross_section_neutron_iron.root", tgraphname2);  
        lead_neutron_elastic_tgraph = loadTGraph("/data/hadronReweight/cross_section_neutron_lead.root", tgraphname2); 
        // Proton
        carbon_proton_elastic_tgraph = loadTGraph("/data/hadronReweight/cross_section_proton_carbon.root", tgraphname2); 
        iron_proton_elastic_tgraph = loadTGraph("/data/hadronReweight/cross_section_proton_iron.root", tgraphname2); 
        lead_proton_elastic_tgraph = loadTGraph("/data/hadronReweight/cross_section_proton_lead.root", tgraphname2);
        // Pion
        carbon_pion_elastic_tgraph = loadTGraph("/data/hadronReweight/cross_section_pion_carbon.root", tgraphname2); 
        iron_pion_elastic_tgraph = loadTGraph("/data/hadronReweight/cross_section_pion_iron.root", tgraphname2); 
        lead_pion_elastic_tgraph = loadTGraph("/data/hadronReweight/cross_section_pion_lead.root", tgraphname2); 
        const char* tgraphname3 = "totalXS";
        // Neutron
        carbon_neutron_total_tgraph = loadTGraph("/data/hadronReweight/cross_section_neutron_carbon.root", tgraphname3); 
        iron_neutron_total_tgraph = loadTGraph("/data/hadronReweight/cross_section_neutron_iron.root", tgraphname3);  
        lead_neutron_total_tgraph = loadTGraph("/data/hadronReweight/cross_section_neutron_lead.root", tgraphname3); 
        // Proton
        carbon_proton_total_tgraph = loadTGraph("/data/hadronReweight/cross_section_proton_carbon.root", tgraphname3); 
        iron_proton_total_tgraph = loadTGraph("/data/hadronReweight/cross_section_proton_iron.root", tgraphname3); 
        lead_proton_total_tgraph = loadTGraph("/data/hadronReweight/cross_section_proton_lead.root", tgraphname3);
        // Pion
        carbon_pion_total_tgraph = loadTGraph("/data/hadronReweight/cross_section_pion_carbon.root", tgraphname3); 
        iron_pion_total_tgraph = loadTGraph("/data/hadronReweight/cross_section_pion_iron.root", tgraphname3); 
        lead_pion_total_tgraph = loadTGraph("/data/hadronReweight/cross_section_pion_lead.root", tgraphname3); 
    }
    else
    {
        std::cout<<"Using the HD XSec"<<std::endl;
        const char* tgraphname1 = "inelXS";
        // Neutron
        carbon_neutron_tgraph = loadTGraph("/data/hadronReweight/cross_section_neutron_carbon_hd.root", tgraphname1); 
        iron_neutron_tgraph = loadTGraph("/data/hadronReweight/cross_section_neutron_iron_hd.root", tgraphname1);  
        lead_neutron_tgraph = loadTGraph("/data/hadronReweight/cross_section_neutron_lead_hd.root", tgraphname1); 
        hydrogen_neutron_inelastic_tgraph = loadTGraph("/data/hadronReweight/cross_section_neutron_hydrogen_hd.root", tgraphname1); 
        // Proton
        carbon_proton_tgraph = loadTGraph("/data/hadronReweight/cross_section_proton_carbon.root", tgraphname1); 
        iron_proton_tgraph = loadTGraph("/data/hadronReweight/cross_section_proton_iron.root", tgraphname1); 
        lead_proton_tgraph = loadTGraph("/data/hadronReweight/cross_section_proton_lead.root", tgraphname1);
        // Pion
        carbon_pion_tgraph = loadTGraph("/data/hadronReweight/cross_section_pion_carbon.root", tgraphname1); 
        iron_pion_tgraph = loadTGraph("/data/hadronReweight/cross_section_pion_iron.root", tgraphname1); 
        lead_pion_tgraph = loadTGraph("/data/hadronReweight/cross_section_pion_lead.root", tgraphname1); 
        const char* tgraphname2 = "elXS";
        // Neutron
        carbon_neutron_elastic_tgraph = loadTGraph("/data/hadronReweight/cross_section_neutron_carbon_hd.root", tgraphname2); 
        iron_neutron_elastic_tgraph = loadTGraph("/data/hadronReweight/cross_section_neutron_iron_hd.root", tgraphname2);  
        lead_neutron_elastic_tgraph = loadTGraph("/data/hadronReweight/cross_section_neutron_lead_hd.root", tgraphname2); 
        hydrogen_neutron_elastic_tgraph = loadTGraph("/data/hadronReweight/cross_section_neutron_hydrogen_hd.root", tgraphname2); 
        // Proton
        carbon_proton_elastic_tgraph = loadTGraph("/data/hadronReweight/cross_section_proton_carbon.root", tgraphname2); 
        iron_proton_elastic_tgraph = loadTGraph("/data/hadronReweight/cross_section_proton_iron.root", tgraphname2); 
        lead_proton_elastic_tgraph = loadTGraph("/data/hadronReweight/cross_section_proton_lead.root", tgraphname2);
        // Pion
        carbon_pion_elastic_tgraph = loadTGraph("/data/hadronReweight/cross_section_pion_carbon.root", tgraphname2); 
        iron_pion_elastic_tgraph = loadTGraph("/data/hadronReweight/cross_section_pion_iron.root", tgraphname2); 
        lead_pion_elastic_tgraph = loadTGraph("/data/hadronReweight/cross_section_pion_lead.root", tgraphname2); 
        const char* tgraphname3 = "totalXS";
        // Neutron
        carbon_neutron_total_tgraph = loadTGraph("/data/hadronReweight/cross_section_neutron_carbon_hd.root", tgraphname3); 
        iron_neutron_total_tgraph = loadTGraph("/data/hadronReweight/cross_section_neutron_iron_hd.root", tgraphname3);  
        lead_neutron_total_tgraph = loadTGraph("/data/hadronReweight/cross_section_neutron_lead_hd.root", tgraphname3); 
        hydrogen_neutron_total_tgraph = loadTGraph("/data/hadronReweight/cross_section_neutron_hydrogen_hd.root", tgraphname3); 
        // Proton
        carbon_proton_total_tgraph = loadTGraph("/data/hadronReweight/cross_section_proton_carbon.root", tgraphname3); 
        iron_proton_total_tgraph = loadTGraph("/data/hadronReweight/cross_section_proton_iron.root", tgraphname3); 
        lead_proton_total_tgraph = loadTGraph("/data/hadronReweight/cross_section_proton_lead.root", tgraphname3);
        // Pion
        carbon_pion_total_tgraph = loadTGraph("/data/hadronReweight/cross_section_pion_carbon.root", tgraphname3); 
        iron_pion_total_tgraph = loadTGraph("/data/hadronReweight/cross_section_pion_iron.root", tgraphname3); 
        lead_pion_total_tgraph = loadTGraph("/data/hadronReweight/cross_section_pion_lead.root", tgraphname3); 
    }
    
    /*std::cout<<"TF1s:"<<std::endl; 
    
    std::cout<<carbon_pion_inelastic.GetName()<<std::endl;
    std::cout<<carbon_proton_inelastic.GetName()<<std::endl;
    std::cout<<carbon_neutron_inelastic.GetName()<<std::endl;
    std::cout<<carbon_pion_elastic.GetName()<<std::endl;
    std::cout<<carbon_proton_elastic.GetName()<<std::endl;
    std::cout<<carbon_neutron_elastic.GetName()<<std::endl;
    std::cout<<carbon_pion_total.GetName()<<std::endl;
    std::cout<<carbon_proton_total.GetName()<<std::endl;
    std::cout<<carbon_neutron_total.GetName()<<std::endl;
  //TF1 carbon_kaon("carbon_kaon", kaon_t, 0.0, 5000.0, 0);
  // Lead
    std::cout<<lead_pion_inelastic.GetName()<<std::endl;
    std::cout<<lead_proton_inelastic.GetName()<<std::endl;
    std::cout<<lead_neutron_inelastic.GetName()<<std::endl;
    std::cout<<lead_pion_elastic.GetName()<<std::endl;
    std::cout<<lead_proton_elastic.GetName()<<std::endl;
    std::cout<<lead_neutron_elastic.GetName()<<std::endl;
    std::cout<<lead_pion_total.GetName()<<std::endl;
    std::cout<<lead_proton_total.GetName()<<std::endl;
    std::cout<<lead_neutron_total.GetName()<<std::endl;
  // Iron
    std::cout<<iron_pion_inelastic.GetName()<<std::endl;
    std::cout<<iron_proton_inelastic.GetName()<<std::endl;
    std::cout<<iron_neutron_inelastic.GetName()<<std::endl;
    std::cout<<iron_pion_elastic.GetName()<<std::endl;
    std::cout<<iron_proton_elastic.GetName()<<std::endl;
    std::cout<<iron_neutron_elastic.GetName()<<std::endl;
    std::cout<<iron_pion_total.GetName()<<std::endl;
    std::cout<<iron_proton_total.GetName()<<std::endl;
    std::cout<<iron_neutron_total.GetName()<<std::endl; */
    return true;
  }
} // end XSecFunctions

/*MnvH1D* MnvHadronReweight::getPlotFromVertError(MnvH1D* hist, const std::string &name, int number)
{
  MnvVertErrorBand* v = hist->GetVertErrorBand(name);
  MnvH1D* h = new MnvH1D(*v->GetHist(number));
  h->AddMissingErrorBandsAndFillWithCV(*hist);
  //MnvH1D* output = new MnvH1D(*hist);
  //output->Multiply(hist, h);
  //output = 
  //delete h;
  return h;
}*/


HadronBranchContainer::HadronBranchContainer(TTree* tree)
{
      fOSF = false;
      fChain = tree;
      if (!tree->GetListOfBranches()->FindObject("truth_hadronReweightNPaths"))
      {
        std::cout<<"Hadron Reweight Error: Missing hadron reweight branches."<<std::endl;
        std::cout<<"Did you forget declareHadronReweightBranches and fillHadronReweightBranches?"<<std::endl;
        std::cout<<"Or are you using an old version? Revert to old version of MnvHadronReweight if you need to use it"<<std::endl;
        std::cout<<"but it's better to respin your ntuples."<<std::endl;
        std::cout<<"Switching to 'default mode' which turns off hadron reweight effectively."<<std::endl;
        defaultmode = true;
        return;
      }
      defaultmode = false;

      fChain->SetBranchStatus("truth_hadronReweightNPaths", 1);
      fChain->SetBranchStatus("truth_hadronReweightNPoints", 1);
      fChain->SetBranchStatus("truth_hadronReweightIntCode_sz", 1);
      fChain->SetBranchStatus("truth_hadronReweightIntCode", 1);
      fChain->SetBranchStatus("truth_hadronReweightNuke_sz", 1);
      fChain->SetBranchStatus("truth_hadronReweightNuke", 1);
      fChain->SetBranchStatus("truth_hadronReweightPDG_sz", 1);
      fChain->SetBranchStatus("truth_hadronReweightPDG", 1);
      fChain->SetBranchStatus("truth_hadronReweightTrackID_sz", 1);
      fChain->SetBranchStatus("truth_hadronReweightTrackID", 1);
      fChain->SetBranchStatus("truth_hadronReweightColumnarDensity_sz", 1);
      fChain->SetBranchStatus("truth_hadronReweightColumnarDensity", 1);
      fChain->SetBranchStatus("truth_hadronReweightFinalE_sz", 1);
      fChain->SetBranchStatus("truth_hadronReweightFinalE", 1);
      fChain->SetBranchStatus("truth_hadronReweightFinalSigmaE_sz", 1);
      fChain->SetBranchStatus("truth_hadronReweightFinalSigmaE", 1);
      fChain->SetBranchStatus("truth_hadronReweightInitialE_sz", 1);
      fChain->SetBranchStatus("truth_hadronReweightInitialE", 1);
      fChain->SetBranchStatus("truth_hadronReweightInitialSigmaE_sz", 1);
      fChain->SetBranchStatus("truth_hadronReweightInitialSigmaE", 1);
      fChain->SetBranchStatus("truth_hadronReweightPosX_sz", 1);
      fChain->SetBranchStatus("truth_hadronReweightPosX", 1);
      fChain->SetBranchStatus("truth_hadronReweightPosY_sz", 1);
      fChain->SetBranchStatus("truth_hadronReweightPosY", 1);
      fChain->SetBranchStatus("truth_hadronReweightPosZ_sz", 1);
      fChain->SetBranchStatus("truth_hadronReweightPosZ", 1);
      fChain->SetBranchStatus("truth_hadronReweightIntCodePerSegment_sz", 1);
      fChain->SetBranchStatus("truth_hadronReweightIntCodePerSegment", 1);

      fChain->SetBranchAddress("truth_hadronReweightNPaths", &truth_hadronReweightNPaths, &b_truth_hadronReweightNPaths);
      fChain->SetBranchAddress("truth_hadronReweightNPoints", &truth_hadronReweightNPoints, &b_truth_hadronReweightNPoints);
      fChain->SetBranchAddress("truth_hadronReweightIntCode_sz", &truth_hadronReweightIntCode_sz, &b_truth_hadronReweightIntCode_sz);
      fChain->SetBranchAddress("truth_hadronReweightIntCode", truth_hadronReweightIntCode, &b_truth_hadronReweightIntCode);
      fChain->SetBranchAddress("truth_hadronReweightNuke_sz", &truth_hadronReweightNuke_sz, &b_truth_hadronReweightNuke_sz);
      fChain->SetBranchAddress("truth_hadronReweightNuke", truth_hadronReweightNuke, &b_truth_hadronReweightNuke);
      fChain->SetBranchAddress("truth_hadronReweightPDG_sz", &truth_hadronReweightPDG_sz, &b_truth_hadronReweightPDG_sz);
      fChain->SetBranchAddress("truth_hadronReweightPDG", truth_hadronReweightPDG, &b_truth_hadronReweightPDG);
      fChain->SetBranchAddress("truth_hadronReweightTrackID_sz", &truth_hadronReweightTrackID_sz, &b_truth_hadronReweightTrackID_sz);
      fChain->SetBranchAddress("truth_hadronReweightTrackID", truth_hadronReweightTrackID, &b_truth_hadronReweightTrackID);
      fChain->SetBranchAddress("truth_hadronReweightColumnarDensity_sz", &truth_hadronReweightColumnarDensity_sz, &b_truth_hadronReweightColumnarDensity_sz);
      fChain->SetBranchAddress("truth_hadronReweightColumnarDensity", truth_hadronReweightColumnarDensity, &b_truth_hadronReweightColumnarDensity);
      fChain->SetBranchAddress("truth_hadronReweightFinalE_sz", &truth_hadronReweightFinalE_sz, &b_truth_hadronReweightFinalE_sz);
      fChain->SetBranchAddress("truth_hadronReweightFinalE", truth_hadronReweightFinalE, &b_truth_hadronReweightFinalE);
      fChain->SetBranchAddress("truth_hadronReweightFinalSigmaE_sz", &truth_hadronReweightFinalSigmaE_sz, &b_truth_hadronReweightFinalSigmaE_sz);
      fChain->SetBranchAddress("truth_hadronReweightFinalSigmaE", truth_hadronReweightFinalSigmaE, &b_truth_hadronReweightFinalSigmaE);
      fChain->SetBranchAddress("truth_hadronReweightInitialE_sz", &truth_hadronReweightInitialE_sz, &b_truth_hadronReweightInitialE_sz);
      fChain->SetBranchAddress("truth_hadronReweightInitialE", truth_hadronReweightInitialE, &b_truth_hadronReweightInitialE);
      fChain->SetBranchAddress("truth_hadronReweightInitialSigmaE_sz", &truth_hadronReweightInitialSigmaE_sz, &b_truth_hadronReweightInitialSigmaE_sz);
      fChain->SetBranchAddress("truth_hadronReweightInitialSigmaE", truth_hadronReweightInitialSigmaE, &b_truth_hadronReweightInitialSigmaE);
      fChain->SetBranchAddress("truth_hadronReweightPosX_sz", &truth_hadronReweightPosX_sz, &b_truth_hadronReweightPosX_sz);
      fChain->SetBranchAddress("truth_hadronReweightPosX", truth_hadronReweightPosX, &b_truth_hadronReweightPosX);
      fChain->SetBranchAddress("truth_hadronReweightPosY_sz", &truth_hadronReweightPosY_sz, &b_truth_hadronReweightPosY_sz);
      fChain->SetBranchAddress("truth_hadronReweightPosY", truth_hadronReweightPosY, &b_truth_hadronReweightPosY);
      fChain->SetBranchAddress("truth_hadronReweightPosZ_sz", &truth_hadronReweightPosZ_sz, &b_truth_hadronReweightPosZ_sz);
      fChain->SetBranchAddress("truth_hadronReweightPosZ", truth_hadronReweightPosZ, &b_truth_hadronReweightPosZ);
      fChain->SetBranchAddress("truth_hadronReweightIntCodePerSegment_sz", &truth_hadronReweightIntCodePerSegment_sz, &b_truth_hadronReweightIntCodePerSegment_sz);
      fChain->SetBranchAddress("truth_hadronReweightIntCodePerSegment", truth_hadronReweightIntCodePerSegment, &b_truth_hadronReweightIntCodePerSegment);

      if( fChain->GetBranch("mc_targetZ")->GetAddress() == nullptr )
      {
        //Suppress this message unless compiling in debug mode.  The -g flag to g++ will turn this line on for example.
        #ifndef NDEBUG
        std::cout<<"Setting up branches for the MHRW.  If you haven't set branch addresses before doing this (ie. MakeClass), this can cause problems for you."<<std::endl;
        #endif //NDEBUG

        fChain->SetBranchStatus("mc_targetZ",1);
        fChain->SetBranchStatus("mc_vtx",1);

        fChain->SetBranchStatus("mc_nFSPart",1);
        fChain->SetBranchStatus("mc_FSPartPDG",1);
        fChain->SetBranchStatus("mc_FSPartE",1);
        fChain->SetBranchStatus("mc_FSPartPx",1);
        fChain->SetBranchStatus("mc_FSPartPy",1);
        fChain->SetBranchStatus("mc_FSPartPz",1);
      
        fChain->SetBranchAddress("mc_targetZ"  , &mc_targetZ , &b_mc_targetZ   );
        fChain->SetBranchAddress("mc_vtx"      , mc_vtx      , &b_mc_vtx   );

        fChain->SetBranchAddress("mc_nFSPart"  , &mc_nFSPart , &b_mc_nFSPart   );
        fChain->SetBranchAddress("mc_FSPartPDG", mc_FSPartPDG, &b_mc_FSPartPDG );
        fChain->SetBranchAddress("mc_FSPartE"  , mc_FSPartE  , &b_mc_FSPartE   );
        fChain->SetBranchAddress("mc_FSPartPx" , mc_FSPartPx , &b_mc_FSPartPx  );
        fChain->SetBranchAddress("mc_FSPartPy" , mc_FSPartPy , &b_mc_FSPartPy  );
        fChain->SetBranchAddress("mc_FSPartPz" , mc_FSPartPz , &b_mc_FSPartPz  );
	      setOSF(false);
      }
      else 
      {
      	setOSF(true);
        mc_targetZ = *((Double_t*)fChain->GetBranch("mc_targetZ")->GetAddress());
        *mc_vtx     = *(Double_t*)fChain->GetLeaf("mc_vtx")->GetValuePointer();

        mc_nFSPart    = *((Int_t*)fChain->GetBranch("mc_nFSPart")->GetAddress());
        *mc_FSPartPDG = *(Double_t*)fChain->GetLeaf("mc_FSPartPDG")->GetValuePointer();
        *mc_FSPartE   = *(Double_t*)fChain->GetLeaf("mc_FSPartE")->GetValuePointer();
        *mc_FSPartPx  = *(Double_t*)fChain->GetLeaf("mc_FSPartPx")->GetValuePointer();
        *mc_FSPartPy  = *(Double_t*)fChain->GetLeaf("mc_FSPartPy")->GetValuePointer();
        *mc_FSPartPz  = *(Double_t*)fChain->GetLeaf("mc_FSPartPz")->GetValuePointer();
      }
}

HadronBranchContainer::~HadronBranchContainer( )
{}

Int_t HadronBranchContainer::GetEntry(Long64_t entry, Int_t getall)
{
  if (!fChain) std::cout<<"Warning: fChain not initialized."<<std::endl;
  double retval = fChain->GetEntry(entry, getall);

  //Sad hack... someone more clever than me can figure out the proper code.
  if(fOSF){
    this->mc_targetZ = *((Double_t*)fChain->GetBranch("mc_targetZ")->GetAddress());
    Double_t *t_mc_vtx = (Double_t*)fChain->GetBranch("mc_vtx")->GetAddress();
    for(int i=0;i<3;i++){
      this->mc_vtx[i] = t_mc_vtx[i];
    }
    this->mc_nFSPart    = *((Int_t*)fChain->GetBranch("mc_nFSPart")->GetAddress());
    //    std::cout << "What I think I have for FSI " << this->mc_nFSPart << std::endl;
    
    Int_t *t_mc_FSPartPDG = (Int_t*)fChain->GetBranch("mc_FSPartPDG")->GetAddress();
    Double_t *t_mc_FSPartE   = (Double_t*)fChain->GetBranch("mc_FSPartE")->GetAddress();
    Double_t *t_mc_FSPartPx  = (Double_t*)fChain->GetBranch("mc_FSPartPx")->GetAddress();
    Double_t *t_mc_FSPartPy  = (Double_t*)fChain->GetBranch("mc_FSPartPy")->GetAddress();
    Double_t *t_mc_FSPartPz  = (Double_t*)fChain->GetBranch("mc_FSPartPz")->GetAddress();
    for(int i=0;i<this->mc_nFSPart;i++){
      //      std::cout << i << "\t" <<t_mc_FSPartPDG[i] << std::endl;
      this->mc_FSPartPDG[i] = t_mc_FSPartPDG[i];
      this->mc_FSPartE[i]   = t_mc_FSPartE[i];
      this->mc_FSPartPx[i]  = t_mc_FSPartPx[i];
      this->mc_FSPartPy[i]  = t_mc_FSPartPy[i];
      this->mc_FSPartPz[i]  = t_mc_FSPartPz[i];
    }
  }

  return retval;
}

//TargetUtils
TargetUtils* MnvHadronReweight::m_TargetUtils;
//================================
// Constructors
//================================
MnvHadronReweight::MnvHadronReweight(TTree* truth, TTree* data, bool hdXSec) :
  fData(NULL), fTruth(NULL), 
  m_directory(Form("%s/data/mhrwKineRenorm",getenv("PLOTUTILSROOT"))),
  m_projectname("Renorm_Kine_Truth"), m_reweightNeutronCV(false),
  m_current_entry_neutronCV(-1),m_current_entry_weight(-1),
  m_current_rwgt_neutronCV(1.0)
{
      if (data != NULL) {
        setDataTree(data);
      }
      
      if (truth != NULL) {
        setTruthTree(truth);
      }
      
      // Load xsec tgraphs
      XSecFunctions::loadTGraphs(hdXSec);
      
      deltascale = 1; 
      
      // We have no fiducial volume defined to start. The user needs to define one
      fiducialType = -1;
      
      defaultmode = false; // For testing
      
      doElasticReweight = true;
      
      minimum_fake_elastic_threshold = 1.0;
      look_for_fake_elastics = false; // Need to turn on explicitly
      ResetFakeElasticsSeed();
      
      elasticReweightForAllParticles = false;
      m_reweightNeutronCV = false;
      m_playlist = "";
      m_kinefilename = "";     
 
      renorm = false;
      renormKine = false;

      initializeWeights(m_current_weights, 1.0); // init to 1

      m_TargetUtils=new TargetUtils();
}

void MnvHadronReweight::TurnOffElasticReweight() {doElasticReweight = false;}
void MnvHadronReweight::TurnOnElasticReweightForParticlesOtherThanNeutron() {elasticReweightForAllParticles = true;}
void MnvHadronReweight::TurnOnFakeElastics() {look_for_fake_elastics = true;
std::cout<<"ERROR: Fake Elastics are no longer needed. Are you sure you want to use them?\nI'm going to throw an error now. Don't use MnvHadronReweight::TurnOnFakeElastics."<<std::endl; throw 66;}
int MnvHadronReweight::GetNFakeElastics() {return nfakeelastics;}
void MnvHadronReweight::ResetFakeElasticsSeed() 
{
      srand (1123581321);
      nfakeelastics = 0;
}
void MnvHadronReweight::useHDXSec(bool useHDXSec) 
{
      XSecFunctions::loadTGraphs(useHDXSec);
}

void MnvHadronReweight::useReweightedNeutronCV()
{
  // Load XSec
  XSecFunctions::new_carbon_neutron_elastic_tgraph = XSecFunctions::loadTGraph("/data/hadronReweight/new_geant4_xsec.root", "elastic_xsec");  
  XSecFunctions::new_carbon_neutron_inelastic_tgraph = XSecFunctions::loadTGraph("/data/hadronReweight/new_geant4_xsec.root", "inelastic_xsec");  
  XSecFunctions::new_carbon_neutron_total_tgraph = XSecFunctions::loadTGraph("/data/hadronReweight/new_geant4_xsec.root", "total_xsec");  
  // Really want the HDXSec
  useHDXSec(true); 
  m_reweightNeutronCV = true;
}

//================================
// Singleton
//================================
MnvHadronReweight* MnvHadronReweight::instance;
MnvHadronReweight* MnvHadronReweight::get() 
{
  if (!instance) 
    instance = new MnvHadronReweight();
  return instance;
}

//================================
// Destructor
//================================
MnvHadronReweight::~MnvHadronReweight( )
{
  // 
}

void MnvHadronReweight::setdefaultmode(bool mode) 
{ 
  // This is just to test. Sets all weights to 1 if true
  defaultmode = mode; 
}

void MnvHadronReweight::setParticle(int pdg, bool treatAntiParticleSame)
{
  // The first arg is which pdg
  // The second arg is whether to treat the anti-particle the same as the particle
  // For example, you want to treat anti-pions same as pions
  particles[pdg] = treatAntiParticleSame;  // TODO This doesn't work at all
}

void MnvHadronReweight::removeParticle(int pdg)
{
  // Removes the pdg from the particle list
  particles.erase(pdg);
}

void MnvHadronReweight::useDefaultParticles()
{
  // Adds pions, protons and neutrons to list
  setParticle(PDGPION, true); // treat anti-pions same as pions
  setParticle(PDGPROTON, false);
  setParticle(PDGNEUTRON, false);
}

void MnvHadronReweight::setTruthTree(TTree* truth) 
{
  fTruthChain = truth;
  fTruthChain->GetEntry(0);
  //fTruthChain->SetMakeClass(1);// TODO was this needed?
  if(fTruth)
  {
    delete fTruth;
    fTruth = NULL;
  }
  fTruth = new HadronBranchContainer(truth);
}

void MnvHadronReweight::setDataTree(TTree* data) 
{
  fChain = data;
  fChain->GetEntry(0);
  //fCurrent = -1;
  //fChain->SetMakeClass(1);
  if(fData)
  {
    delete fData;
    fData = NULL;
  }
  fData = new HadronBranchContainer(data); 
  //fData = NULL;
}

TTree* MnvHadronReweight::makeTree(const char* filename, const char* treename)
{
  return makeTree(filename, false, treename);
}

TTree* MnvHadronReweight::makeTree(const char* filename, bool usetruth, const char* treename)
{
  // Makes a ttree that you can friend with another ttree. This makes testing simpler.
  
  TFile *file=new TFile(filename,"recreate");
  std::cout<<"Saving in "<<filename<<std::endl;
  TTree* out = new TTree(treename, treename);
  
  
  std::map<int, double> weightup;
  std::map<int, double> weightdown;
  std::map<int, bool> eventhas;
  for (std::map<int, bool>::iterator it=particles.begin(); it!=particles.end(); ++it) 
  {
    const int pdg = it->first;
    weightup[pdg] = 0;
    weightdown[pdg] = 0;
    eventhas[pdg] = false;
    out->Branch(Form("weight_up_%d", pdg), &weightup[pdg], Form("weight_up_%d/D", pdg));
    out->Branch(Form("weight_down_%d", pdg), &weightdown[pdg], Form("weight_down_%d/D", pdg));
    out->Branch(Form("event_has_%d", pdg), &eventhas[pdg], Form("event_has_%d/O", pdg));
  } // end for (std::map<int, bool>::iterator it=particles.begin(); it!=particles.end(); ++it) 
  
  // Loop over all entries
  Long64_t entries = (usetruth) ? fTruthChain->GetEntries() : fChain->GetEntries();
  for (Long64_t i = 0; i < entries; i++) 
  {
    out->GetEntry(i);

    InelXSecWeights weights = getWeights(i, usetruth);
    // Iterate over all particles
    for (std::map<int, bool>::iterator it=particles.begin(); it!=particles.end(); ++it) 
    {
      const int pdg = it->first;
      weightup[pdg] = weights.weightUp[pdg];
      weightdown[pdg] = weights.weightDown[pdg];
      eventhas[pdg] = weights.eventHas[pdg];
    }
    out->Fill();
  }
  
  file->Write();
  file->Close();
  delete file;
  
  return out;
} 

std::string MnvHadronReweight::makefullfilename(const std::string filename, const std::string projectname)
{
  // Version numbering lets me recreate the weights without 
  // having people remember to delete their cache files.
  std::string suffix;
  suffix = "_v10"; 
  if (XSecFunctions::useHDXSec) suffix += "_HD";
  if (doElasticReweight) suffix += "_withElastics"; 
  if (m_reweightNeutronCV) suffix += "_RwgtNeutronCV"; 
  
  std::string root = ".root";
  return filename + projectname + suffix + root;
}

bool MnvHadronReweight::tryLoadingFromFile(const char* filename, const char* projectname)
{
  if (filename == NULL) {
    std::cout<<"MnvHadronReweight::tryLoadingFromFile: filename was null."<<std::endl;
    return false;
  }
  if (projectname == NULL) {
    std::cout<<"MnvHadronReweight::tryLoadingFromFile: projectname was null."<<std::endl;
    return false;
  }
  std::string fullfilenamestring = makefullfilename(filename, projectname);
  const char* fullfilename = fullfilenamestring.c_str();
  
  // open root file if possible, return false if not
  // Get ttree, return false if not there
  // Get pdg, renorm up, renorm down tuples
  // Put in renorm map
  
  // Make a file
  TFile file(fullfilename);
  if (file.IsZombie()) 
  {
    std::cout<<"Unable to load "<<fullfilename<<std::endl;
    return false;
  }
  //set projectname 
  setProjectName(projectname);
  
  // Load ttree
  TTree* foo = (TTree*) file.Get("renorm");
  // Check if the ttree's there
  if (foo == NULL || foo->GetEntries() == 0) 
  {
    std::cout<<"TTree is missing in "<<fullfilename<<std::endl;
    return false;
  }
  
  int npdgs;
  TBranch        *b_npdgs;
  foo->SetBranchAddress("npdgs", &npdgs, &b_npdgs);
  foo->GetEntry(0);
  
  const int length = npdgs;
  
  // Holds the number of entries
  int nEntries;
  // Holds the pdg
  int pdg[length];
  // Holds the renorm amount for up
  double renormAmountUp[length];  
  double renormAmountDown[length];  
  TBranch        *b_nEntries;
  TBranch        *b_pdg;
  TBranch        *b_renormAmountUp; 
  TBranch        *b_renormAmountDown; 
  foo->SetBranchAddress("nEntries", &nEntries, &b_nEntries);
  foo->SetBranchAddress("pdg", &pdg, &b_pdg);
  foo->SetBranchAddress("renormAmountUp", &renormAmountUp, &b_renormAmountUp);
  foo->SetBranchAddress("renormAmountDown", &renormAmountDown, &b_renormAmountDown);
  foo->GetEntry(0);
  // Get all the info
  for (int i = 0; i < npdgs; i++) 
  {
    renormWeights.weightUp[pdg[i]] = renormAmountUp[i];
    renormWeights.weightDown[pdg[i]] = renormAmountDown[i];
  }
  
  return true; 
}

bool MnvHadronReweight::saveRenormFactorsToFile(const char* filename, const char* projectname)
{
  if (filename == NULL) {
    std::cout<<"MnvHadronReweight::saveRenormFactorsToFile: filename was null."<<std::endl;
    return false;
  }
  if (projectname == NULL) {
    std::cout<<"MnvHadronReweight::saveRenormFactorsToFile: projectname was null."<<std::endl;
    return false;
  }
  
  std::string fullfilenamestring = makefullfilename(filename, projectname);
  const char* fullfilename = fullfilenamestring.c_str();
  
  // Make a file
  TFile file(fullfilename, "RECREATE");
  if (file.IsZombie() || !file.IsWritable()) 
  {
    std::cout<<"Unable to write in "<<fullfilename<<std::endl;
    return false;
  }
  setProjectName(projectname);
  
  // We want to save on a per-projectname basis
  TTree foo("renorm", "renorm");
  foo.GetEntry(0);
  
  int length;
  if (!m_reweightNeutronCV) length = particles.size();
  else length = particles.size() + 1;
  // Holds the number of entries
  int nEntries;
  int npdgs = length;
  // Holds the pdg
  int pdg[length];
  // Holds the renorm amount for up
  double renormAmountUp[length];  
  double renormAmountDown[length];  
  foo.Branch("nEntries", &nEntries, "nEntries/I");
  foo.Branch("npdgs", &npdgs, "npdgs/I");
  foo.Branch("pdg", &pdg, "pdg[npdgs]/I");
  foo.Branch("renormAmountUp", &renormAmountUp, "renormAmountUp[npdgs]/D");
  foo.Branch("renormAmountDown", &renormAmountDown, "renormAmountDown[npdgs]/D");
  nEntries = fChain->GetEntries();
  int i = 0;
  for (std::map<int, bool>::iterator it=particles.begin(); it!=particles.end(); ++it) 
  {
    pdg[i] = it->first;
    renormAmountUp[i] = renormWeights.weightUp[pdg[i]];
    renormAmountDown[i] = renormWeights.weightDown[pdg[i]];
    i++;
  }
  if (m_reweightNeutronCV) 
  {
    pdg[i] = PDGNEUTRONCV;
    std::cout<<"Also saving neutron CV. PDG name = "<<pdg[i]<<std::endl;
    renormAmountUp[i] = renormWeights.weightUp[pdg[i]];
    renormAmountDown[i] = renormWeights.weightDown[pdg[i]];
    i++;
  }
  foo.Fill();
  foo.Write();
  file.Close();

  return true;
}
      
// Note: this uses TBranch technology. It's not doing too much harm.  May want
// to re-write in terms of chains or universes.
void MnvHadronReweight::getRenormFactors(const char* filename, const char* projectname, TTree* truth)
{
  // Are you trying to use two renorm factors at once?
  // Shoot out a message saying that's probably not a great idea
  if(renormKine) 
  {
    std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
    std::cout<<"MnvHadronReweighter WARNING: it seems that you're trying to use both the renormalization scale"<<std::endl;
    std::cout<<"  as well as the kinematic renormalizations.  While this is allowed, the default files assume "<<std::endl;
    std::cout<<"  the use of one and not the other."                                                           <<std::endl;
    std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
  }

  // TODO get renorm factors correctly
  if (tryLoadingFromFile(filename, projectname)) 
  {
    renorm = true;
    std::cout<<"    ***Loaded renorm factors from file:*** "<<std::endl;
  } // end if (tryLoadingFromFile(filename)) 
  else 
  {
    if (particles.size() == 0) 
    {
      std::cout<<"MnvHadronReweight error: There are no particles to calc the renom values to."<<std::endl;
      std::cout<<"Please specify at least one particle using MnvHadronReweight::setParticle(int pdg, bool treatAntiParticleSame)"<<std::endl;
      std::cout<<"Or use MnvHadronReweight::useDefaultParticles() to set pion, proton and neutron."<<std::endl;
      throw 66;
    }
    std::cout<<"    ***Creating renorm factors (this could take a while):*** "<<std::endl;

    if (truth == NULL) 
    {
      ; // TODO The TTree* truth argument should be removed altogether. It's
        // not a big ask to have the user set fTruth outside of this function
        // (else we throw). And usetruth is hardcoded false anyway.
    }
    
    HadronBranchContainer *containerToUse;
    bool usetruth = false; // We decided against using truth
    if (usetruth)
    {
      std::cout<<"Calculating renorm factors based on truth tree."<<std::endl;
      containerToUse = fTruth;
    }
    else
    {
      std::cout<<"Calculating renorm factors based on data and not truth tree."<<std::endl;
      containerToUse = fData; 
    }
    
    // Check if in default mode
    if (containerToUse->defaultmode) 
    {
      // Return a renorm amount with initial values of 1.0
      std::cout<<"Detected default mode. Returning default renorm values (1.0)"<<std::endl;
      initializeWeights(renormWeights, 1.0);
      return;
    } // end if (containerToUse->defaultmode) 

    // Get a container to hold the total wgt
    std::cout<<"Initializing sum container"<<std::endl;
    InelXSecWeights renormW;
    initializeWeights(renormW, 0.0);
    
    // Initialize a map to hold the actual numbers of entries
    std::cout<<"Initializing map of entries"<<std::endl;
    std::map<int, Long64_t> actual_entries;
    std::map<int, Long64_t> missing_entries;
    for (std::map<int, bool>::iterator it=particles.begin(); it!=particles.end(); ++it) 
    {
      const int pdg = it->first;
      actual_entries[pdg] = 0;
      missing_entries[pdg] = 0;
    } // end for (std::map<int, bool>::iterator it=particles.begin(); it!=particles.end(); ++it) 
    if (m_reweightNeutronCV)
    {
      actual_entries[PDGNEUTRONCV] = 0;
      missing_entries[PDGNEUTRONCV] = 0;
    }
    
    // Set the renormalization to false temporarily
    bool oldrenorm = renorm;
    renorm = false;
    
    // Loop over all entries in the tree
    std::cout<<"Starting loop"<<std::endl;
    for (Long64_t i = 0; i < containerToUse->fChain->GetEntriesFast(); i++)
    {
      // https://root.cern.ch/phpBB3/viewtopic.php?t=14946     

      Long64_t local_entry = containerToUse->fChain->LoadTree(i);
      containerToUse->fChain->GetEntry(i);
      // Long64_t local_entry = containerToUse->fChain->GetEntry(i);
      if (local_entry < 0) break;
      else if(local_entry < 1){
        std::cout<<"Bad Entry : " << i <<std::endl;              
        continue;
      }
      
      // std::cout<<"Entry Loaded : " << i <<std::endl;      
      if(i % 10000 == 0) std::cout<<(double(i) / 1000)<<"k "<<std::flush;
      if(i < 100) std::cout<<i<<" "<<std::flush;

      InelXSecWeights weights = getWeights( containerToUse );

      // std::cout<<"Got the Weights " <<std::endl;      

      // Loop over all particles
      for (std::map<int, bool>::iterator it=particles.begin(); it!=particles.end(); ++it)
      {
        // TODO how to do antiparticles?
        
        // Get the pdg
        const int pdg = it->first;

        //const bool treatAntiSame = it->second;
        // Check if this event has this type of particle
        if (weights.eventHas[pdg]) 
        {

          // Get the pdg
          const int pdg = it->first;
          // Add up the reweighted weights
          renormW.weightUp[pdg] += weights.weightUp[pdg];
          renormW.weightDown[pdg] += weights.weightDown[pdg];
          // Add up the actual entries
          actual_entries[pdg] += 1;
        } // end if (weights.eventHas[pdg])
        else missing_entries[pdg] += 1;
      } // end for (std::map<int, bool>::iterator it=particles.begin(); it!=particles.end(); ++it)
      
      // std::cout<<"Filled Norms " <<std::endl;      


      // Now do CV Neutron reweighting
      if (m_reweightNeutronCV)
      {
        // std::cout<<"Entering reweightNeutronCV" <<std::endl;      
        double wgt = reweightNeutronCV(i, usetruth);
        renormW.weightUp[PDGNEUTRONCV] += wgt;
        renormW.weightDown[PDGNEUTRONCV] += wgt;
        actual_entries[PDGNEUTRONCV] += 1;
      }
      // std::cout<<"End of sequence" <<std::endl;      
      
    } // end for (Long64_t i = 0; i < entries; i++)
    std::cout<<std::endl;
    std::cout<<"Finished loop"<<std::endl;

    // Reset the old renorm value
    //renorm = oldrenorm; 

    // Never understood why you'd build renorm factors and then proceed not to use them
    // Set so the user will use them
    renorm = true;
    
    std::cout<<"    ***Created renorm factors:*** "<<std::endl;
    for (std::map<int, bool>::iterator it=particles.begin(); it!=particles.end(); ++it) 
    {
      // Get the pdg
      const int pdg = it->first;
      double up = renormW.weightUp[pdg];
      double down = renormW.weightDown[pdg];
      if (up == 0 || down == 0 || actual_entries[pdg] == 0) 
      {
        std::cout<<"Error: Particle sum is zero. This could only happen if your sample doesn't contain this pdg or you're in default mode"<<std::endl;
        std::cout<<pdg<<" sum up:      "<<up<<std::endl;
        std::cout<<pdg<<" sum down:    "<<down<<std::endl;
        std::cout<<pdg<<" sum cv:      "<<actual_entries[pdg]<<std::endl;
        std::cout<<pdg<<" sum missing: "<<missing_entries[pdg]<<std::endl;
        std::cout<<"Setting a default value of zero for the renorm amount."<<std::endl;
        renormW.weightUp[pdg] = 0; // Set default values
        renormW.weightDown[pdg] = 0;
      }
      // Print some info
      std::cout<<pdg<<" sum up:      "<<up<<std::endl;
      std::cout<<pdg<<" sum down:    "<<down<<std::endl;
      std::cout<<pdg<<" sum cv:      "<<actual_entries[pdg]<<std::endl;
      std::cout<<pdg<<" sum missing: "<<missing_entries[pdg]<<std::endl;
      // Calculate the renorm amount needed
      renormW.weightUp[pdg] = actual_entries[pdg] / up;
      renormW.weightDown[pdg] = actual_entries[pdg] / down;
      // Print it
      std::cout<<pdg<<" pdg:      "<<pdg<<std::endl;
      std::cout<<pdg<<" up:    "<<renormW.weightUp[pdg]<<std::endl;
      std::cout<<pdg<<" down:  "<<renormW.weightDown[pdg]<<std::endl;
    } // end for (std::map<int, bool>::iterator it=particles.begin(); it!=particles.end(); ++it) 
    if (m_reweightNeutronCV)
    {
      std::cout<<"CV Neutron Renorm:"<<std::endl;
      int pdg = PDGNEUTRONCV;
      // Print some info
      double up = renormW.weightUp[pdg];
      double down = renormW.weightDown[pdg];
      std::cout<<" sum new CV:      "<<up<<std::endl;
      std::cout<<" sum old CV:      "<<actual_entries[pdg]<<std::endl;
      // Calculate the renorm amount needed
      renormW.weightUp[pdg] = actual_entries[pdg] / up;
      renormW.weightDown[pdg] = actual_entries[pdg] / down;
      // Print it
      std::cout<<" renorm factor:    "<<renormW.weightUp[pdg]<<std::endl;
    }
    
    // Set the final value
    renormWeights = renormW;
    
    // Reset the fake elastics
    ResetFakeElasticsSeed();
    
    // Finally save the values to the cache
    if (!saveRenormFactorsToFile(filename, projectname)) 
      std::cout<<"Failed saving renorm factors to cache file."<<std::endl;
    else std::cout<<"Saved renorm factors to cache file successfully"<<std::endl;
  }
}

bool MnvHadronReweight::tryLoadingKineFile()
{
  //I'm turning the HDXSec on permanently
  useHDXSec(true);
  std::string fullfilenamestring = makekinefilename();

  if( strcmp( fullfilenamestring.c_str(), m_kinefilename.c_str() ) == 0 ) 
  {  
    std::cout<<"Already loaded renorm kine from file"<<std::endl;
    return true;
  }
  m_kinefilename = fullfilenamestring;
  const char* fullfilename = fullfilenamestring.c_str();
  
  // open root file if possible, return false if not
  // Get ttree, return false if not there
  // Get pdg, renorm up, renorm down tuples
  // Put in renorm map
  
  // Make a file
  TFile file(fullfilename);
  if (file.IsZombie()) 
  {
    std::cout<<"Unable to load "<<fullfilename<<std::endl;
    return false;
  }
  std::cout<<"Loading file "<<fullfilename<<std::endl;

  //Passive targets, tracker
  for( int iTar = 0; iTar < mhrw_nTars; ++iTar )
  {
    std::string tar_name = getTargetNucleiName( mhrw_targets[iTar], mhrw_nuclei[iTar] );

    for( uint iShift = 0; iShift < 3; ++iShift )
    {
      for( auto &particle : particles ) 
      {
        int pdg = particle.first;

        std::string hist_name = Form("kine_rw_%s_p_theta_%d_%s",tar_name.c_str(),pdg,mhrw_shift_names[iShift]);
        TH2D* tmp_h = (TH2D*)file.Get(hist_name.c_str());
        if( tmp_h == NULL ) std::cout<<"Cannot load "<<hist_name<<" from kine file"<<std::endl;

        kineRenormWeights[mhrw_targets[iTar]][mhrw_nuclei[iTar]][pdg][iShift] = tmp_h;
      } 
      if (m_reweightNeutronCV)
      {
        std::string hist_name = Form("kine_rw_%s_p_theta_%d_%s",tar_name.c_str(),PDGNEUTRONCV,mhrw_shift_names[iShift]);
        TH2D* tmp_h = (TH2D*)file.Get(hist_name.c_str());
        if( tmp_h == NULL ) std::cout<<"Cannot load "<<hist_name<<" from kine file"<<std::endl;

        kineRenormWeights[mhrw_targets[iTar]][mhrw_nuclei[iTar]][PDGNEUTRONCV][iShift] = tmp_h;
      }
    }
  }
  //Scintillator in nuclear targets
  for( int iTar = 1; iTar <= mhrw_nNTR_Scint; ++iTar )
  {
    for( uint iShift = 0; iShift < 3; ++iShift )
    {
      for( auto &particle : particles ) 
      {
        int pdg = particle.first;
 
        std::string hist_name = Form("kine_rw_NTR_Scint%d_p_theta_%d_%s",iTar,pdg,mhrw_shift_names[iShift]); 
        TH2D* tmp_h = (TH2D*)file.Get(hist_name.c_str());
        if( tmp_h == NULL ) std::cout<<"Cannot load "<<hist_name<<" from kine file"<<std::endl;
 
        kineRenormWeights[100+iTar][6][pdg][iShift] = tmp_h;
      } 
      if (m_reweightNeutronCV)
      {
        std::string hist_name = Form("kine_rw_NTR_Scint%d_p_theta_%d_%s",iTar,PDGNEUTRONCV,mhrw_shift_names[iShift]); 
        TH2D* tmp_h = (TH2D*)file.Get(hist_name.c_str());
        if( tmp_h == NULL ) std::cout<<"Cannot load "<<hist_name<<" from kine file"<<std::endl;

        kineRenormWeights[100+iTar][6][PDGNEUTRONCV][iShift] = tmp_h;
      }
    }
  }

  //OD, ECAL/HCAL
  for( uint iShift = 0; iShift < 3; ++iShift )
  {
    for( auto &particle : particles ) 
    {
      int pdg = particle.first;
 
      std::string hist_name = Form("kine_rw_OD_p_theta_%d_%s",pdg,mhrw_shift_names[iShift]); 
      TH2D* tmp_h = (TH2D*)file.Get(hist_name.c_str());
      if( tmp_h == NULL ) std::cout<<"Cannot load "<<hist_name<<" from kine file"<<std::endl;
 
      kineRenormWeights[200][0][pdg][iShift] = tmp_h;

      ///////////////////////////////////////////////////////////////////////////
       
      hist_name = Form("kine_rw_EHCAL_p_theta_%d_%s",pdg,mhrw_shift_names[iShift]); 
      tmp_h = (TH2D*)file.Get(hist_name.c_str());
      if( tmp_h == NULL ) std::cout<<"Cannot load "<<hist_name<<" from kine file"<<std::endl;
 
      kineRenormWeights[300][0][pdg][iShift] = tmp_h;
    } 
    if (m_reweightNeutronCV)
    {
      std::string hist_name = Form("kine_rw_OD_p_theta_%d_%s",PDGNEUTRONCV,mhrw_shift_names[iShift]); 
      TH2D* tmp_h = (TH2D*)file.Get(hist_name.c_str());
      if( tmp_h == NULL ) std::cout<<"Cannot load "<<hist_name<<" from kine file"<<std::endl;

      kineRenormWeights[200][0][PDGNEUTRONCV][iShift] = tmp_h;

      ///////////////////////////////////////////////////////////////////////////

      hist_name = Form("kine_rw_EHCAL_p_theta_%d_%s",PDGNEUTRONCV,mhrw_shift_names[iShift]); 
      tmp_h = (TH2D*)file.Get(hist_name.c_str());
      if( tmp_h == NULL ) std::cout<<"Cannot load "<<hist_name<<" from kine file"<<std::endl;

      kineRenormWeights[300][0][PDGNEUTRONCV][iShift] = tmp_h;
    }
  }

  for( auto& tarmap : kineRenormWeights )
  {
    for( auto& nucmap : tarmap.second )
    {
      for( auto& partmap : nucmap.second )
      {
        for( auto& shiftmap : partmap.second )
        {
          shiftmap.second->SetDirectory(0);
          if(false) std::cout<<tarmap.first<<" "<<nucmap.first<<" "<<partmap.first<<" "<<shiftmap.first<<" "<<shiftmap.second->GetName()<<std::endl;
        }
      }
    }
  }

  return true; 
}

bool MnvHadronReweight::saveKineRenormFactorsToFile()
{
  std::string fullfilenamestring = makekinefilename();
  const char* fullfilename = fullfilenamestring.c_str();
  
  // Make a file
  TFile file(fullfilename, "RECREATE");
  if (file.IsZombie() || !file.IsWritable()) 
  {
    std::cout<<"Unable to write in "<<fullfilename<<std::endl;
    return false;
  }
  
  file.cd();
  m_kineRenormHists.Write();
  file.Close();

  return true;
}
      
std::string MnvHadronReweight::makekinefilename()
{
  std::string suffix;
  suffix = !m_playlist.empty() ? "_Minerva"+m_playlist : "";
  suffix += "_v10"; 
  if (XSecFunctions::useHDXSec) suffix += "_HD";
  if (doElasticReweight) suffix += "_withElastics"; 
  if (m_reweightNeutronCV) suffix += "_RwgtNeutronCV"; //Definitely needed
 
  if( fiducialType == 1 ) suffix += Form("_z%d-%d_a%d",(int)fiducialMinZ,(int)fiducialMaxZ,(int)fiducialApothem);
 
  std::string root = ".root";
  return m_directory + "/" + m_projectname + suffix + root;
}
      
std::string MnvHadronReweight::getTargetNucleiName( int tar, int nuc )
{
  std::string ret = "";
  if( tar == 6 )     ret = "tracker";
  else if( tar == 7) ret = "water"; 
  else               ret = Form("tar%d_nuc%d",tar, nuc);

  return ret;
}

std::map<int, LeadingFSParticle> MnvHadronReweight::getLeadingFSParticles( HadronBranchContainer* b )
{
         //pdg            idx  energy
  std::map<int, std::pair<int, double > > leadingPartIdxE;

  for( int iPart = 0; iPart < b->mc_nFSPart; ++iPart )
  {
    const int pdg = b->mc_FSPartPDG[iPart];
    bool bPart     = particles.find(pdg) != particles.end();
    bool bAntiPart = particles.find( getAntiPDG(pdg) ) != particles.end();

    // Is the particle or antiparticle in particles vector?
    if( !bPart && !bAntiPart ) continue;

    // If it is the antiparticle, do we care?
    if( bAntiPart && !particles[getAntiPDG(pdg)] ) continue;

    //Make the pdg we store this match those in particles vector
    const int part_pdg = bPart ? pdg : getAntiPDG(pdg);

    if( leadingPartIdxE.find(part_pdg) == leadingPartIdxE.end() )
    {
      leadingPartIdxE[part_pdg].first  = iPart;
      leadingPartIdxE[part_pdg].second = b->mc_FSPartE[iPart];
    }
    else
    {
      if( b->mc_FSPartE[iPart] > leadingPartIdxE[part_pdg].second )
      {
        leadingPartIdxE[part_pdg].first  = iPart;
        leadingPartIdxE[part_pdg].second = b->mc_FSPartE[iPart];
      }

    } 
  }

  std::map<int, LeadingFSParticle> FSParts;

  //Loop through map
  //Make TVector3, make magnitude and angle, import into FSParts
  for( auto &lPart : leadingPartIdxE )
  {
    const int pdg = lPart.first;
    const int fs_idx = lPart.second.first;    

    TVector3 p3FS( b->mc_FSPartPx[fs_idx], b->mc_FSPartPy[fs_idx], b->mc_FSPartPz[fs_idx] );
    p3FS.RotateX(MinervaUnits::numi_beam_angle_rad); //Do I actually need this if it's internal?

    FSParts[pdg].p     = p3FS.Mag();
    FSParts[pdg].theta = p3FS.Theta()*180/3.14159;
    //std::cout<<pdg<<" "<<FSParts[pdg].p<<" "<<FSParts[pdg].theta<<std::endl;
  }
  return FSParts;
}

void MnvHadronReweight::getTruthKineRenorm()
{
  // Are you trying to use two renorm factors at once?
  // Shoot out a message saying that's probably not a great idea
  if(renorm) 
  {
    std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
    std::cout<<"MnvHadronReweighter WARNING: it seems that you're trying to use both the renormalization scale"<<std::endl;
    std::cout<<"  as well as the kinematic renormalizations.  While this is allowed, the default files assume "<<std::endl;
    std::cout<<"  the use of one and not the other."                                                           <<std::endl;
    std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
  }
  if (tryLoadingKineFile()) 
  {
    renormKine = true;
    std::cout<<"    ***Loaded renorm factors from file:*** "<<std::endl;
    return;
  } // end if (tryLoadingFromFile(directory)) 

  if (particles.size() == 0) 
  {
    std::cout<<"MnvHadronReweight error: There are no particles to calc the kine renorm values to."<<std::endl;
    std::cout<<"Please specify at least one particle using MnvHadronReweight::setParticle(int pdg, bool treatAntiParticleSame)"<<std::endl;
    std::cout<<"Or use MnvHadronReweight::useDefaultParticles() to set pion, proton and neutron."<<std::endl;
    throw 66;
  }
  std::cout<<"    ***Creating kinematic renorm factors (this could take a while):*** "<<std::endl;
  
  //Check and warn the user if the name isn't "Truth"
  if( fData == NULL || std::string( fData->fChain->GetName() ).find("Truth") == std::string::npos )
  {
    std::cout<<"MnvHadronReweight error: you are trying to run over a tree that isn't labelled 'Truth'"<<std::endl;
    std::cout<<"  to create the kinematic renorm. The kinematic renorm weights assume truth tree info."<<std::endl;
    std::cout<<"  Please check your inputs."<<std::endl;
    if( fData == NULL )
    {
      std::cout<<"Note, for MnvHadronReweight(TTree* truth, TTree* data, bool hdXSec), the <truth> tree is"<<std::endl;
      std::cout<<"  depreciated when creating renorm weights.  Please put the truth tree in the <data> input"<<std::endl;
    }
    throw 67;
  }

  //Since the xsec effects the error bands, I'm removing this
  //Turn on neutron CV reweighting to fill the right histograms
  //bool old_rwgtNeutronCV = m_reweightNeutronCV;
  //useReweightedNeutronCV(); 

  HadronBranchContainer *containerToUse;
  containerToUse = fData;
  
  bool debug = false;
  bool usetruth = false; // This makes sure we don't use the *truth input
  // Check if in default mode
  if (containerToUse->defaultmode) 
  {
    // Return a renorm amount with initial values of 1.0
    std::cout<<"Detected default mode. Returning default kine renorm values (1.0)"<<std::endl;
    initKineRenormWeights(kineRenormWeights, 1.0);
    return;
  } // end if (containerToUse->defaultmode) 

  // Get a container to hold the total wgt
  std::cout<<"Initializing sum container"<<std::endl;

  TruthKineRenormWeights kineRenormW;
  initKineRenormWeights(kineRenormW, 0.0);

  // Set the renormalization to false temporarily
  renormKine = false;

  // Loop over all entries in the tree
  std::cout<<"Starting loop"<<std::endl;
  for (Long64_t i = 0; i < containerToUse->fChain->GetEntriesFast(); i++)
  {
    // https://root.cern.ch/phpBB3/viewtopic.php?t=14946     

    Long64_t local_entry = containerToUse->fChain->LoadTree(i);
    containerToUse->fChain->GetEntry(i);
    if (local_entry < 0) break;
    else if(local_entry < 1){
      std::cout<<"Bad Entry : " << i <<std::endl;              
      continue;
    }
  
    if(debug) std::cout<<containerToUse->mc_nFSPart<<std::endl;

    if(debug && i == 10000 ) break; 
    if(i % 10000 == 0) std::cout<<(double(i) / 1000)<<"k "<<std::flush;
    if(i < 100) std::cout<<i<<" "<<std::flush;

    //Get Targets
    std::pair<int, int> tarnuc = getTargetNuclei( containerToUse );

    int target = tarnuc.first;
    int nuclei = tarnuc.second;

    if( target < 0 ) 
    {  
      if(debug) std::cout<<"Can't find a target"<<std::endl;
      continue;
    }

    if(debug) std::cout<<target<<" "<<nuclei<<" "<<std::flush; 
    //Get weights, leading FS particles, where is the interaction vertex?
    InelXSecWeights weights = getWeights( containerToUse );
 
    std::map<int, LeadingFSParticle> fsParticles = getLeadingFSParticles( containerToUse );

    if(debug) std::cout<<"PDG "<<std::flush; 
    // Loop over all particles
    for (std::map<int, bool>::iterator it=particles.begin(); it!=particles.end(); ++it)
    {
      // Get the pdg.  This pdg will be written out (so antiparticles will not)
      const int pdg = it->first;

      // Check if this event has this type of particle
      if( weights.eventHas[pdg] || (it->second && weights.eventHas[getAntiPDG(pdg)]) ) 
      {
        if(debug) std::cout<<pdg<<" "<<std::flush; 
        kineRenormW[target][nuclei][pdg][UP]  ->Fill( fsParticles[pdg].p, fsParticles[pdg].theta, weights.weightUp[pdg] );
        kineRenormW[target][nuclei][pdg][DOWN]->Fill( fsParticles[pdg].p, fsParticles[pdg].theta, weights.weightDown[pdg] );
        kineRenormW[target][nuclei][pdg][CV]  ->Fill( fsParticles[pdg].p, fsParticles[pdg].theta, 1 );
      } // end if (weights.eventHas[pdg])
    } // end for (std::map<int, bool>::iterator it=particles.begin(); it!=particles.end(); ++it)
    
    // Now do CV Neutron reweighting
    //Check to see if there is a neutron in this event
    if (m_reweightNeutronCV && weights.eventHas[PDGNEUTRON] )
    {
      if(debug) std::cout<<"Entering reweightNeutronCV" <<std::endl;      

      double wgt = reweightNeutronCV(i, usetruth);
      kineRenormW[target][nuclei][PDGNEUTRONCV][UP]  ->Fill( fsParticles[PDGNEUTRON].p, fsParticles[PDGNEUTRON].theta, wgt );
      kineRenormW[target][nuclei][PDGNEUTRONCV][DOWN]->Fill( fsParticles[PDGNEUTRON].p, fsParticles[PDGNEUTRON].theta, wgt );
      kineRenormW[target][nuclei][PDGNEUTRONCV][CV]  ->Fill( fsParticles[PDGNEUTRON].p, fsParticles[PDGNEUTRON].theta, 1 );
    }
    if(debug) std::cout<<"End of sequence" <<std::endl;      
    
  } // end for (Long64_t i = 0; i < entries; i++)
  std::cout<<std::endl;
  std::cout<<"Finished loop"<<std::endl;

  //Set renormKine to true
  renormKine = true;
   
  std::cout<<"    ***Created renorm factors:*** "<<std::endl;

  for( auto& tarmap : kineRenormW )
  {
    int tar = tarmap.first;
    for( auto& nucmap : tarmap.second )
    {
      int nuc = nucmap.first;
      for( auto& partmap : nucmap.second ) //neutron CV should be included here
      {
        int pdg = partmap.first;
        
        std::cout<<"TAR "<<tar<<" NUC  "<<nuc<<" PDG "<<pdg<<std::endl;
    
        //Integral of pdg.  Check if one is 0
        double up   = partmap.second[UP]->Integral();
        double down = partmap.second[DOWN]->Integral();
        double cv   = partmap.second[CV]->Integral();
        if (up == 0 || down == 0 || cv == 0) 
        {
          std::cout<<"Error: integral of kinematic particle plots are zero. This could only happen if your sample doesn't contain this pdg or you're in default mode"<<std::endl;
        }
        // Print some info
        std::cout<<pdg<<" integral up:      "<<up<<std::endl;
        std::cout<<pdg<<" integral down:    "<<down<<std::endl;
        std::cout<<pdg<<" integral cv:      "<<cv<<std::endl;
      
        std::cout<<pdg<<" dividing and creating weights" <<std::endl;
        
        //Divide the shifts, so we can just look up the reweight later
        partmap.second[UP]->Divide(partmap.second[CV]); 
        partmap.second[DOWN]->Divide(partmap.second[CV]); 
      }
    }
  }
  // Set the final value
  kineRenormWeights = kineRenormW;
 
  // Reset the fake elastics
  ResetFakeElasticsSeed();

  //Since the xsec effects the error bands, I'm removing this
  //Set the neutron CV reweight  options back to the original 
  //setReweightedNeutronCV( old_rwgtNeutronCV );

  // Finally save the values to the cache
  if (!saveKineRenormFactorsToFile()) std::cout<<"Failed saving kine renorm factors to cache file."<<std::endl;
  else std::cout<<"Saved kine renorm factors to cache file successfully"<<std::endl;
}

void MnvHadronReweight::initKineRenormWeights( TruthKineRenormWeights& tkrw, int initialValue )
{
  const int nPBins = 80;
  double p_bins[nPBins+1];
  // 0-6000, 75 width
  for( int p = 0; p <= nPBins; ++p ) p_bins[p] = 75*p;

  const int nThetaBins = 11; 
  double theta_bins[nThetaBins+1] = { 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 90.0, 120.0, 150.0, 180.0};

  //Passive targets, tracker
  for( int iTar = 0; iTar < mhrw_nTars; ++iTar )
  {
    std::string tar_name = getTargetNucleiName( mhrw_targets[iTar], mhrw_nuclei[iTar] );

    for( uint iShift = 0; iShift < 3; ++iShift )
    {
      for( auto &particle : particles ) 
      {
        int pdg = particle.first;
        tkrw[mhrw_targets[iTar]][mhrw_nuclei[iTar]][pdg][iShift] = new TH2D(Form("kine_rw_%s_p_theta_%d_%s",tar_name.c_str(),pdg,mhrw_shift_names[iShift]),";p (MeV);#theta (deg)",nPBins,p_bins,nThetaBins,theta_bins); 
        //Makes things easier to write out
        m_kineRenormHists.Add(tkrw[mhrw_targets[iTar]][mhrw_nuclei[iTar]][pdg][iShift]);

        //Initialize
        if( initialValue >= 1e-6 )//If we don't set a initial value, assume we're filling these, so use the default of 0
        {
          for( int iBinX = 0; iBinX <= nPBins; ++iBinX )
          {
            for( int iBinY = 0; iBinY <= nThetaBins; ++iBinY )
            {
              tkrw[mhrw_targets[iTar]][mhrw_nuclei[iTar]][pdg][iShift]->SetBinContent( iBinX, iBinY, initialValue );
            } 
          } 
        } //End initialization
      } //End of particles
      
      //Neutron CV
      if (m_reweightNeutronCV)
      {
        tkrw[mhrw_targets[iTar]][mhrw_nuclei[iTar]][PDGNEUTRONCV][iShift] = new TH2D(Form("kine_rw_%s_p_theta_%d_%s",tar_name.c_str(),PDGNEUTRONCV,mhrw_shift_names[iShift]),";p (MeV);#theta (deg)",nPBins,p_bins,nThetaBins,theta_bins); 
        //Makes things easier to write out
        m_kineRenormHists.Add(tkrw[mhrw_targets[iTar]][mhrw_nuclei[iTar]][PDGNEUTRONCV][iShift]);

        //Initialize
        if( initialValue >= 1e-6 )//If we don't set a initial value, assume we're filling these, so use the default of 0
        {
          for( int iBinX = 0; iBinX <= nPBins; ++iBinX )
          {
            for( int iBinY = 0; iBinY <= nThetaBins; ++iBinY )
            {
              tkrw[mhrw_targets[iTar]][mhrw_nuclei[iTar]][PDGNEUTRONCV][iShift]->SetBinContent( iBinX, iBinY, initialValue );
            } 
          } 
        } //End initialization
      }//End of neutron CV

    }//End of shifts
  }//End of targets

  //Scintillator in nuclear targets
  for( int iTar = 1; iTar <= mhrw_nNTR_Scint; ++iTar )
  {
    for( uint iShift = 0; iShift < 3; ++iShift )
    {
      for( auto &particle : particles ) 
      {
        int pdg = particle.first;
        tkrw[100+iTar][6][pdg][iShift] = new TH2D(Form("kine_rw_NTR_Scint%d_p_theta_%d_%s",iTar,pdg,mhrw_shift_names[iShift]),";p (MeV);#theta (deg)",nPBins,p_bins,nThetaBins,theta_bins); 
        //Makes things easier to write out
        m_kineRenormHists.Add(tkrw[100+iTar][6][pdg][iShift]);

        //Initialize
        if( initialValue >= 1e-6 )//If we don't set a initial value, assume we're filling these, so use the default of 0
        {
          for( int iBinX = 0; iBinX <= nPBins; ++iBinX )
          {
            for( int iBinY = 0; iBinY <= nThetaBins; ++iBinY )
            {
              tkrw[100+iTar][6][pdg][iShift]->SetBinContent( iBinX, iBinY, initialValue );
            } 
          } 
        } //End initialization
      } //End of particles

      //Neutron CV
      if (m_reweightNeutronCV)
      {
        tkrw[100+iTar][6][PDGNEUTRONCV][iShift] = new TH2D(Form("kine_rw_NTR_Scint%d_p_theta_%d_%s",iTar,PDGNEUTRONCV,mhrw_shift_names[iShift]),";p (MeV);#theta (deg)",nPBins,p_bins,nThetaBins,theta_bins); 
        //Makes things easier to write out
        m_kineRenormHists.Add(tkrw[100+iTar][6][PDGNEUTRONCV][iShift]);

        //Initialize
        if( initialValue >= 1e-6 )//If we don't set a initial value, assume we're filling these, so use the default of 0
        {
          for( int iBinX = 0; iBinX <= nPBins; ++iBinX )
          {
            for( int iBinY = 0; iBinY <= nThetaBins; ++iBinY )
            {
              tkrw[100+iTar][6][PDGNEUTRONCV][iShift]->SetBinContent( iBinX, iBinY, initialValue );
            } 
          } 
        } //End initialization
      }//End of neutron CV

    }//End of shifts
  }//End of targets

  //OD, ECAL/HCAL
  for( uint iShift = 0; iShift < 3; ++iShift )
  {
    for( auto &particle : particles ) 
    {
      int pdg = particle.first;
      tkrw[200][0][pdg][iShift] = new TH2D(Form("kine_rw_OD_p_theta_%d_%s",pdg,mhrw_shift_names[iShift]),";p (MeV);#theta (deg)",nPBins,p_bins,nThetaBins,theta_bins); 
      tkrw[300][0][pdg][iShift] = new TH2D(Form("kine_rw_EHCAL_p_theta_%d_%s",pdg,mhrw_shift_names[iShift]),";p (MeV);#theta (deg)",nPBins,p_bins,nThetaBins,theta_bins); 
      m_kineRenormHists.Add(tkrw[200][0][pdg][iShift]);
      m_kineRenormHists.Add(tkrw[300][0][pdg][iShift]);

      //Initialize
      if( initialValue >= 1e-6 )//If we don't set a initial value, assume we're filling these, so use the default of 0
      {
        for( int iBinX = 0; iBinX <= nPBins; ++iBinX )
        {
          for( int iBinY = 0; iBinY <= nThetaBins; ++iBinY )
          {
            tkrw[200][0][pdg][iShift]->SetBinContent( iBinX, iBinY, initialValue );
            tkrw[300][0][pdg][iShift]->SetBinContent( iBinX, iBinY, initialValue );
          } 
        } 
      } //End initialization
    } //End of particles

    //Neutron CV
    if (m_reweightNeutronCV)
    {
      tkrw[200][0][PDGNEUTRONCV][iShift] = new TH2D(Form("kine_rw_OD_p_theta_%d_%s",PDGNEUTRONCV,mhrw_shift_names[iShift]),";p (MeV);#theta (deg)",nPBins,p_bins,nThetaBins,theta_bins); 
      tkrw[300][0][PDGNEUTRONCV][iShift] = new TH2D(Form("kine_rw_EHCAL_p_theta_%d_%s",PDGNEUTRONCV,mhrw_shift_names[iShift]),";p (MeV);#theta (deg)",nPBins,p_bins,nThetaBins,theta_bins); 
      m_kineRenormHists.Add(tkrw[200][0][PDGNEUTRONCV][iShift]);
      m_kineRenormHists.Add(tkrw[300][0][PDGNEUTRONCV][iShift]);

      //Initialize
      if( initialValue >= 1e-6 )//If we don't set a initial value, assume we're filling these, so use the default of 0
      {
        for( int iBinX = 0; iBinX <= nPBins; ++iBinX )
        {
          for( int iBinY = 0; iBinY <= nThetaBins; ++iBinY )
          {
            tkrw[200][0][PDGNEUTRONCV][iShift]->SetBinContent( iBinX, iBinY, initialValue );
            tkrw[300][0][PDGNEUTRONCV][iShift]->SetBinContent( iBinX, iBinY, initialValue );
          } 
        } 
      } //End initialization
    }//End of neutron CV

  }//End of shifts
}

void MnvHadronReweight::initializeWeights(InelXSecWeights& weights, double initialValue)
{
  // Loop over particles
  for (std::map<int, bool>::iterator it=particles.begin(); it!=particles.end(); ++it) {
    int pdg = it->first;
    // Initialize weights
    weights.weightUp[pdg] = initialValue;
    weights.weightDown[pdg] = initialValue;
    weights.eventHas[pdg] = false;

    // Check whether to treat the antiparticle as the same particle.
    // Like pi+ and pi-
    if (it->second) {
      int antiPDG = getAntiPDG(pdg);
      weights.weightUp[antiPDG] = initialValue;
      weights.weightDown[antiPDG] = initialValue;
      weights.eventHas[antiPDG] = false;
    }
  }
  
  if (m_reweightNeutronCV)
  {
    weights.weightUp[PDGNEUTRONCV] = initialValue;
    weights.weightDown[PDGNEUTRONCV] = initialValue;
  }
}

double MnvHadronReweight::getMass(const int pdg) 
{
  // Gets the masses in MeV.
  // Assuming anti-pdg is same mass
  switch (getNormalPDG(pdg)) 
  {
    case PDGPION:  return pion_mass;
    case PDGPROTON: return proton_mass;
    case PDGKAON:  return kaon_mass;
    case PDGNEUTRON: return neutron_mass;
    default: std::cout<<"Do not recognize PDG "<<pdg<<". Please add this pdg to the 'MnvHadronReweight::getMass' function."<<std::endl; throw 66;
  }
}

double MnvHadronReweight::getCorrectedDensity(const int nuke, const double density)
{
  // density is 
  // g4material->GetDensity() / (Gaudi::Units::g / Gaudi::Units::cm3) (units of g/cm^3)
  // times path_length / Gaudi::Units::cm; (units of cm)
  // So density has units of g/cm^2
  // I need molecules / cm^2
  // This is not the best way to do this. Ideally this would be calc'd when
  // saving the density, so you can get a difference between air and carbon.
  // But density of nitrogen ~ carbon so it's not too different.
  double density_correction = 0;
  switch (nuke) // The correction is 1/mass(nuke)
  {
    case -6: density_correction = 4.626e22; break;
    case  6: density_correction = 5.01398e22; break;
    case 26: density_correction = 1.0784e22; break;
    case 82: density_correction = 2.906e21; break;
    case 13: density_correction = 2.2319486e22; break; // aluminium
    case  8: density_correction = 3.3428e22; break; // Water
    default: std::cout<<"The number of molecules per g for this nuke is not known. nuke="<<nuke<<std::endl; throw 66; 
  }
  return density * density_correction;
}

bool MnvHadronReweight::getWeightWithElastic(const int nuke, const int pdg, const int intCode, const double density, const double Ei, const double Ef, const double delta_inelastic, const double delta_elastic, double& up, double& down)
{
  // Check for trivial case here
  if (density == 0) return false; // Weight = 1.0 in all cases
  // Get the mass in MeV
  const double mass = getMass(pdg);

  // Get the kinetic E in MeV.
  double Ti = Ei - mass;
  double Tf = Ef - mass;
  if (Ti < 0 || Tf < 0) 
  { 
    if (Ef == 0) return false; // This doesn't make sense but can happen bc of genie
    if (Ti <= 0 && Tf <= 0) return false; // Don't care.
    if (Ti > 0 && Tf < 0) Tf = 0.0; // Round Tf to zero
  }
  
  // Changes density from g/cm^2 to molecules/cm^2
  const double newdensity = getCorrectedDensity(nuke, density);
  
  double rhoX = newdensity;
  //double sigma_total_geant = getAverageCrossSection2(pdg, nuke, Ti, Tf, xsectype_total);
  double sigma_elastic_geant = getAverageCrossSection2(pdg, nuke, Ti, Tf, xsectype_elastic);
  double sigma_inelastic_geant = getAverageCrossSection2(pdg, nuke, Ti, Tf, xsectype_inelastic);
  double sigma_total_geant = sigma_elastic_geant + sigma_inelastic_geant;
  if (sigma_total_geant <= 0 || sigma_elastic_geant <= 0 || sigma_inelastic_geant <= 0) return false;
  
  double sigma_total_data_up = sigma_elastic_geant * (1 + delta_elastic) + sigma_inelastic_geant * (1 + delta_inelastic);
  double sigma_total_data_down = sigma_elastic_geant * (1 - delta_elastic) + sigma_inelastic_geant * (1 - delta_inelastic);
  
  if (intCode != 0 && intCode != 2 && intCode != 3) // Inelastic_Interacted, intcode=1,4
  {
    double sigma_partial_geant = sigma_inelastic_geant;
    double sigma_partial_data_up = sigma_partial_geant * (1 + delta_inelastic);
    if (sigma_partial_data_up < 0) sigma_partial_data_up = 0;
    double sigma_partial_data_down = sigma_partial_geant * (1 - delta_inelastic);
    if (sigma_partial_data_down < 0) sigma_partial_data_down = 0;
    if (sigma_total_data_up <= 0) up = 0;
    else up = doWeightCalcInteracted(rhoX, sigma_partial_data_up, sigma_partial_geant, sigma_total_data_up, sigma_total_geant);
    if (sigma_total_data_down <= 0) down = 0;
    else down = doWeightCalcInteracted(rhoX, sigma_partial_data_down, sigma_partial_geant, sigma_total_data_down, sigma_total_geant);
  } // end if (intCode != 0 && intCode != 2)
  else if (intCode == 3) // Elastic, intcode=3
  {
    double sigma_partial_geant = sigma_elastic_geant;
    double sigma_partial_data_up = sigma_partial_geant * (1 + delta_elastic);
    if (sigma_partial_data_up < 0) sigma_partial_data_up = 0;
    double sigma_partial_data_down = sigma_partial_geant * (1 - delta_elastic);
    if (sigma_partial_data_down < 0) sigma_partial_data_down = 0;
    if (sigma_total_data_up <= 0) up = 0;
    else up = doWeightCalcInteracted(rhoX, sigma_partial_data_up, sigma_partial_geant, sigma_total_data_up, sigma_total_geant);
    if (sigma_total_data_down <= 0) down = 0;
    else down = doWeightCalcInteracted(rhoX, sigma_partial_data_down, sigma_partial_geant, sigma_total_data_down, sigma_total_geant);
    
  }
  else // Inelastic_NotInteracted, intcode=0,2
  {
    up = doWeightCalcNoninteracted(rhoX, sigma_total_data_up, sigma_total_geant);
    down = doWeightCalcNoninteracted(rhoX, sigma_total_data_down, sigma_total_geant);
  }
  
  //std::cout<<"Up="<<up<<", Down="<<down<<", rhoX="<<rhoX<<", newdensity="<<newdensity<<", sigma_total_geant="<<sigma_total_geant<<", sigma_elastic_geant="<<sigma_elastic_geant<<", sigma_inelastic_geant="<<sigma_inelastic_geant<<", sigma_total_data_up="<<sigma_total_data_up<<", sigma_total_data_down="<<sigma_total_data_down<<", Ti="<<Ti<<", Tf="<<Tf<<", nuke="<<nuke<<", pdg="<<pdg<<", int_code="<<intCode<<", mass="<<mass<<", delta_inelastic="<<delta_inelastic<<", delta_elastic="<<delta_elastic<<", sigma_elastic_data_up="<<sigma_elastic_geant * (1 + delta_elastic)<<", sigma_elastic_data_down="<<sigma_elastic_geant * (1 - delta_elastic)<<", sigma_total_data_up==sigma_total_geant="<<(sigma_total_geant == sigma_total_data_up)<<std::endl; // This code is simply there to print values for manual calculation
  
  
  if (std::isnan(up) || std::isnan(down))
  {
    std::cout<<"Got NaN. Up="<<up<<", Down="<<down<<", rhoX="<<rhoX<<", newdensity="<<newdensity<<", sigma_total_geant="<<sigma_total_geant<<", sigma_elastic_geant="<<sigma_elastic_geant<<", sigma_inelastic_geant="<<sigma_inelastic_geant<<", Ti="<<Ti<<", Tf="<<Tf<<", nuke="<<nuke<<", pdg="<<pdg<<", int_code="<<intCode<<", mass="<<mass<<", delta_inelastic="<<delta_inelastic<<", delta_elastic="<<delta_elastic<<std::endl;
    throw 66;
  }
  return true;
}

bool MnvHadronReweight::getWeightToNewNeutronXSec(int nuke, int intCode, double density, double Ei, double Ef, double& wgt)
{
  // Check for trivial case here
  if (density == 0) return false; // Weight = 1.0 in all cases
  if (nuke != -6) return false; // Don't support anything but scint right now
  // Get the mass in MeV
  const double mass = getMass(PDGNEUTRON);

  // Get the kinetic E in MeV.
  double Ti = Ei - mass;
  double Tf = Ef - mass;
  if (Ti < 0 || Tf < 0) 
  { 
    if (Ef == 0) return false; // This doesn't make sense but can happen bc of genie
    if (Ti <= 0 && Tf <= 0) return false; // Don't care.
    if (Ti > 0 && Tf < 0) Tf = 0.0; // Round Tf to zero
  }
  
  // Changes density from g/cm^2 to molecules/cm^2
  const double newdensity = getCorrectedDensity(nuke, density);
  
  double rhoX = newdensity;
  double sigma_elastic_geant = getAverageCrossSection2(PDGNEUTRON, nuke, Ti, Tf, xsectype_elastic, true); // isCVCalc tells it that it's this calc and it wants the old XSec
  double sigma_inelastic_geant = getAverageCrossSection2(PDGNEUTRON, nuke, Ti, Tf, xsectype_inelastic, true);
  double sigma_total_geant = sigma_elastic_geant + sigma_inelastic_geant;
  if (sigma_total_geant <= 0 || sigma_elastic_geant <= 0 || sigma_inelastic_geant <= 0) return false;
  
  double sigma_elastic_data_new = getAverageNewNeutronCrossSection(nuke, Ti, Tf, xsectype_elastic);
  double sigma_inelastic_data_new = getAverageNewNeutronCrossSection(nuke, Ti, Tf, xsectype_inelastic);
  double sigma_total_data_new = sigma_elastic_data_new + sigma_inelastic_data_new;
  if (sigma_total_data_new <= 0 || sigma_elastic_data_new <= 0 || sigma_inelastic_data_new <= 0) return false;
  
  if (intCode != 0 && intCode != 2 && intCode != 3) // Inelastic_Interacted, intcode=1,4
  {
    wgt = doWeightCalcInteracted(rhoX, sigma_inelastic_data_new, sigma_inelastic_geant, sigma_total_data_new, sigma_total_geant);
  } // end if (intCode != 0 && intCode != 2)
  else if (intCode == 3) // Elastic, intcode=3
  {
    wgt = doWeightCalcInteracted(rhoX, sigma_elastic_data_new, sigma_elastic_geant, sigma_total_data_new, sigma_total_geant);
    
  }
  else // Inelastic_NotInteracted, intcode=0,2
  {
    wgt = doWeightCalcNoninteracted(rhoX, sigma_total_data_new, sigma_total_geant);
  }
  
  
  if (std::isnan(wgt))
  {
    std::cout<<"Got NaN. Wgt="<<wgt<<", rhoX="<<rhoX<<", newdensity="<<newdensity<<", sigma_total_geant="<<sigma_total_geant<<", sigma_elastic_geant="<<sigma_elastic_geant<<", sigma_inelastic_geant="<<sigma_inelastic_geant<<", Ti="<<Ti<<", Tf="<<Tf<<", nuke="<<nuke<<", int_code="<<intCode<<", mass="<<mass<<std::endl;
    throw 66;
  }
  return true;
}

bool MnvHadronReweight::getWeight(const int nuke, const int pdg, const int intCode, const double density, const double Ei, const double Ef, const double delta, double& up, double& down)
{
  // Check for trivial case here
  if (density == 0) return false; // Weight = 1.0 in all cases
  // Get the mass in MeV
  const double mass = getMass(pdg);

  // Get the kinetic E in MeV.
  double Ti = Ei - mass;
  double Tf = Ef - mass;
  if (Ti < 0 || Tf < 0) 
  { 
    if (Ef == 0) return false; // This doesn't make sense but can happen bc of genie
    if (Ti <= 0 && Tf <= 0) return false; // Don't care.
    if (Ti > 0 && Tf < 0) Tf = 0.0; // Round Tf to zero
  }
  
  // Changes density from g/cm^2 to molecules/cm^2
  const double newdensity = getCorrectedDensity(nuke, density);
  
  double rhoX = newdensity;
  double sigma_inelastic_geant = getAverageCrossSection2(pdg, nuke, Ti, Tf, xsectype_inelastic);
  if (sigma_inelastic_geant <= 0) return false;
  
  double sigma_inelastic_data_up = sigma_inelastic_geant * (1 + delta);
  double sigma_inelastic_data_down = sigma_inelastic_geant * (1 - delta);
  
  if (intCode != 0 && intCode != 2 && intCode != 3) // Inelastic_Interacted
  {
    up = doWeightCalcInelasticOnly(rhoX, sigma_inelastic_data_up, sigma_inelastic_geant);
    down = doWeightCalcInelasticOnly(rhoX, sigma_inelastic_data_down, sigma_inelastic_geant);
  } // end if (intCode != 0 && intCode != 2)
  else // Inelastic_NotInteracted
  {
    up = doWeightCalcNoninteracted(rhoX, sigma_inelastic_data_up, sigma_inelastic_geant);
    down = doWeightCalcNoninteracted(rhoX, sigma_inelastic_data_down, sigma_inelastic_geant);
  }
  
  
  //std::cout<<"Up="<<up<<", Down="<<down<<", rhoX="<<rhoX<<", newdensity="<<newdensity<<", sigma_total_geant="<<sigma_total_geant<<", sigma_elastic_geant="<<sigma_elastic_geant<<", sigma_inelastic_geant="<<sigma_inelastic_geant<<", Ti="<<Ti<<", Tf="<<Tf<<", nuke="<<nuke<<", pdg="<<pdg<<", int_code="<<intCode<<", mass="<<mass<<", delta_inelastic="<<delta_inelastic<<", delta_elastic="<<delta_elastic<<std::endl; // TODO Remove. JDK
  
  
  if (std::isnan(up) || std::isnan(down))
  {
    std::cout<<"Got NaN in getWeight. Up="<<up<<", Down="<<down<<", rhoX="<<rhoX<<", newdensity="<<newdensity<<", sigma_inelastic_geant="<<sigma_inelastic_geant<<", Ti="<<Ti<<", Tf="<<Tf<<", nuke="<<nuke<<", pdg="<<pdg<<", int_code="<<intCode<<", mass="<<mass<<", delta="<<delta<<std::endl;
    throw 66;
  }
  return true;
}

double MnvHadronReweight::doWeightCalcInteracted(double rhoX, double sigma_partial_data, double sigma_partial_geant, double sigma_total_data, double sigma_total_geant)
{
  double denom = (1.0 - exp(-1.0*rhoX*sigma_total_geant));
  if (denom <= 0) return 0.0;
  double nom = (1.0 - exp(-1.0*rhoX*sigma_total_data));
  double a = sigma_partial_data / sigma_total_data;
  double b = sigma_partial_geant / sigma_total_geant;
  return (nom / denom) * a / b;
}

double MnvHadronReweight::doWeightCalcInelasticOnly(double rhoX, double sigma_inelastic_data, double sigma_inelastic_geant)
{
  double denom = (1.0 - exp(-1.0*rhoX*sigma_inelastic_geant));
  if (denom <= 0) return 0.0;
  double nom = (1.0 - exp(-1.0*rhoX*sigma_inelastic_data));
  return (nom / denom);
}

double MnvHadronReweight::doWeightCalcNoninteracted(double rhoX, double sigma_total_data, double sigma_total_geant)
{
  return exp(-1.0*rhoX*(sigma_total_data - sigma_total_geant));
}

void MnvHadronReweight::setReadoutVolume( std::string volname )
{
  if( strcmp( volname.c_str(), "Tracker" ) == 0 || strcmp( volname.c_str(), "tracker" ) == 0 )
  {
    setBasicFiducial( NSFDefaults::TrackerFace, NSFDefaults::TrackerBack, NSFDefaults::StandardApothem );
    std::cout<<"Setting readout volume for Tracker: minZ "<<fiducialMinZ<<" maxZ "<<fiducialMaxZ<<" "<<fiducialApothem<<std::endl;
  } 
  else if( strcmp( volname.c_str(), "Nuke" ) == 0 || strcmp( volname.c_str(), "nuke" ) == 0 )
  {
    setBasicFiducial( NSFDefaults::nuclearTargetZFace, NSFDefaults::TrackerBack, NSFDefaults::StandardApothem );
    std::cout<<"Setting readout volume for Nuke: minZ "<<fiducialMinZ<<" maxZ "<<fiducialMaxZ<<" "<<fiducialApothem<<std::endl;
  }
  else
  {
    if( volname.size() > 0 ) 
    {
      std::cout<<volname<<" is not at standard readout volume!"<<std::endl;
      std::cout<<" Current options are \"Tracker\" or \"Nuke\""<<std::endl;
    }
    else
    {
      std::cout<<"No readout volume specified"<<std::endl;
      std::cout<<"Set the volume with setBasicFiducial( double minZ, double maxZ, double apothem )"<<std::endl;
    }
  }
}

void MnvHadronReweight::setBasicFiducial(double minZ, double maxZ, double apothem)
{
  // minZ, maxZ, and apothem in mm
  // Sets this basic type of fiducial volume cut
  fiducialType = 1;
  fiducialMinZ = minZ;
  fiducialMaxZ = maxZ;
  fiducialApothem = apothem;
}

void MnvHadronReweight::getBasicFiducial(double& minZ, double& maxZ, double& apothem)
{
  minZ = fiducialMinZ;
  maxZ = fiducialMaxZ;
  apothem = fiducialApothem;
}

void MnvHadronReweight::setFiducialVolumeType(int type)
{
  if (type == 1)
  {
    // It makes no sense to use this function because the max/min z and apothem haven't been set
    std::cout<<"MnvHadronReweight::setFiducialVolumeType: For type 1 you should set the parameters using setBasicFiducial(double minZ, double maxZ, double apothem)"<<std::endl;
    throw 66;
  }
  fiducialType = type;
}

// Branch Version
bool MnvHadronReweight::inFiducial(HadronBranchContainer* b, int iPrimTraj) 
{
  const double X = b->truth_hadronReweightPosX[iPrimTraj];
  const double Y = b->truth_hadronReweightPosY[iPrimTraj];
  const double Z = b->truth_hadronReweightPosZ[iPrimTraj];
  return inFiducial(X, Y, Z);
}

// Generic calculator
bool MnvHadronReweight::inFiducial(const double X, const double Y,
                                   const double Z) {
  // Check which method of fiducial volume checking to use.
  // The default is 
  if (fiducialType == 1)
  {
    if (Z < fiducialMinZ || Z > fiducialMaxZ) return false;
    return isInHexagon(X, Y, fiducialApothem);
  } // end if (fiducialType == 1)
  else if (fiducialType == -1) {
    std::cout<<"MnvHadronReweight::inFiducial: No fiducial volume was set"<<std::endl;
    throw 66;
  }
  else 
  {
    std::cout<<"MnvHadronReweight::inFiducial: Fiducial volume type "<<fiducialType<<" is not defined. Throwing error."<<std::endl;
    throw 66;
  }
}

std::pair<int, int> MnvHadronReweight::getTargetNuclei( HadronBranchContainer* b ) 
{
  int nuc = b->mc_targetZ;

  const double vtx_x = b->mc_vtx[0];
  const double vtx_y = b->mc_vtx[1];
  const double vtx_z = b->mc_vtx[2];
  //  std::cout<<vtx_x<<", "<<vtx_y<<", "<<vtx_z << "\t" << b->mc_nFSPart<<std::endl;

  int tar = -1;

  bool debug = false;
  if(debug) std::cout<<" "<<vtx_z<<" "<<nuc<<" "<<std::flush;
  //assuming apothem = 850

  //Check to see if it in hexagon.  we have a separate class
  if( !isInHexagon( vtx_x, vtx_y, 850. ) )
  {
    if(debug) std::cout<<"Not in hexagon "<<std::flush;
    if(debug && !inFiducial( vtx_x, vtx_y, vtx_z) ) std::cout<<"Not in fid vol "<<std::flush;

    tar = 200;
    nuc = 0;
    std::pair<int, int> tar_nuc = std::make_pair(tar,nuc);
    return tar_nuc;
  }

  if( vtx_z < PlotUtils::TargetProp::Tracker::Face )
  {
    //if(      m_TargetUtils->InTarget1MC( vtx_x, vtx_y, vtx_z, nuc ) ) tar = 1;
    //else if( m_TargetUtils->InTarget2MC( vtx_x, vtx_y, vtx_z, nuc ) ) tar = 2;
    //else if( m_TargetUtils->InTarget3MC( vtx_x, vtx_y, vtx_z, nuc ) ) tar = 3;
    //else if( m_TargetUtils->InTarget4MC( vtx_x, vtx_y, vtx_z, nuc ) ) tar = 4;
    //else if( m_TargetUtils->InTarget5MC( vtx_x, vtx_y, vtx_z, nuc ) ) tar = 5;
    //else if( m_TargetUtils->InWaterTargetMC( vtx_x, vtx_y, vtx_z, nuc ) ){ tar = 7; nuc = 8; }//Water is O or H, but we can only have 1 integer
    if(      m_TargetUtils->InIron1VolMC(   vtx_x, vtx_y, vtx_z ) ){ tar = 1; nuc = 26; }
    else if( m_TargetUtils->InLead1VolMC(   vtx_x, vtx_y, vtx_z ) ){ tar = 1; nuc = 82; }
    else if( m_TargetUtils->InIron2VolMC(   vtx_x, vtx_y, vtx_z ) ){ tar = 2; nuc = 26; }
    else if( m_TargetUtils->InLead2VolMC(   vtx_x, vtx_y, vtx_z ) ){ tar = 2; nuc = 82; }
    else if( m_TargetUtils->InCarbon3VolMC( vtx_x, vtx_y, vtx_z ) ){ tar = 3; nuc = 6 ; }
    else if( m_TargetUtils->InIron3VolMC(   vtx_x, vtx_y, vtx_z ) ){ tar = 3; nuc = 26; }
    else if( m_TargetUtils->InLead3VolMC(   vtx_x, vtx_y, vtx_z ) ){ tar = 3; nuc = 82; }
    else if( m_TargetUtils->InLead4VolMC(   vtx_x, vtx_y, vtx_z ) ){ tar = 4; nuc = 82; }
    else if( m_TargetUtils->InIron5VolMC(   vtx_x, vtx_y, vtx_z ) ){ tar = 5; nuc = 26; }
    else if( m_TargetUtils->InLead5VolMC(   vtx_x, vtx_y, vtx_z ) ){ tar = 5; nuc = 82; }
    else if( m_TargetUtils->InWaterTargetVolMC( vtx_x, vtx_y, vtx_z ) ){ tar = 7; nuc = 8; }//Water is O or H, but we can only have 1 integer
    if( tar > 0 ) 
    {
      std::pair<int, int> tar_nuc = std::make_pair(tar,nuc);
      return tar_nuc;
    }
    if(debug) std::cout<<"Active "<<std::flush;
    //Which active region are we in?
    if( vtx_z < ( m_TargetUtils->GetTarget1CenterZMC() - PlotUtils::TargetProp::ThicknessMC::Tgt1::Pb/2 ) )      { tar = 101; nuc = 6; }
    else if( ( m_TargetUtils->GetTarget1CenterZMC() + PlotUtils::TargetProp::ThicknessMC::Tgt1::Pb/2 ) < vtx_z &&
             vtx_z < ( m_TargetUtils->GetTarget2CenterZMC() - PlotUtils::TargetProp::ThicknessMC::Tgt2::Pb/2 ) ) { tar = 102; nuc = 6; }
    else if( ( m_TargetUtils->GetTarget2CenterZMC() + PlotUtils::TargetProp::ThicknessMC::Tgt2::Pb/2 ) < vtx_z &&
             vtx_z < ( m_TargetUtils->GetTarget3CenterZMC() - PlotUtils::TargetProp::ThicknessMC::Tgt3::Pb/2 ) ) { tar = 103; nuc = 6; }
    else if( ( m_TargetUtils->GetTarget3CenterZMC() + PlotUtils::TargetProp::ThicknessMC::Tgt3::Pb/2 ) < vtx_z &&
             vtx_z < PlotUtils::TargetProp::WaterTarget::Face )                                                  { tar = 104; nuc = 6; }
    else if( PlotUtils::TargetProp::WaterTarget::Back  < vtx_z &&
             vtx_z < ( m_TargetUtils->GetTarget4CenterZMC() - PlotUtils::TargetProp::ThicknessMC::Tgt4::Pb/2 ) ) { tar = 105; nuc = 6; }
    else if( ( m_TargetUtils->GetTarget4CenterZMC() + PlotUtils::TargetProp::ThicknessMC::Tgt4::Pb/2 ) < vtx_z &&
             vtx_z < ( m_TargetUtils->GetTarget5CenterZMC() - PlotUtils::TargetProp::ThicknessMC::Tgt5::Pb/2 ) ) { tar = 106; nuc = 6; }
    else if( ( m_TargetUtils->GetTarget5CenterZMC() + PlotUtils::TargetProp::ThicknessMC::Tgt5::Pb/2 ) < vtx_z &&
             vtx_z < PlotUtils::TargetProp::Tracker::Face  )                                                     { tar = 107; nuc = 6; }
  }
  else if( PlotUtils::TargetProp::Tracker::Face < vtx_z && vtx_z < PlotUtils::TargetProp::Tracker::Back )
  {
    if( m_TargetUtils->InTracker( vtx_x, vtx_y, vtx_z ) ) { tar = 6; nuc = 6; }
    if(debug && !inFiducial( vtx_x, vtx_y, vtx_z) ) std::cout<<"Not in tracker fid vol "<<std::flush;
  }
  else//ECAL, HCAL
  {
    tar = 300; nuc = 0;  
  } 
 
  std::pair<int, int> tar_nuc = std::make_pair(tar,nuc);
  return tar_nuc;
}

bool MnvHadronReweight::isInHexagon( double x, double y, double apothem )
{
  // This code is taken from NukeCCQE code
  double lenOfSide = apothem*(2/sqrt(3)); 
  double slope     = (lenOfSide/2.0)/apothem;
  double xp        = fabs(x);
  double yp        = fabs(y);

  if( (xp*xp + yp*yp) < apothem*apothem )             return true;
  else if( xp <= apothem && yp*yp < lenOfSide/2.0 )   return true; 
  else if( xp <= apothem && yp < lenOfSide-xp*slope ) return true;

  return false;
}

double MnvHadronReweight::NeutronReweightAmount(double neutronKE)
{
  // This function was created by Kevin way back and represents the uncertainty as a function of KE.
  return std::max(0.0,std::min(150.0,-150 + neutronKE)/1500.) + std::min(0.25,0.1 + (3*std::max(0.0,50 - neutronKE))/500.);
}

void MnvHadronReweight::setDeltaScale(double newscale) { deltascale = newscale; }

// Branch version
double MnvHadronReweight::getDelta(HadronBranchContainer* b, int iPrimTraj) 
{
  const int index = b->truth_hadronReweightTrackID[iPrimTraj];
  const int pdg = getNormalPDG(b->truth_hadronReweightPDG[index]);
  const int target = b->truth_hadronReweightNuke[iPrimTraj];
  const double ke = b->truth_hadronReweightFinalE[iPrimTraj] - getMass(pdg);
  return MnvHadronReweight::getDelta(pdg, target, ke);
}

// Generic calculator
double MnvHadronReweight::getDelta(const int pdg, const int target, const double ke) 
{
  if (m_reweightNeutronCV && pdg == PDGNEUTRON && target == -6) {
    if (ke > 200) return 0.02;
    if (ke > 100) return 0.03;
    if (ke > 25)
      return 0.15;  // Rik moves to 0.15 // 3% uncertainty above 25 MeV.
                    // Hydrogen is still not known well
    if (ke > 10) return 0.20;  // Rik moves to 0.20 10% uncertainty above 10 MeV
    return 0.25;               // 25% uncertainty below 10 MeV
  }
  if (target == -6 || target == 6) {
    if (pdg == PDGPION || pdg == PDGPROTON)
      return 0.1;
    else if (pdg == PDGNEUTRON) {
      return NeutronReweightAmount(ke);
    }
  } else if (target == 26) {
    return 0.15;  // 15%
  } else if (target == 82) {
    return 0.20;  // 20%
  } else if (target == 13)
    return 0;  // aluminium
  else if (target == 8)
    return 0.15;  // Water, don't have uncertainty so using number of iron
  // This shouldn't happen, must have forgotten to include numbers above
  std::cout << "MnvHadronReweight This shouldn't happen!! :( " << pdg << ", "
            << target << std::endl;
  throw 66;
  return 0;
}

double MnvHadronReweight::getDeltaElastic(const int pdg, const int nuke, const double ke)
{
  /*
  Max errors from docdb 16284:
    Carbon > 80 MeV: 20%
    Carbon > 30 MeV: 40%
    Carbon > 10 MeV: 25%
    Carbon > 8 MeV: 40%
    Carbon > 5 MeV: 70%
    Carbon other: 40%
    Iron > 400 MeV: 20%
    Iron > 160 MeV: 25%
    Iron other: 15%
    Lead > 110 MeV: 15%
    Lead > 30 MeV: 25%
    Lead other: 10%
  */
  
  // Special case for reweighting to new geant4. See email from 10/26/2017
  if (m_reweightNeutronCV && pdg == PDGNEUTRON && nuke == -6) 
  {
    if (ke > 25) return 0.03; // 3% uncertainty above 25 MeV. Hydrogen is still not known well
    if (ke > 10) return 0.10; // 10% uncertainty above 10 MeV
    return 0.25; // 25% uncertainty below 10 MeV
  }
  
  double foo = pdg + ke; // Remove warnings for now
  foo += 1;
  if (pdg == PDGNEUTRON)
  {
    if (nuke == 82) // Lead
    {
      if (ke > 110) return 0.15;
      else if (ke > 30) return 0.25;
      else return 0.10;
    }
    else if (nuke == 26) // Iron
    {
      if (ke > 400) return 0.20;
      else if (ke > 160) return 0.25;
      else return 0.15;
    }
    else // Scint/carbon
    {
      if (ke > 80) return 0.20;
      else if (ke > 30) return 0.40;
      else if (ke > 10) return 0.25;
      else if (ke > 8) return 0.40;
      else if (ke > 5) return 0.70;
      else return 0.40;
    }
  }
  // TODO add water
  return 0.8; // Defaults to 80% error. TODO pick proton/pion delta
  
}

double MnvHadronReweight::ezrand()
{
  float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
  return r;
}

void MnvHadronReweight::changeLowEElasticXSec(double multiplier)
{
  XSecFunctions::changeLowEElasticXSec = multiplier;
}

InelXSecWeights MnvHadronReweight::getWeights( HadronBranchContainer* b )
{
  // Will hold the weights
  InelXSecWeights weights;
  initializeWeights(weights, 1.0); // init to 1
  // If in default mode, return a bunch of 1.0s
  if (b->defaultmode || defaultmode) return weights;
  
  // std::cout << "Checking b->truth_hadronReweightNPoints";// << std::endl;

  // Double check that the number of points is within the size of the array
  // Or else you'll get a seg fault
  if (b->truth_hadronReweightNPoints > AL) 
  {
    std::cout<<"MnvHadronReweight::getWeights:FATAL: Not enough room to load branches: "<<b->truth_hadronReweightNPoints<<" < Array length = "<<AL<<std::endl;
    std::cout<<"MnvHadronReweight::getWeights:FATAL: Go inside MnvHadronReweight and increase the size of the const variable AL"<<std::endl;
    throw 66;
  } // end if (b->truth_hadronReweightNPoints > AL)
  
  // std::cout << " Good. Entering iPrimTraj loop..." << std::endl;  

  int prevElasticTrackIndex = -1;
  
  // Check that all the hadron elements are the right size:
  // std::cout << "     ReweightNPoints = " << b->truth_hadronReweightNPoints << std::endl;
  // std::cout << "  ColumnarDensity_sz = " << b->truth_hadronReweightColumnarDensity_sz << std::endl;
  // std::cout << "          TrackID_sz = " << b->truth_hadronReweightTrackID_sz << std::endl;
  // std::cout << "IntCodePerSegment_sz = " << b->truth_hadronReweightIntCodePerSegment_sz << std::endl;
  // std::cout << "             Nuke_sz = " << b->truth_hadronReweightNuke_sz << std::endl;
  // std::cout << "           FinalE_sz = " << b->truth_hadronReweightFinalE_sz << std::endl;
  // std::cout << "         InitialE_sz = " << b->truth_hadronReweightInitialE_sz << std::endl;

  // Loop over all the points
  for ( int iPrimTraj = 0; iPrimTraj < b->truth_hadronReweightNPoints; ++iPrimTraj )
  {
    // std::cout << "Traj: " << iPrimTraj << ":"<< b->truth_hadronReweightNPoints << std::endl;

    // std::cout << "  ColumnarDensity[" << iPrimTraj << "] = " << b->truth_hadronReweightColumnarDensity[iPrimTraj] << std::endl;
    // std::cout << "          TrackID[" << iPrimTraj << "] = " << b->truth_hadronReweightTrackID[iPrimTraj] << std::endl;
    // std::cout << "IntCodePerSegment[" << iPrimTraj << "] = " << b->truth_hadronReweightIntCodePerSegment[iPrimTraj] << std::endl;
    // std::cout << "             Nuke[" << iPrimTraj << "] = " << b->truth_hadronReweightNuke[iPrimTraj] << std::endl;
    // std::cout << "           FinalE[" << iPrimTraj << "] = " << b->truth_hadronReweightFinalE[iPrimTraj] << std::endl;
    // std::cout << "         InitialE[" << iPrimTraj << "] = " << b->truth_hadronReweightInitialE[iPrimTraj] << std::endl;

    // Get the density and make sure it's not zero
    const double density     = b->truth_hadronReweightColumnarDensity[iPrimTraj]; 
    if ( density == 0 ) continue;
    
    // Get the index of the PDG and Interation code of the track
    // which is saved on a track-by-track basis instead of point-by-point
    const int index = b->truth_hadronReweightTrackID[iPrimTraj];

    // std::cout << "index == " << index << std::endl;

    if(index < 0 || index >= b->truth_hadronReweightPDG_sz)
      continue;

    // This makes it act like inelastic, where we only weight based on the first part of the track.
    //if (index == prevElasticTrackIndex) continue;
    
    // Now check if this segment is where the interaction occurs.
    // To do this, check the next track ID. 
    // If it is different, then this is the last segment of the track.
    // If it is the same, then there is another segment in this track.
    // If iPrimTraj + 1 < b->truth_hadronReweightNPoints, then it's 
    // the same answer as if the track IDs were different.
    /*bool isLastSegmentOfTrack;
    if (iPrimTraj + 1 >= b->truth_hadronReweightNPoints) isLastSegmentOfTrack = true;
    else if (index != b->truth_hadronReweightTrackID[iPrimTraj+1]) isLastSegmentOfTrack = true;
    else isLastSegmentOfTrack = false; // Next track segment ID == this ID */
    
    // Now get the interaction code
    //const int intCode    = b->truth_hadronReweightIntCode[index];
    /*int intCodeTemp    = b->truth_hadronReweightIntCode[index];  
    if (!isLastSegmentOfTrack) intCodeTemp = 0; // Change the int code based on this
    const int intCode    = intCodeTemp;*/
    
    // We no longer need to look for the last segment of the track
    // because we are saving int codes per segment. This is done because
    // elastics happen during trajectories instead of at the end only like
    // inelastics, stopping and decays.
    const int intCode    = b->truth_hadronReweightIntCodePerSegment[iPrimTraj];
  
    // std::cout << "b->truth_hadronReweightPDG[" << index << "] = " << b->truth_hadronReweightPDG[index] << std::endl;
    // std::cout<<"index="<<index<<", isLastSegment="<<isLastSegmentOfTrack<<", intCode="<<intCode<<", true_int_code="<<b->truth_hadronReweightIntCode[index]<<", iPrimTraj="<<iPrimTraj<<std::endl; // TODO remove. JDK
    
    // Get the pdg and make sure we care about that pdg
    const int pdg = b->truth_hadronReweightPDG[index];

    const int antiPDG = getAntiPDG(pdg);
    bool treatAntiTemp = false;
    // Check if particle is in a particle we care about
    if (particles.find(pdg) == particles.end()) 
    {
      // Check if anti-particle is something we care about
      if (particles.find(antiPDG) == particles.end()) continue;
      // Now we have to check if the anti-pdg is treated the same as the pdg
      else if (!particles[antiPDG]) continue;
      else treatAntiTemp = particles[antiPDG];
    } // end if (particles.find(pdg) == particles.end()) 
    else treatAntiTemp = particles[pdg];
    const bool treatAnti = treatAntiTemp; // how we want to treat the antiparticles
    
    // Get the material (scintillator=-6, pure carbon [target 3]=6, lead=82, steel=26, aluminum=13)
    const int nuke       = b->truth_hadronReweightNuke[iPrimTraj];  
    // Get the final and initial energy
    const double Ef      = b->truth_hadronReweightFinalE[iPrimTraj];  
    const double Ei      = b->truth_hadronReweightInitialE[iPrimTraj];  
    
    double minE = 0.0;
    if (abs(pdg) == 211) minE = 141.8;
    if (abs(pdg) == 2212) minE = 943.6;
    if (abs(pdg) == 2112) minE = 951.5;
    // skip low E traj, also ignored in ana utils code
    if (Ef < minE && Ei < minE) continue; 

    // Mark that this event has this PDG
    weights.eventHas[pdg] = true;
    // Mark that the antiparticle has this antiPDG, if needed
    if (treatAnti) weights.eventHas[antiPDG] = true;
    
    // Check the fiducial volume
    // Condider what to do what to do with non-fiducial
    // and what to do with the partial bit when it exits the fiducial.
    if (!inFiducial(b, iPrimTraj)) continue;
    
    // Now initialize the weights
    double up = 1.0;
    double down = 1.0;
    
    // Gets the amount of variation to use for this segment
    // It's based on the pdg, target type and whatever someone chooses
    // For example, the neutron is also based on KE
    double delta_inelastic = getDelta(b, iPrimTraj) * deltascale;
    bool success;
    if ((pdg == PDGNEUTRON && doElasticReweight) || elasticReweightForAllParticles) 
    {
      // Skip any elastics after the first elastic in traj 
      // bc we don't trust that the time matching isn't adding the same
      // elastic event to multiple segments.
      if (intCode == 3 && index == prevElasticTrackIndex) continue; 
      // Get the elastic delta
      const double mass = getMass(pdg);
      double ke = Ei - mass;
      double delta_elastic = getDeltaElastic(pdg, nuke, ke); 
      
      // Check if there was an elastic interaction that wasn't saved.
      // To do this, look at the final E of this segment and the E of the next segment.
      // If they are more different than minimum_missing_elastic_threshold
      // then there was an elastic interaction that wasn't saved.
      // Only check for missing elastic if
      // a) It's not the last segment of the traj.
      // b) You've asked to look for missing elastics
      /*if (!isLastSegmentOfTrack && look_for_missing_elastics && iPrimTraj+1 < b->truth_hadronReweightNPoints && index != b->truth_hadronReweightTrackID[iPrimTraj+1])
      {
        double nextE = b->truth_hadronReweightInitialE[iPrimTraj+1];
        if (nextE > 0)
        {
            double diff = Ef - nextE;
            double ratio = diff / nextE;
            if (ratio > minimum_missing_elastic_threshold) fakeelastic = true;
        }
      }*/
      
      bool fakeelastic = false;
      if (look_for_fake_elastics && false)
      {
        // Finding missing elastics is difficult. So we'll randomly create ones in segments instead
        double fake_segments_density = getCorrectedDensity(nuke, density);
        double Ti = Ei - mass;
        double Tf = Ef - mass;
        double elastic_xsec = getAverageCrossSection2(pdg, nuke, Ti, Tf, xsectype_elastic);
        double total_xsec = getAverageCrossSection2(pdg, nuke, Ti, Tf, xsectype_total);
        double prob;
        //if (total_xsec <= 0.0) prob = 0;
        //else 
        prob = (1 - exp(-fake_segments_density * total_xsec)) * (elastic_xsec / total_xsec);
        double throw_ = ezrand();
        if (throw_ < prob) fakeelastic = true; 
        //std::cout<<"Throw: "<<throw_<<", prob: "<<prob<<", density: "<<density<<", Ti="<<Ti<<", Tf="<<Tf<<", xsec_total: "<<total_xsec<<", elastic_xsec: "<<elastic_xsec;
        //std::cout<<", (elastic_xsec / total_xsec)="<<(elastic_xsec / total_xsec);
      }
      
      if (fakeelastic && false) // This is off
      {
        nfakeelastics += 1;
        double lerp = ezrand();
        // This has a fake elastic
        double up1 = 1.0;
        double down1 = 1.0;
        bool success1 = getWeightWithElastic(nuke, pdg, 3, density * lerp, Ei, Ei*(1-lerp)+Ef*lerp, delta_inelastic, delta_elastic, up1, down1);
        // Now the saved interaction
        double up2 = 1.0;
        double down2 = 1.0;
        bool success2 = getWeightWithElastic(nuke, pdg, intCode, density * (1-lerp), Ei*lerp+Ef*(1-lerp), Ef, delta_inelastic, delta_elastic, up2, down2);
        success = success1 && success2;
        up = up1 * up2;
        down = down1 * down2;
        //std::cout<<"Up1: "<<up1<<", down1: "<<down1<<", up2: "<<up1<<", down2: "<<down2<<", up: "<<up<<", down: "<<down<<", success: "<<success;
      }
      else
      {
        // Get the actual weights
        success = getWeightWithElastic(nuke, pdg, intCode, density, Ei, Ef, delta_inelastic, delta_elastic, up, down);
      }
      //std::cout<<std::endl;
      
    }
    else
    {
      // Default to previous weights where elastic isn't treated
      success = getWeight(nuke, pdg, intCode, density, Ei, Ef, delta_inelastic, up, down);
    }
    // Check if it was successful
    if (success) 
    {
      // Set the weight
      weights.weightUp[pdg] *= up;
      // Set the weight of the antiparticle, if needed
      if (treatAnti) weights.weightUp[antiPDG] *= up;
      // Set the weight
      weights.weightDown[pdg] *= down;
      // Set the weight of the antiparticle, if needed
      if (treatAnti) weights.weightDown[antiPDG] *= down;
      
      // This makes it act like the inelastic trajectories, where we only measure the first part.
      // Or this can be used to skip all elastics after the first elastic. The assumption is
      // the time matching alg doesn't work well enough.
      // TODO revisit
      if (intCode == 3) prevElasticTrackIndex = index; 
      
    } // end if (success)
  } // end for ( int iPrimTraj = 0; iPrimTraj < b->truth_hadronReweightNPoints; ++iPrimTraj )

  // Renorm if needed
  if ( renorm )
  {
    // Loop over particles
    for (std::map<int, bool>::iterator it=particles.begin(); it!=particles.end(); ++it)
    {
      // TODO check the logic here. If we treat antiparticles right
      
      // Get the pdg
      int pdg = it->first;
      // Skip this pdg if this event doesn't have that particle
      if (! weights.eventHas[pdg]) continue;
      // Renormalize the weights
      weights.weightUp[pdg] *= renormWeights.weightUp[pdg];
      weights.weightDown[pdg] *= renormWeights.weightDown[pdg];
      // Check whether to treat the antiparticle as the same particle.
      // Like pi+ and pi-
      if (it->second) 
      {
        // Get the anti particle PDG
        int antiPDG = getAntiPDG(pdg);
        // We want to skip renorm of the anti-pdg if it's also in the particle list
        // to avoid renormalizing twice
        if (particles.find(antiPDG) == particles.end()) 
        {
          // Renormalize
          weights.weightUp[antiPDG] *= renormWeights.weightUp[antiPDG];
          weights.weightDown[antiPDG] *= renormWeights.weightDown[antiPDG];
        } // end if (particles.find(antiPDG) == particles.end()) 
      } // end if (it->second) 
    } // end for (std::map<char,int>::iterator it=particles.begin(); it!=particles.end(); ++it)
  } // end if (renorm)
  
  //Kinematics renormalization
  //  std::cout<<"myrenorm status " << renormKine<<std::endl;
  if( renormKine )
  {
    //Get Targets
    std::pair<int, int> tarnuc = getTargetNuclei( b );

    int target = tarnuc.first;
    int nuclei = tarnuc.second;
    //std::cout<<nuclei<<std::endl;
    //    std::cout << "My target and nuclei are " << target << "\t" << nuclei << std::endl;
    if( target > 0 ) 
    {  
      std::map<int, LeadingFSParticle> fsParticles = getLeadingFSParticles( b );
      //      std::cout<< "number of FS particles = " << fsParticles.size()<<std::endl;

      for ( auto & fsParticle : fsParticles )
      {
        int pdg = fsParticle.first;
	//	std::cout << "MY PDG " << pdg << std::endl;
        if( !weights.eventHas[pdg] ) continue;        

        int binNum = kineRenormWeights[target][nuclei][pdg][CV]->FindBin( fsParticle.second.p, fsParticle.second.theta );
	//	std::cout << "MY entries " << kineRenormWeights[target][nuclei][pdg][CV]->GetBinContent(binNum) << std::endl;
        if( kineRenormWeights[target][nuclei][pdg][CV]->GetBinContent(binNum) > minKineEntries ) //reweight only on bin with enough events to make a godo reweight
        {
          weights.weightUp[pdg]   /= kineRenormWeights[target][nuclei][pdg][UP]->GetBinContent(binNum);
          weights.weightDown[pdg] /= kineRenormWeights[target][nuclei][pdg][DOWN]->GetBinContent(binNum);
        }
      }
    }
  }

  // Setting the current weights which can be retrieved later
  m_current_weights = weights;
  
  /*double largenum = 100;
  if (weights.weightUp[211] > largenum || weights.weightDown[211] > largenum || weights.weightUp[2212] > largenum || weights.weightDown[2212] > largenum || weights.weightUp[2112] > largenum || weights.weightDown[2112] > largenum)
  {
    std::cout<<"Entry="<<entry<<std::endl;
    std::cout<<"MnvHadronReweight::getWeights:WARNING: Weights big: Pion, up="<<weights.weightUp[211]<<", down="<<weights.weightDown[211]<<std::endl;
    std::cout<<"MnvHadronReweight::getWeights:WARNING: Weights big: Proton, up="<<weights.weightUp[2212]<<", down="<<weights.weightDown[2212]<<std::endl;
    std::cout<<"MnvHadronReweight::getWeights:WARNING: Weights big: Neutron, up="<<weights.weightUp[2112]<<", down="<<weights.weightDown[2112]<<std::endl;
  }*/ // Debugging code

  return weights;
}

// Branch Version
InelXSecWeights MnvHadronReweight::getWeights(Long64_t entry, bool usetruth)
{
  //Check if weights have been calculated for this entry
  if( entry == m_current_entry_weight ) return m_current_weights;
  else m_current_entry_weight = entry;
  //  std::cout << entry << "\t" << m_current_entry_weight << std::endl;
  HadronBranchContainer* b;
  if (usetruth) b = fTruth;
  else b = fData;

  // The above is a pointer to the same TChain in the loop... Why is this 
  // called again?
  // Get the right entry
  b->GetEntry(entry);


  Int_t local_entry = b->GetEntry(entry);
  // std::cout << "getWeights: local_entry = " << local_entry << std::endl;
  // if(local_entry < 1){
  //   std::cout<<"Bad Entry : " <<std::endl;
  // }

  return getWeights( b );
}

double MnvHadronReweight::reweightNeutronCV( HadronBranchContainer *b )
{
  if (!m_reweightNeutronCV) return 1.0;
  double out = 1.0;

  // If in default mode, return a bunch of 1.0s
  if (b->defaultmode || defaultmode) return 1.0;
  
  // Double check that the number of points is within the size of the array
  // Or else you'll get a seg fault
  if (b->truth_hadronReweightNPoints > AL) 
  {
    std::cout<<"MnvHadronReweight::reweightNeutronCV:FATAL: Not enough room to load branches: "<<b->truth_hadronReweightNPoints<<" < Array length = "<<AL<<std::endl;
    std::cout<<"MnvHadronReweight::reweightNeutronCV:FATAL: Go inside MnvHadronReweight and increase the size of the const variable AL"<<std::endl;
    throw 66;
  } // end if (b->truth_hadronReweightNPoints > AL)
  
  int prevElasticTrackIndex = -1;
  
  // Loop over all the points
  for ( int iPrimTraj = 0; iPrimTraj < b->truth_hadronReweightNPoints; ++iPrimTraj )
  {
    
    // Get the density and make sure it's not zero
    const double density     = b->truth_hadronReweightColumnarDensity[iPrimTraj]; 
    if ( density == 0 ) continue;
    
    // Get the index of the PDG and Interation code of the track
    // which is saved on a track-by-track basis instead of point-by-point
    const int index = b->truth_hadronReweightTrackID[iPrimTraj];

    
    // We no longer need to look for the last segment of the track
    // because we are saving int codes per segment. This is done because
    // elastics happen during trajectories instead of at the end only like
    // inelastics, stopping and decays.
    const int intCode    = b->truth_hadronReweightIntCodePerSegment[iPrimTraj];
    
    // Get the pdg and make sure we care about that pdg
    const int pdg = b->truth_hadronReweightPDG[index];
    if (pdg != PDGNEUTRON) continue; // Only neutron in this function

    // Get the material (scintillator=-6, pure carbon [target 3]=6, lead=82, steel=26, aluminum=13)
    const int nuke       = b->truth_hadronReweightNuke[iPrimTraj]; 
    if (nuke != -6) continue; // Only scintillator in this function 
    // Get the final and initial energy
    const double Ef      = b->truth_hadronReweightFinalE[iPrimTraj];  
    const double Ei      = b->truth_hadronReweightInitialE[iPrimTraj];  
    
    double minE = 0.0;
    if (abs(pdg) == 2112) minE = 951.5;
    // skip low E traj, also ignored in ana utils code
    if (Ef < minE && Ei < minE) continue; 

    // Check the fiducial volume
    // Condider what to do what to do with non-fiducial
    // and what to do with the partial bit when it exits the fiducial.
    if (!inFiducial(b, iPrimTraj)) continue;
    
    // Now initialize the weights
    double wgt = 1.0;
    
    // Gets the amount of variation to use for this segment
    // It's based on the pdg, target type and whatever someone chooses
    // For example, the neutron is also based on KE
    bool success;
    
    // Skip any elastics after the first elastic in traj 
    // bc we don't trust that the time matching isn't adding the same
    // elastic event to multiple segments.
    if (intCode == 3 && index == prevElasticTrackIndex) continue; 
      
    // Get the actual weights
    success = getWeightToNewNeutronXSec(nuke, intCode, density, Ei, Ef, wgt);
      
    // Check if it was successful
    if (success) 
    {
      // Set the weight
      out *= wgt;
      
      // This makes it act like the inelastic trajectories, where we only measure the first part.
      // Or this can be used to skip all elastics after the first elastic. The assumption is
      // the time matching alg doesn't work well enough.
      // TODO revisit
      if (intCode == 3) prevElasticTrackIndex = index; 
      
    } // end if (success)
  } // end for ( int iPrimTraj = 0; iPrimTraj < b->truth_hadronReweightNPoints; ++iPrimTraj )

  // Renorm if needed
  if ( renorm )
  {
    out *= renormWeights.weightUp[PDGNEUTRONCV];
  } // end if (renorm)

  //Kinematics renormalization
  if( renormKine )
  {
    std::map<int, LeadingFSParticle> fsParticles = getLeadingFSParticles( b );

    if( fsParticles.find(PDGNEUTRON) != fsParticles.end() ) //There's a neutron here
    {
      //Get Targets
      std::pair<int, int> tarnuc = getTargetNuclei( b );

      int target = tarnuc.first;
      int nuclei = tarnuc.second;

      if( target > 0 ) 
      {  
        int binNum = kineRenormWeights[target][nuclei][PDGNEUTRONCV][CV]->FindBin( fsParticles[PDGNEUTRON].p, fsParticles[PDGNEUTRON].theta );
        if( kineRenormWeights[target][nuclei][PDGNEUTRONCV][CV]->GetBinContent(binNum) > minKineEntries ) //reweight only on bin with enough events to make a godo reweight
        {
          out /= kineRenormWeights[target][nuclei][PDGNEUTRONCV][UP]->GetBinContent(binNum);
        }
      }
    }
  }
  
  // Setting the current reweight which can be retrieved later
  m_current_rwgt_neutronCV = out;
  return out;
}

double MnvHadronReweight::reweightNeutronCV(Long64_t entry, bool usetruth)
{
  if (!m_reweightNeutronCV) return 1.0;
  
  //Check if weights have been calculated for this entry
  if( entry == m_current_entry_neutronCV ) return m_current_rwgt_neutronCV;
  else m_current_entry_neutronCV = entry;

  // Get the branch containing all event data
  HadronBranchContainer* b;
  if (usetruth) b = fTruth;
  else b = fData;
  // Get the right entry
  b->GetEntry(entry);
 
  return reweightNeutronCV( b );
}

inline double safeEval(TF1 function, double Ti, double Tf)
{

  if (Ti == Tf) return function.Eval(Ti);
  else return function.Integral(Ti, Tf) / (Tf - Ti);
}

bool MnvHadronReweight::getXSecFunc(int pdg, int nuke, MnvHadronReweight::XSecType xsectype, TF1& function)
{
  if (pdg == -211) pdg = 211;
  if (nuke == 6 || nuke == -6)
  { 
    if (xsectype == xsectype_inelastic)
    {
      if (pdg == PDGPION) function = XSecFunctions::carbon_pion_inelastic;
      else if (pdg == PDGPROTON) function = XSecFunctions::carbon_proton_inelastic;
      else if (pdg == PDGNEUTRON) 
      {
        if (nuke == -6) function = XSecFunctions::scint_neutron_inelastic;
        else function = XSecFunctions::carbon_neutron_inelastic;
      }
      else if (pdg == PDGKAON) return false;
      else return false;
    }
    else if (xsectype == xsectype_elastic)
    {
      if (pdg == PDGPION) function = XSecFunctions::carbon_pion_elastic;
      else if (pdg == PDGPROTON) function = XSecFunctions::carbon_proton_elastic;
      else if (pdg == PDGNEUTRON) 
      {
        if (nuke == -6) function = XSecFunctions::scint_neutron_elastic;
        else function = XSecFunctions::carbon_neutron_elastic;
      }
      else if (pdg == PDGKAON) return false;
      else return false;
    }
    else if (xsectype == xsectype_total)
    {
      if (pdg == PDGPION) function = XSecFunctions::carbon_pion_total;
      else if (pdg == PDGPROTON) function = XSecFunctions::carbon_proton_total;
      else if (pdg == PDGNEUTRON)
      {
        if (nuke == -6) function = XSecFunctions::scint_neutron_total;
        else function = XSecFunctions::carbon_neutron_total;
      }
      else if (pdg == PDGKAON) return false;
      else return false;
    }
  }
  // TODO add water
  else if (nuke == 26)
  {
    if (xsectype == xsectype_inelastic)
    {
      if (pdg == PDGPION) function = XSecFunctions::iron_pion_inelastic;
      else if (pdg == PDGPROTON) function = XSecFunctions::iron_proton_inelastic;
      else if (pdg == PDGNEUTRON) function = XSecFunctions::iron_neutron_inelastic;
      else if (pdg == PDGKAON) return false;
      else return false;
    }
    else if (xsectype == xsectype_elastic)
    {
      if (pdg == PDGPION) function = XSecFunctions::iron_pion_elastic;
      else if (pdg == PDGPROTON) function = XSecFunctions::iron_proton_elastic;
      else if (pdg == PDGNEUTRON) function = XSecFunctions::iron_neutron_elastic;
      else if (pdg == PDGKAON) return false;
      else return false;
    }
    else if (xsectype == xsectype_total)
    {
      if (pdg == PDGPION) function = XSecFunctions::iron_pion_total;
      else if (pdg == PDGPROTON) function = XSecFunctions::iron_proton_total;
      else if (pdg == PDGNEUTRON) function = XSecFunctions::iron_neutron_total;
      else if (pdg == PDGKAON) return false;
      else return false;
    }
  }
  else if (nuke == 82)
  {
    if (xsectype == xsectype_inelastic)
    {
      if (pdg == PDGPION) function = XSecFunctions::lead_pion_inelastic;
      else if (pdg == PDGPROTON) function = XSecFunctions::lead_proton_inelastic;
      else if (pdg == PDGNEUTRON) function = XSecFunctions::lead_neutron_inelastic;
      else if (pdg == PDGKAON) return false;
      else return false;
    }
    else if (xsectype == xsectype_elastic)
    {
      if (pdg == PDGPION) function = XSecFunctions::lead_pion_elastic;
      else if (pdg == PDGPROTON) function = XSecFunctions::lead_proton_elastic;
      else if (pdg == PDGNEUTRON) function = XSecFunctions::lead_neutron_elastic;
      else if (pdg == PDGKAON) return false;
      else return false;
    }
    else if (xsectype == xsectype_total)
    {
      if (pdg == PDGPION) function = XSecFunctions::lead_pion_total;
      else if (pdg == PDGPROTON) function = XSecFunctions::lead_proton_total;
      else if (pdg == PDGNEUTRON) function = XSecFunctions::lead_neutron_total;
      else if (pdg == PDGKAON) return false;
      else return false;
    }
  }
  else return false;
  return true;
}

double MnvHadronReweight::getAverageCrossSection2(int pdg, int nuke, double Ti, double Tf, XSecType xsectype, bool isCVCalc)
{

  if (m_reweightNeutronCV && pdg == PDGNEUTRON && nuke == -6 && !isCVCalc) return getAverageNewNeutronCrossSection(nuke, Ti, Tf, xsectype);
  double integral = 0.0;
  TF1 func;
  if (getXSecFunc(pdg, nuke, xsectype, func)) integral = safeEval(func, Ti, Tf);
  return integral;
}

double MnvHadronReweight::getAverageNewNeutronCrossSection(int nuke, double Ti, double Tf, XSecType xsectype)
{
  double integral = 0.0;
  if (nuke != -6) return integral; // Not scint, not supported
  TF1 func;
  if (xsectype == xsectype_elastic) func = XSecFunctions::new_carbon_neutron_elastic;
  else if (xsectype == xsectype_inelastic) func = XSecFunctions::new_carbon_neutron_inelastic;
  else if (xsectype == xsectype_total) func = XSecFunctions::new_carbon_neutron_total;
  integral = safeEval(func, Ti, Tf);
  return integral; 
}


#endif // MNVHADRONREWEIGHT_CXX
