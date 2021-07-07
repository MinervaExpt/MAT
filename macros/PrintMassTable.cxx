#include <PlotUtils/TargetUtils.h>

#include <TString.h>

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <iomanip>


using namespace std;
using namespace PlotUtils;



//==================================
// User options
//==================================
//tracker fiducial size
const int N_FID_PLANES = 12 * 9; //12 planes, 9 faux targets (modules 27-80)

//formatting
const int PRECISION = 2; //N decimal places

//do you want a tabular latex table?
const bool tabular = true;
const string colSep = tabular ? " & " : ", ";
const string colEnd = tabular ? " \\\\ \n" : "\n";



//declare functions
void next( ostream &fout ) { fout << endl << endl; };
void PrintTracker( ostream &fout, int nplanes, bool isMC );
void PrintTrackerElement( ostream &fout, const std::string& elementName, int Z, bool isMC);
void PrintPassives( ostream &fout, bool isMC );
int PrintMassTable();


// Help with using passive targets
typedef pair<int,int> MatPair; //targetID->targetZ
typedef vector<MatPair> MatPairs;
MatPairs GetStdMatPairs()
{
  MatPairs mats;
  mats.push_back( MatPair(1,26) );
  mats.push_back( MatPair(1,82) );
  mats.push_back( MatPair(2,26) );
  mats.push_back( MatPair(2,82) );
  mats.push_back( MatPair(3, 6) );
  mats.push_back( MatPair(3,26) );
  mats.push_back( MatPair(3,82) );
  mats.push_back( MatPair(4,82) );
  mats.push_back( MatPair(5,26) );
  mats.push_back( MatPair(5,82) );
  return mats;
}


//============================================
//============================================
//============================================

void PrintTracker( ostream &fout, int nplanes, bool isMC )
{
  TargetUtils t = TargetUtils::Get();
  string header = "Areal Mass of a Plane (g/cm^2), Fiducial Area (cm^2), Mass of a Plane (kg), Protons per Plane, Neutrons per Plane, Nucleons per Plane";
  if( tabular )
  {
    header = "\\begin{table}\n";
    header += "\\begin{center}\n";
    header += "\\begin{tabular}{cccccc}\n";
    header += "\\toprule\n";
    header += "Areal Mass ($\\frac{g}{cm^{2}}$)& Fiducial Area ($cm^{2}$)& Mass (kg)&N Protons&N Neutrons& N Nucleons \\\\ \n";
    header += "\\midrule\n";
  }

  string caption( Form("Contents of one plane in %s", ( isMC ? " MC " : " Data " ) ) );

  using namespace TargetProp::Scint;

  const double onePlane_mass = t.GetTrackerMass(1, isMC) / (1000.);
  const double onePlane_prot = t.GetTrackerNProtons(1, isMC);
  const double onePlane_neut = t.GetTrackerNNeutrons(1, isMC);
  const double onePlane_nucl = t.GetTrackerNNucleons(1, isMC);
  const double dens = isMC ? areal_density_MC : areal_density;


  if( ! tabular )
    fout << caption << endl;
  fout << header << endl;
  fout
    << setprecision(PRECISION)
    << scientific
    << dens * (10.*10.) << colSep
    << t.GetHexArea( ) / (10.*10.) << colSep
    << onePlane_mass << colSep
    << onePlane_prot << colSep
    << onePlane_neut << colSep
    << onePlane_nucl 
    << colEnd;

  if( tabular )
  {
    string footer = "\\bottomrule\n";
    footer += "\\end{tabular}\n";
    footer += "\\end{center}\n";
    footer += Form("\\caption{%s}\n", caption.c_str() );
    footer += Form("\\label{tab:plane_contents_%s}\n", ( isMC ? "MC" : "Data" ) );
    footer += "\\end{table}\n";
    fout << footer;
  }

  next(fout);

  caption = string( Form("Tracker plane composition in %s", ( isMC ? " MC " : " Data " ) ) );
  header = "Element, \% of Plane Mass, Total Mass (kg), Atoms Contributed, Protons Contributed, Neutrons Contributed, Nucleons Contributed";
  if( tabular )
  {
    header = "\\begin{table}\n";
    header += "\\begin{center}\n";
    header += "\\begin{tabular}{c|cccccc}\n";
    header += "\\toprule\n";
    header += "Element& \\\% of Plane's Mass& Mass (kg)&N Atoms&N Protons&N Neutrons&N Nucleons \\\\ \n";
    header += "\\midrule\n";
  }

  if( ! tabular )
    fout << caption << endl;
  fout << header << endl;
  PrintTrackerElement( fout, "Carbon", 6, isMC );
  PrintTrackerElement( fout, "Hydrogen", 1, isMC  );
  PrintTrackerElement( fout, "Oxygen", 8, isMC  );
  PrintTrackerElement( fout, "Aluminum", 13, isMC  );
  PrintTrackerElement( fout, "Titanium", 22, isMC  );
  PrintTrackerElement( fout, "Silicon", 14, isMC  );
  PrintTrackerElement( fout, "Chlorine", 17, isMC  );

  if( tabular )
  {
    string footer = "\\bottomrule\n";
    footer += "\\end{tabular}\n";
    footer += "\\end{center}\n";
    footer += Form("\\caption{%s}\n", caption.c_str() );
    footer += Form("\\label{tab:plane_composition_%s}\n", ( isMC ? "MC" : "Data" ) );
    footer += "\\end{table}\n";
    fout << footer;
  }


}

void PrintTrackerElement( ostream &fout, const std::string& elementName, int Z, bool isMC )
{
  TargetUtils t = TargetUtils::Get();

  const double onePlane_massTotal = t.GetTrackerMass(1, isMC);
  const double massFrac = t.GetTrackerElementMassFraction( Z, isMC );
  const double onePlane_mass = massFrac * onePlane_massTotal;
  const double onePlane_atoms = onePlane_mass*t.GetTrackerElementAtomsPerGram( Z );
  const double onePlane_prot  = onePlane_atoms * Z;
  const double onePlane_neut  = onePlane_atoms * t.GetTrackerElementA( Z );
  const double onePlane_nucl = onePlane_prot + onePlane_neut;


  fout
    << elementName << colSep
    << setprecision(PRECISION)
    << fixed
    << massFrac*100. << colSep
    << scientific
    << onePlane_mass /1000. << colSep
    << onePlane_atoms << colSep
    << onePlane_prot << colSep
    << onePlane_neut << colSep
    << onePlane_nucl
    << colEnd;
}

//============================================
void PrintPassives( ostream &fout, bool isMC )
{

  TargetUtils &t = TargetUtils::Get();
  MatPairs mats = GetStdMatPairs();

  string header = "TargetID, TargetZ, Fiducial Area (cm^2), Areal Mass (g/cm^2), Mass (kg), N Protons, N Neutrons, N Nucleons";
  if( tabular )
  {
    header = "\\begin{table}\n";
    header += "\\begin{center}\n";
    header += "\\begin{tabular}{cc|ccc|ccc}\n";
    header += "\\toprule\n";
    header += "TargetID& TargetZ& Fiducial Area ($cm^{2}$)& Areal Mass ($\\frac{g}{cm^{2}}$)& Mass (kg)& N Protons& N Neutrons& N Nucleons \\\\ \n";
    header += "\\midrule\n";
  }

  string caption( Form("Passive targets in %s", ( isMC ? " MC " : " Data " ) ) );

  if( !tabular )
    fout << caption << endl;

  fout << header << endl;

  for( MatPairs::iterator mat = mats.begin(); mat != mats.end(); ++mat )
  {
    int targetID = mat->first;
    int targetZ  = mat->second;

    fout
      << targetID << colSep
      << targetZ << colSep
      << scientific
      << setprecision(PRECISION)
      << t.GetPassiveTargetArea( targetID, targetZ )  / (10.*10.) << colSep
      << t.GetPassiveTargetArealDensity( targetID, targetZ, isMC )  * (10.*10.) << colSep
      << t.GetPassiveTargetMass( targetID, targetZ, isMC ) / 1000. << colSep
      << t.GetPassiveTargetNProtons( targetID, targetZ, isMC ) << colSep
      << t.GetPassiveTargetNNeutrons( targetID, targetZ, isMC ) << colSep
      << t.GetPassiveTargetNNucleons( targetID, targetZ, isMC ) 
      << fixed << colEnd;
  }

  if( tabular )
  {
    string footer = "\\bottomrule\n";
    footer += "\\end{tabular}\n";
    footer += "\\end{center}\n";
    footer += Form("\\caption{%s}\n", caption.c_str() );
    footer += Form("\\label{tab:passive_targets_%s}\n", ( isMC ? "MC" : "Data" ) );
    footer += "\\end{table}\n";
    fout << footer;
  }

}

int PrintMassTable()
{
  //or open a file
  ostream& fout = cout;

  PrintTracker( fout, N_FID_PLANES, false );
  next(fout);
  PrintTracker( fout, N_FID_PLANES, true );

  next(fout);
  PrintPassives( fout, false ); //data

  next(fout);
  PrintPassives( fout, true ); //MC

  return 0;
}

int main()
{
  return PrintMassTable();
}
