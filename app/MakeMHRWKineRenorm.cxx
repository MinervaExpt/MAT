#include "TApplication.h"
#include "PlotUtils/NSFDefaults.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/GeantHadronSystematics.h"
#include "Cintex/Cintex.h"
#include <getopt.h>

int main( int argc, char **argv )
{
  ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);

  //REQUIRED
  std::string filelist     = "";

  //optional
  double minZ              = NSFDefaults::TrackerFace;
  double maxZ              = NSFDefaults::TrackerBack;
  double apothem           = NSFDefaults::StandardApothem;
  bool bElastics           = true;
  std::string playlist     = ""; 
  std::string project_name = "Renorm_Kine_Truth";
  std::string output_dir   = Form("%s/data/mhrwKineRenorm", getenv("MATFLUXANDWEIGHTFILES"));

  const char* const short_options = "f:m:M:a:eP:p:o:h";
  static struct option long_options[]=
  {
    {"filelist",     required_argument, nullptr, 'f'},
    {"minZ",         required_argument, nullptr, 'm'},
    {"maxZ",         required_argument, nullptr, 'M'},
    {"apothem",      required_argument, nullptr, 'a'},
    {"elasticsOff",  no_argument,       nullptr, 'e'},
    {"playlist",     required_argument, nullptr, 'P'},
    {"project_name", required_argument, nullptr, 'p'},
    {"output_dir",   required_argument, nullptr, 'o'},
    {"help",         no_argument,       nullptr, 'h'},
    {nullptr,        0,                 nullptr, 0}
  };

  int cc;
  while ((cc = getopt_long(argc, argv, short_options, long_options, nullptr)) != -1) {
    switch (cc)
    {
      case 'f': filelist     = std::string(optarg); break;
      case 'm': minZ         = atof(optarg);        break;
      case 'M': maxZ         = atof(optarg);        break;
      case 'a': apothem      = atof(optarg);        break;
      case 'e': bElastics    = false;               break;
      case 'P': playlist     = std::string(optarg); break;
      case 'p': project_name = std::string(optarg); break;
      case 'o': output_dir   = std::string(optarg); break;
      case '?':
      case 'h':
      default:
        std::cout << "|********************************* Options ************************************************************|" << std::endl
                  << "| This program creates the kinetic renormalizations for the MnvHadronReweighter.                       |" << std::endl
                  << "| Please see Aaron Bercellie's talk (docdb 28556)                                                      |" << std::endl
                  << "| These must be MC files, preferably full detector processing                                          |" << std::endl
                  << "| This will use the Truth trees                                                                        |" << std::endl
                  << "|                                                                                                      |" << std::endl
                  << "| Options:                                                                                             |" << std::endl
                  << "|   -f, --filelist     : (Required) Text file with list of anatuples                                   |" << std::endl
                  << "|   -m, --minZ         : (Optional) Mininmum Z of readout volume (Default: NSFDefaults::TrackerFace)   |" << std::endl
                  << "|   -M, --maxZ         : (Optional) Maximum Z of readout volume  (Default: NSFDefaults::TrackerBack)   |" << std::endl
                  << "|   -a, --apothem      : (Optional) Apothem of readout volume (Default: NSFDefaults::StandardApothem)  |" << std::endl
                  << "|   -e, --elasticsOff  : (Optional) Turn elastic cross sections off (Default: On)                      |" << std::endl
                  << "|   -P, --playlist     : (Optional) Playlist name to add to output file (Default: "" )                 |" << std::endl
                  << "|   -p,   project_name : (Optional) Output file name (Default: Renorm_Kine_Truth)                      |" << std::endl
                  << "|   -o, --output_dir   : (Optional) Output directories (Default: $MATFLUXANDWEIGHTFILES/data/mhrwKineRenorm)   |" << std::endl
                  << "|   -h, --help         : Print this.                                                                   |" << std::endl
                  << "|******************************************************************************************************|" << std::endl;
         return 0;
    }
  }
  if( filelist.empty() ) { std::cout<<"Need a list of files (-f, --filelist)"<<std::endl; return 1; }

  //Make chain from file
  PlotUtils::ChainWrapper* chain = makeChainWrapperPtr( filelist, "Truth");
                                                                                                          //Turning neutron CV on as it doesn't really affect the renormalization   
  //Set up MnvHadronReweighter instance, run the renormalization                                            it just fills some additional histograms
  MnvHadronReweight mhrw = weight_hadron<PlotUtils::ChainWrapper*>( chain, playlist, minZ, maxZ, apothem, true, bElastics, project_name, output_dir);
  mhrw.getTruthKineRenorm();
  return 0;
}
