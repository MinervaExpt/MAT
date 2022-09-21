#ifndef MNV_MnvPlotter_cxx
#define MNV_MnvPlotter_cxx 1

#include "PlotUtils/MnvPlotter.h"
#include "PlotUtils/HistogramUtils.h" //for IsAutoAxisLimit
#include "PlotUtils/MnvColors.h"

#include "TROOT.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLegend.h"
#include "TLine.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TList.h"

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <assert.h>
#include <math.h>


//set this to 1 to turn on garbage collection by deletion of tmpObjects
#define DO_GARBAGE_COLLECTION 0
// Some weirdness copied from TColor::SetPalette()
// #define fgPalettesList TColor__PalettesList()

using namespace std;
using namespace PlotUtils;

// Some weirdness copied from TColor::SetPalette()
// namespace {
//    static TArrayD& TColor__PalettesList() {
//       static TArrayD globalPalettesList(0);
//       return globalPalettesList;
//    }
// }

/*#####################################################

  Below is an example of how to use it :
  ... to be updated ...
  }
########################################################*/

//=========================================================
// Constructors and Destructors
//=========================================================
MnvPlotter::MnvPlotter()
{
#if DO_GARBAGE_COLLECTION
    gROOT->GetListOfCleanups()->Add( &fTmpObjects );
#endif
    ApplyStyle( kDefaultStyle );
}

MnvPlotter::MnvPlotter( PlotUtils::t_PlotStyle style )
{
#if DO_GARBAGE_COLLECTION
    gROOT->GetListOfCleanups()->Add( &fTmpObjects );
#endif
    ApplyStyle( style );
}


MnvPlotter::~MnvPlotter()
{
#if DO_GARBAGE_COLLECTION

    //! Delete all heap based objects in tmpObjects array
    CleanTmp();

#endif
}


void MnvPlotter::CleanTmp()
{
#if DO_GARBAGE_COLLECTION

    //! Get rid of the NULL entries
    fTmpObjects.Compress();

    //Error("MnvPlotter::CleanTmp", Form( "Number of objects in tmp before clean: %d", fTmpObjects.GetEntries() ) );

    //! Call Delete to get rid of the tmp objects
    fTmpObjects.Delete();


    //Error("MnvPlotter::CleanTmp", Form( "Number of objects in tmp after clean: %d", fTmpObjects.GetEntries() ) );

#endif
}

void MnvPlotter::AddToTmp( TObject* obj )
{
    //note: build this everytime to avoid obj unused compiler warning
    if ( 0 == obj )
    {
        Warning( "MnvPlotter::AddToTmp", "Attempting to add NULL object to garbage.");
        return;
    }

#if DO_GARBAGE_COLLECTION

    //Error("MnvPlotter::AddToTmp", Form( "Adding object '%s' at index %d.", obj->GetName(), fTmpObjects.GetEntries() ) );

    //! Add to array if the object isn't in the array already
    if ( ! fTmpObjects.FindObject( obj ) )
        fTmpObjects.AddLast( obj );
#endif
}

//=========================================================
// utility to set up basic root environment
//=========================================================
void MnvPlotter::SetRootEnv()
{
    gStyle->SetPalette(palette_style);

    // Canvas Styles
    gStyle->SetCanvasDefW(900);
    gStyle->SetCanvasDefH(750);
    gStyle->SetOptStat(0000);
    gStyle->SetOptFit(0000);
    gStyle->SetOptTitle(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadTopMargin(0.09);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetHistLineWidth(2);

    // Axis Styles
    gStyle->SetHistMinimumZero( hist_min_zero );
    gStyle->SetTitleOffset( axis_title_offset_x, "X" );
    gStyle->SetTitleSize( axis_title_size_x, "X" );
    gStyle->SetTitleFont( axis_title_font_x, "X" );
    gStyle->SetTitleOffset( axis_title_offset_y, "Y" );
    gStyle->SetTitleSize( axis_title_size_y, "Y" );
    gStyle->SetTitleFont( axis_title_font_y, "Y" );
    gStyle->SetTitleOffset( axis_title_offset_z, "Z" );
    gStyle->SetTitleSize( axis_title_size_z, "Z" );
    gStyle->SetTitleFont( axis_title_font_z, "Z" );
    gStyle->SetLabelFont( axis_label_font, "XYZ" );
    gStyle->SetLabelSize( axis_label_size, "XYZ" );
    TGaxis::SetMaxDigits(axis_max_digits);
    gStyle->SetPadGridX( axis_draw_grid_x );
    gStyle->SetPadGridY( axis_draw_grid_y );

    // Marker Styles
    gStyle->SetMarkerStyle(data_marker);
    gStyle->SetMarkerSize(data_marker_size);
    gStyle->SetMarkerColor(data_color);

    gStyle->SetEndErrorSize(2);
    gStyle->SetErrorX(0.5);
}

//================================================================
// set the variables of this namespace to defaults and apply style
//================================================================
void MnvPlotter::ApplyStyle( PlotUtils::t_PlotStyle style /* = kDefaultStyle */ )
{
  n_color_contours = 999;

  if ( style == kDefaultStyle ) {
    //-- chi2 calculation
    chi2_use_overflow_err = false;

    //-- do you want to draw histograms normalized to bin width
    draw_normalized_to_bin_width = true;

    //-- rainbow
    palette_style = 1;

    //-- marker settings
    data_marker = 20;
    ratio_marker = 20;
    data_marker_size = 1.0;
    ratio_marker_size = 1.0;

    //-- line settings
    data_line_width = 1;
    data_line_style = 1;
    mc_line_width = 3;
    mc_line_style = 1;
    ratio_line_width = 3;

    //-- cut arrow settings
    arrow_line_width = 4;
    arrow_line_style = 1;
    arrow_line_color = 1;
    arrow_size = 0.01;
    arrow_type = "|>";

    //-- color settings
    data_color = 1;
    mc_color   = 2;
    mc_error_color = kRed-10;//45;
    mc_error_style = 1001;
    ratio_color = 1;

    mc_bkgd_color = 14;
    mc_bkgd_width = 1;
    mc_bkgd_line_color = 1;
    mc_bkgd_style = 3005;

    data_bkgd_color = 12;//gray
    data_bkgd_style = 24;//circle
    data_bkgd_size   = 1.;

    //-- correlation
    draw_corr_max1     = false;
    draw_corr_red_blue = true;

    //-- title settings
    title_font    = 62;
    title_size = 0.06;

    //-- axis options
    hist_min_zero     = true;
    axis_draw_grid_x = false;
    axis_draw_grid_y = false;
    axis_max_digits   = 3;
    axis_title_font_x = 62;
    axis_title_font_y = 62;
    axis_title_font_z = 62;
    axis_title_offset_x = 1.15;
    axis_title_offset_y = 1.2;
    axis_title_offset_z = .75;
    axis_title_size_x = 0.06;
    axis_title_size_y = 0.06;
    axis_title_size_z = 0.06;
    axis_minimum      = MnvHist::AutoAxisLimit;
    axis_maximum      = MnvHist::AutoAxisLimit;
    axis_maximum_group= MnvHist::AutoAxisLimit; //0.5;

    //-- axis label options
    axis_label_font = 42;
    axis_label_size = 0.05;

    //-- margins
    extra_top_margin    = -.02; //negative means go closer to edge
    extra_bottom_margin = 0.;
    extra_left_margin   = 0.;
    extra_right_margin  = -0.50;

    //-- layout
    headroom = 1.5;  //old 1.5
    footroom = 1.25;  //old 1

    //-- legend
    height_nspaces_per_hist = 2.;
    width_xspace_per_letter = .5;
    legend_border_size      = 0;
    legend_fill_color       = -1;
    legend_text_size        = .035;
    legend_offset_x         = 0.;
    legend_offset_y         = 0.;
    legend_n_columns        = 1;
    legend_text_font        = 62;

    //-- define good colors for general use
    //-- used in particular for DrawErrorSummary
    good_colors.clear();
    good_colors = MnvColors::GetColors(MnvColors::kOkabeItoPalette);
    good_colors.insert(good_colors.begin(), kBlack);

    //-- define colors of the standard errors
    error_color_map.clear();
    error_color_map["Flux_Tertiary"]  = kYellow-3;
    error_color_map["Flux_BeamFocus"] = kOrange+2;
    error_color_map["Flux_NA49"]      = kRed+2;
    error_color_map["GENIE"]          = kGreen+2;
    error_color_map["Normalization"]  = kTeal+2;
    error_color_map["Muon_Energy"]    = kPink+2;
    error_color_map["Hadronic_Energy"]= kMagenta+2;
    error_color_map["BackgroundFit"]  = kBlue+2;

    //-- stat error labeling
    stat_error_name = "Statistical";

    print_formats.clear();
    print_formats.push_back( "png" );
    print_formats.push_back( "eps" );
    print_formats.push_back( "C" );

    print_topdir = "";
  }
  else if ( style == kCompactStyle ) {
    // Start from defaults
    ApplyStyle( kDefaultStyle );

    //tweak margins
    extra_top_margin = -.035; //go slightly closer to top of pad

    mc_bkgd_color = 46;
    mc_bkgd_line_color = 46;

    data_bkgd_color = 12;//gray
    data_bkgd_style = 24;//circle

    //legend entries are closer
    height_nspaces_per_hist = 1.2;
    width_xspace_per_letter = .4;
    legend_text_size        = .03;
  }
  else if ( style == kCCNuPionIncStyle) {
    ApplyStyle( kCompactStyle );

    //-- do you want to draw histograms normalized to bin width
    draw_normalized_to_bin_width = false;

    legend_n_columns        = 1;

    //-- define colors of the standard errors
    error_color_map.clear();
    error_summary_group_map.clear();

    //Systematic color scheme
    error_color_map["Flux"]        = kYellow-3;
    error_color_map["Interaction Model"] = kGreen+2;
    error_color_map["Michel"] = kPink+2;
    error_color_map["Normalization"] = kCyan+2;
    error_color_map["Detector Model"]  = kOrange+2;
    error_color_map["Energy Response"] = kRed+2;
    error_color_map["Angle & Vertex"]  = kViolet+2;
    error_color_map["Other"]  = kMagenta+2;

    vector<string> genieGroup;
    genieGroup.push_back("GENIE_AhtBY");
    genieGroup.push_back("GENIE_BhtBY");
    genieGroup.push_back("GENIE_CCQEPauliSupViaFK");
    genieGroup.push_back("GENIE_CV1uBY");
    genieGroup.push_back("GENIE_CV2uBY");
    genieGroup.push_back("GENIE_EtaNCEL");
    genieGroup.push_back("GENIE_MaCCQE");
    genieGroup.push_back("GENIE_MaCCQEshape");
    genieGroup.push_back("GENIE_MaNCEL");
    genieGroup.push_back("GENIE_MaRES");
    genieGroup.push_back("GENIE_MvRES");
    genieGroup.push_back("GENIE_NormCCQE");
    genieGroup.push_back("GENIE_NormCCRES");
    genieGroup.push_back("GENIE_NormDISCC");
    genieGroup.push_back("GENIE_NormNCRES");
    genieGroup.push_back("GENIE_Rvn1pi");
    genieGroup.push_back("GENIE_Rvn2pi");
    genieGroup.push_back("GENIE_Rvp1pi");
    genieGroup.push_back("GENIE_Rvp2pi");
    genieGroup.push_back("GENIE_VecFFCCQEshape");

    genieGroup.push_back("GENIE_AGKYxF1pi");
    genieGroup.push_back("GENIE_FrAbs_N");
    genieGroup.push_back("GENIE_FrAbs_pi");
    genieGroup.push_back("GENIE_FrCEx_N");
    genieGroup.push_back("GENIE_FrCEx_pi");
    genieGroup.push_back("GENIE_FrElas_N");
    genieGroup.push_back("GENIE_FrElas_pi");
    genieGroup.push_back("GENIE_FrInel_N");
    genieGroup.push_back("GENIE_FrInel_pi");
    genieGroup.push_back("GENIE_FrPiProd_N");
    genieGroup.push_back("GENIE_FrPiProd_pi");
    genieGroup.push_back("GENIE_MFP_N");
    genieGroup.push_back("GENIE_MFP_pi");
    genieGroup.push_back("GENIE_RDecBR1gamma");
    genieGroup.push_back("GENIE_Theta_Delta2Npi");
    genieGroup.push_back("GENIE_EFNUCR");
    genieGroup.push_back("GENIE_FZONE");
    genieGroup.push_back("GENIE_AKGY");
    error_summary_group_map["Interaction Model"] = genieGroup;

    vector<string> fluxGroup;
    fluxGroup.push_back("Flux_Tertiary");
    fluxGroup.push_back("Flux_BeamFocus");
    fluxGroup.push_back("Flux_NA49");
    error_summary_group_map["Flux"] = fluxGroup;

    vector<string> michelGroup;
    michelGroup.push_back("Michel_Eff");
    michelGroup.push_back("Michel_Bkg");
    michelGroup.push_back("Michel_AlternateCuts");
    error_summary_group_map["Michel"] = michelGroup;

    vector<string> geantGroup;
    geantGroup.push_back("Geant4_TotalPionInelastic");
    geantGroup.push_back("Geant4_TotalProtonInelastic");
    geantGroup.push_back("Geant4_ComponentPionInelastic");
    error_summary_group_map["Detector Model"] = geantGroup;

    vector<string> normGroup;
    normGroup.push_back("Normalization");
    normGroup.push_back("Mass_Scale");
    normGroup.push_back("POT_Scale");
    normGroup.push_back("Proton_PID_Eff");
    error_summary_group_map["Normalization"] = normGroup;

    vector<string> anaGroup;
    anaGroup.push_back("BackgroundFit");
    anaGroup.push_back("Unfolding");
    anaGroup.push_back("Unfolding_Theta");
    error_summary_group_map["Other"] = anaGroup;

    vector<string> angleGroup;
    angleGroup.push_back("Track_Angle");
    angleGroup.push_back("Beam_Angle");
    angleGroup.push_back("Vertex");
    error_summary_group_map["Angle & Vertex"] = angleGroup;

    vector<string> detGroup;
    detGroup.push_back("MINOS_Muon_Energy");
    detGroup.push_back("BetheBloch");
    detGroup.push_back("Mass");
    detGroup.push_back("Birks");
    detGroup.push_back("Particle_Response");
    detGroup.push_back("Res_Calorimetry");
    error_summary_group_map["Energy Response"] = detGroup;
  }
  else if ( kNukeCCStyle == style ) {
    // Start from compact
    ApplyStyle( kCompactStyle );

    data_marker_size = ratio_marker_size = data_bkgd_size = 1.75;
    data_line_width = 5; //make these thicker now that everyone has 4k UHD monitors...
    mc_line_width = 5; //make these thicker now that everyone has 4k UHD monitors...
    mc_bkgd_style = 3002; //change the hatched histogram
    data_bkgd_style = 22;//triangle
    gStyle->SetEndErrorSize(0.0); //DO NOT WANT
    legend_text_size = .04;

    axis_draw_grid_x = true;
    axis_draw_grid_y = true;
    legend_border_size      = 1;
    legend_fill_color       = 10;

    //change the BG color scheme
    //-- define colors of the standard errors
    error_color_map.clear();
    error_summary_group_map.clear();
    /*
    //verbose errors
    error_color_map["BG Scale"]       = kViolet+2;
    error_color_map["BG MC Stat."]    = kRed;
    error_color_map["Eff. MC Stat."]  = kYellow+2;
    error_color_map["Flux"]           = kOrange+2;
    error_color_map["XSec Models"]    = kGreen+2;
    error_color_map["FSI Models"]     = kCyan+2;
    error_color_map["Normalization"]  = kTeal+2;
    error_color_map["Vertex Rec."]    = kAzure+2;
    error_color_map["Muon Energy Rec."]    = kPink+2;
    error_color_map["Muon Angle Rec."]     = kRed+2;
    error_color_map["Hadronic Energy Rec."]= kMagenta+2;
    */
    error_color_map["Scint. BG"] = kMagenta;
    error_color_map["BG Scale"]       = kMagenta;
    error_color_map["BG MC Stat."]    = kTeal+2;
    error_color_map["Interaction Models"] = kGreen+2;
    error_color_map["FSI Models"]  = kCyan+2;
    error_color_map["Flux"]        = kYellow-3;
    error_color_map["Detector Res."] = kRed;
    error_color_map["Flux + Mass"]  = kOrange+2;
    error_color_map["MC Stats."] = kPink+2;;
    error_color_map["Sidebands"] = kViolet+2;

    /*vector<string> scintBGGroup;
      scintBGGroup.push_back("BG Scale");
      scintBGGroup.push_back("BG MC Stat.");
      error_summary_group_map["Scint. BG"] = scintBGGroup;*/

    vector<string> mcStats;
    mcStats.push_back("BG MC Stat.");
    mcStats.push_back("Eff Stat.");
    error_summary_group_map["MC Stats."] = mcStats;

    //stat error labeling
    stat_error_name = "Statistical";

    vector<string> detResGroup;
    detResGroup.push_back("Vertex Rec.");
    detResGroup.push_back("Muon Energy Rec.");
    detResGroup.push_back("Muon Angle Rec.");
    detResGroup.push_back("Hadronic Energy Rec.");
    detResGroup.push_back("Birks' Parameter");
    detResGroup.push_back( "Muon_Energy_MINOS" );
    detResGroup.push_back( "Muon_Energy_MINERvA" );
    detResGroup.push_back( "Muon_Energy_Resolution" );
    detResGroup.push_back( "BeamAngleX" );
    detResGroup.push_back( "BeamAngleY" );
    error_summary_group_map["Detector Res."] = detResGroup;

    vector<string> otherGroup;
    otherGroup.push_back("Flux");
    otherGroup.push_back("Normalization");
    otherGroup.push_back("C Mass");
    otherGroup.push_back("Fe Mass");
    otherGroup.push_back("Pb Mass");
    otherGroup.push_back("CH Mass");
    error_summary_group_map["Flux + Mass"] = otherGroup;

    vector<string> xsecErrGroup, fsiErrGroup;
    xsecErrGroup.push_back("AhtBY");
    xsecErrGroup.push_back("BhtBY");
    xsecErrGroup.push_back("CCQEPauliSupViaKF");
    xsecErrGroup.push_back("CV1uBY");
    xsecErrGroup.push_back("CV2uBY");
    xsecErrGroup.push_back("EtaNCEL");
    xsecErrGroup.push_back("MaCCQE");
    xsecErrGroup.push_back("MaNCEL");
    xsecErrGroup.push_back("MaRES");
    xsecErrGroup.push_back("MvRES");
    xsecErrGroup.push_back("NormCCQE");
    xsecErrGroup.push_back("NormCCRES");
    xsecErrGroup.push_back("NormDISCC");
    xsecErrGroup.push_back("NormNCRES");
    xsecErrGroup.push_back("Rvn1pi");
    xsecErrGroup.push_back("Rvn2pi");
    xsecErrGroup.push_back("Rvp1pi");
    xsecErrGroup.push_back("Rvp2pi");
    xsecErrGroup.push_back("VecFFCCQEshape");
    xsecErrGroup.push_back("RPA-Model");
    xsecErrGroup.push_back("Non Resonant Pion");
    xsecErrGroup.push_back("2p2h-Model"); //
    //xsecErrGroup.push_back("Eff Stat.");

    fsiErrGroup.push_back("AGKYxF1pi");
    fsiErrGroup.push_back("FrAbs_N");
    fsiErrGroup.push_back("FrAbs_pi");
    fsiErrGroup.push_back("FrCEx_N");
    fsiErrGroup.push_back("FrCEx_pi");
    fsiErrGroup.push_back("FrElas_N");
    fsiErrGroup.push_back("FrElas_pi");
    fsiErrGroup.push_back("FrInel_N");
    fsiErrGroup.push_back("FrInel_pi");
    fsiErrGroup.push_back("FrPiProd_N");
    fsiErrGroup.push_back("FrPiProd_pi");
    fsiErrGroup.push_back("MFP_N");
    fsiErrGroup.push_back("MFP_pi");
    fsiErrGroup.push_back("RDecBR1gamma");
    fsiErrGroup.push_back("Theta_Delta2Npi");
    fsiErrGroup.push_back("Nuclear Radius");
    fsiErrGroup.push_back("Formation Time");
    fsiErrGroup.push_back("AGKY Model");

    error_summary_group_map["Interaction Models"] = xsecErrGroup;
    error_summary_group_map["FSI Models"] = fsiErrGroup;

    vector<string> sb;
    sb.push_back("Plastic_SB");
    sb.push_back("Phys_SB");
    error_summary_group_map["Sidebands"] = sb;

    print_formats.clear();
    print_formats.push_back( "png" );
    //print_formats.push_back( "eps" );
    //print_formats.push_back( "pdf" );
    print_formats.push_back( "C" );
    //print_formats.push_back( "root" );
  }
  else if ( style == kNukeCCPrintStyle) {
    //same as NukeCC, only make lines a little thinner for printing
    ApplyStyle(kNukeCCStyle);
    data_line_width = 3;
    mc_line_width = 4; //make these thicker now that everyone has 4k UHD monitors...
  }
  else if ( style == kCCCohStyle) {
    ApplyStyle( kCompactStyle );

    legend_n_columns        = 2;
    axis_maximum_group      =0.2;
    height_nspaces_per_hist = 1.;
    legend_text_size        = .035;
    hist_min_zero = false;
    MnvHist::IsAutoAxisLimit(axis_minimum);

    //-- define colors of the standard errors
    error_color_map.clear();
    error_summary_group_map.clear();

    //Systematic color scheme
    error_color_map["Flux"]        = kYellow-3;
    error_color_map["Interaction Model"] = kGreen+2;
    error_color_map["MINOS Matching Eff"] = kBlue+2;
    error_color_map["Detector Model"]  = kOrange+2;
    error_color_map["Energy Response"] = kRed+2;
    error_color_map["Sideband Model"]  = kViolet+2;
    error_color_map["Vertex Energy"]  = kMagenta+2;

    vector<string> genieGroup;
    genieGroup.push_back( "GENIE_NormCCQE" );
    genieGroup.push_back( "GENIE_MaCCQEshape" );
    genieGroup.push_back( "GENIE_VecFFCCQEshape" );
    genieGroup.push_back( "GENIE_CCQEPauliSupViaFK" );
    genieGroup.push_back( "GENIE_MaNCEL" );
    genieGroup.push_back( "GENIE_EtaNCEL" );
    genieGroup.push_back( "GENIE_MaRES" );
    genieGroup.push_back( "GENIE_MvRES" );
    genieGroup.push_back( "GENIE_NormCCRES" );
    genieGroup.push_back( "GENIE_NormNCRES" );
    genieGroup.push_back( "GENIE_RDecBR1gamma" );
    genieGroup.push_back( "GENIE_Theta_Delta2Npi" );
    genieGroup.push_back( "GENIE_Rvn1pi" );
    genieGroup.push_back( "GENIE_Rvn2pi" );
    genieGroup.push_back( "GENIE_Rvp1pi" );
    genieGroup.push_back( "GENIE_Rvp2pi" );
    genieGroup.push_back( "GENIE_AhtBY" );
    genieGroup.push_back( "GENIE_BhtBY" );
    genieGroup.push_back( "GENIE_CV1uBY" );
    genieGroup.push_back( "GENIE_CV2uBY" );
    genieGroup.push_back( "GENIE_NormDISCC" );
    genieGroup.push_back( "GENIE_FrAbs_N" );
    genieGroup.push_back( "GENIE_FrAbs_pi" );
    genieGroup.push_back( "GENIE_FrCEx_N" );
    genieGroup.push_back( "GENIE_FrCEx_pi" );
    genieGroup.push_back( "GENIE_FrElas_N" );
    genieGroup.push_back( "GENIE_FrElas_pi" );
    genieGroup.push_back( "GENIE_FrInel_N" );
    genieGroup.push_back( "GENIE_FrInel_pi" );
    genieGroup.push_back( "GENIE_FrPiProd_N" );
    genieGroup.push_back( "GENIE_FrPiProd_pi" );
    genieGroup.push_back( "GENIE_MFP_N" );
    genieGroup.push_back( "GENIE_MFP_pi" );
    genieGroup.push_back( "GENIE_AGKYxF1pi" );
    error_summary_group_map["Interaction Model"] = genieGroup;

    vector<string> fluxGroup;
    fluxGroup.push_back("Flux_Tertiary");
    fluxGroup.push_back("Flux_BeamFocus");
    fluxGroup.push_back("Flux_NA49");
    fluxGroup.push_back("Flux");
    error_summary_group_map["Flux"] = fluxGroup;

    vector<string> detGroup;
    detGroup.push_back("Pion Inel. XSec");
    detGroup.push_back("Proton Inel. XSec");
    detGroup.push_back("Neutron Pathlength");
    detGroup.push_back("BeamAngleXZ");
    detGroup.push_back("BeamAngleYZ");
    //detGroup.push_back("NuMI Beam Angle");
    error_summary_group_map["Detector Model"] = detGroup;

    vector<string> normGroup;
    normGroup.push_back("Norm. Corrections");
    error_summary_group_map["Tracking Eff"] = normGroup;

    vector<string> angleGroup;
    angleGroup.push_back("Sideband Angle");
    error_summary_group_map["Sideband Model"] = angleGroup;

    vector<string> eRespGroup;
    eRespGroup.push_back("Muon Energy");
    eRespGroup.push_back("Hadron Response");
    error_summary_group_map["Energy Response"] = eRespGroup;

    vector<string> vertexEGroup;
    vertexEGroup.push_back("Vertex Energy");
    error_summary_group_map["Vertex Energy"] =vertexEGroup;
  }
  else if ( style == kCCQENuStyle) {
    ApplyStyle( kCompactStyle );
    //-- define colors of the standard errors
    error_color_map.clear();
    error_summary_group_map.clear();

    error_color_map["Flux"]                    = kViolet+6;
    error_color_map["Recoil Reconstruction"]   = kOrange+2;
    error_color_map["Cross Section Models"]         = kMagenta;
    error_color_map["FSI Models"]              = kRed;
    error_color_map["Muon Reconstruction"]     = kOrange-3;
    error_color_map["Others"]                  = kGreen+3;
    error_color_map["Low Recoil Fits"]         = kRed+3;

    std::vector< string > flux;
    flux.push_back("Flux");
    error_summary_group_map["Flux"] = flux;

    std::vector< string > xsection;
    xsection.push_back( "GENIE_CCQEPauliSupViaKF" );
    xsection.push_back( "GENIE_EtaNCEL" );
    //xsection.push_back( "GENIE_MaCCQE" );
    xsection.push_back( "GENIE_MaCCQEshape" );
    xsection.push_back( "GENIE_MaNCEL" );
    xsection.push_back( "GENIE_MaRES" );
    xsection.push_back( "GENIE_MvRES" );
    xsection.push_back( "GENIE_NormCCQE" );
    xsection.push_back( "GENIE_NormCCRES" );
    xsection.push_back( "GENIE_NormDISCC" );
    xsection.push_back( "GENIE_NormNCRES" );
    xsection.push_back( "GENIE_Rvn1pi" );
    xsection.push_back( "GENIE_Rvn2pi" );
    xsection.push_back( "GENIE_Rvp1pi" );
    xsection.push_back( "GENIE_Rvp2pi" );
    xsection.push_back( "GENIE_VecFFCCQEshape" );
    xsection.push_back( "GENIE_AhtBY" );
    xsection.push_back( "GENIE_BhtBY" );
    xsection.push_back( "GENIE_CV1uBY" );
    xsection.push_back( "GENIE_CV2uBY" );
    xsection.push_back( "RPA_HighQ2" );
    xsection.push_back( "RPA_LowQ2" );
    error_summary_group_map["Cross Section Models"] = xsection;

    std::vector< string > fsi;
    fsi.push_back( "GENIE_AGKYxF1pi" );
    fsi.push_back( "GENIE_FrAbs_N" );
    fsi.push_back( "GENIE_FrAbs_pi" );
    fsi.push_back( "GENIE_FrCEx_N" );
    fsi.push_back( "GENIE_FrCEx_pi" );
    fsi.push_back( "GENIE_FrElas_N" );
    fsi.push_back( "GENIE_FrElas_pi" );
    fsi.push_back( "GENIE_FrInel_N" );
    fsi.push_back( "GENIE_FrInel_pi" );
    fsi.push_back( "GENIE_FrPiProd_N" );
    fsi.push_back( "GENIE_FrPiProd_pi" );
    fsi.push_back( "GENIE_MFP_N" );
    fsi.push_back( "GENIE_MFP_pi" );
    fsi.push_back( "GENIE_RDecBR1gamma" );
    fsi.push_back( "GENIE_Theta_Delta2Npi" );
    error_summary_group_map["FSI Models"] = fsi;

    std::vector< string > lowrecoilfit;
    lowrecoilfit.push_back("Low_Recoil_2p2h_Tune");
    error_summary_group_map["Low Recoil Fits"] = lowrecoilfit;

    std::vector< string > muon;
    //    muon.push_back( "Muon_Energy" );
    muon.push_back( "Muon_Energy_MINOS" );
    muon.push_back( "Muon_Energy_MINERvA" );
    //    muon.push_back( "Muon_Theta" );
    muon.push_back( "BeamAngleX" );
    muon.push_back( "BeamAngleY" );
    muon.push_back( "Muon_Energy_Resolution" );
    error_summary_group_map["Muon Reconstruction"] = muon;

    std::vector< string > other;
    //    other.push_back( "Normalization_Factors" );
    other.push_back( "MINOS_Reconstruction_Efficiency" );
    other.push_back( "Michel_Efficiency");
    other.push_back( "Target_Mass" );
    other.push_back( "GENIE_nonreweightable" );
    other.push_back( "Muon_Response" );
    other.push_back( "Proton_Response" );
    other.push_back( "Low_Neutron_Response" );
    other.push_back( "Mid_Neutron_Response" );
    other.push_back( "High_Neutron_Response" );
    other.push_back( "Pion_Response" );
    other.push_back( "EM_Response" );
    other.push_back( "Other_Response" );
    other.push_back( "Crosstalk" );
    other.push_back( "MEU_Recoil" );
    other.push_back( "Binding_Energy" );
    other.push_back( "Reweight_Pion" );
    other.push_back( "Reweight_Proton" );
    other.push_back( "Reweight_Neutron" );
    other.push_back( "Bethe_Bloch" );
    other.push_back( "MEU_Proton" );
    other.push_back( "Mass_Model_Proton" );
    other.push_back( "Birks_Response_Proton" );
    other.push_back( "Proton_TrackEff" );
    other.push_back( "Proton_Angle" );
    error_summary_group_map["Others"] = other;
  }

  else if ( style == kCCQENuInclusiveStyle) {
    ApplyStyle( kCompactStyle );
    //-- define colors of the standard errors
    error_color_map.clear();
    error_summary_group_map.clear();

    error_color_map["Flux"]                    = kViolet+6;
    error_color_map["Hadrons"]   = kCyan+2;
    error_color_map["Models"]         = kRed;
    error_color_map["Normalization"]              = kBlue;
    error_color_map["Muon Reconstruction"]     = kOrange-3;


    std::vector< string > flux;
    flux.push_back("Flux");
    error_summary_group_map["Flux"] = flux;

    std::vector< string > xsection;
    xsection.push_back( "GENIE_CCQEPauliSupViaKF" );
    xsection.push_back( "GENIE_EtaNCEL" );
    //xsection.push_back( "GENIE_MaCCQE" );
    xsection.push_back( "GENIE_MaCCQEshape" );
    xsection.push_back( "GENIE_MaNCEL" );
    xsection.push_back( "GENIE_MaRES" );
    xsection.push_back( "GENIE_MvRES" );
    xsection.push_back( "GENIE_NormCCQE" );
    xsection.push_back( "GENIE_NormCCRES" );
    xsection.push_back( "GENIE_NormDISCC" );
    xsection.push_back( "GENIE_NormNCRES" );
    xsection.push_back( "GENIE_Rvn1pi" );
    xsection.push_back( "GENIE_Rvn2pi" );
    xsection.push_back( "GENIE_Rvp1pi" );
    xsection.push_back( "GENIE_Rvp2pi" );
    xsection.push_back( "GENIE_VecFFCCQEshape" );
    xsection.push_back( "GENIE_AhtBY" );
    xsection.push_back( "GENIE_BhtBY" );
    xsection.push_back( "GENIE_CV1uBY" );
    xsection.push_back( "GENIE_CV2uBY" );
    xsection.push_back( "RPA_HighQ2" );
    xsection.push_back( "RPA_LowQ2" );
    xsection.push_back( "GENIE_AGKYxF1pi" );
    xsection.push_back( "GENIE_FrAbs_N" );
    xsection.push_back( "GENIE_FrAbs_pi" );
    xsection.push_back( "GENIE_FrCEx_N" );
    xsection.push_back( "GENIE_FrCEx_pi" );
    xsection.push_back( "GENIE_FrElas_N" );
    xsection.push_back( "GENIE_FrElas_pi" );
    xsection.push_back( "GENIE_FrInel_N" );
    xsection.push_back( "GENIE_FrInel_pi" );
    xsection.push_back( "GENIE_FrPiProd_N" );
    xsection.push_back( "GENIE_FrPiProd_pi" );
    xsection.push_back( "GENIE_MFP_N" );
    xsection.push_back( "GENIE_MFP_pi" );
    xsection.push_back( "GENIE_RDecBR1gamma" );
    xsection.push_back( "GENIE_Theta_Delta2Npi" );
    xsection.push_back("Low_Recoil_2p2h_Tune");
    error_summary_group_map["Models"] = xsection;

    std::vector< string > muon;
    muon.push_back( "Muon_Energy_MINOS" );
    muon.push_back( "Muon_Energy_MINERvA" );
    muon.push_back( "BeamAngleX" );
    muon.push_back( "BeamAngleY" );
    muon.push_back( "Muon_Energy_Resolution" );
    error_summary_group_map["Muon Reconstruction"] = muon;

    std::vector< string > normalization;
    normalization.push_back("MINOS_Matching_Efficiency");
    normalization.push_back( "Target_Mass" );
    error_summary_group_map["Normalization"] = normalization;

    std::vector< string > other;
    other.push_back( "Muon_Response" );
    other.push_back( "Proton_Response" );
    other.push_back( "Low_Neutron_Response" );
    other.push_back( "Mid_Neutron_Response" );
    other.push_back( "High_Neutron_Response" );
    other.push_back( "Pion_Response" );
    other.push_back( "EM_Response" );
    other.push_back( "Other_Response" );
    other.push_back( "Crosstalk" );
    other.push_back( "MEU_Recoil" );
    other.push_back( "Binding_Energy" );
    other.push_back( "Reweight_Pion" );
    other.push_back( "Reweight_Proton" );
    other.push_back( "Reweight_Neutron" );
    other.push_back( "Proton_Angle" );
    error_summary_group_map["Hadrons"] = other;
  }

     else if ( style == kCCQEAntiNuStyle) {

     ApplyStyle( kCompactStyle );

      error_color_map.clear();
    error_summary_group_map.clear();

    error_color_map["Flux"]                    = kViolet+6;
    error_color_map["Recoil Reconstruction"]   = kOrange+2; 
    error_color_map["Cross Section Models"]         = kMagenta;
    error_color_map["FSI Models"]              = kRed; 
    error_color_map["Muon Reconstruction"]     = kOrange-3; 
    error_color_map["Other"]                  = kGreen+3; 
    error_color_map["Low Recoil Fits"]         = kRed+3;
    error_color_map["GEANT4"]                = kBlue;
    error_color_map["Background Subtraction"] = kGreen;

    std::vector< string > flux;
    flux.push_back("Flux");
    error_summary_group_map["Flux"] = flux;

    std::vector< string > xsection;
    xsection.push_back( "GENIE_CCQEPauliSupViaKF" );
    xsection.push_back( "GENIE_EtaNCEL" );
    xsection.push_back( "GENIE_MaCCQE" );
    //xsection.push_back( "GENIE_MaCCQEshape" );
    xsection.push_back( "GENIE_MaNCEL" );
    xsection.push_back( "GENIE_MaRES" );
    xsection.push_back( "GENIE_MvRES" );
    xsection.push_back( "GENIE_NormCCQE" );
    xsection.push_back( "GENIE_NormCCRES" );
    xsection.push_back( "GENIE_NormDISCC" );
    xsection.push_back( "GENIE_NormNCRES" );
    xsection.push_back( "GENIE_Rvn1pi" );
    xsection.push_back( "GENIE_Rvn2pi" );
    xsection.push_back( "GENIE_Rvp1pi" );
    xsection.push_back( "GENIE_Rvp2pi" );
    xsection.push_back( "GENIE_VecFFCCQEshape" );
    xsection.push_back( "GENIE_AhtBY" ); 
    xsection.push_back( "GENIE_BhtBY" ); 
    xsection.push_back( "GENIE_CV1uBY" ); 
    xsection.push_back( "GENIE_CV2uBY" ); 
    xsection.push_back( "RPA_HighQ2" );
    xsection.push_back( "RPA_LowQ2" );
    error_summary_group_map["Cross Section Models"] = xsection;

    std::vector< string > fsi;
    fsi.push_back( "GENIE_AGKYxF1pi" );
    fsi.push_back( "GENIE_FrAbs_N" );
    fsi.push_back( "GENIE_FrAbs_pi" );
    fsi.push_back( "GENIE_FrCEx_N" );
    fsi.push_back( "GENIE_FrCEx_pi" );
    fsi.push_back( "GENIE_FrElas_N" );
    fsi.push_back( "GENIE_FrElas_pi" );
    fsi.push_back( "GENIE_FrInel_N" );
    fsi.push_back( "GENIE_FrInel_pi" );
    fsi.push_back( "GENIE_FrPiProd_N" );
    fsi.push_back( "GENIE_FrPiProd_pi" );
    fsi.push_back( "GENIE_MFP_N" );
    fsi.push_back( "GENIE_MFP_pi" );
    fsi.push_back( "GENIE_RDecBR1gamma" );
    fsi.push_back( "GENIE_Theta_Delta2Npi" );
    error_summary_group_map["FSI Models"] = fsi;

    std::vector<string> sigfrac;
    sigfrac.push_back("SignalFraction_00");
    sigfrac.push_back("SignalFraction_01");
    sigfrac.push_back("SignalFraction_02");
    sigfrac.push_back("SignalFraction_03");
    sigfrac.push_back("SignalFraction_04");
    sigfrac.push_back("SignalFraction_05");
    sigfrac.push_back("SignalFraction_06");
    sigfrac.push_back("SignalFraction_07");
    sigfrac.push_back("SignalFraction_08");
    sigfrac.push_back("SignalFraction_09");
    sigfrac.push_back("SignalFraction_10");
    sigfrac.push_back("SignalFraction_11");
    sigfrac.push_back("SignalFraction_12");
    sigfrac.push_back("SignalFraction_13");
    error_summary_group_map["Background Subtraction"] = sigfrac;
    
    
    std::vector< string > lowrecoilfit;
    lowrecoilfit.push_back("Low_Recoil_2p2h_Tune");
    error_summary_group_map["Low Recoil Fits"] = lowrecoilfit;

    std::vector< string > muon;
    //    muon.push_back( "Muon_Energy" );
    muon.push_back( "Muon_Energy_MINOS" );
    muon.push_back( "Muon_Energy_MINERvA" );
    //    muon.push_back( "Muon_Theta" );
    muon.push_back( "BeamAngleX" );
    muon.push_back( "BeamAngleY" );
    muon.push_back( "Muon_Energy_Resolution" );
    error_summary_group_map["Muon Reconstruction"] = muon;
    
    std::vector< string > normalization;
    normalization.push_back("MINOS_Reconstruction_Efficiency");
    error_summary_group_map["Normalization"] = normalization;
    //normalization.push_back( "Target_Mass" );

    std::vector< string > recoil_reconstruct;
    //    other.push_back( "Normalization_Factors" );
    //other.push_back( "Michel_Efficiency");
    //other.push_back( "Target_Mass" );
    //other.push_back( "GENIE_nonreweightable" );
    //other.push_back( "Muon_Response" );
    recoil_reconstruct.push_back( "response_em" );
    recoil_reconstruct.push_back( "response_meson" );
    recoil_reconstruct.push_back( "response_proton" );
    recoil_reconstruct.push_back( "response_other" );
    //other.push_back("TrueProtonKECut");
    //other.push_back( "Pion_Response" );
    //other.push_back( "EM_Response" );
    //other.push_back( "Other_Response" );
    //other.push_back( "Crosstalk" );
    //other.push_back( "MEU_Recoil" );
    //other.push_back( "Binding_Energy" );
    //other.push_back( "Reweight_Pion" ); 
    //other.push_back( "Reweight_Proton" ); 
    //other.push_back( "Reweight_Neutron" ); 
    //other.push_back( "Bethe_Bloch" ); 
    //other.push_back( "MEU_Proton" ); 
    //other.push_back( "Mass_Model_Proton" ); 
    //other.push_back( "Birks_Response_Proton" ); 
    //other.push_back( "Proton_TrackEff" );
    error_summary_group_map["Recoil Reconstruction"] = recoil_reconstruct;


    std::vector<string>geant4;
    geant4.push_back("GEANT_Neutron");
    geant4.push_back("GEANT_Proton");
    geant4.push_back("GEANT_Pion");
    error_summary_group_map["GEANT"] = geant4;
  
    std::vector<string>other;
    other.push_back("TrueProtonKECut");
    other.push_back( "GENIE_nonreweightable" );
    other.push_back("TrueProtonKECut");
    other.push_back("RecoProtonKECut");
    other.push_back("bethe_bloch");
    other.push_back("birks_response_proton");
    other.push_back("mass_model_proton");
    other.push_back("meu_proton");
    error_summary_group_map["Other"] = other;


  }
else if(style == kCCInclusiveHeliumStyle){
  ApplyStyle( kCompactStyle );
  //-- define colors of the standard errors
  draw_normalized_to_bin_width = false; //Bin width norm inside plotting function because we have substraction thus we bin_width normalize pot scale than substract  


  error_color_map.clear();
  error_summary_group_map.clear();

  error_color_map["Flux"]                    = kViolet+6;
  //error_color_map["Recoil Reconstruction"]   = kOrange+2;
  error_color_map["Beam Angle"]   = kOrange+2;
  error_color_map["GENIE Cross Section"]         = kMagenta;
  error_color_map["GENIE FSI"]              = kRed;
  error_color_map["2p2h RPA Mvn1Tune"]         = kRed+3;
  error_color_map["Muon Reconstruction"]     = kOrange-3;
  error_color_map["Others"]                  = kGreen+3;


  std::vector< string > flux;
  flux.push_back("Flux");
  flux.push_back("Flux_Tertiary");
  flux.push_back("Flux_BeamFocus");
  flux.push_back("Flux_NA49");
  error_summary_group_map["Flux"] = flux;

  std::vector< string > xsection;
  xsection.push_back( "GENIE_CCQEPauliSupViaKF" );
  xsection.push_back( "GENIE_EtaNCEL" );
  xsection.push_back( "GENIE_MaCCQE" );
  xsection.push_back( "GENIE_MaCCQEshape" );
  xsection.push_back( "GENIE_MaNCEL" );
  xsection.push_back( "GENIE_MaRES" );
  xsection.push_back( "GENIE_MvRES" );
  xsection.push_back( "GENIE_NormCCQE" );
  xsection.push_back( "GENIE_NormCCRES" );
  xsection.push_back( "GENIE_NormDISCC" );
  xsection.push_back( "GENIE_NormNCRES" );
  xsection.push_back( "GENIE_Rvn1pi" );
  xsection.push_back( "GENIE_Rvn2pi" );
  xsection.push_back( "GENIE_Rvp1pi" );
  xsection.push_back( "GENIE_Rvp2pi" );
  xsection.push_back( "GENIE_VecFFCCQEshape" );
  xsection.push_back( "GENIE_AhtBY" );
  xsection.push_back( "GENIE_BhtBY" );
  xsection.push_back( "GENIE_CV1uBY" );
  xsection.push_back( "GENIE_CV2uBY" );
  //xsection.push_back( "RPA_HighQ2" );
  //xsection.push_back( "RPA_LowQ2" );
  error_summary_group_map["GENIE Cross Section"] = xsection;

  std::vector< string > fsi;
  fsi.push_back( "GENIE_AGKYxF1pi" );
  fsi.push_back( "GENIE_FrAbs_N" );
  fsi.push_back( "GENIE_FrAbs_pi" );
  fsi.push_back( "GENIE_FrCEx_N" );
  fsi.push_back( "GENIE_FrCEx_pi" );
  fsi.push_back( "GENIE_FrElas_N" );
  fsi.push_back( "GENIE_FrElas_pi" );
  fsi.push_back( "GENIE_FrInel_N" );
  fsi.push_back( "GENIE_FrInel_pi" );
  fsi.push_back( "GENIE_FrPiProd_N" );
  fsi.push_back( "GENIE_FrPiProd_pi" );
  fsi.push_back( "GENIE_MFP_N" );
  fsi.push_back( "GENIE_MFP_pi" );
  fsi.push_back( "GENIE_RDecBR1gamma" );
  fsi.push_back( "GENIE_Theta_Delta2Npi" );
  error_summary_group_map["GENIE FSI"] = fsi;

  std::vector< string > lowrecoilfit;
  lowrecoilfit.push_back("Low_Recoil_2p2h_Tune");
  lowrecoilfit.push_back("RPA_LowQ2");
  lowrecoilfit.push_back("RPA_HighQ2");
  error_summary_group_map["2p2h RPA Mvn1Tune"] = lowrecoilfit;

  std::vector< string > muon;
  //    muon.push_back( "Muon_Energy" );
  muon.push_back( "Muon_Energy_MINOS" );
  muon.push_back( "Muon_Energy_MINERvA" );
  //    muon.push_back( "Muon_Theta" );
  muon.push_back("MINOS_Reconstruction_Efficiency");
  muon.push_back( "Muon_Energy_Resolution" );
  error_summary_group_map["Muon Reconstruction"] = muon;

  std::vector< string > beam_angle;
  beam_angle.push_back( "BeamAngleX" );
  beam_angle.push_back( "BeamAngleY" );
  error_summary_group_map["Beam Angle"] = beam_angle;

  std::vector< string > other;
  //    other.push_back( "Normalization_Factors" );
  //other.push_back( "MINOS_Reconstruction_Efficiency" );
  //other.push_back( "Michel_Efficiency");
  other.push_back( "Target_Mass" );
  other.push_back( "GENIE_nonreweightable" );
  other.push_back( "Muon_Response" );
  other.push_back( "Proton_Response" );
  other.push_back( "Low_Neutron_Response" );
  other.push_back( "Mid_Neutron_Response" );
  other.push_back( "High_Neutron_Response" );
  other.push_back( "Pion_Response" );
  other.push_back( "EM_Response" );
  other.push_back( "Other_Response" );
  other.push_back( "Crosstalk" );
  //other.push_back( "MEU_Recoil" );
  other.push_back( "Binding_Energy" );
  other.push_back( "Reweight_Pion" );
  other.push_back( "Reweight_Proton" );
  other.push_back( "Reweight_Neutron" );
  other.push_back( "Bethe_Bloch" );
  other.push_back( "MEU_Proton" );
  other.push_back( "Mass_Model_Proton" );
  other.push_back( "Birks_Response_Proton" );
  other.push_back( "Proton_TrackEff" );
  error_summary_group_map["Others"] = other;




}//end of helium style
  else if ( style == kIMDStyle) {
    ApplyStyle( kCompactStyle );
    //-- define colors of the standard errors
    error_color_map.clear();
    error_summary_group_map.clear();

    error_color_map["Flux"]                    = kViolet+6;
    error_color_map["Cross Section Models"]         = kMagenta;
    error_color_map["FSI Models"]              = kRed;
    error_color_map["Muon Reconstruction"]     = kOrange-3;
    error_color_map["Others"]                  = kGreen+3;

    std::vector< string > flux;
    flux.push_back("Flux");
    error_summary_group_map["Flux"] = flux;

    std::vector< string > xsection;
    xsection.push_back( "GENIE_CCQEPauliSupViaKF" );
    xsection.push_back( "GENIE_EtaNCEL" );
    //xsection.push_back( "GENIE_MaCCQE" );
    xsection.push_back( "GENIE_MaCCQEshape" );
    xsection.push_back( "GENIE_MaNCEL" );
    xsection.push_back( "GENIE_MaRES" );
    xsection.push_back( "GENIE_MvRES" );
    xsection.push_back( "GENIE_NormCCQE" );
    xsection.push_back( "GENIE_NormCCRES" );
    xsection.push_back( "GENIE_NormDISCC" );
    xsection.push_back( "GENIE_NormNCRES" );
    xsection.push_back( "GENIE_Rvn1pi" );
    xsection.push_back( "GENIE_Rvn2pi" );
    xsection.push_back( "GENIE_Rvp1pi" );
    xsection.push_back( "GENIE_Rvp2pi" );
    xsection.push_back( "GENIE_VecFFCCQEshape" );
    xsection.push_back( "GENIE_AhtBY" );
    xsection.push_back( "GENIE_BhtBY" );
    xsection.push_back( "GENIE_CV1uBY" );
    xsection.push_back( "GENIE_CV2uBY" );
    xsection.push_back( "RPA_HighQ2" );
    xsection.push_back( "RPA_LowQ2" );
    xsection.push_back( "ZExpansion" );
    error_summary_group_map["Cross Section Models"] = xsection;

    std::vector< string > fsi;
    fsi.push_back( "GENIE_AGKYxF1pi" );
    fsi.push_back( "GENIE_FrAbs_N" );
    fsi.push_back( "GENIE_FrAbs_pi" );
    fsi.push_back( "GENIE_FrCEx_N" );
    fsi.push_back( "GENIE_FrCEx_pi" );
    fsi.push_back( "GENIE_FrElas_N" );
    fsi.push_back( "GENIE_FrElas_pi" );
    fsi.push_back( "GENIE_FrInel_N" );
    fsi.push_back( "GENIE_FrInel_pi" );
    fsi.push_back( "GENIE_FrPiProd_N" );
    fsi.push_back( "GENIE_FrPiProd_pi" );
    fsi.push_back( "GENIE_MFP_N" );
    fsi.push_back( "GENIE_MFP_pi" );
    fsi.push_back( "GENIE_RDecBR1gamma" );
    fsi.push_back( "GENIE_Theta_Delta2Npi" );
    error_summary_group_map["FSI Models"] = fsi;

    std::vector< string > muon;
    //    muon.push_back( "Muon_Energy" );
    muon.push_back( "Muon_Energy_MINOS" );
    muon.push_back( "Muon_Energy_MINERvA" );
    //    muon.push_back( "Muon_Theta" );
    muon.push_back( "BeamAngleX" );
    muon.push_back( "BeamAngleY" );
    muon.push_back( "Muon_Energy_Resolution" );
    error_summary_group_map["Muon Reconstruction"] = muon;

    std::vector< string > other;
    other.push_back("Low_Recoil_2p2h_Tune");
    //    other.push_back( "Normalization_Factors" );
    other.push_back( "MINOS_Reconstruction_Efficiency" );
    other.push_back( "Michel_Efficiency");
    other.push_back( "Target_Mass" );
    other.push_back( "GENIE_nonreweightable" );
    other.push_back( "Muon_Response" );
    other.push_back( "Proton_Response" );
    other.push_back( "Low_Neutron_Response" );
    other.push_back( "Mid_Neutron_Response" );
    other.push_back( "High_Neutron_Response" );
    other.push_back( "Pion_Response" );
    other.push_back( "EM_Response" );
    other.push_back( "Other_Response" );
    other.push_back( "Crosstalk" );
    other.push_back( "MEU_Recoil" );
    other.push_back( "Binding_Energy" );
    other.push_back( "Reweight_Pion" );
    other.push_back( "Reweight_Proton" );
    other.push_back( "Reweight_Neutron" );
    other.push_back( "Bethe_Bloch" );
    other.push_back( "MEU_Proton" );
    other.push_back( "Mass_Model_Proton" );
    other.push_back( "Birks_Response_Proton" );
    other.push_back( "Proton_TrackEff" );
    other.push_back( "Proton_Angle" );
    error_summary_group_map["Others"] = other;
  }

  else if ( style == kCCQENuTransverseStyle) {
    ApplyStyle( kCompactStyle );
    //-- define colors of the standard errors
    error_color_map.clear();
    error_summary_group_map.clear();

    error_color_map["Flux"]                    = kViolet+6;
    error_color_map["Recoil Reconstruction"]   = kOrange+2;
    error_color_map["Cross Section Models"]         = kMagenta;
    error_color_map["FSI Models"]              = kRed;
    error_color_map["Muon Reconstruction"]     = kOrange-3;
    error_color_map["Proton Reconstruction"]     = kAzure+8;
    error_color_map["Others"]                  = kGreen+3;
    error_color_map["Low Recoil Fits"]         = kRed+3;

    std::vector< string > flux;
    flux.push_back("Flux");
    error_summary_group_map["Flux"] = flux;

    std::vector< string > xsection;
    xsection.push_back( "GENIE_CCQEPauliSupViaKF" );
    xsection.push_back( "GENIE_EtaNCEL" );
    //xsection.push_back( "GENIE_MaCCQE" );
    xsection.push_back( "GENIE_MaCCQEshape" );
    xsection.push_back( "GENIE_MaNCEL" );
    xsection.push_back( "GENIE_MaRES" );
    xsection.push_back( "GENIE_MvRES" );
    xsection.push_back( "GENIE_NormCCQE" );
    xsection.push_back( "GENIE_NormCCRES" );
    xsection.push_back( "GENIE_NormDISCC" );
    xsection.push_back( "GENIE_NormNCRES" );
    xsection.push_back( "GENIE_Rvn1pi" );
    xsection.push_back( "GENIE_Rvn2pi" );
    xsection.push_back( "GENIE_Rvp1pi" );
    xsection.push_back( "GENIE_Rvp2pi" );
    xsection.push_back( "GENIE_VecFFCCQEshape" );
    xsection.push_back( "GENIE_AhtBY" );
    xsection.push_back( "GENIE_BhtBY" );
    xsection.push_back( "GENIE_CV1uBY" );
    xsection.push_back( "GENIE_CV2uBY" );
    xsection.push_back( "RPA_HighQ2" );
    xsection.push_back( "RPA_LowQ2" );
    error_summary_group_map["Cross Section Models"] = xsection;

    std::vector< string > fsi;
    fsi.push_back( "GENIE_AGKYxF1pi" );
    fsi.push_back( "GENIE_FrAbs_N" );
    fsi.push_back( "GENIE_FrAbs_pi" );
    fsi.push_back( "GENIE_FrCEx_N" );
    fsi.push_back( "GENIE_FrCEx_pi" );
    fsi.push_back( "GENIE_FrElas_N" );
    fsi.push_back( "GENIE_FrElas_pi" );
    fsi.push_back( "GENIE_FrInel_N" );
    fsi.push_back( "GENIE_FrInel_pi" );
    fsi.push_back( "GENIE_FrPiProd_N" );
    fsi.push_back( "GENIE_FrPiProd_pi" );
    fsi.push_back( "GENIE_MFP_N" );
    fsi.push_back( "GENIE_MFP_pi" );
    fsi.push_back( "GENIE_RDecBR1gamma" );
    fsi.push_back( "GENIE_Theta_Delta2Npi" );
    error_summary_group_map["FSI Models"] = fsi;

    std::vector< string > lowrecoilfit;
    lowrecoilfit.push_back("Low_Recoil_2p2h_Tune");
    error_summary_group_map["Low Recoil Fits"] = lowrecoilfit;

    std::vector< string > muon;
    //    muon.push_back( "Muon_Energy" );
    muon.push_back( "Muon_Energy_MINOS" );
    muon.push_back( "Muon_Energy_MINERvA" );
    //    muon.push_back( "Muon_Theta" );
    muon.push_back( "BeamAngleX" );
    muon.push_back( "BeamAngleY" );
    muon.push_back( "Muon_Energy_Resolution" );
    error_summary_group_map["Muon Reconstruction"] = muon;

    std::vector< string > proton;
    proton.push_back( "Bethe_Bloch" );
    proton.push_back( "MEU_Proton" );
    proton.push_back( "Mass_Model_Proton" );
    proton.push_back( "Birks_Response_Proton" );
    error_summary_group_map["Proton Reconstruction"] = proton;

    std::vector< string > other;
    //    other.push_back( "Normalization_Factors" );
    other.push_back( "MINOS_Reconstruction_Efficiency" );
    other.push_back( "Michel_Efficiency");
    other.push_back( "Target_Mass" );
    other.push_back( "GENIE_nonreweightable" );
    other.push_back( "Muon_Response" );
    other.push_back( "Proton_Response" );
    other.push_back( "Low_Neutron_Response" );
    other.push_back( "Mid_Neutron_Response" );
    other.push_back( "High_Neutron_Response" );
    other.push_back( "Pion_Response" );
    other.push_back( "EM_Response" );
    other.push_back( "Other_Response" );
    other.push_back( "Crosstalk" );
    other.push_back( "MEU_Recoil" );
    other.push_back( "Binding_Energy" );
    other.push_back( "Reweight_Pion" );
    other.push_back( "Reweight_Proton" );
    other.push_back( "Reweight_Neutron" );
    other.push_back( "Proton_TrackEff" );
    other.push_back( "Proton_Angle" );
    error_summary_group_map["Others"] = other;
  }
  
  else if ( style == kCCPi0AnaStyle ) {
    // Apply compact style first
    ApplyStyle(kCompactStyle);
    
    // Keep this as default unless otherwise
    draw_normalized_to_bin_width = true;
    
    // MC and data histogram settings
    mc_line_width     = 3;
    data_line_width   = 3;
    data_marker_size  = 1.3;
    ratio_marker_size = 1.3;
    
    // Axis options
//     hist_min_zero     = true;
    axis_draw_grid_x    = true;
    axis_draw_grid_y    = true;
    axis_max_digits     = 3;
    axis_title_font_x   = 62;
    axis_title_font_y   = 62;
    axis_title_font_z   = 62;
    axis_title_offset_x = 1.15;
    axis_title_offset_y = 1.15;
    axis_title_offset_z = 0.9;
    axis_title_size_x   = 0.05;
    axis_title_size_y   = 0.05;
    axis_title_size_z   = 0.05;
    axis_label_font     = 42;
    axis_label_size     = 0.045;
//     axis_minimum      = MnvHist::AutoAxisLimit;
//     axis_maximum      = MnvHist::AutoAxisLimit;
//     axis_maximum_group= MnvHist::AutoAxisLimit; //0.5;
    
    // Legend settings
    height_nspaces_per_hist = 1.0;
    width_xspace_per_letter = 0.4;
    legend_border_size      = 1;
    legend_fill_color       = 10;
    legend_text_size        = 0.03;
    legend_offset_x         = 0.0;
    legend_offset_y         = 0.0;
    legend_n_columns        = 1;
    legend_text_font        = 62;
    
    // Extra margin
    extra_right_margin = -1.0;
    
    // Systematics map
    error_summary_group_map.clear();
    
    std::vector<std::string> flux;
    flux.push_back("Flux");
    error_summary_group_map["Neutrino Flux"] = flux;
    
    std::vector<std::string> genie_interaction_model;
    genie_interaction_model.push_back("GENIE_AhtBY");
    genie_interaction_model.push_back("GENIE_BhtBY");
    genie_interaction_model.push_back("GENIE_CCQEPauliSupViaKF");
    genie_interaction_model.push_back("GENIE_CV1uBY");
    genie_interaction_model.push_back("GENIE_CV2uBY");
    genie_interaction_model.push_back("GENIE_MaCCQE");
    genie_interaction_model.push_back("GENIE_MaNCEL");
    genie_interaction_model.push_back("GENIE_MaRES");
    genie_interaction_model.push_back("GENIE_MvRES");
    genie_interaction_model.push_back("GENIE_EtaNCEL");
    genie_interaction_model.push_back("GENIE_NormDISCC");
    genie_interaction_model.push_back("GENIE_NormNCRES");
    genie_interaction_model.push_back("GENIE_Rvn1pi");
    genie_interaction_model.push_back("GENIE_Rvn2pi");
    genie_interaction_model.push_back("GENIE_Rvp1pi");
    genie_interaction_model.push_back("GENIE_Rvp2pi");
    genie_interaction_model.push_back("GENIE_VecFFCCQEshape");
    genie_interaction_model.push_back("GENIE_D2_MaRES");
    genie_interaction_model.push_back("GENIE_D2_NormCCRES");
    genie_interaction_model.push_back("GENIE_EP_MvRES");
    error_summary_group_map["GENIE Interaction Models"] = genie_interaction_model;
    
    std::vector<std::string> genie_nucleon_fsi;
    genie_nucleon_fsi.push_back("GENIE_FrAbs_N");
    genie_nucleon_fsi.push_back("GENIE_FrCEx_N");
    genie_nucleon_fsi.push_back("GENIE_FrElas_N");
    genie_nucleon_fsi.push_back("GENIE_FrInel_N");
    genie_nucleon_fsi.push_back("GENIE_FrPiProd_N");
    genie_nucleon_fsi.push_back("GENIE_MFP_N");
    error_summary_group_map["GENIE Nucleon FSI"] = genie_nucleon_fsi;
    
    std::vector<std::string> genie_pion_fsi;
    genie_pion_fsi.push_back("GENIE_AGKYxF1pi");
    genie_pion_fsi.push_back("GENIE_FrAbs_pi");
    genie_pion_fsi.push_back("GENIE_FrCEx_pi");
    genie_pion_fsi.push_back("GENIE_FrElas_pi");
    genie_pion_fsi.push_back("GENIE_FrPiProd_pi");
    genie_pion_fsi.push_back("GENIE_MFP_pi");
    genie_pion_fsi.push_back("GENIE_RDecBR1gamma");
    genie_pion_fsi.push_back("GENIE_Theta_Delta2Npi");
    error_summary_group_map["GENIE Pion FSI"] = genie_pion_fsi;
    
    std::vector<std::string> mnvgenie;
    mnvgenie.push_back("Low_Recoil_2p2h_Tune");
    mnvgenie.push_back("RPA_HighQ2");
    mnvgenie.push_back("RPA_LowQ2");
    mnvgenie.push_back("LowQ2Pi");
    error_summary_group_map["MnvGENIE Tune"] = mnvgenie;
    
    std::vector<std::string> muon;
    muon.push_back("Muon_Energy_MINERvA");
    muon.push_back("Muon_Energy_MINOS");
    muon.push_back("Muon_Energy_Resolution");
    muon.push_back("MuonAngleXResolution");
    muon.push_back("MuonAngleYResolution");
    muon.push_back("BeamAngleX");
    muon.push_back("BeamAngleY");
    muon.push_back("MINOS_Reconstruction_Efficiency");
    error_summary_group_map["Muon Reconstruction"] = muon;
    
    std::vector<std::string> detector;
    detector.push_back("GEANT_Neutron");
    detector.push_back("GEANT_Pion");
    detector.push_back("GEANT_Proton");
    detector.push_back("MichelEfficiency");
    error_summary_group_map["Detector Model"] = detector;
    
    std::vector<std::string> other;
    other.push_back("Target_Mass_CH");
    other.push_back("Target_Mass_C");
    other.push_back("Target_Mass_H2O");
    other.push_back("Target_Mass_Fe");
    other.push_back("Target_Mass_Pb");
    error_summary_group_map["Other"] = other;
    
    // Systematics color scheme
    error_color_map.clear();
    error_color_map["Neutrino Flux"]            = kRed+2;
    error_color_map["GENIE Interaction Models"] = kGreen+1;
    error_color_map["GENIE Nucleon FSI"]        = kBlue+2;
    error_color_map["GENIE Pion FSI"]           = kMagenta+1;
    error_color_map["MnvGENIE Tune"]            = kOrange+2;
    error_color_map["Muon Reconstruction"]      = kCyan+2;
    error_color_map["Detector Model"]           = kViolet+6;
    error_color_map["Other"]                    = kPink+2;
  }

  else {
    Error( "ApplyStyle", "This plot style is not recognized.  Using the default: kDefaultStyle" );
    ApplyStyle( kDefaultStyle );
  }

  //set the environment to apply styles
  SetRootEnv();
}

//================================================================
// set the number of columns to use in legends
//================================================================
void MnvPlotter::SetLegendNColumns( int user_legend_n_columns )
{
  legend_n_columns = user_legend_n_columns;
}

//=========================================================
// Set the global color palette to hot/cold style
//=========================================================
void MnvPlotter::SetCorrelationPalette() const
{
    // A colour palette that goes blue->white->red, useful for
    // correlation matrices
    const int NRGBs = 3;
    static bool initialized=false;
    static int* colors=new int[n_color_contours];

    if (!initialized) {
        gStyle->SetNumberContours(n_color_contours);
        Double_t stops[NRGBs] = { 0.00, 0.50, 1.00};
        Double_t red[NRGBs]   = { 0.00, 1.00, 1.00};
        Double_t green[NRGBs] = { 0.00, 1.00, 0.00};
        Double_t blue[NRGBs]  = { 1.00, 1.00, 0.00};
        int colmin=TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, n_color_contours);
        for (uint i=0; i<n_color_contours; ++i) colors[i]=colmin+i;

        initialized=true;
    }
    gStyle->SetNumberContours(n_color_contours);
    gStyle->SetPalette(n_color_contours, colors);
}

void MnvPlotter::SetWhiteRainbowPalette() const
{
    const int NRGBs = 7;
    static bool initialized=false;
    static int* colors=new int[n_color_contours];

    if (!initialized) {
        gStyle->SetNumberContours(n_color_contours);
        Double_t stops[NRGBs] = { 0.00, 0.05, 0.23, 0.45, 0.60, 0.85, 1.00 };
        Double_t red[NRGBs]   = { 1.00, 0.00, 0.00, 0.00, 1.00, 1.00, 0.33 };
        Double_t green[NRGBs] = { 1.00, 1.00, 0.30, 0.40, 1.00, 0.00, 0.00 };
        Double_t blue[NRGBs]  = { 1.00, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00 };
        int colmin=TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, n_color_contours);
        for (uint i=0; i<n_color_contours; ++i) colors[i]=colmin+i;

        initialized=true;
    }
    gStyle->SetNumberContours(n_color_contours);
    gStyle->SetPalette(n_color_contours, colors);
}

void MnvPlotter::SetRedHeatPalette() const
{
    const int NRGBs = 9;
    static bool initialized=false;
    static int* colors=new int[n_color_contours];

    if (!initialized) {
        // White -> red
        Double_t stops[NRGBs] = { 0.00, 0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1.000};
        Double_t red[NRGBs]   = { 1.00, 1.00, 0.99, 0.99, 0.98, 0.94, 0.80, 0.65, 0.40 };
        Double_t green[NRGBs] = { 0.96, 0.88, 0.73, 0.57, 0.42, 0.23, 0.09, 0.06, 0.00 };
        Double_t blue[NRGBs]  = { 0.94, 0.82, 0.63, 0.45, 0.29, 0.17, 0.11, 0.08, 0.05 };
        int colmin=TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, n_color_contours);
        for (uint i=0; i<n_color_contours; ++i) colors[i]=colmin+i;

        initialized=true;
    }
    gStyle->SetNumberContours(n_color_contours);
    gStyle->SetPalette(n_color_contours, colors);
}

void MnvPlotter::SetBlackbodyPalette() const
{
    // Available as gStyle->SetPalette(56) in sufficiently recent versions of ROOT, but not ours
    const int nRGBs = 5;
    const int NCont = 99;
    static bool initialized=false;
    static int colors[99];

    if (!initialized) {
        gStyle->SetNumberContours(NCont);
        Double_t stops[nRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00};
        Double_t red[nRGBs] = { 0.00, 0.50, 1.00, 1.00, 1.00};
        Double_t green[nRGBs] = { 0.00, 0.00, 0.55, 1.00, 1.00};
        Double_t blue[nRGBs] = { 0.00, 0.00, 0.00, 0.00, 1.00};

        int colmin=TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, NCont);
        for (int i=0; i<NCont; ++i) colors[i]=colmin+i;

        initialized=true;
    }
    gStyle->SetNumberContours(NCont);
    gStyle->SetPalette(NCont, colors);
}

//==========================================================================================================
// Set the global color palette to any ROOT 6 option (see: https://root.cern.ch/doc/master/classTColor.html)
//
// Code copied wholesale from ROOT 6's TColor
//==========================================================================================================
void MnvPlotter::SetROOT6Palette(Int_t ncolors)
{
   Int_t i;

   TArrayD fgPalettesList(0);
   static Int_t paletteType = 0;

   // High quality palettes (255 levels)
   if (ncolors>50) {

      if (!fgPalettesList.fN) fgPalettesList.Set(62);        // Right now 62 high quality palettes
      Int_t Idx = (Int_t)fgPalettesList.fArray[ncolors-51];  // High quality palettes indices start at 51

      // This high quality palette has already been created. Reuse it.
      if (Idx > 0) {
         if (paletteType == ncolors) return; // The current palette is already this one.

         const int NCont=255;
         int colors[NCont];
         for (i=0;i<NCont;i++) colors[i] = Idx+i;
         gStyle->SetNumberContours(NCont);
         gStyle->SetPalette(NCont, colors);

         return;
      }

      TColor::InitializeColors();
      Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};

      switch (ncolors) {
      // Deep Sea
      case 51:
         {
            Double_t red[9]   = {  0./255.,  9./255., 13./255., 17./255., 24./255.,  32./255.,  27./255.,  25./255.,  29./255.};
            Double_t green[9] = {  0./255.,  0./255.,  0./255.,  2./255., 37./255.,  74./255., 113./255., 160./255., 221./255.};
            Double_t blue[9]  = { 28./255., 42./255., 59./255., 78./255., 98./255., 129./255., 154./255., 184./255., 221./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Grey Scale
      case 52:
         {
            Double_t red[9]   = { 0./255., 32./255., 64./255., 96./255., 128./255., 160./255., 192./255., 224./255., 255./255.};
            Double_t green[9] = { 0./255., 32./255., 64./255., 96./255., 128./255., 160./255., 192./255., 224./255., 255./255.};
            Double_t blue[9]  = { 0./255., 32./255., 64./255., 96./255., 128./255., 160./255., 192./255., 224./255., 255./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Dark Body Radiator
      case 53:
         {
            Double_t red[9]   = { 0./255., 45./255., 99./255., 156./255., 212./255., 230./255., 237./255., 234./255., 242./255.};
            Double_t green[9] = { 0./255.,  0./255.,  0./255.,  45./255., 101./255., 168./255., 238./255., 238./255., 243./255.};
            Double_t blue[9]  = { 0./255.,  1./255.,  1./255.,   3./255.,   9./255.,   8./255.,  11./255.,  95./255., 230./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Two-color hue (dark blue through neutral gray to bright yellow)
      case 54:
         {
            Double_t red[9]   = {  0./255.,  22./255., 44./255., 68./255., 93./255., 124./255., 160./255., 192./255., 237./255.};
            Double_t green[9] = {  0./255.,  16./255., 41./255., 67./255., 93./255., 125./255., 162./255., 194./255., 241./255.};
            Double_t blue[9]  = { 97./255., 100./255., 99./255., 99./255., 93./255.,  68./255.,  44./255.,  26./255.,  74./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Rain Bow
      case 55:
         {
            Double_t red[9]   = {  0./255.,   5./255.,  15./255.,  35./255., 102./255., 196./255., 208./255., 199./255., 110./255.};
            Double_t green[9] = {  0./255.,  48./255., 124./255., 192./255., 206./255., 226./255.,  97./255.,  16./255.,   0./255.};
            Double_t blue[9]  = { 99./255., 142./255., 198./255., 201./255.,  90./255.,  22./255.,  13./255.,   8./255.,   2./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Inverted Dark Body Radiator
      case 56:
         {
            Double_t red[9]   = { 242./255., 234./255., 237./255., 230./255., 212./255., 156./255., 99./255., 45./255., 0./255.};
            Double_t green[9] = { 243./255., 238./255., 238./255., 168./255., 101./255.,  45./255.,  0./255.,  0./255., 0./255.};
            Double_t blue[9]  = { 230./255.,  95./255.,  11./255.,   8./255.,   9./255.,   3./255.,  1./255.,  1./255., 0./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Bird
      case 57:
         {
            Double_t red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
            Double_t green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
            Double_t blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Cubehelix
      case 58:
         {
            Double_t red[9]   = { 0.0000, 0.0956, 0.0098, 0.2124, 0.6905, 0.9242, 0.7914, 0.7596, 1.0000};
            Double_t green[9] = { 0.0000, 0.1147, 0.3616, 0.5041, 0.4577, 0.4691, 0.6905, 0.9237, 1.0000};
            Double_t blue[9]  = { 0.0000, 0.2669, 0.3121, 0.1318, 0.2236, 0.6741, 0.9882, 0.9593, 1.0000};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Green Red Violet
      case 59:
         {
            Double_t red[9]   = {13./255., 23./255., 25./255., 63./255., 76./255., 104./255., 137./255., 161./255., 206./255.};
            Double_t green[9] = {95./255., 67./255., 37./255., 21./255.,  0./255.,  12./255.,  35./255.,  52./255.,  79./255.};
            Double_t blue[9]  = { 4./255.,  3./255.,  2./255.,  6./255., 11./255.,  22./255.,  49./255.,  98./255., 208./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Blue Red Yellow
      case 60:
         {
            Double_t red[9]   = {0./255.,  61./255.,  89./255., 122./255., 143./255., 160./255., 185./255., 204./255., 231./255.};
            Double_t green[9] = {0./255.,   0./255.,   0./255.,   0./255.,  14./255.,  37./255.,  72./255., 132./255., 235./255.};
            Double_t blue[9]  = {0./255., 140./255., 224./255., 144./255.,   4./255.,   5./255.,   6./255.,   9./255.,  13./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Ocean
      case 61:
         {
            Double_t red[9]   = { 14./255.,  7./255.,  2./255.,  0./255.,  5./255.,  11./255.,  55./255., 131./255., 229./255.};
            Double_t green[9] = {105./255., 56./255., 26./255.,  1./255., 42./255.,  74./255., 131./255., 171./255., 229./255.};
            Double_t blue[9]  = {  2./255., 21./255., 35./255., 60./255., 92./255., 113./255., 160./255., 185./255., 229./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Color Printable On Grey
      case 62:
         {
            Double_t red[9]   = { 0./255.,   0./255.,   0./255.,  70./255., 148./255., 231./255., 235./255., 237./255., 244./255.};
            Double_t green[9] = { 0./255.,   0./255.,   0./255.,   0./255.,   0./255.,  69./255.,  67./255., 216./255., 244./255.};
            Double_t blue[9]  = { 0./255., 102./255., 228./255., 231./255., 177./255., 124./255., 137./255.,  20./255., 244./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Alpine
      case 63:
         {
            Double_t red[9]   = { 50./255., 56./255., 63./255., 68./255.,  93./255., 121./255., 165./255., 192./255., 241./255.};
            Double_t green[9] = { 66./255., 81./255., 91./255., 96./255., 111./255., 128./255., 155./255., 189./255., 241./255.};
            Double_t blue[9]  = { 97./255., 91./255., 75./255., 65./255.,  77./255., 103./255., 143./255., 167./255., 217./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Aquamarine
      case 64:
         {
            Double_t red[9]   = { 145./255., 166./255., 167./255., 156./255., 131./255., 114./255., 101./255., 112./255., 132./255.};
            Double_t green[9] = { 158./255., 178./255., 179./255., 181./255., 163./255., 154./255., 144./255., 152./255., 159./255.};
            Double_t blue[9]  = { 190./255., 199./255., 201./255., 192./255., 176./255., 169./255., 160./255., 166./255., 190./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Army
      case 65:
         {
            Double_t red[9]   = { 93./255.,   91./255.,  99./255., 108./255., 130./255., 125./255., 132./255., 155./255., 174./255.};
            Double_t green[9] = { 126./255., 124./255., 128./255., 129./255., 131./255., 121./255., 119./255., 153./255., 173./255.};
            Double_t blue[9]  = { 103./255.,  94./255.,  87./255.,  85./255.,  80./255.,  85./255., 107./255., 120./255., 146./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Atlantic
      case 66:
         {
            Double_t red[9]   = { 24./255., 40./255., 69./255.,  90./255., 104./255., 114./255., 120./255., 132./255., 103./255.};
            Double_t green[9] = { 29./255., 52./255., 94./255., 127./255., 150./255., 162./255., 159./255., 151./255., 101./255.};
            Double_t blue[9]  = { 29./255., 52./255., 96./255., 132./255., 162./255., 181./255., 184./255., 186./255., 131./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Aurora
      case 67:
         {
            Double_t red[9]   = { 46./255., 38./255., 61./255., 92./255., 113./255., 121./255., 132./255., 150./255., 191./255.};
            Double_t green[9] = { 46./255., 36./255., 40./255., 69./255., 110./255., 135./255., 131./255.,  92./255.,  34./255.};
            Double_t blue[9]  = { 46./255., 80./255., 74./255., 70./255.,  81./255., 105./255., 165./255., 211./255., 225./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Avocado
      case 68:
         {
            Double_t red[9]   = { 0./255.,  4./255., 12./255.,  30./255.,  52./255., 101./255., 142./255., 190./255., 237./255.};
            Double_t green[9] = { 0./255., 40./255., 86./255., 121./255., 140./255., 172./255., 187./255., 213./255., 240./255.};
            Double_t blue[9]  = { 0./255.,  9./255., 14./255.,  18./255.,  21./255.,  23./255.,  27./255.,  35./255., 101./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Beach
      case 69:
         {
            Double_t red[9]   = { 198./255., 206./255., 206./255., 211./255., 198./255., 181./255., 161./255., 171./255., 244./255.};
            Double_t green[9] = { 103./255., 133./255., 150./255., 172./255., 178./255., 174./255., 163./255., 175./255., 244./255.};
            Double_t blue[9]  = {  49./255.,  54./255.,  55./255.,  66./255.,  91./255., 130./255., 184./255., 224./255., 244./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Black Body
      case 70:
         {
            Double_t red[9]   = { 243./255., 243./255., 240./255., 240./255., 241./255., 239./255., 186./255., 151./255., 129./255.};
            Double_t green[9] = {   0./255.,  46./255.,  99./255., 149./255., 194./255., 220./255., 183./255., 166./255., 147./255.};
            Double_t blue[9]  = {   6./255.,   8./255.,  36./255.,  91./255., 169./255., 235./255., 246./255., 240./255., 233./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Blue Green Yellow
      case 71:
         {
            Double_t red[9]   = { 22./255., 19./255.,  19./255.,  25./255.,  35./255.,  53./255.,  88./255., 139./255., 210./255.};
            Double_t green[9] = {  0./255., 32./255.,  69./255., 108./255., 135./255., 159./255., 183./255., 198./255., 215./255.};
            Double_t blue[9]  = { 77./255., 96./255., 110./255., 116./255., 110./255., 100./255.,  90./255.,  78./255.,  70./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Brown Cyan
      case 72:
         {
            Double_t red[9]   = { 68./255., 116./255., 165./255., 182./255., 189./255., 180./255., 145./255., 111./255.,  71./255.};
            Double_t green[9] = { 37./255.,  82./255., 135./255., 178./255., 204./255., 225./255., 221./255., 202./255., 147./255.};
            Double_t blue[9]  = { 16./255.,  55./255., 105./255., 147./255., 196./255., 226./255., 232./255., 224./255., 178./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // CMYK
      case 73:
         {
            Double_t red[9]   = {  61./255.,  99./255., 136./255., 181./255., 213./255., 225./255., 198./255., 136./255., 24./255.};
            Double_t green[9] = { 149./255., 140./255.,  96./255.,  83./255., 132./255., 178./255., 190./255., 135./255., 22./255.};
            Double_t blue[9]  = { 214./255., 203./255., 168./255., 135./255., 110./255., 100./255., 111./255., 113./255., 22./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Candy
      case 74:
         {
            Double_t red[9]   = { 76./255., 120./255., 156./255., 183./255., 197./255., 180./255., 162./255., 154./255., 140./255.};
            Double_t green[9] = { 34./255.,  35./255.,  42./255.,  69./255., 102./255., 137./255., 164./255., 188./255., 197./255.};
            Double_t blue[9]  = { 64./255.,  69./255.,  78./255., 105./255., 142./255., 177./255., 205./255., 217./255., 198./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Cherry
      case 75:
         {
            Double_t red[9]   = { 37./255., 102./255., 157./255., 188./255., 196./255., 214./255., 223./255., 235./255., 251./255.};
            Double_t green[9] = { 37./255.,  29./255.,  25./255.,  37./255.,  67./255.,  91./255., 132./255., 185./255., 251./255.};
            Double_t blue[9]  = { 37./255.,  32./255.,  33./255.,  45./255.,  66./255.,  98./255., 137./255., 187./255., 251./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Coffee
      case 76:
         {
            Double_t red[9]   = { 79./255., 100./255., 119./255., 137./255., 153./255., 172./255., 192./255., 205./255., 250./255.};
            Double_t green[9] = { 63./255.,  79./255.,  93./255., 103./255., 115./255., 135./255., 167./255., 196./255., 250./255.};
            Double_t blue[9]  = { 51./255.,  59./255.,  66./255.,  61./255.,  62./255.,  70./255., 110./255., 160./255., 250./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Dark Rain Bow
      case 77:
         {
            Double_t red[9]   = {  43./255.,  44./255., 50./255.,  66./255., 125./255., 172./255., 178./255., 155./255., 157./255.};
            Double_t green[9] = {  63./255.,  63./255., 85./255., 101./255., 138./255., 163./255., 122./255.,  51./255.,  39./255.};
            Double_t blue[9]  = { 121./255., 101./255., 58./255.,  44./255.,  47./255.,  55./255.,  57./255.,  44./255.,  43./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Dark Terrain
      case 78:
         {
            Double_t red[9]   = {  0./255., 41./255., 62./255., 79./255., 90./255., 87./255., 99./255., 140./255., 228./255.};
            Double_t green[9] = {  0./255., 57./255., 81./255., 93./255., 85./255., 70./255., 71./255., 125./255., 228./255.};
            Double_t blue[9]  = { 95./255., 91./255., 91./255., 82./255., 60./255., 43./255., 44./255., 112./255., 228./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Fall
      case 79:
         {
            Double_t red[9]   = { 49./255., 59./255., 72./255., 88./255., 114./255., 141./255., 176./255., 205./255., 222./255.};
            Double_t green[9] = { 78./255., 72./255., 66./255., 57./255.,  59./255.,  75./255., 106./255., 142./255., 173./255.};
            Double_t blue[9]  = { 78./255., 55./255., 46./255., 40./255.,  39./255.,  39./255.,  40./255.,  41./255.,  47./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Fruit Punch
      case 80:
         {
            Double_t red[9]   = { 243./255., 222./255., 201./255., 185./255., 165./255., 158./255., 166./255., 187./255., 219./255.};
            Double_t green[9] = {  94./255., 108./255., 132./255., 135./255., 125./255.,  96./255.,  68./255.,  51./255.,  61./255.};
            Double_t blue[9]  = {   7./255.,  9./255.,   12./255.,  19./255.,  45./255.,  89./255., 118./255., 146./255., 118./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Fuchsia
      case 81:
         {
            Double_t red[9]   = { 19./255., 44./255., 74./255., 105./255., 137./255., 166./255., 194./255., 206./255., 220./255.};
            Double_t green[9] = { 19./255., 28./255., 40./255.,  55./255.,  82./255., 110./255., 159./255., 181./255., 220./255.};
            Double_t blue[9]  = { 19./255., 42./255., 68./255.,  96./255., 129./255., 157./255., 188./255., 203./255., 220./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Grey Yellow
      case 82:
         {
            Double_t red[9]   = { 33./255., 44./255., 70./255.,  99./255., 140./255., 165./255., 199./255., 211./255., 216./255.};
            Double_t green[9] = { 38./255., 50./255., 76./255., 105./255., 140./255., 165./255., 191./255., 189./255., 167./255.};
            Double_t blue[9]  = { 55./255., 67./255., 97./255., 124./255., 140./255., 166./255., 163./255., 129./255.,  52./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Green Brown Terrain
      case 83:
         {
            Double_t red[9]   = { 0./255., 33./255., 73./255., 124./255., 136./255., 152./255., 159./255., 171./255., 223./255.};
            Double_t green[9] = { 0./255., 43./255., 92./255., 124./255., 134./255., 126./255., 121./255., 144./255., 223./255.};
            Double_t blue[9]  = { 0./255., 43./255., 68./255.,  76./255.,  73./255.,  64./255.,  72./255., 114./255., 223./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Green Pink
      case 84:
         {
            Double_t red[9]   = {  5./255.,  18./255.,  45./255., 124./255., 193./255., 223./255., 205./255., 128./255., 49./255.};
            Double_t green[9] = { 48./255., 134./255., 207./255., 230./255., 193./255., 113./255.,  28./255.,   0./255.,  7./255.};
            Double_t blue[9]  = {  6./255.,  15./255.,  41./255., 121./255., 193./255., 226./255., 208./255., 130./255., 49./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Island
      case 85:
         {
            Double_t red[9]   = { 180./255., 106./255., 104./255., 135./255., 164./255., 188./255., 189./255., 165./255., 144./255.};
            Double_t green[9] = {  72./255., 126./255., 154./255., 184./255., 198./255., 207./255., 205./255., 190./255., 179./255.};
            Double_t blue[9]  = {  41./255., 120./255., 158./255., 188./255., 194./255., 181./255., 145./255., 100./255.,  62./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Lake
      case 86:
         {
            Double_t red[9]   = {  57./255.,  72./255.,  94./255., 117./255., 136./255., 154./255., 174./255., 192./255., 215./255.};
            Double_t green[9] = {   0./255.,  33./255.,  68./255., 109./255., 140./255., 171./255., 192./255., 196./255., 209./255.};
            Double_t blue[9]  = { 116./255., 137./255., 173./255., 201./255., 200./255., 201./255., 203./255., 190./255., 187./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Light Temperature
      case 87:
         {
            Double_t red[9]   = {  31./255.,  71./255., 123./255., 160./255., 210./255., 222./255., 214./255., 199./255., 183./255.};
            Double_t green[9] = {  40./255., 117./255., 171./255., 211./255., 231./255., 220./255., 190./255., 132./255.,  65./255.};
            Double_t blue[9]  = { 234./255., 214./255., 228./255., 222./255., 210./255., 160./255., 105./255.,  60./255.,  34./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Light Terrain
      case 88:
         {
            Double_t red[9]   = { 123./255., 108./255., 109./255., 126./255., 154./255., 172./255., 188./255., 196./255., 218./255.};
            Double_t green[9] = { 184./255., 138./255., 130./255., 133./255., 154./255., 175./255., 188./255., 196./255., 218./255.};
            Double_t blue[9]  = { 208./255., 130./255., 109./255.,  99./255., 110./255., 122./255., 150./255., 171./255., 218./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Mint
      case 89:
         {
            Double_t red[9]   = { 105./255., 106./255., 122./255., 143./255., 159./255., 172./255., 176./255., 181./255., 207./255.};
            Double_t green[9] = { 252./255., 197./255., 194./255., 187./255., 174./255., 162./255., 153./255., 136./255., 125./255.};
            Double_t blue[9]  = { 146./255., 133./255., 144./255., 155./255., 163./255., 167./255., 166./255., 162./255., 174./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Neon
      case 90:
         {
            Double_t red[9]   = { 171./255., 141./255., 145./255., 152./255., 154./255., 159./255., 163./255., 158./255., 177./255.};
            Double_t green[9] = { 236./255., 143./255., 100./255.,  63./255.,  53./255.,  55./255.,  44./255.,  31./255.,   6./255.};
            Double_t blue[9]  = {  59./255.,  48./255.,  46./255.,  44./255.,  42./255.,  54./255.,  82./255., 112./255., 179./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Pastel
      case 91:
         {
            Double_t red[9]   = { 180./255., 190./255., 209./255., 223./255., 204./255., 228./255., 205./255., 152./255.,  91./255.};
            Double_t green[9] = {  93./255., 125./255., 147./255., 172./255., 181./255., 224./255., 233./255., 198./255., 158./255.};
            Double_t blue[9]  = { 236./255., 218./255., 160./255., 133./255., 114./255., 132./255., 162./255., 220./255., 218./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Pearl
      case 92:
         {
            Double_t red[9]   = { 225./255., 183./255., 162./255., 135./255., 115./255., 111./255., 119./255., 145./255., 211./255.};
            Double_t green[9] = { 205./255., 177./255., 166./255., 135./255., 124./255., 117./255., 117./255., 132./255., 172./255.};
            Double_t blue[9]  = { 186./255., 165./255., 155./255., 135./255., 126./255., 130./255., 150./255., 178./255., 226./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Pigeon
      case 93:
         {
            Double_t red[9]   = { 39./255., 43./255., 59./255., 63./255., 80./255., 116./255., 153./255., 177./255., 223./255.};
            Double_t green[9] = { 39./255., 43./255., 59./255., 74./255., 91./255., 114./255., 139./255., 165./255., 223./255.};
            Double_t blue[9]  = { 39./255., 50./255., 59./255., 70./255., 85./255., 115./255., 151./255., 176./255., 223./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Plum
      case 94:
         {
            Double_t red[9]   = { 0./255., 38./255., 60./255., 76./255., 84./255., 89./255., 101./255., 128./255., 204./255.};
            Double_t green[9] = { 0./255., 10./255., 15./255., 23./255., 35./255., 57./255.,  83./255., 123./255., 199./255.};
            Double_t blue[9]  = { 0./255., 11./255., 22./255., 40./255., 63./255., 86./255.,  97./255.,  94./255.,  85./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Red Blue
      case 95:
         {
            Double_t red[9]   = { 94./255., 112./255., 141./255., 165./255., 167./255., 140./255.,  91./255.,  49./255.,  27./255.};
            Double_t green[9] = { 27./255.,  46./255.,  88./255., 135./255., 166./255., 161./255., 135./255.,  97./255.,  58./255.};
            Double_t blue[9]  = { 42./255.,  52./255.,  81./255., 106./255., 139./255., 158./255., 155./255., 137./255., 116./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Rose
      case 96:
         {
            Double_t red[9]   = { 30./255., 49./255., 79./255., 117./255., 135./255., 151./255., 146./255., 138./255., 147./255.};
            Double_t green[9] = { 63./255., 60./255., 72./255.,  90./255.,  94./255.,  94./255.,  68./255.,  46./255.,  16./255.};
            Double_t blue[9]  = { 18./255., 28./255., 41./255.,  56./255.,  62./255.,  63./255.,  50./255.,  36./255.,  21./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Rust
      case 97:
         {
            Double_t red[9]   = {  0./255., 30./255., 63./255., 101./255., 143./255., 152./255., 169./255., 187./255., 230./255.};
            Double_t green[9] = {  0./255., 14./255., 28./255.,  42./255.,  58./255.,  61./255.,  67./255.,  74./255.,  91./255.};
            Double_t blue[9]  = { 39./255., 26./255., 21./255.,  18./255.,  15./255.,  14./255.,  14./255.,  13./255.,  13./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Sandy Terrain
      case 98:
         {
            Double_t red[9]   = { 149./255., 140./255., 164./255., 179./255., 182./255., 181./255., 131./255., 87./255., 61./255.};
            Double_t green[9] = {  62./255.,  70./255., 107./255., 136./255., 144./255., 138./255., 117./255., 87./255., 74./255.};
            Double_t blue[9]  = {  40./255.,  38./255.,  45./255.,  49./255.,  49./255.,  49./255.,  38./255., 32./255., 34./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Sienna
      case 99:
         {
            Double_t red[9]   = { 99./255., 112./255., 148./255., 165./255., 179./255., 182./255., 183./255., 183./255., 208./255.};
            Double_t green[9] = { 39./255.,  40./255.,  57./255.,  79./255., 104./255., 127./255., 148./255., 161./255., 198./255.};
            Double_t blue[9]  = { 15./255.,  16./255.,  18./255.,  33./255.,  51./255.,  79./255., 103./255., 129./255., 177./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Solar
      case 100:
         {
            Double_t red[9]   = { 99./255., 116./255., 154./255., 174./255., 200./255., 196./255., 201./255., 201./255., 230./255.};
            Double_t green[9] = {  0./255.,   0./255.,   8./255.,  32./255.,  58./255.,  83./255., 119./255., 136./255., 173./255.};
            Double_t blue[9]  = {  5./255.,   6./255.,   7./255.,   9./255.,   9./255.,  14./255.,  17./255.,  19./255.,  24./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // South West
      case 101:
         {
            Double_t red[9]   = { 82./255., 106./255., 126./255., 141./255., 155./255., 163./255., 142./255., 107./255.,  66./255.};
            Double_t green[9] = { 62./255.,  44./255.,  69./255., 107./255., 135./255., 152./255., 149./255., 132./255., 119./255.};
            Double_t blue[9]  = { 39./255.,  25./255.,  31./255.,  60./255.,  73./255.,  68./255.,  49./255.,  72./255., 188./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Starry Night
      case 102:
         {
            Double_t red[9]   = { 18./255., 29./255., 44./255.,  72./255., 116./255., 158./255., 184./255., 208./255., 221./255.};
            Double_t green[9] = { 27./255., 46./255., 71./255., 105./255., 146./255., 177./255., 189./255., 190./255., 183./255.};
            Double_t blue[9]  = { 39./255., 55./255., 80./255., 108./255., 130./255., 133./255., 124./255., 100./255.,  76./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Sunset
      case 103:
         {
            Double_t red[9]   = { 0./255., 48./255., 119./255., 173./255., 212./255., 224./255., 228./255., 228./255., 245./255.};
            Double_t green[9] = { 0./255., 13./255.,  30./255.,  47./255.,  79./255., 127./255., 167./255., 205./255., 245./255.};
            Double_t blue[9]  = { 0./255., 68./255.,  75./255.,  43./255.,  16./255.,  22./255.,  55./255., 128./255., 245./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Temperature Map
      case 104:
         {
            Double_t red[9]   = {  34./255.,  70./255., 129./255., 187./255., 225./255., 226./255., 216./255., 193./255., 179./255.};
            Double_t green[9] = {  48./255.,  91./255., 147./255., 194./255., 226./255., 229./255., 196./255., 110./255.,  12./255.};
            Double_t blue[9]  = { 234./255., 212./255., 216./255., 224./255., 206./255., 110./255.,  53./255.,  40./255.,  29./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Thermometer
      case 105:
         {
            Double_t red[9]   = {  30./255.,  55./255., 103./255., 147./255., 174./255., 203./255., 188./255., 151./255., 105./255.};
            Double_t green[9] = {   0./255.,  65./255., 138./255., 182./255., 187./255., 175./255., 121./255.,  53./255.,   9./255.};
            Double_t blue[9]  = { 191./255., 202./255., 212./255., 208./255., 171./255., 140./255.,  97./255.,  57./255.,  30./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Valentine
      case 106:
         {
            Double_t red[9]   = { 112./255., 97./255., 113./255., 125./255., 138./255., 159./255., 178./255., 188./255., 225./255.};
            Double_t green[9] = {  16./255., 17./255.,  24./255.,  37./255.,  56./255.,  81./255., 110./255., 136./255., 189./255.};
            Double_t blue[9]  = {  38./255., 35./255.,  46./255.,  59./255.,  78./255., 103./255., 130./255., 152./255., 201./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Visible Spectrum
      case 107:
         {
            Double_t red[9]   = { 18./255.,  72./255.,   5./255.,  23./255.,  29./255., 201./255., 200./255., 98./255., 29./255.};
            Double_t green[9] = {  0./255.,   0./255.,  43./255., 167./255., 211./255., 117./255.,   0./255.,  0./255.,  0./255.};
            Double_t blue[9]  = { 51./255., 203./255., 177./255.,  26./255.,  10./255.,   9./255.,   8./255.,  3./255.,  0./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Water Melon
      case 108:
         {
            Double_t red[9]   = { 19./255., 42./255., 64./255.,  88./255., 118./255., 147./255., 175./255., 187./255., 205./255.};
            Double_t green[9] = { 19./255., 55./255., 89./255., 125./255., 154./255., 169./255., 161./255., 129./255.,  70./255.};
            Double_t blue[9]  = { 19./255., 32./255., 47./255.,  70./255., 100./255., 128./255., 145./255., 130./255.,  75./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Cool
      case 109:
         {
            Double_t red[9]   = {  33./255.,  31./255.,  42./255.,  68./255.,  86./255., 111./255., 141./255., 172./255., 227./255.};
            Double_t green[9] = { 255./255., 175./255., 145./255., 106./255.,  88./255.,  55./255.,  15./255.,   0./255.,   0./255.};
            Double_t blue[9]  = { 255./255., 205./255., 202./255., 203./255., 208./255., 205./255., 203./255., 206./255., 231./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Copper
      case 110:
         {
            Double_t red[9]   = { 0./255., 25./255., 50./255., 79./255., 110./255., 145./255., 181./255., 201./255., 254./255.};
            Double_t green[9] = { 0./255., 16./255., 30./255., 46./255.,  63./255.,  82./255., 101./255., 124./255., 179./255.};
            Double_t blue[9]  = { 0./255., 12./255., 21./255., 29./255.,  39./255.,  49./255.,  61./255.,  74./255., 103./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Gist Earth
      case 111:
         {
            Double_t red[9]   = { 0./255., 13./255.,  30./255.,  44./255.,  72./255., 120./255., 156./255., 200./255., 247./255.};
            Double_t green[9] = { 0./255., 36./255.,  84./255., 117./255., 141./255., 153./255., 151./255., 158./255., 247./255.};
            Double_t blue[9]  = { 0./255., 94./255., 100./255.,  82./255.,  56./255.,  66./255.,  76./255., 131./255., 247./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Viridis
      case 112:
         {
            Double_t red[9]   = { 26./255., 51./255.,  43./255.,  33./255.,  28./255.,  35./255.,  74./255., 144./255., 246./255.};
            Double_t green[9] = {  9./255., 24./255.,  55./255.,  87./255., 118./255., 150./255., 180./255., 200./255., 222./255.};
            Double_t blue[9]  = { 30./255., 96./255., 112./255., 114./255., 112./255., 101./255.,  72./255.,  35./255.,   0./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;

      // Invert Sunset
      case 113:
         {
            // Double_t red[9]   = { 0./255., 48./255., 119./255., 173./255., 212./255., 224./255., 228./255., 228./255., 245./255.};
            // Double_t green[9] = { 0./255., 13./255.,  30./255.,  47./255.,  79./255., 127./255., 167./255., 205./255., 245./255.};
            // Double_t blue[9]  = { 0./255., 68./255.,  75./255.,  43./255.,  16./255.,  22./255.,  55./255., 128./255., 245./255.};

            Double_t red[9]   = { 0./255., 48./255., 119./255., 173./255., 212./255., 224./255., 228./255., 228./255., 245./255.};
            Double_t green[9] = { 0./255., 13./255.,  30./255.,  47./255.,  79./255., 127./255., 167./255., 205./255., 245./255.};
            Double_t blue[9]  = { 0./255., 68./255.,  75./255.,  43./255.,  16./255.,  22./255.,  55./255., 128./255., 245./255.};

	    std::reverse(red,red+9);
	    std::reverse(green,green+9);
	    std::reverse(blue,blue+9);
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
         }
         break;


      default:
         std::cerr << "Unknown palette number " << ncolors << std::endl;
         return;
      }
      paletteType = ncolors;
      
      if (Idx>0) fgPalettesList.fArray[paletteType-51] = (Double_t)Idx;
      else       fgPalettesList.fArray[paletteType-51] = 0.;

      return;
   }
 
   // // set user defined palette
   // if (colors)  {
   //    fgPalette.Set(ncolors);
   //    for (i=0;i<ncolors;i++) fgPalette.fArray[i] = colors[i];
   // } else {
   //    fgPalette.Set(TMath::Min(50,ncolors));
   //    for (i=0;i<TMath::Min(50,ncolors);i++) fgPalette.fArray[i] = palette[i];
   // }

   paletteType = 3;
}


//=========================================================
//  Apply the MnvPlotter axis style to all histograms in a TDirectory
//=========================================================
bool MnvPlotter::ApplyAxisStyle(
        TDirectory *td,
        bool centerXTitle /*= true*/,
        bool centerYTitle /*= false*/,
        bool centerZTitle /*= false*/
        ) const
{
    //! Complain and return false if no TDirectory was given
    if ( !td )
    {
        Error( "ApplyAxisStyle", "You gave me a NULL TDirectory." );
        return false;
    }

    //! Loop over all keys in the TDirectory
    TIter nextkey( td->GetListOfKeys() );
    TKey *key = 0;
    while ( (key = (TKey*)nextkey()))
    {
        TObject *obj = key->ReadObj();

        //! Skip the object if it doesn't inherit from a TH1
        TH1 *h = dynamic_cast<TH1*>( obj );
        if ( ! h )
            continue;

        ApplyAxisStyle( h, centerXTitle, centerYTitle, centerZTitle );
    }
    return true;
}

bool MnvPlotter::ApplyAxisStyle(
        TH1 *h,
        bool centerXTitle /*= true*/,
        bool centerYTitle /*= false*/,
        bool centerZTitle /*= false*/
        ) const
{

    //!Set the X axis
    h->GetXaxis()->SetNdivisions(509);
    h->GetXaxis()->CenterTitle( centerXTitle );
    h->GetXaxis()->SetTitleOffset( axis_title_offset_x );
    h->GetXaxis()->SetTitleSize( axis_title_size_x );
    h->GetXaxis()->SetTitleFont( axis_title_font_x );

    h->GetXaxis()->SetLabelFont( axis_label_font );
    h->GetXaxis()->SetLabelSize( axis_label_size );

    //!Set the Y axis
    h->GetYaxis()->CenterTitle( centerYTitle );
    h->GetYaxis()->SetTitleOffset( axis_title_offset_y );
    h->GetYaxis()->SetTitleSize( axis_title_size_y );
    h->GetYaxis()->SetTitleFont( axis_title_font_y );

    h->GetYaxis()->SetLabelFont( axis_label_font );
    h->GetYaxis()->SetLabelSize( axis_label_size );

    //! Set the Z axis
    if ( h->GetZaxis() != NULL )
    {
        h->GetZaxis()->CenterTitle( centerZTitle );
        h->GetZaxis()->SetTitleOffset( axis_title_offset_z );
        h->GetZaxis()->SetTitleSize( axis_title_size_z );
        h->GetZaxis()->SetTitleFont( axis_title_font_z );

        h->GetZaxis()->SetLabelFont( axis_label_font );
        h->GetZaxis()->SetLabelSize( axis_label_size );
    }

    return true;
}

bool MnvPlotter::ApplyAxisStyle(
        bool centerXTitle /*= true*/,
        bool centerYTitle /*= false*/,
        bool centerZTitle /*= false*/
        ) const
{
    return ApplyAxisStyle( gDirectory, centerXTitle, centerYTitle, centerZTitle );
}


void MnvPlotter::UseAutoAxisLimits()
{
    axis_minimum = axis_maximum = MnvHist::AutoAxisLimit;
}

//=========================================================


//=========================================================
// easily add Latex labels on plots
//=========================================================

void MnvPlotter::AddPlotLabel(
        const char* label,
        const double x,
        const double y,
        const double size /*= 0.05*/,
        const int color /*= 1*/,
        const int font /*= 62*/,
        const int align /*= 22*/,
        const double angle /*= 0*/
        )
{
    TLatex *latex = new TLatex( x, y, label );
    AddToTmp( latex );

    latex->SetNDC();
    latex->SetTextSize(size);
    latex->SetTextColor(color);
    latex->SetTextFont(font);
    latex->SetTextAlign(align);
    latex->SetTextAngle(angle);
    latex->Draw();
}

void MnvPlotter::AddHistoTitle(
        const char* title,
        double titleSize,
        int titleFont
        )
{
    AddPlotLabel(title, 0.5, 1-gStyle->GetPadTopMargin()-extra_top_margin, titleSize, 1, titleFont, 21);
}




size_t MnvPlotter::GetLegendEntrySize(
        const std::string& title
        )  const
{

    unsigned int raw_size = title.size();

    // Disabled this

    if (0) {
        size_t new_pos = title.find('#');
        while ( new_pos != std::string::npos && raw_size>1 ) {
            raw_size -= 2; // do not count # variable and the following character towards size (Latex)
            new_pos = title.find('#',new_pos);
        }
    }

    //Remove characters: {  }  ^  _   from size count
    //
    for (unsigned int i=0; i!= title.size(); ++i) {
        char c = title[i];
        if (raw_size>0) {
            if ( c=='}' || c=='{') {
                --raw_size;
            }
            if (i < title.size()-1 && (c=='^'||c=='_') && (char)title[i+1]=='{') {
                --raw_size;
            }
        }
    }

    return (size_t)raw_size;
}




void MnvPlotter::DecodePosition(
        const std::string& opts,
        double size,
        int &align,
        double &xLabel,
        double &yLabel
        ) const
{
    //! if opts is a mixture of two strings like TR-TC, then use the average of those two positions and align of the first
    size_t dashLoc = opts.find("-");
    if ( dashLoc != string::npos )
    {
        const string opts1 = opts.substr(0, dashLoc);
        int align1;
        double x1, y1;
        DecodePosition( opts1, size, align1, x1, y1 );

        const string opts2 = opts.substr(dashLoc+1);
        int align2;
        double x2, y2;
        DecodePosition( opts2, size, align2, x2, y2 );

        align = align1;
        xLabel = ( x1 + x2 ) / 2.;
        yLabel = ( y1 + y2 ) / 2.;
        return;
    }

    const double xLeft  = gStyle->GetPadLeftMargin() + 0.03;
    const double xCenter = .5;
    const double xRight = 1 - gStyle->GetPadRightMargin() - 0.025;

    const double yBottom = gStyle->GetPadBottomMargin() + size/2.;
    const double yCenter = .5;
    const double yTop    = 1 - gStyle->GetPadTopMargin() - size/2.;

    if ( opts == "TC" || opts == "") {
        // default is TC (top center)
        align = 23;
        xLabel = xCenter;
        yLabel = yTop;
    }else if ( opts == "TR" ) {
        align = 33;
        xLabel = xRight;
        yLabel = yTop;
    }else if ( opts == "TL" ) {
        align = 13;
        xLabel = xLeft;
        yLabel = yTop;
    }else if ( opts == "BR" ) {
        align = 31;
        xLabel = xRight;
        yLabel = yBottom;
    }else if ( opts == "BL" ) {
        align = 11;
        xLabel = xLeft;
        yLabel = yBottom;
    }else if ( opts == "BC" ) {
        align = 21;
        xLabel = xCenter;
        yLabel = yBottom;
    }else if ( opts == "L" ) {
        align = 12;
        xLabel = xLeft;
        yLabel = yCenter;
    }else if ( opts == "C" ) {
        align = 22;
        xLabel = xCenter;
        yLabel = yCenter;
    }else if ( opts == "R" ) {
        align = 32;
        xLabel = xRight;
        yLabel = yCenter;
    }else{
        Warning("DecodePosition", Form("Position option '%s' is not valid,  No values have been set.", opts.c_str() ) );
    }

}

size_t MnvPlotter::GetLegendWidthInLetters( const std::vector<std::string>& titles ) const
{
    //find the longest title in each column
    size_t longestTitle[legend_n_columns];
    for ( size_t i = 0; i != legend_n_columns; ++i )
        longestTitle[i] = 0;

    for ( unsigned int i = 0; i != titles.size(); ++i )
    {
        unsigned int col = i % legend_n_columns;
        longestTitle[col] = std::max( GetLegendEntrySize( titles[i] ), longestTitle[col] );
    }
    size_t sumLongestTitles = 0;
    for ( size_t i = 0; i != legend_n_columns; ++i )
        sumLongestTitles += longestTitle[i];

    return sumLongestTitles;
}



void MnvPlotter::DecodeLegendPosition(
        double &x1,
        double &y1,
        double &x2,
        double &y2,
        const std::string& opts,
        const size_t nHists,
        const size_t longestTitleSize /*=14*/,
        const double textSize /* = legend_text_size */
        ) const
{
    //! start by using DecodePosition
    int align = 22;
    DecodePosition( opts, textSize, align, x1, y1 );

    //! how tall does the legend need to be for this number of histograms?
    const double yspace = height_nspaces_per_hist * textSize * nHists / float( legend_n_columns );

    //! how wide does the legend need to be based on the longest title size?
    const double xspace_marker = .06;  //~space taken up by marker
    const double xspace = ((float)legend_n_columns)*xspace_marker + longestTitleSize*width_xspace_per_letter*textSize; //x extent necessary

    //! set vertical position
    if ( align % 10 == 2 ) //center should be at current y_legend
        y1 -= yspace / 2.;
    else if ( align % 10 == 3 ) //right edge should be at current y_legend
        y1 -= yspace;

    //! set horizontal position
    if ( align / 10 == 2 ) //center should be at current x_legend
        x1 -= xspace / 2.;
    else if ( align / 10 == 3 ) //top edge should be at current x_legend
        x1 -= xspace;

    x1 += legend_offset_x;
    y1 += legend_offset_y;

    x2 = x1 + xspace;
    y2 = y1 + yspace;

}


void MnvPlotter::WritePreliminary(
        const std::string& opts,
        double size,
        double yOffset,
        double xOffset,
        bool isWorkInProgress /*=false*/
        )
{
    int align;
    double xLabel, yLabel;
    DecodePosition( opts, size, align, xLabel, yLabel );
    yLabel += yOffset;
    xLabel += xOffset;
    if (isWorkInProgress) AddPlotLabel( "MINER#nuA Work In Progress", xLabel, yLabel, size, kRed+1, 112, align );
    else AddPlotLabel( "MINER#nuA Preliminary", xLabel, yLabel, size, 4, 112, align );
}

void MnvPlotter::WritePreliminary(
        double x,
        double y,
        double size,
        bool isWorkInProgress /*=false*/
        )
{
    if (isWorkInProgress) AddPlotLabel( "MINER#nuA Work In Progress", x, y, size, kRed+1, 112, 22 );
    else AddPlotLabel( "MINER#nuA Preliminary", x, y, size, 4, 112, 22 );
}

void MnvPlotter::WriteNorm(
        const char *norm,
        const std::string& opts,
        double size,
        double yOffset,
        double xOffset,
        double pot )
{
    int align;
    double xLabel, yLabel;
    DecodePosition( opts, size, align, xLabel, yLabel );
    yLabel += yOffset;
    xLabel += xOffset;
    AddPlotLabel( norm, xLabel, yLabel, size, 9, 42, align );
    if ( pot != 0 )
        AddPlotLabel( Form("%4.2e POT", pot), xLabel, yLabel-size, 0.8*size, 1, 42, align );
}

void MnvPlotter::WriteNorm(
        const char *norm,
        double x,
        double y,
        double size,
        double pot
        )
{
    AddPlotLabel( norm, x, y, size, 9, 42, 22 );
    if ( pot != 0 )
        AddPlotLabel( Form("%4.2e POT", pot), x, y-size, 0.8*size, 1, 42, 22 );
}

void MnvPlotter::AddPOTNormBox(
        const double dataPOT,
        const double mcPOT,
        const double xLeft,
        const double yTop,
        const double size /*= .03*/
        )
{
    WriteNorm( "POT-Normalized",                 xLeft, yTop, size );
    WriteNorm( Form( "Data POT: %.2E", dataPOT), xLeft, yTop - size, size );
    WriteNorm( Form( "MC POT: %.2E", mcPOT),     xLeft, yTop - 2*size, size );
}

void MnvPlotter::AddAreaNormBox(
        const double dataScale,
        const double mcScale,
        const double xLeft,
        const double yTop,
        double size /*= 0.03*/
        )
{
    WriteNorm( "Area-Normalized",                               xLeft, yTop, size );
    WriteNorm( Form( "MC Scale = %.2E ", dataScale / mcScale ), xLeft, yTop - size, size );
    WriteNorm( Form( "  = %.2E/%.2E", dataScale, mcScale ),     xLeft, yTop - 2.*size, size );
}

void MnvPlotter::AddChi2Label(
        const TH1* dataHist,
        const TH1* mcHist,
        const std::string& opts
        )
{
    AddChi2Label( dataHist, mcHist, 1.0, opts );
}

void MnvPlotter::AddChi2Label(
        const MnvH1D* dataHist,
        const MnvH1D* mcHist,
        const std::string& opts,
        const bool useDataErrorMatrix,
        const bool useOnlyShapeErrors
        )
{
    AddChi2Label( dataHist, mcHist, 1.0, opts, 0.04, 0.0, useDataErrorMatrix, useOnlyShapeErrors );
}

void MnvPlotter::AddChi2Label(
        const TH1* dataHist,
        const TH1* mcHist,
        Double_t mcScale,
        const std::string& opts,
        double size,
        double yOffset
        )
{
    int align;
    double xLabel, yLabel;
    DecodePosition( opts, size, align, xLabel, yLabel );
    yLabel += yOffset;

    Int_t ndf;
    Double_t chi2 = Chi2DataMC( dataHist, mcHist, ndf, mcScale );

    char *words = Form("#chi^{2}/ndf = %3.2f/%d = %3.2f", chi2, ndf, chi2/(Double_t)ndf);
    AddPlotLabel( words, xLabel, yLabel, size, 1, 62, align );

}

void MnvPlotter::AddChi2Label(
        const MnvH1D* dataHist,
        const MnvH1D* mcHist,
        Double_t mcScale,
        const std::string& opts,
        double size,
        double yOffset,
        const bool useDataErrorMatrix,
        const bool useOnlyShapeErrors,
        const std::string& pre_tag /*=""*/
        )
{
    int align;
    double xLabel, yLabel;
    DecodePosition( opts, size, align, xLabel, yLabel );
    yLabel += yOffset;

    Int_t ndf;
    Double_t chi2 = Chi2DataMC( dataHist, mcHist, ndf, mcScale, useDataErrorMatrix, useOnlyShapeErrors );

    char *words = Form("%s #chi^{2}/ndf = %3.2f/%d = %3.2f", pre_tag.c_str(), chi2, ndf, chi2/(Double_t)ndf);
    AddPlotLabel( words, xLabel, yLabel, size, 1, 62, align );
}

void MnvPlotter::AddCutArrow(
        const double cut_location,
        const double y1,
        const double y2,
        const double arrow_length,
        const std::string& arrow_direction
        ) const
{

    double arrow_tip = arrow_length;
    if (arrow_direction == "L") {
        arrow_tip *= -1.0;
    }
    else if (arrow_direction != "R") {
        std::cout<<"Do not understand supplied arrow direction.  L (left) and R (right) only!"<<std::endl;
    }

    TLine line;
    line.SetLineWidth(arrow_line_width);
    line.SetLineStyle(arrow_line_style);
    line.SetLineColor(arrow_line_color);
    line.DrawLine(cut_location,y1,cut_location,y2);

    TArrow arrow;
    arrow.SetLineWidth(arrow_line_width);
    arrow.SetLineStyle(arrow_line_style);
    arrow.SetLineColor(arrow_line_color);
    arrow.DrawArrow(cut_location,y2,cut_location+arrow_tip,y2,arrow_size,arrow_type.c_str());
}

//==========================================================
// Get a distinct line style
//==========================================================
void MnvPlotter::ApplyNextLineStyle(
        TH1* h,
        bool startOver /* = false */,
        bool changeStyle /* = true */
        )const
{
    static unsigned int nseen = 0;
    if ( startOver)
        nseen = 0;

    if ( ! h )
        return;

    //! Pick the next good color and apply it to the line
    h->SetLineColor( good_colors[ nseen % good_colors.size() ] );

    //! If changing the line style, apply a new style
    if ( changeStyle )
        h->SetLineStyle( nseen % 10 + 1 );

    h->SetLineWidth( mc_line_width );
    ++nseen;
}
//==========================================================
// draw MC data points on current Pad - MnvH1D overload.
//==========================================================

void MnvPlotter::DrawDataMC(
        const TH1* dataHist,
        const TH1* mcHist,
        const Double_t mcScale, /*= 1.0*/
        const std::string& legPos, /*= "L"*/
        const bool useHistTitles /*=false*/
        )
{
    if (!gPad)
      throw std::runtime_error("MnvPlotter requires a TCanvas. Please make one first.");

    TH1D *tmpData = (TH1D*)dataHist->Clone("tmpData");
    TH1D *tmpMC = (TH1D*)mcHist->Clone("tmpMC");

    tmpMC->Scale(mcScale);

    //respect max/min setting the user may have used
    if ( MnvHist::IsAutoAxisLimit( axis_minimum ) )
    {
        if ( hist_min_zero && !gPad->GetLogy() )
            tmpMC->SetMinimum( 0. );
        else
            tmpMC->SetMinimum( footroom * std::min( tmpData->GetMinimum(), tmpMC->GetMinimum() ) );
    }
    else
        tmpMC->SetMinimum( axis_minimum );

    if ( MnvHist::IsAutoAxisLimit( axis_maximum ) )
        tmpMC->SetMaximum( headroom * std::max( tmpData->GetMaximum(), tmpMC->GetMaximum() ) );
    else
        tmpMC->SetMaximum( axis_maximum );

    tmpData->SetMarkerStyle(data_marker);
    tmpData->SetMarkerSize(data_marker_size);
    tmpData->SetMarkerColor(data_color);
    tmpData->SetLineWidth(data_line_width);
    tmpData->SetLineStyle(data_line_style);
    tmpData->SetLineColor(data_color);

    tmpMC->GetXaxis()->SetTitleFont(axis_title_font_x);
    tmpMC->GetYaxis()->SetTitleFont(axis_title_font_y);
    tmpMC->GetXaxis()->SetTitleSize(axis_title_size_x);
    tmpMC->GetYaxis()->SetTitleSize(axis_title_size_y);
    tmpMC->GetXaxis()->SetLabelFont(axis_label_font);
    tmpMC->GetYaxis()->SetLabelFont(axis_label_font);
    tmpMC->GetXaxis()->SetLabelSize(axis_label_size);
    tmpMC->GetYaxis()->SetLabelSize(axis_label_size);
    tmpMC->GetXaxis()->CenterTitle(kTRUE);

    tmpMC->SetLineColor(mc_color);
    tmpMC->SetLineWidth(mc_line_width);
    tmpMC->SetLineStyle(mc_line_style);

    tmpMC->Draw("HIST");
    tmpData->DrawCopy("E1X0SAME");


    if ( legPos != "N" )
    {
        const string data_name    = (useHistTitles ? tmpData->GetTitle() : "Data" );
        const string mc_name      = (useHistTitles ? tmpMC->GetTitle()     : "Simulation" );
        double x1, y1, x2, y2;
        //figure out where to put the legend
        vector<string> titles;
        titles.push_back( data_name );
        titles.push_back( mc_name );
        size_t legendWidth = GetLegendWidthInLetters( titles );
        DecodeLegendPosition( x1, y1, x2, y2, legPos, titles.size(), legendWidth );

        TLegend *leg = new TLegend(x1, y1, x2, y2);
        leg->SetNColumns( legend_n_columns );
        static int legN = 0;
        leg->SetName( Form("This Legend %d", legN++ ) );
        AddToTmp( leg );

        leg->SetBorderSize( legend_border_size );
        leg->SetFillColor( legend_fill_color );
        if ( legend_fill_color < 0 )
            leg->SetFillStyle(0);
        leg->SetTextSize( legend_text_size);
        leg->SetTextFont( legend_text_font );
        leg->AddEntry(tmpData, data_name.c_str(), "lep");
        leg->AddEntry(tmpMC, mc_name.c_str(), "l");

        leg->Draw();
    }

    gPad->RedrawAxis();
    gPad->Update();

}


//==========================================================
// draw MC with error band + data points on current Pad - MnvH1D overload.
//==========================================================

void MnvPlotter::DrawDataMCWithErrorBand(
        const MnvH1D* dataHist,
        const MnvH1D* mcHist,
        const Double_t mcScale /*= 1.0*/,
        const std::string& legPos /*= "L"*/,
        const bool useHistTitles /*=false*/,
        const MnvH1D* bkgdHist /*= NULL*/,
        const MnvH1D* dataBkgdHist /*= NULL*/,
        const bool covAreaNormalize/*= false, Area Normalize considerations for covariance matrix*/,
        const bool statPlusSys /* = false */,
        const bool isSmooth /* = false */)
{

    TH1Ptr tmpData;
    TH1Ptr tmpDataStat;
    if (statPlusSys)
    {
        tmpData = TH1Ptr( (TH1*)dataHist->GetCVHistoWithError(true, covAreaNormalize).Clone( "mnv_tmp_data_statplussys" ) );
        tmpDataStat = TH1Ptr( (TH1*)dataHist->GetCVHistoWithStatError().Clone( "mnv_tmp_data_statonly" ) );
    }
    else
        tmpData = TH1Ptr( (TH1*)dataHist->GetCVHistoWithError(true, covAreaNormalize).Clone( "mnv_tmp_data" ) );

    TH1Ptr tmpMC   = TH1Ptr( (TH1*)mcHist  ->GetCVHistoWithError(true, covAreaNormalize).Clone( "mnv_tmp_mc"   ) );
    if (covAreaNormalize && tmpData->Integral() > 0)
        tmpMC->Scale(tmpData->Integral() / tmpMC->Integral());

    if ( draw_normalized_to_bin_width )
    {
        if ( dataHist->GetNormBinWidth() > 0 )
        {
            tmpData->Scale( dataHist->GetNormBinWidth(), "width" );
            if (tmpDataStat != NULL)
                tmpDataStat->Scale( dataHist->GetNormBinWidth(), "width" );
        }

        if ( mcHist->GetNormBinWidth() > 0 )
            tmpMC->Scale( mcHist  ->GetNormBinWidth(), "width" );
    }

    if ( bkgdHist != NULL )
    {
        TH1Ptr tmpBK = TH1Ptr( (TH1*)bkgdHist->GetCVHistoWithError(true, covAreaNormalize).Clone( Form("mnv_tmp_bkgd_%d", __LINE__ ) ) );
        if ( draw_normalized_to_bin_width && bkgdHist->GetNormBinWidth() > 0  )
            tmpBK->Scale( bkgdHist->GetNormBinWidth(), "width" );

        if ( dataBkgdHist != NULL )
        {
            TH1Ptr tmpDataBK = TH1Ptr( (TH1*)dataBkgdHist->GetCVHistoWithError(true, covAreaNormalize).Clone( "mnv_tmp_data_bkgd" ) );
            if ( draw_normalized_to_bin_width  && dataBkgdHist->GetNormBinWidth() > 0 )
                tmpDataBK->Scale( dataBkgdHist->GetNormBinWidth(), "width" );

            DrawDataMCWithErrorBand( tmpData, tmpMC, mcScale, legPos, useHistTitles, tmpBK, tmpDataBK );
        } // if (dataBkgdHist != NULL)
        else
            DrawDataMCWithErrorBand( tmpData, tmpMC, mcScale, legPos, useHistTitles, tmpBK );
    } // if ( bkgdHist != NULL )
    else
        DrawDataMCWithErrorBand( tmpData, tmpMC, mcScale, legPos, useHistTitles, NULL, NULL, isSmooth);

    if (statPlusSys)
    {
        tmpDataStat->SetMarkerStyle(data_marker);
        tmpDataStat->SetMarkerSize(data_marker_size);
        tmpDataStat->SetMarkerColor(data_color);
        tmpDataStat->SetLineWidth(data_line_width);
        tmpDataStat->SetLineStyle(data_line_style);
        tmpDataStat->SetLineColor(data_color);
        tmpDataStat->DrawCopy("SAME E1 X0");
    }

}

//==========================================================
// draw MC with error band + data points on current Pad
//==========================================================

void MnvPlotter::DrawDataMCWithErrorBand(
        const TH1* dataHist,
        const TH1 *mcHist,
        const Double_t mcScale     /*= 1.0*/,
        const std::string& legPos /*= "L"*/,
        const bool useHistTitles    /*=false */,
        const TH1* bkgdHist        /*= NULL*/,
        const TH1* dataBkgdHist   /*= NULL*/,
        const bool isSmooth /* =false */)
{
    if (!gPad)
      throw std::runtime_error("MnvPlotter requires a TCanvas. Please make one first.");

    TH1 * tmpData = dynamic_cast<TH1*>(dataHist->Clone( "tmp_data" ));
    TH1 * tmpMC = dynamic_cast<TH1*>(mcHist->Clone( "tmpMC" ));

    TH1 * tmpBkgd = NULL;
    if ( bkgdHist )
    {
        tmpBkgd = dynamic_cast<TH1*>(bkgdHist->Clone( "bktmp" ) );
        tmpBkgd->SetFillColor(mc_bkgd_color);
        tmpBkgd->SetFillStyle(mc_bkgd_style);
        tmpBkgd->SetLineColor(mc_bkgd_line_color);
        tmpBkgd->SetLineWidth(mc_bkgd_width);
    }

    TH1 * tmpDataBkgd = NULL;
    if ( dataBkgdHist )
    {
        tmpDataBkgd = dynamic_cast<TH1*>(dataBkgdHist->Clone( "databktmp" ) );
        tmpDataBkgd->SetMarkerColor(data_bkgd_color);
        tmpDataBkgd->SetMarkerStyle(data_bkgd_style);
        tmpDataBkgd->SetMarkerSize(data_bkgd_size);
    }


    if ( tmpData->GetSumw2N() == 0 )
        tmpData->Sumw2();

    if ( tmpMC->GetSumw2N() == 0 )
        tmpMC->Sumw2();

    //respect max/min setting the user may have used
    if ( MnvHist::IsAutoAxisLimit( axis_minimum ) )
    {
        if ( hist_min_zero && !gPad->GetLogy() )
            tmpMC->SetMinimum( 0. );
        else
            tmpMC->SetMinimum( footroom * std::min( tmpData->GetMinimum(), mcScale*tmpMC->GetMinimum() ) );
    }
    else
        tmpMC->SetMinimum( axis_minimum );

    if ( MnvHist::IsAutoAxisLimit( axis_maximum ) )
        tmpMC->SetMaximum( headroom * std::max( tmpData->GetMaximum(), mcScale*tmpMC->GetMaximum() ) );
    else
        tmpMC->SetMaximum( axis_maximum );

    {
        //temporarily set the min/max so DrawDataMCWithErrorBand doesn't change anything
        const double oldMin = axis_minimum;
        const double oldMax = axis_maximum;
        axis_minimum = tmpMC->GetMinimum();
        axis_maximum = tmpMC->GetMaximum();

        if (bkgdHist)
            DrawMCWithErrorBand( tmpMC, mcScale, tmpBkgd, isSmooth );
        else
            DrawMCWithErrorBand( tmpMC, mcScale, NULL, isSmooth );

        //put back the axis min/max
        axis_minimum = oldMin;
        axis_maximum = oldMax;
    }

    tmpData->SetMarkerStyle(data_marker);
    tmpData->SetMarkerSize(data_marker_size);
    tmpData->SetMarkerColor(data_color);
    tmpData->SetLineWidth(data_line_width);
    tmpData->SetLineStyle(data_line_style);
    tmpData->SetLineColor(data_color);
    tmpData->DrawCopy("SAME E1 X0");
    if ( dataBkgdHist )
        tmpDataBkgd->DrawCopy("SAME HIST P");

    if ( legPos != "N" )
    {
        const string data_name    = (useHistTitles ? tmpData->GetTitle() : "Data" );
        const string mc_name      = (useHistTitles ? tmpMC->GetTitle()     : "Simulation" );
        const string data_bg_name = (bkgdHist && useHistTitles ? tmpDataBkgd->GetTitle() : "Data Background" );
        const string bg_name      = (bkgdHist     && useHistTitles ? tmpBkgd->GetTitle()     : "Sim. Background" );

        //figure out where to put the legend
        double x1, y1, x2, y2;
        vector<string> titles;
        titles.push_back( data_name );
        titles.push_back( mc_name );
        if ( dataBkgdHist )
            titles.push_back( data_bg_name );
        if ( bkgdHist )
            titles.push_back( bg_name );
        size_t legendWidth = GetLegendWidthInLetters( titles );
        DecodeLegendPosition( x1, y1, x2, y2, legPos, titles.size(), legendWidth );

        tmpMC->SetLineColor(mc_color);
        tmpMC->SetLineWidth(mc_line_width);
        tmpMC->SetLineStyle(mc_line_style);
        tmpMC->SetFillColor(mc_error_color);
        tmpMC->SetFillStyle(mc_error_style);

        TLegend *leg = new TLegend(x1, y1, x2, y2);
        leg->SetNColumns( legend_n_columns );
        static int legN = 0;
        leg->SetName( Form("This Legend %d", legN++ ) );
        AddToTmp( leg );

        leg->SetBorderSize( legend_border_size );
        leg->SetFillColor( legend_fill_color );
        if ( legend_fill_color < 0 )
            leg->SetFillStyle(0);
        leg->SetTextSize( legend_text_size);
        leg->SetTextFont( legend_text_font );
        leg->AddEntry(tmpData, data_name.c_str(), "lep");

        // If the MC error is smaller than 1% in every bin, you can't see it so take it off the legend?
        bool has_error = false;
        for ( int b = 1; b <= tmpMC->GetNbinsX(); ++b ) {
            if ( tmpMC->GetBinContent(b) > 0.0 && tmpMC->GetBinError(b) / tmpMC->GetBinContent(b) > 0.01 ) has_error = true;
        }
        if ( has_error ) leg->AddEntry(tmpMC, mc_name.c_str(), "fl");
        else            leg->AddEntry(tmpMC, mc_name.c_str(), "l");

        if ( dataBkgdHist )
            leg->AddEntry(tmpDataBkgd, data_bg_name.c_str(), "p" );
        if ( bkgdHist )
            leg->AddEntry(tmpBkgd, bg_name.c_str(), "f" );

        leg->Draw();
    }

    gPad->RedrawAxis();
    gPad->Update();
}


//=========================================================
// draw data/MC ratio in current pad
//=========================================================

void MnvPlotter::DrawMCWithErrorBand(
        const TH1 *mcHist,
        const Double_t mcScale /*= 1.0*/,
        const TH1* bkgdHist /*= NULL*/,
        const bool isSmooth /* = false */)
{
    if (!gPad)
      throw std::runtime_error("MnvPlotter requires a TCanvas. Please make one first.");

    TH1Ptr tmpMC( (TH1*)mcHist->Clone( "tmpMC" ) );
    tmpMC->Scale(mcScale);

    //respect max/min setting the user may have used
    if ( MnvHist::IsAutoAxisLimit( axis_minimum ) )
    {
        if ( hist_min_zero && !gPad->GetLogy() )
            tmpMC->SetMinimum( 0. );
        else
            tmpMC->SetMinimum( footroom * mcScale*tmpMC->GetMinimum() );
    }
    else
        tmpMC->SetMinimum( axis_minimum );

    if ( MnvHist::IsAutoAxisLimit( axis_maximum ) ) {
      //tmpMC->SetMaximum( headroom * mcScale*tmpMC->GetMaximum() );
      //better maximum hist setting. Xianguo
      const Double_t hmax = tmpMC->GetBinContent(tmpMC->GetMaximumBin())*1.5;
      tmpMC->SetMaximum(hmax);
    }
    else
        tmpMC->SetMaximum( axis_maximum );

    tmpMC->GetXaxis()->SetTitleFont(axis_title_font_x);
    tmpMC->GetYaxis()->SetTitleFont(axis_title_font_y);
    tmpMC->GetXaxis()->SetTitleSize(axis_title_size_x);
    tmpMC->GetYaxis()->SetTitleSize(axis_title_size_y);
    tmpMC->GetXaxis()->SetTitleOffset(axis_title_offset_x);
    tmpMC->GetYaxis()->SetTitleOffset(axis_title_offset_y);
    tmpMC->GetXaxis()->SetLabelFont(axis_label_font);
    tmpMC->GetYaxis()->SetLabelFont(axis_label_font);
    tmpMC->GetXaxis()->SetLabelSize(axis_label_size);
    tmpMC->GetYaxis()->SetLabelSize(axis_label_size);
    tmpMC->GetXaxis()->CenterTitle(kTRUE);

    //draw the error part of the MC histogram
    tmpMC->SetFillColor(mc_error_color);
    tmpMC->SetFillStyle(mc_error_style);
    tmpMC->SetMarkerStyle(0);
    // No errors for smooth curve
    if (!isSmooth) {
        tmpMC->DrawCopy("E2");
        tmpMC->Draw("SAME AXIS");
    }
    //draw the line part of the MC histogram
    tmpMC->SetFillColor(0);
    tmpMC->SetLineColor(mc_color);
    tmpMC->SetLineStyle(mc_line_style);
    tmpMC->SetLineWidth(mc_line_width);
    if (isSmooth) {
        tmpMC->Smooth(2);
        tmpMC->DrawCopy("HISTC");
        tmpMC->Draw("SAME AXIS");
    }else{
        tmpMC->DrawCopy("SAME HIST");
    }

    if ( bkgdHist != NULL )
    {
        TH1Ptr bkTmp( (TH1*)bkgdHist->Clone( "tmp_bkgd_hist") );

        bkTmp->Scale(mcScale);
        bkTmp->SetFillColor(mc_bkgd_color);
        bkTmp->SetFillStyle(mc_bkgd_style);
        bkTmp->SetLineColor(mc_bkgd_line_color);
        bkTmp->SetLineWidth(mc_bkgd_width);
        bkTmp->DrawCopy("SAME HIST");
    }

    gPad->Update();
}


bool MnvPlotter::AddSysError( TH1D *h, const TH1D *hErr, const bool sumSquares /*= true*/) const
{
    if ( h->GetNbinsX() != hErr->GetNbinsX() ) {
        Error("MnvPlotter::AddSysError", "The histogram of errors must have the same number of bins as the histgram you are adding errors to");
        return false;
    }

    for ( Int_t i = 0; i <= h->GetNbinsX(); i++ ) {
        double err = hErr->GetBinContent(i) + h->GetBinError(i);
        if ( sumSquares )
            err = sqrt( pow(err,2) + pow(h->GetBinError(i),2) );
        h->SetBinError( i, err );
    }

    return true;
}


void MnvPlotter::DrawErrorBand(
        const TH1 *h,
        const int color)
{
    if (!gPad)
      throw std::runtime_error("MnvPlotter requires a TCanvas. Please make one first.");

    TH1Ptr tmpErrUp(   (TH1D*)h->Clone("tmp_err_up"  ) );
    TH1Ptr tmpErrDown( (TH1D*)h->Clone("tmp_err_down") );

    for ( Int_t i = 1; i <= h->GetNbinsX(); i++ )
    {
        double err = h->GetBinError(i);
        tmpErrUp  ->SetBinContent(i,h->GetBinContent(i) + err );
        tmpErrDown->SetBinContent(i,h->GetBinContent(i) - err );
    }

    tmpErrUp->SetFillColor(color);
    tmpErrUp->SetFillStyle(3001);
    tmpErrUp->SetLineColor(color);
    tmpErrUp->SetLineWidth(1);

    tmpErrDown->SetFillColor(10);
    tmpErrDown->SetLineColor(color);
    tmpErrDown->SetLineWidth(1);

    tmpErrUp  ->DrawCopy("A HIST SAME");
    tmpErrDown->DrawCopy("SAME HIST");

    gPad->Update();
}


//==========================================================
// draw MC with error band + data points on current Pad
//==========================================================

void MnvPlotter::DrawDataMCRatio(
        const TH1* dataHist,
        const TH1 *mcHist,
        const Double_t mcScale,    /*= 1.0*/
        const bool     drawOneLine, /*= true*/
        const double   plotMin,     /*=-1 (auto)*/
        const double   plotMax,     /*=-1 (auto)*/
        const char*    yaxisLabel   /*="Data / MC"*/
        )
{
    if (!gPad)
      throw std::runtime_error("MnvPlotter requires a TCanvas. Please make one first.");

    TH1Ptr tmpData( (TH1*)dataHist->Clone("tmp_data") );
    TH1Ptr tmpMC(   (TH1*)mcHist  ->Clone("tmp_MC"  ) );

    if ( tmpData->GetSumw2N() == 0 )
        tmpData->Sumw2();
    if ( tmpMC->GetSumw2N() == 0 )
        tmpMC->Sumw2();

    tmpMC->Scale(mcScale);

    TH1Ptr tmpRatio( (TH1*)tmpMC->Clone( "tmp_ratio" ) );

    tmpRatio->Divide(tmpData, tmpMC);

    tmpRatio->GetYaxis()->SetTitle( yaxisLabel );
    tmpRatio->GetXaxis()->SetTitleFont(axis_title_font_x);
    tmpRatio->GetYaxis()->SetTitleFont(axis_title_font_y);
    tmpRatio->GetXaxis()->SetTitleSize(axis_title_size_x);
    tmpRatio->GetYaxis()->SetTitleSize(axis_title_size_y);
    tmpRatio->GetXaxis()->SetLabelFont(axis_label_font);
    tmpRatio->GetYaxis()->SetLabelFont(axis_label_font);
    tmpRatio->GetXaxis()->SetLabelSize(axis_label_size);
    tmpRatio->GetYaxis()->SetLabelSize(axis_label_size);
    tmpRatio->GetXaxis()->CenterTitle(kTRUE);
    tmpRatio->SetMarkerStyle(ratio_marker);
    tmpRatio->SetMarkerSize(ratio_marker_size);
    tmpRatio->SetLineWidth(ratio_line_width);
    tmpRatio->SetLineColor(ratio_color);
    tmpRatio->GetXaxis()->SetNdivisions(509);

    //!set the max bin based on the taller of the ratio or error
    double setMax = headroom*std::max(1.0, tmpRatio->GetBinContent(tmpRatio->GetMaximumBin()) + tmpRatio->GetBinError(tmpRatio->GetMaximumBin()) );

    tmpRatio->SetMaximum( setMax );
    if ( plotMin != -1 )
        tmpRatio->SetMinimum( plotMin );
    if ( plotMax != -1 )
        tmpRatio->SetMaximum( plotMax );


    tmpRatio->DrawCopy("X0");

    if ( drawOneLine )
    {
        const TAxis *axis = tmpRatio->GetXaxis();
        double lowX  = axis->GetBinLowEdge( axis->GetFirst() );
        double highX = axis->GetBinUpEdge(  axis->GetLast() );

        TLine line;
        line.SetLineStyle(2);
        line.SetLineWidth(3);
        line.SetLineColor(36);
        line.DrawLine(lowX, 1., highX, 1.); //creates a new line which is owned by gPad
    }

    gPad->RedrawAxis();
    gPad->Update();
}

//=====================================================================
// draw MC with error band + data points on current Pad with stat error
// draw sys error band as an envelope around 1
//=====================================================================
void MnvPlotter::DrawDataMCRatio(
        const MnvH1D* dataHist,
        const MnvH1D *mcHist,
        const Double_t mcScale,     /*= 1.0*/
        const bool     drawSysLines, /*= true*/
        const bool     drawOneLine,  /*= true*/
        const double   plotMin,      /*=-1 (auto)*/
        const double   plotMax,      /*=-1 (auto)*/
        const char*    yaxisLabel,    /*="Data / MC"*/
        const bool     covAreaNormalize /*= false*/
        )
{
    if (!gPad)
      throw std::runtime_error("MnvPlotter requires a TCanvas. Please make one first.");

    //!make clones
    TH1Ptr tmpData( (TH1*)dataHist->Clone("tmp_data") );
    TH1Ptr tmpMC  ( (TH1*)mcHist  ->Clone("tmp_mc"  ) );

    if ( draw_normalized_to_bin_width )
    {
        if ( dataHist->GetNormBinWidth() > 0 )
            tmpData->Scale( dataHist->GetNormBinWidth(), "width" );

        if ( mcHist->GetNormBinWidth() > 0 )
            tmpMC  ->Scale( mcHist  ->GetNormBinWidth(), "width" );
    }


    //scale MC to data
    tmpMC->Scale(mcScale);

    //get the ratio
    TH1Ptr tmpRatio( (TH1*)tmpMC->Clone("tmp_ratio") );
    tmpRatio->Divide(tmpData, tmpMC);

    //get the total systematic error as a fraction of the CV
    TH1D sysErr = mcHist->GetTotalError(false, true, covAreaNormalize); //no stat error, get as fraction

    //!set the max bin based on the taller of the ratio or error
    double setMax = headroom * max(
            tmpRatio->GetBinContent( tmpRatio->GetMaximumBin() ),
            1. + sysErr.GetBinContent( sysErr.GetMaximumBin() )
            );

    tmpRatio->SetMaximum( setMax );
    if ( plotMin != -1 )
        tmpRatio->SetMinimum( plotMin );
    if ( plotMax != -1 )
        tmpRatio->SetMaximum( plotMax );


    //apply style
    tmpRatio->GetYaxis()->SetTitle( yaxisLabel );
    tmpRatio->GetXaxis()->SetTitleFont(axis_title_font_x);
    tmpRatio->GetYaxis()->SetTitleFont(axis_title_font_y);
    tmpRatio->GetXaxis()->SetTitleSize(axis_title_size_x);
    tmpRatio->GetYaxis()->SetTitleSize(axis_title_size_y);
    tmpRatio->GetXaxis()->SetLabelFont(axis_label_font);
    tmpRatio->GetYaxis()->SetLabelFont(axis_label_font);
    tmpRatio->GetXaxis()->SetLabelSize(axis_label_size);
    tmpRatio->GetYaxis()->SetLabelSize(axis_label_size);
    tmpRatio->GetXaxis()->CenterTitle(kTRUE);
    tmpRatio->SetMarkerStyle(ratio_marker);
    tmpRatio->SetMarkerSize(ratio_marker_size);
    tmpRatio->SetLineWidth(ratio_line_width);
    tmpRatio->SetLineColor(ratio_color);
    tmpRatio->GetXaxis()->SetNdivisions(509);
    tmpRatio->DrawCopy("X0");

    //add sys error lines if desired
    if ( drawSysLines )
    {
        //get the total sys error.  then fill a histogram of total negative sys error for the lower lines
        sysErr.SetMaximum( setMax );

        //content of sys error histogram should be 1 (data=MC)
        //error should be the error
        int firstNonZeroBin = 0 ;
        int lastNonZeroBin = 0 ;
        for ( int ibin = 1; ibin <= sysErr.GetNbinsX(); ++ibin )
        {
            if ( firstNonZeroBin == 0 && sysErr.GetBinContent(ibin) > 0 )
                firstNonZeroBin = ibin;
            if ( sysErr.GetBinContent(ibin) > 0 )
                lastNonZeroBin = ibin;

            double err = sysErr.GetBinContent(ibin);
            sysErr.SetBinContent( ibin, 1. );
            sysErr.SetBinError( ibin, err );
        }

        if ( lastNonZeroBin < 0 ) //just in case
            lastNonZeroBin = sysErr.GetNbinsX();

        if ( firstNonZeroBin != lastNonZeroBin ) {

            sysErr.GetXaxis()->SetRange( firstNonZeroBin, lastNonZeroBin );
            sysErr.SetFillColor( mc_error_color );
            sysErr.SetFillStyle( mc_error_style );
            sysErr.SetLineColor( mc_error_color );
            sysErr.SetLineWidth( mc_line_width );
            sysErr.SetMarkerStyle( 0 );
            sysErr.DrawCopy("E2 same ][");
            tmpRatio->DrawCopy("SAME AXIS");//get tickmarks back
            tmpRatio->DrawCopy("SAME X0"); //need to draw the ratio again
        }
    }

    //draw line at 1. if desired
    if ( drawOneLine )
    {
        const TAxis *axis = tmpRatio->GetXaxis();
        double lowX  = axis->GetBinLowEdge( axis->GetFirst() );
        double highX = axis->GetBinUpEdge(  axis->GetLast() );

        TLine line;
        line.SetLineStyle(2);
        line.SetLineWidth(3);
        line.SetLineColor(36);
        line.DrawLine(lowX, 1., highX, 1.); //creates a new line which is owned by gPad
    }

    gPad->RedrawAxis();
    gPad->Update();
}

//==============================================================
// calculate the chi2 between two MnvH2Ds
//==============================================================
Double_t MnvPlotter::Chi2DataMC(
        const MnvH2D* dataHist,
        const MnvH2D* mcHist,
        Int_t& ndf,
        const Double_t mcScale,
        const bool useDataErrorMatrix,
        const bool useOnlyShapeErrors,
	const bool useModelStat,
	TMatrixD *Chi2ByBin
        )
{
    ndf=0; // This will be filled with the number of degrees of freedom

    //We get the number of bins and make sure it's the same in data and MC
    if ( dataHist->GetNbinsX() != mcHist->GetNbinsX() )
    {
        Error("MnvPlotter::Chi2DataMC", "The number of x bins from Data and MC 2D histograms differ. Returning -1.");
        return -1.;
    }

    if ( dataHist->GetNbinsY() != mcHist->GetNbinsY() )
    {
        Error("MnvPlotter::Chi2DataMC", "The number of y bins from Data and MC 2D histograms differ. Returning -1.");
        return -1.;
    }

    //only consider the plotted range (no underflow and overflow unless specified)
    const int lowBin  = chi2_use_overflow_err?0:1; // Either they both use underflow, or neither of them does
    int nbinsx=mcHist->GetNbinsX();
    int nbinsy=mcHist->GetNbinsY();
    const int highBinX = nbinsx + (chi2_use_overflow_err?1:0);
    const int highBinY = nbinsy + (chi2_use_overflow_err?1:0);

    //Scaling MC to Data
    MnvH2D* tmpMCHist = (MnvH2D*) mcHist->Clone("tmpMCHist");
    tmpMCHist->Scale(mcScale);

    // Defining Error Matrix dimensions
    // Either use the requested range or the full error matrix with under/overflow
    Int_t Nbins = (highBinX-lowBin+1)* (highBinY-lowBin+1); // from low to high inclusive in each dimension

    //get the covariance matrix
    TMatrixD covMatrix(Nbins, Nbins);
    {
        const Int_t NbinsTotal = (nbinsx + 2)*(nbinsy + 2); //+2 for under/overflow
        TMatrixD covMatrixTmp( NbinsTotal, NbinsTotal ); // covariance matrix including under- and over-flow
        const bool includeStatError = true;
        const bool errorAsFraction  = false;

        // Use Error Matrix from Data or MC?
        if ( useDataErrorMatrix )
        {
            covMatrixTmp = dataHist->GetTotalErrorMatrix( includeStatError, errorAsFraction, useOnlyShapeErrors);
            if (useModelStat)covMatrixTmp += tmpMCHist->GetStatErrorMatrix( );
        }
        else
        {
            covMatrixTmp = tmpMCHist->GetTotalErrorMatrix(  includeStatError, errorAsFraction, useOnlyShapeErrors);
            if (useModelStat)covMatrixTmp += dataHist->GetStatErrorMatrix( );
        }

        //select only covariance elements in the histogram range
        // unless using the contributions to covariance from overflow
        if ( chi2_use_overflow_err )
        {
            covMatrix = covMatrixTmp;
        }
        else // this is going to be a nuisance for 2 dimensions
        {
            // We need to weed out anything that is underflow or overflow in either the x or the y dimension. Ugh.
            for ( int i = 0; i != Nbins; ++i )
            {
                for ( int j = 0; j != Nbins; ++j )
                {
                    // I think this is right... it checks out on a small sample
                    int old_i_bin = (i/nbinsx + 1)* (nbinsx +2) +1 + i%nbinsx;
                    int old_j_bin = (j/nbinsx + 1)* (nbinsx +2) +1 + j%nbinsx;
                    covMatrix[i][j] = covMatrixTmp[old_i_bin][old_j_bin];
                }
            }
        }

    }// end of block to get covariance
    //Now, we invert the covariance error Matrix and store the result in "errorMatrix"
    // Note: TDecompSVD can handle singular matrices
    covMatrix*=1e80;//ROOT can't seem to handle small entries with lots of zeros? Suggested scaling the histogram and then rescaling the inverted matrix

    TDecompSVD error(covMatrix);
    TMatrixD errorMatrix(covMatrix);
    if ( ! error.Invert( errorMatrix ) )
    {
        Warning("MnvPlotter::Chi2DataMC", "Cannot invert total covariance matrix. You could use statistical errors only for Chi2 calculation. But it isn't implemented yet.");
        return -1;
    }
    errorMatrix*=1e80;//ROOT can't seem to handle small entries with lots of zeros? Suggested scaling the histogram and then rescaling the inverted matrix
    if (Chi2ByBin) {
      //cout << "Resizing the matrix you gave me for the chi2. Orignal size " << Chi2ByBin->GetNrows()<<endl;
      Chi2ByBin->ResizeTo(Nbins,Nbins);
      //cout << "New size row " << Chi2ByBin->GetNrows()<<endl;
    }
    //Calculating chi2
    ndf = 0;
    Double_t chi2 = 0.;
    //cout << "Calculating 2D chi2" << endl;
    for ( int i = 0; i != Nbins; ++i )
    {
        int hist_i_bin = chi2_use_overflow_err?i:((i/nbinsx + 1)* (nbinsx +2) +1 + i%nbinsx); // Translate to the histogram bin, if we aren't using the overflow errors, meaning the covariance matrix is smaller than the histogram
        const Double_t x_data_i = dataHist->GetBinContent(hist_i_bin);
        const Double_t x_mc_i   = tmpMCHist->GetBinContent(hist_i_bin);
        for ( int j = 0; j != Nbins; ++j )
        {
            // Each element of the inverted covariance matrix corresponds to a pair of data and MC
            int hist_j_bin = chi2_use_overflow_err?j:((j/nbinsx + 1)* (nbinsx +2) +1 + j%nbinsx);
            const Double_t x_data_j = dataHist->GetBinContent(hist_j_bin);
            const Double_t x_mc_j   = tmpMCHist->GetBinContent(hist_j_bin);
            const double chi2_ij = (x_data_i - x_mc_i) * errorMatrix[i][j] * (x_data_j - x_mc_j);

	    if (Chi2ByBin) Chi2ByBin[0][i][j]=chi2_ij;//Dunno why this work... TBD
            chi2 += chi2_ij;
        }
        ++ndf;
    }
 //   cout << "Done calculating 2D chi2" << endl;
    // if this is a shape comparison, subtract one degree of freedom
    if (useOnlyShapeErrors)
        --ndf;

    delete tmpMCHist;

    return chi2;
}

Double_t MnvPlotter::Chi2DataMC(
                                const MnvH2D* dataHist,
                                const MnvH2D* mcHist,
                                const Double_t mcScale,
                                const bool useDataErrorMatrix,
                                const bool useOnlyShapeErrors,
                                const bool useModelStat,
                                TMatrixD *Chi2ByBin
                                )
{
    int ndf=0; // This will be filled with the number of degrees of freedom

    //We get the number of bins and make sure it's the same in data and MC
    if ( dataHist->GetNbinsX() != mcHist->GetNbinsX() )
    {
        Error("MnvPlotter::Chi2DataMC", "The number of x bins from Data and MC 2D histograms differ. Returning -1.");
        return -1.;
    }

    if ( dataHist->GetNbinsY() != mcHist->GetNbinsY() )
    {
        Error("MnvPlotter::Chi2DataMC", "The number of y bins from Data and MC 2D histograms differ. Returning -1.");
        return -1.;
    }

    //only consider the plotted range (no underflow and overflow unless specified)
    const int lowBin  = chi2_use_overflow_err?0:1; // Either they both use underflow, or neither of them does
    int nbinsx=mcHist->GetNbinsX();
    int nbinsy=mcHist->GetNbinsY();
    const int highBinX = nbinsx + (chi2_use_overflow_err?1:0);
    const int highBinY = nbinsy + (chi2_use_overflow_err?1:0);

    //Scaling MC to Data
    MnvH2D* tmpMCHist = (MnvH2D*) mcHist->Clone("tmpMCHist");
    tmpMCHist->Scale(mcScale);

    // Defining Error Matrix dimensions
    // Either use the requested range or the full error matrix with under/overflow
    Int_t Nbins = (highBinX-lowBin+1)* (highBinY-lowBin+1); // from low to high inclusive in each dimension

    //get the covariance matrix
    TMatrixD covMatrix(Nbins, Nbins);
    {
        const Int_t NbinsTotal = (nbinsx + 2)*(nbinsy + 2); //+2 for under/overflow
        TMatrixD covMatrixTmp( NbinsTotal, NbinsTotal ); // covariance matrix including under- and over-flow
        const bool includeStatError = true;
        const bool errorAsFraction  = false;

        // Use Error Matrix from Data or MC?
        if ( useDataErrorMatrix )
        {
            covMatrixTmp = dataHist->GetTotalErrorMatrix( includeStatError, errorAsFraction, useOnlyShapeErrors);
            if (useModelStat)covMatrixTmp += tmpMCHist->GetStatErrorMatrix( );
        }
        else
        {
            covMatrixTmp = tmpMCHist->GetTotalErrorMatrix(  includeStatError, errorAsFraction, useOnlyShapeErrors);
            if (useModelStat)covMatrixTmp += dataHist->GetStatErrorMatrix( );
        }

        //select only covariance elements in the histogram range
        // unless using the contributions to covariance from overflow
        if ( chi2_use_overflow_err )
        {
            covMatrix = covMatrixTmp;
        }
        else // this is going to be a nuisance for 2 dimensions
        {
            // We need to weed out anything that is underflow or overflow in either the x or the y dimension. Ugh.
            for ( int i = 0; i != Nbins; ++i )
            {
                for ( int j = 0; j != Nbins; ++j )
                {
                    // I think this is right... it checks out on a small sample
                    int old_i_bin = (i/nbinsx + 1)* (nbinsx +2) +1 + i%nbinsx;
                    int old_j_bin = (j/nbinsx + 1)* (nbinsx +2) +1 + j%nbinsx;
                    covMatrix[i][j] = covMatrixTmp[old_i_bin][old_j_bin];
                }
            }
        }

    }// end of block to get covariance


    //Now, we invert the covariance error Matrix and store the result in "errorMatrix"
    // Note: TDecompSVD can handle singular matrices
    covMatrix*=1e80;//ROOT can't seem to handle small entries with lots of zeros? Suggested scaling the histogram and then rescaling the inverted matrix

    TDecompSVD error(covMatrix);
    TMatrixD errorMatrix(covMatrix);
    if ( ! error.Invert( errorMatrix ) )
    {
        Warning("MnvPlotter::Chi2DataMC", "Cannot invert total covariance matrix. You could use statistical errors only for Chi2 calculation. But it isn't implemented yet.");
        return -1;
    }
    errorMatrix*=1e80;//ROOT can't seem to handle small entries with lots of zeros? Suggested scaling the histogram and then rescaling the inverted matrix
    if (Chi2ByBin) {
        //cout << "Resizing the matrix you gave me for the chi2. Orignal size " << Chi2ByBin->GetNrows()<<endl;
        Chi2ByBin->ResizeTo(Nbins,Nbins);
        //cout << "New size row " << Chi2ByBin->GetNrows()<<endl;
    }
    //Calculating chi2
    ndf = 0;
    Double_t chi2 = 0.;
    //cout << "Calculating 2D chi2" << endl;
    for ( int i = 0; i != Nbins; ++i )
    {
        int hist_i_bin = chi2_use_overflow_err?i:((i/nbinsx + 1)* (nbinsx +2) +1 + i%nbinsx); // Translate to the histogram bin, if we aren't using the overflow errors, meaning the covariance matrix is smaller than the histogram
        const Double_t x_data_i = dataHist->GetBinContent(hist_i_bin);
        const Double_t x_mc_i   = tmpMCHist->GetBinContent(hist_i_bin);
        for ( int j = 0; j != Nbins; ++j )
        {
            // Each element of the inverted covariance matrix corresponds to a pair of data and MC
            int hist_j_bin = chi2_use_overflow_err?j:((j/nbinsx + 1)* (nbinsx +2) +1 + j%nbinsx);
            const Double_t x_data_j = dataHist->GetBinContent(hist_j_bin);
            const Double_t x_mc_j   = tmpMCHist->GetBinContent(hist_j_bin);
            const double chi2_ij = (x_data_i - x_mc_i) * errorMatrix[i][j] * (x_data_j - x_mc_j);

            if (Chi2ByBin) Chi2ByBin[0][i][j]=chi2_ij;//Dunno why this work... TBD
            chi2 += chi2_ij;
        }
        ++ndf;
    }
    //   cout << "Done calculating 2D chi2" << endl;
    // if this is a shape comparison, subtract one degree of freedom
    if (useOnlyShapeErrors)
        --ndf;

    delete tmpMCHist;

    return chi2;
}




//==============================================================
// calculate the chi2 between two histograms
//==============================================================

Double_t MnvPlotter::Chi2DataMC(
        const TH1* dataHist,
        const TH1 *mcHist,
        int& ndf,
        const Double_t mcScale /*= 1.0*/
        )
{
    TH1Ptr tmpData( (TH1D*)dataHist->Clone("tmp_data_chi2") );
    TH1Ptr tmpMC  ( (TH1D*)mcHist  ->Clone("tmp_mc_chi2"  ) );

    if ( tmpData->GetSumw2N() == 0 )
        tmpData->Sumw2();
    if ( tmpMC->GetSumw2N() == 0 )
        tmpMC->Sumw2();

    tmpMC->Scale(mcScale);

    Double_t chi2 = 0;
    ndf = 0;

    const int lowBin  = tmpMC->GetXaxis()->GetFirst();
    const int highBin = tmpMC->GetXaxis()->GetLast();

    for ( int i = lowBin; i <= highBin; ++i )
    {
        if ( tmpData->GetBinError(i)*tmpData->GetBinError(i) + tmpMC->GetBinError(i)*tmpMC->GetBinError(i) > 0 )
        {
            chi2 += (tmpData->GetBinContent(i) - tmpMC->GetBinContent(i))
                *(tmpData->GetBinContent(i) - tmpMC->GetBinContent(i))
                /(tmpData->GetBinError(i)*tmpData->GetBinError(i) + tmpMC->GetBinError(i)*tmpMC->GetBinError(i));
            ++ndf;
        }
    }

    return chi2;

}

Double_t MnvPlotter::Chi2DataMC(
                                const TH1* dataHist,
                                const TH1 *mcHist,
                                const Double_t mcScale /*= 1.0*/
)
{
    TH1Ptr tmpData( (TH1D*)dataHist->Clone("tmp_data_chi2") );
    TH1Ptr tmpMC  ( (TH1D*)mcHist  ->Clone("tmp_mc_chi2"  ) );

    if ( tmpData->GetSumw2N() == 0 )
        tmpData->Sumw2();
    if ( tmpMC->GetSumw2N() == 0 )
        tmpMC->Sumw2();

    tmpMC->Scale(mcScale);

    Double_t chi2 = 0;
    int ndf = 0;

    const int lowBin  = tmpMC->GetXaxis()->GetFirst();
    const int highBin = tmpMC->GetXaxis()->GetLast();

    for ( int i = lowBin; i <= highBin; ++i )
    {
        if ( tmpData->GetBinError(i)*tmpData->GetBinError(i) + tmpMC->GetBinError(i)*tmpMC->GetBinError(i) > 0 )
        {
            chi2 += (tmpData->GetBinContent(i) - tmpMC->GetBinContent(i))
            *(tmpData->GetBinContent(i) - tmpMC->GetBinContent(i))
            /(tmpData->GetBinError(i)*tmpData->GetBinError(i) + tmpMC->GetBinError(i)*tmpMC->GetBinError(i));
            ++ndf;
        }
    }

    return chi2;

}


Double_t MnvPlotter::Chi2DataMC(
        const MnvH1D* dataHist,
        const MnvH1D* mcHist,
        int& ndf,
        const Double_t mcScale,
        const bool useDataErrorMatrix,
        const bool useOnlyShapeErrors,
        const bool useModelStat,
        TMatrixD *Chi2ByBin)
{
    //We get the number of bins and make sure it's compatible with the NxN matrix given
    if ( dataHist->GetNbinsX() != mcHist->GetNbinsX() )
    {
        Error("MnvPlotter::Chi2DataMC", "The number of bins from Data and MC histograms differ. Returning -1.");
        return -1.;
    }

    //only consider the plotted range
    const int lowBin  = chi2_use_overflow_err ? 0 : mcHist->GetXaxis()->GetFirst();
    const int highBin = chi2_use_overflow_err ? mcHist->GetNbinsX()+1 : mcHist->GetXaxis()->GetLast();

    //Scaling MC to Data
    MnvH1D* tmpMCHist = (MnvH1D*) mcHist->Clone("tmpMCHist");
    tmpMCHist->Scale(mcScale);

    // Defining Error Matrix dimensions
    // Either use the requested range or the full error matrix with under/overflow
    Int_t Nbins = highBin - lowBin + 1;

    //get the covariance matrix
    TMatrixD covMatrix(Nbins, Nbins);
    {
        const Int_t NbinsTotal = tmpMCHist->GetNbinsX() + 2; //+2 for under/overflow
        TMatrixD covMatrixTmp( NbinsTotal, NbinsTotal );
        const bool includeStatError = true;
        const bool errorAsFraction  = false;

        // Use Error Matrix from Data or MC?
        if ( useDataErrorMatrix )
        {
            covMatrixTmp = dataHist->GetTotalErrorMatrix( includeStatError, errorAsFraction, useOnlyShapeErrors);
            if (useModelStat) covMatrixTmp += tmpMCHist->GetStatErrorMatrix( );
        }
        else
        {
            covMatrixTmp = tmpMCHist->GetTotalErrorMatrix(  includeStatError, errorAsFraction, useOnlyShapeErrors);
            if (useModelStat) covMatrixTmp += dataHist->GetStatErrorMatrix( );
        }

        //select only covariance elements in the histogram range
        // unless using the contributions to covariance from overflow
        if ( chi2_use_overflow_err )
        {
            covMatrix = covMatrixTmp;
        }
        else
        {
            for ( int i = 0; i != Nbins; ++i )
            {
                for ( int j = 0; j != Nbins; ++j )
                {
                    covMatrix[i][j] = covMatrixTmp[i+lowBin][j+lowBin];
                }
            }
        }

    }// end of block to get covarance

    //Now, we invert the covariance error Matrix and store the result in "errorMatrix"
    // Note: TDecompSVD can handle singular matrices

     //Now, we invert the covariance error Matrix and store the result in "errorMatrix"
    // Note: TDecompSVD can handle singular matrices

    covMatrix*=1e80;//ROOT can't seem to handle small entries with lots of zeros? Suggested scaling the histogram and then rescaling the inverted matrix

    TDecompSVD error(covMatrix);
    TMatrixD errorMatrix(covMatrix);
    if ( ! error.Invert( errorMatrix ) )
    {
        Warning("MnvPlotter::Chi2DataMC", "Cannot invert total covariance matrix.  Using statistical errors only for Chi2 calculation.");
        return Chi2DataMC( (TH1D*)dataHist, (TH1D*)mcHist, ndf, mcScale );
    }
    errorMatrix*=1e80;//ROOT can't seem to handle small entries with lots of zeros? Suggested scaling the histogram and then rescaling the inverted matrix

    if (Chi2ByBin) {
      //cout << "Resizing the matrix you gave me for the chi2. Orignal size " << Chi2ByBin->GetNrows()<<endl;
      Chi2ByBin->ResizeTo(Nbins,Nbins);
      //cout << "New size row " << Chi2ByBin->GetNrows()<<endl;
    }

    //Calculating chi2
    ndf = 0;
    Double_t chi2 = 0.;
    // under/overflow bins not taken into account in the chi2 calculation
    for (int i=lowBin; i<=highBin ; ++i)
    {
        const int iErrBin = i - lowBin;
        const Double_t x_data_i = dataHist->GetBinContent(i);
        const Double_t x_mc_i   = tmpMCHist->GetBinContent(i);

        for (int j=lowBin; j<=highBin; ++j)
        {
            const int jErrBin = j - lowBin;
            const Double_t x_data_j = dataHist->GetBinContent(j);
            const Double_t x_mc_j   = tmpMCHist->GetBinContent(j);
            const double chi2_ij = (x_data_i-x_mc_i!=0 && x_data_j-x_mc_j !=0) ?(x_data_i - x_mc_i) * errorMatrix[iErrBin][jErrBin] * (x_data_j - x_mc_j): 0;
            if (Chi2ByBin) Chi2ByBin[0][iErrBin][jErrBin]=chi2_ij;//Dunno why this work... TBD
            chi2 += chi2_ij;
        }
        ++ndf; // Is this the right way to calcualte ndf?
    }

    // if this is a shape comparison, subtract one degree of freedom?
    if (useOnlyShapeErrors)
        --ndf;

    delete tmpMCHist;

    return chi2;
}

Double_t MnvPlotter::Chi2DataMC(
                                const MnvH1D* dataHist,
                                const MnvH1D* mcHist,
                                const Double_t mcScale,
                                const bool useDataErrorMatrix,
                                const bool useOnlyShapeErrors,
                                const bool useModelStat,
                                TMatrixD *Chi2ByBin)
{
    //We get the number of bins and make sure it's compatible with the NxN matrix given
    if ( dataHist->GetNbinsX() != mcHist->GetNbinsX() )
    {
        Error("MnvPlotter::Chi2DataMC", "The number of bins from Data and MC histograms differ. Returning -1.");
        return -1.;
    }

    //only consider the plotted range
    const int lowBin  = chi2_use_overflow_err ? 0 : mcHist->GetXaxis()->GetFirst();
    const int highBin = chi2_use_overflow_err ? mcHist->GetNbinsX()+1 : mcHist->GetXaxis()->GetLast();

    //Scaling MC to Data
    MnvH1D* tmpMCHist = (MnvH1D*) mcHist->Clone("tmpMCHist");
    tmpMCHist->Scale(mcScale);

    // Defining Error Matrix dimensions
    // Either use the requested range or the full error matrix with under/overflow
    Int_t Nbins = highBin - lowBin + 1;
    if ( chi2_use_overflow_err )
        Nbins = tmpMCHist->GetNbinsX()+2;

    //get the covariance matrix
    TMatrixD covMatrix(Nbins, Nbins);
    {
        const Int_t NbinsTotal = tmpMCHist->GetNbinsX() + 2; //+2 for under/overflow
        TMatrixD covMatrixTmp( NbinsTotal, NbinsTotal );
        const bool includeStatError = true;
        const bool errorAsFraction  = false;

        // Use Error Matrix from Data or MC?
        if ( useDataErrorMatrix )
        {
            covMatrixTmp = dataHist->GetTotalErrorMatrix( includeStatError, errorAsFraction, useOnlyShapeErrors);
            if (useModelStat) covMatrixTmp += tmpMCHist->GetStatErrorMatrix( );
        }
        else
        {
            covMatrixTmp = tmpMCHist->GetTotalErrorMatrix(  includeStatError, errorAsFraction, useOnlyShapeErrors);
            if (useModelStat) covMatrixTmp += dataHist->GetStatErrorMatrix( );
        }

        //select only covariance elements in the histogram range
        // unless using the contributions to covariance from overflow
        if ( chi2_use_overflow_err )
        {
            covMatrix = covMatrixTmp;
        }
        else
        {
            for ( int i = 0; i != Nbins; ++i )
            {
                for ( int j = 0; j != Nbins; ++j )
                {
                    covMatrix[i][j] = covMatrixTmp[i+lowBin][j+lowBin];
                }
            }
        }

    }// end of block to get covarance

    //Now, we invert the covariance error Matrix and store the result in "errorMatrix"
    // Note: TDecompSVD can handle singular matrices

    //Now, we invert the covariance error Matrix and store the result in "errorMatrix"
    // Note: TDecompSVD can handle singular matrices

    covMatrix*=1e80;//ROOT can't seem to handle small entries with lots of zeros? Suggested scaling the histogram and then rescaling the inverted matrix

    TDecompSVD error(covMatrix);
    TMatrixD errorMatrix(covMatrix);
    if ( ! error.Invert( errorMatrix ) )
    {
        Warning("MnvPlotter::Chi2DataMC", "Cannot invert total covariance matrix.  Using statistical errors only for Chi2 calculation.");
        return Chi2DataMC( (TH1D*)dataHist, (TH1D*)mcHist,  mcScale );
    }
    errorMatrix*=1e80;//ROOT can't seem to handle small entries with lots of zeros? Suggested scaling the histogram and then rescaling the inverted matrix

    if (Chi2ByBin) {
        //cout << "Resizing the matrix you gave me for the chi2. Orignal size " << Chi2ByBin->GetNrows()<<endl;
        Chi2ByBin->ResizeTo(Nbins,Nbins);
        //cout << "New size row " << Chi2ByBin->GetNrows()<<endl;
    }

    //Calculating chi2
    int ndf = 0;
    Double_t chi2 = 0.;
    // under/overflow bins not taken into account in the chi2 calculation
    for (int i=lowBin; i<=highBin ; ++i)
    {
        const int iErrBin = i - lowBin;
        const Double_t x_data_i = dataHist->GetBinContent(i);
        const Double_t x_mc_i   = tmpMCHist->GetBinContent(i);

        for (int j=lowBin; j<=highBin; ++j)
        {
            const int jErrBin = j - lowBin;
            const Double_t x_data_j = dataHist->GetBinContent(j);
            const Double_t x_mc_j   = tmpMCHist->GetBinContent(j);
            const double chi2_ij = (x_data_i - x_mc_i) * errorMatrix[iErrBin][jErrBin] * (x_data_j - x_mc_j);
            if (Chi2ByBin) Chi2ByBin[0][iErrBin][jErrBin]=chi2_ij;//Dunno why this work... TBD
            chi2 += chi2_ij;
        }
        ++ndf; // Is this the right way to calcualte ndf?
    }

    // if this is a shape comparison, subtract one degree of freedom?
    if (useOnlyShapeErrors)
        --ndf;

    delete tmpMCHist;

    return chi2;
}


Double_t MnvPlotter::Chi2MCMC(
        const MnvH1D* histA,
        const MnvH1D* histB,
        int& ndf,
        const Double_t bScale /*= 1.0 */,
        const bool useOnlyShapeErrors /*= false*/
        )
{
    ndf = 0;

    //We get the number of bins and make sure it's compatible with the NxN matrix given
    if ( histA->GetNbinsX() != histB->GetNbinsX() )
    {
        Error("MnvPlotter::Chi2MCMC", "The number of bins from histA and histB histograms differ. Returning -1.");
        return -1.;
    }


    // Make a copy of B
    MnvH1D bCopy(*histB);
    bCopy.SetName( "tmpBCopy" );
    bCopy.SetDirectory(0); //destroy when out of scope

    //take the ratio B/A and apply bScale
    bCopy.Divide( histB, histA, bScale );

    //create hist with all bin values of 1 with no error
    MnvH1D oneHist( *dynamic_cast<const TH1D*>(histA) );
    for ( int iBin = 0; iBin <= histA->GetNbinsX()+1; ++iBin )
    {
        oneHist.SetBinContent( iBin, 1. );
        oneHist.SetBinError(   iBin, 0. );
    }

    bCopy.GetXaxis()->SetRange( histA->GetXaxis()->GetFirst(), histA->GetXaxis()->GetLast() );

    // compare ratio to line of 1
    const double tmpScale = 1.;
    const bool useDataSys = false;
    const double chi2 = Chi2DataMC( &oneHist, &bCopy, ndf, tmpScale, useDataSys, useOnlyShapeErrors );

    return chi2;
}



void MnvPlotter::DrawDataMCVariations(
        const MnvH1D* dataHist,
        const TObjArray* mcHists,
        const Double_t mcScale     /*= 1.0   */,
        const std::string& legPos /*= "TR"  */,
        const bool dataAsPoints     /*= true  */,
        const bool allSolidLines    /*= false */,
        const bool leaveStyleAlone  /*= false */,
        const bool covAreaNormalize /*=false*/
        )
{
    if (!gPad)
      throw std::runtime_error("MnvPlotter requires a TCanvas. Please make one first.");

    unsigned int nHists = mcHists->GetEntries();

    vector<string> titles;
    titles.push_back( dataHist->GetTitle() );

    // Find the maximum among all histograms
    // be sure to consider a logY axis.
    // todo: when logY, we really want to know the minimum non-zero value
    double maxmax = dataHist->GetMaximum();
    double minmin = dataHist->GetMinimum();
    if ( gPad->GetLogy() && minmin < 1.E-8 )
        minmin = -1.; //-1 for not set if the y-axis is log scale

    for ( unsigned int i = 0; i != nHists; ++i )
    {
        MnvH1D *mnv = (MnvH1D*)mcHists->At(i);
        maxmax = std::max( maxmax, mcScale * mnv->GetMaximum() );
        if ( 1.E-8 < mnv->GetMinimum() )
            minmin = std::min( minmin, mcScale * mnv->GetMinimum() );
        titles.push_back( mnv->GetTitle() );
    }
    if ( gPad->GetLogy() )
        minmin = std::max( minmin, 1.E-8 );

    size_t titleWidth = GetLegendWidthInLetters( titles );
    double x1, x2, y1, y2;
    DecodeLegendPosition( x1, y1, x2, y2, legPos, nHists+1, titleWidth, legend_text_size );

    TLegend *leg  = new TLegend(x1, y1, x2, y2);
    leg->SetNColumns( legend_n_columns );
    static int legB_n = 0;
    leg->SetName( Form("legend B - %d", legB_n++) );
    AddToTmp( leg );

    leg->SetBorderSize( legend_border_size );
    leg->SetFillColor( legend_fill_color );
    if ( legend_fill_color < 0 )
        leg->SetFillStyle(0);
    leg->SetTextSize( legend_text_size );
    leg->SetTextFont( legend_text_font );

    if ( draw_normalized_to_bin_width && covAreaNormalize )
        Warning("MnvPlotter::DrawDataMCVariations","Area normalized covariance matrix may be incorrect because it is taken after bin normalization.  See DrawDataMCWithErrorBand(MnvH1Ds) for an example what to do");

    //create as clone because it gets added to the leged?
    TH1* tmpData(0);
    if ( draw_normalized_to_bin_width )
        tmpData = (TH1*)dataHist->GetBinNormalizedCopy().GetCVHistoWithError(true, covAreaNormalize).Clone( Form("tmpData_%d", __LINE__) );
    else
        tmpData = (TH1*)dataHist->GetCVHistoWithError(true, covAreaNormalize).Clone( Form("tmpData_%d", __LINE__) );
    AddToTmp( tmpData );

    //! Use ApplyNextLineStyle on data always, to make sure it resets and uses black
    if ( ! leaveStyleAlone )
    {
        tmpData->SetMarkerStyle(data_marker);
        tmpData->SetMarkerSize(data_marker_size);
        ApplyNextLineStyle( tmpData, true, !allSolidLines );
        tmpData->GetXaxis()->SetTitleFont(axis_title_font_x);
        tmpData->GetYaxis()->SetTitleFont(axis_title_font_y);
        tmpData->GetXaxis()->SetTitleSize(axis_title_size_x);
        tmpData->GetYaxis()->SetTitleSize(axis_title_size_y);
        tmpData->GetXaxis()->SetLabelFont(axis_label_font);
        tmpData->GetYaxis()->SetLabelFont(axis_label_font);
        tmpData->GetXaxis()->SetLabelSize(axis_label_size);
        tmpData->GetYaxis()->SetLabelSize(axis_label_size);
        tmpData->GetXaxis()->CenterTitle(kTRUE);
    }
    //Add data to the legend first
    if ( dataAsPoints )
        leg->AddEntry(tmpData,tmpData->GetTitle(),"ple");
    else
        leg->AddEntry(tmpData,tmpData->GetTitle(),"l");



    //! Loop over the MC histograms.  Apply style, add to legend and draw.
    for ( unsigned int i = 0; i < nHists; i++ )
    {
        const MnvH1D *mnvMC = dynamic_cast<const MnvH1D*>( mcHists->At(i) );
        if ( !mnvMC )
        {
            Error( "DrawDataMCVariations", "Could not cast one of the MC histograms to MnvH1D.  Your draw is incomplete." );
            continue;
        }


        TH1 *hst(0);
        if ( draw_normalized_to_bin_width )
            hst = (TH1*)mnvMC->GetBinNormalizedCopy().GetCVHistoWithError(true, covAreaNormalize).Clone( Form("tmp_MCHist_%d_%d", i, __LINE__) );
        else
            hst = (TH1*)mnvMC->GetCVHistoWithError(true, covAreaNormalize).Clone( Form("tmp_MCHist_%d_%d", i, __LINE__) );
        AddToTmp( hst );

        if ( ! leaveStyleAlone )
        {
            ApplyNextLineStyle( hst, false, !allSolidLines );
            hst->SetMarkerStyle(0);
            hst->GetXaxis()->SetTitleFont(axis_title_font_x);
            hst->GetYaxis()->SetTitleFont(axis_title_font_y);
            hst->GetXaxis()->SetTitleSize(axis_title_size_x);
            hst->GetYaxis()->SetTitleSize(axis_title_size_y);
            hst->GetXaxis()->SetLabelFont(axis_label_font);
            hst->GetYaxis()->SetLabelFont(axis_label_font);
            hst->GetXaxis()->SetLabelSize(axis_label_size);
            hst->GetYaxis()->SetLabelSize(axis_label_size);
            hst->GetXaxis()->CenterTitle(kTRUE);


        }

        leg->AddEntry( hst, mnvMC->GetTitle(), "l" );
        hst->Scale( mcScale );
        if (i == 0) {
            //respect max/min setting the user may have used
            if ( MnvHist::IsAutoAxisLimit( axis_minimum ) )
            {
                if ( hist_min_zero && !gPad->GetLogy() )
                    hst->SetMinimum( 0. );
                else
                    hst->SetMinimum( footroom * minmin );
            }

            else
                hst->SetMinimum( axis_minimum );

            if ( MnvHist::IsAutoAxisLimit( axis_maximum ) )
            {
                if ( gPad->GetLogy() )
                    hst->SetMaximum( pow(headroom * maxmax,10.) );
                else
                    hst->SetMaximum( headroom * maxmax );
            }

            else
                hst->SetMaximum( axis_maximum );

            hst->Draw("HIST");
        }

        else
            hst->Draw("HIST SAME");
    }

    if ( dataAsPoints )
    {
        tmpData->Draw( "X0 E1 SAME" );
    }
    else
    {

        tmpData->Draw( "HIST SAME" );
    }

    if ( legPos != "N" )
        leg->Draw();

    gPad->Update();
}





//================================================================
// adds a legend to the plot
//================================================================
void MnvPlotter::AddPlotLegend(
        const std::vector< TH1* >       & hists,
        const std::vector< std::string >& names,
        const std::vector< std::string >& opts,
        const std::string& legPos /*="R"*/
        )
{
    double x1, y1, x2, y2;
    size_t legendWidth = GetLegendWidthInLetters( names );
    DecodeLegendPosition( x1, y1, x2, y2, legPos, names.size(), legendWidth );
    AddPlotLegend( hists, names, opts, x1, y1, x2-x1, y2-y1 );

}

void MnvPlotter::AddPlotLegend(
        const std::vector< TH1* >       & hists,
        const std::vector< std::string >& names,
        const std::vector< std::string >& opts ,
        const double x,
        const double y,
        const double x_width /*= 0.25*/,
        const double y_width /*= 0.15*/,
        const double textSize /*= 0.02*/
        )
{
    if (!gPad)
      throw std::runtime_error("MnvPlotter requires a TCanvas. Please make one first.");

    if ( hists.size() == 0 )
    {
        Error("MnvPlotter::AddPlotLegend","You didn't give any histograms to  AddPlotLegend.");
        return;
    }

    if ( hists.size() != names.size() )
    {
        Error("MnvPlotter::AddPlotLegend",Form( "You gave AddPlotLegend a different number of hists (%d) and names (%d)", (int)hists.size(), (int)names.size() ) );
        return;
    }

    if ( hists.size() != opts.size() )
    {
        Error("MnvPlotter::AddPlotLegend",Form( "You gave AddPlotLegend a different number of hists (%d) and opts (%d)", (int)hists.size(), (int)opts.size() ) );
        return;
    }


    TLegend *leg = new TLegend(x ,y ,x + x_width,y + y_width);
    leg->SetNColumns( legend_n_columns );
    static int legN = 0;
    leg->SetName( Form("add_leg_N = %d", legN++) );
    AddToTmp( leg );

    leg->SetBorderSize( legend_border_size );
    leg->SetFillColor( legend_fill_color );
    if ( legend_fill_color < 0 )
        leg->SetFillStyle(0);
    leg->SetTextSize( textSize );
    leg->SetTextFont( legend_text_font );
    for ( unsigned int i = 0; i < hists.size(); i++ )
        leg->AddEntry(hists[i],names[i].c_str(),opts[i].c_str());
    leg->Draw();
    gPad->Update();
}


//=================================================
// print a canvas in multiple formats
//=================================================
void MnvPlotter::MultiPrint(
        TCanvas *c,
        const std::string& name
        ) const
{
    if ( print_formats.size() == 0 )
    {
        Error("MultiPrint", "Cannot use default print formats because there are none!  Nothing will be printed.");
        return;
    }

    const string printName = name.size() ? name : std::string(c->GetName());
    string typeStr = print_formats.front();
    for ( vector<string>::const_iterator i = ++print_formats.begin(); i != print_formats.end(); ++i )
        typeStr += "," + *i;

    MultiPrint( c, printName, typeStr );
}

//supply a comma-separated list of formats you want to print
void MnvPlotter::MultiPrint(
        TCanvas *c,
        const std::string& name,
        const std::string& typeStr
        ) const
{
    std::vector<std::string> types;
    size_t i = 0;
    size_t j = typeStr.find(',');
    while( j!=std::string::npos ) {
        types.push_back( typeStr.substr(i,j-i) );
        i = ++j;
        j = typeStr.find(',',i);
    }
    if ( j == std::string::npos )
        types.push_back( typeStr.substr(i, typeStr.size()) );

    //we don't need an info statement here...
    const int oldVerbosity = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kWarning;

    for ( vector<string>::const_iterator itType = types.begin(); itType != types.end(); ++itType )
    {
        if ( print_topdir.empty() )
            c->SaveAs( Form("%s.%s", name.c_str(), itType->c_str()), itType->c_str() );
        else
            c->SaveAs( Form("%s/%s.%s", print_topdir.c_str(), name.c_str(), itType->c_str()), itType->c_str() );
    }

    gErrorIgnoreLevel = oldVerbosity;
}



//=========================================================
// get the mean of a histogram in some restricted range

double MnvPlotter::GetHistoMean(
        const TH1 * hist,
        const int minBin /*= -1*/,
        const int maxBin /*= -1*/
        ) const
{
    double err = 0.0;
    return GetHistoMean( hist, err, minBin, maxBin );
}


double MnvPlotter::GetHistoMean(
        const TH1 * hist,
        double & err,
        const int i_minBin /*= -1*/,
        const int i_maxBin /*= -1*/
        ) const
{
    double mean = 0.;
    err = 0.;

    int maxBin = i_maxBin, minBin = i_minBin;

    if ( i_minBin == -1 )
        minBin = 1;
    if ( i_maxBin == -1 )
        maxBin = hist->GetNbinsX();

    for ( int i = minBin; i < maxBin; i++ )
    {
        if ( hist->GetBinError(i) == 0.0 )
            continue;
        mean += hist->GetBinContent(i) / pow( hist->GetBinError(i), 2.0 );
        err += 1.0 / pow( hist->GetBinError(i), 2.0 );
    }

    if ( err > 0. )
    {
        mean /= err;
        err = sqrt( 1.0 / err );
    }
    else
    {
        err = 0.;
        mean = 0.;
    }

    return mean;
}


//====================================================
// Draw a hexagon
//====================================================
//all in mm
void MnvPlotter::DrawHex(
        const double apothem,
        const double xCenter /*= 0.0*/,
        const double yCenter /*= 0.0*/
        ) const
{
    /*
       Common Values
       fiducial volume    apothem = 850
       outter edge of ID: apothem = 1070   (width = 2140 )
       outter edge of OD: apothem = 1727.2 (width = 3454.4)

*/
    if (!gPad)
      throw std::runtime_error("MnvPlotter requires a TCanvas. Please make one first.");

    double w = apothem; // distance from center to nearest edge == apothem
    double h = apothem * 2.0 / sqrt(3); // distance from center to farthest point = length of a side = apothem / cos(pi/6)

    TLine line;
    line.SetLineWidth(2);
    line.DrawLine( xCenter    , yCenter + h    , xCenter + w, yCenter + h/2.0 );//top to upper right
    line.DrawLine( xCenter + w, yCenter + h/2.0, xCenter + w, yCenter - h/2.0 ); //upper right to lower right
    line.DrawLine( xCenter + w, yCenter - h/2.0, xCenter    , yCenter - h     ); //lower right to bottom
    line.DrawLine( xCenter    , yCenter - h    , xCenter - w, yCenter - h/2.0 );//bottom to lower left
    line.DrawLine( xCenter - w, yCenter - h/2.0, xCenter - w, yCenter + h/2.0 ); //lower left to upper left
    line.DrawLine( xCenter - w, yCenter + h/2.0, xCenter    , yCenter + h     ); //upper left to top

    gPad->Update();
}

//////////////////////////////////////////////
// Reverse the axis
// usually TObject is a histogram
//
//////////////////////////////////////////////
void MnvPlotter::ReverseXAxis( TH1 *h )
{
    if (!gPad)
      throw std::runtime_error("MnvPlotter requires a TCanvas. Please make one first.");

    if ( !h )
    {
        Warning("MnvPlotter::ReverseXAxis"," I can't reverse and axis if you pass me a NULL object.");
        return;
    }

    // Remove the current axis
    h->GetXaxis()->SetLabelOffset(999);
    h->GetXaxis()->SetTickLength(0);

    // Redraw the new axis
    gPad->Update();
    TGaxis *newaxis = new TGaxis(gPad->GetUxmax(),
            gPad->GetUymin(),
            gPad->GetUxmin(),
            gPad->GetUymin()-.0001,
            h->GetXaxis()->GetXmin(),
            h->GetXaxis()->GetXmax(),
            510,"+");
    AddToTmp( newaxis );


    newaxis->SetLabelOffset(0.03);
    newaxis->SetLabelSize(0.03);
    newaxis->Draw();
}

void MnvPlotter::ReverseXAxis( TH2 *h )
{
    if (!gPad)
      throw std::runtime_error("MnvPlotter requires a TCanvas. Please make one first.");

    if ( !h )
    {
        Warning("MnvPlotter::ReverseXAxis"," I can't reverse and axis if you pass me a NULL object.");
        return;
    }

    // Remove the current axis
    h->GetXaxis()->SetLabelOffset(999);
    h->GetXaxis()->SetTickLength(0);

    // Redraw the new axis
    gPad->Update();
    TGaxis *newaxis = new TGaxis(gPad->GetUxmax(),
            gPad->GetUymin(),
            gPad->GetUxmin(),
            gPad->GetUymin()-.0001,
            h->GetXaxis()->GetXmin(),
            h->GetXaxis()->GetXmax(),
            510,"+");
    AddToTmp( newaxis );


    newaxis->SetLabelOffset(0.03);
    newaxis->SetLabelSize(0.03);
    newaxis->Draw();
}

//===========================================
// Find the useful range in bins
//===========================================
bool MnvPlotter::GetBinRangeWithMinimumBinContent(
        const TH1 *h,
        int &lowbin,
        int &highbin,
        double minContentLow /*= 0.*/,
        double minContentHigh /*= 0.*/
        ) const
{
    lowbin = -1;
    highbin = h->GetNbinsX() + 2;
    for ( int ibin = 0; ibin <= h->GetNbinsX(); ++ibin )
    {
        if ( lowbin == -1 && h->GetBinContent(ibin) > minContentLow )
            lowbin = ibin;
        if ( h->GetBinContent(ibin) > minContentHigh )
            highbin = ibin;
    }
    return ( lowbin != -1 && highbin !=  h->GetNbinsX() + 2 );
}

bool MnvPlotter::GetNonZeroBinRange(
        const TH1 *h,
        int &lowbin,
        int &highbin
        ) const
{
    return GetBinRangeWithMinimumBinContent( h, lowbin, highbin, 0., 0. );
}


bool MnvPlotter::GetHistsRange( const TObjArray& hists, double& minmin, double& maxmax ) const
{
    maxmax = minmin = 0.;

    if ( hists.IsEmpty() )
    {
        Warning("MnvPlotter::GetHistsMax", "You gave me an empty array.  You get 0 in return.");
        return false;
    }

    const TH1 *firstHist = dynamic_cast<const TH1*>( hists.At(0) );
    if (!firstHist)
    {
        Error("MnvPlotter::GetHistsMax", "Objects in the array must be castable to TH1*");
        return false;
    }
    minmin = firstHist->GetBinContent( firstHist->GetMinimumBin() );
    maxmax = firstHist->GetBinContent( firstHist->GetMaximumBin() );


    for ( Int_t i = 1, nhists=hists.GetEntriesFast(); i != nhists; ++i )
    {
        const TH1 *hist = dynamic_cast<const TH1*>( hists.At(i) );
        if (!hist)
        {
            Error("MnvPlotter::GetHistsMax", "Objects in the array must be castable to TH1*");
            minmin = maxmax = 0.;
            return false;
        }
        minmin = std::min( minmin, hist->GetBinContent( hist->GetMinimumBin() ) );
        maxmax = std::min( maxmax, hist->GetBinContent( hist->GetMaximumBin() ) );
    }

    return true;
}



//=================================================
// Create bins suitable for plotting in log10 scale
//=================================================
bool MnvPlotter::SetLogBins(
        const int nbins,
        const double min,
        const double max,
        double *bins
        )
{
    if ( min == 0  || max == 0 ) {
        Error("SetLogBins","Min and Max bins must be positive.");
        return false;
    }
    const double lmin =  log10(min);
    const double lmax =  log10(max);
    const double width = (double)(lmax - lmin) / nbins;
    bins[0] = min;
    for ( int i=1; i <= nbins; ++i )
        bins[i] = min + pow(10, lmin + i*width);


    return true;
}

//=======================================
// Get histograms of errors using user defined error groups
//=======================================
std::vector<TH1*> MnvPlotter::GetSysErrorGroupHists(
        MnvH1D* h,
        const bool doFractional/* = false */,
        const bool covAreaNormalize/* = false*/,
        const double ignoreThreshold /* = 0.00001 */
        ) const
{
    //the return vector
    vector<TH1*> hists;

    if ( ! h )
    {
        Error("MnvPlotter::GetSysErrorGroupHists", "You passed me a NULL MnvH1D.  Nothing to do.");
        return hists;
    }


    map<string,TH1D*> errGroupHists;

    vector<string> vertNames = h->GetVertErrorBandNames();
    for ( unsigned int i = 0; i != vertNames.size(); ++i )
    {
        const MnvVertErrorBand *errBand = h->GetVertErrorBand(vertNames[i]);
        TH1D *hErr = (TH1D*)errBand->GetErrorBand( doFractional, covAreaNormalize ).Clone(Form("tmp_vert_error%d_%d", i, __LINE__) );
        hErr->SetTitle( vertNames[i].c_str() );

        //is this histogram part of a group?
        bool inGroup = false;
        for (ErrorSummaryGroupMap::const_iterator itGroup = error_summary_group_map.begin(); itGroup != error_summary_group_map.end(); ++itGroup )
        {
            const string& errName = itGroup->first;
            const vector<string>& histNames = itGroup->second;

            //if this histogram is in the group add it to the group
            if ( find( histNames.begin(), histNames.end(), vertNames[i]) != histNames.end() )
            {
                map<string,TH1D*>::iterator itGroupHist = errGroupHists.find(errName);

                //if this group has no histogram yet, add it
                //otherwise, add in quadrature
                if ( errGroupHists.end() == itGroupHist )
                    errGroupHists[errName] = hErr;
                else
                {
                    MnvHist::AddInQuadrature( itGroupHist->second, hErr );
                    delete hErr;
                }

                inGroup = true;
                break;
            }
        }
        //if we added to the group, move on
        if ( inGroup )
            continue;

        //delet and move on if it is too small
        if ( 0 < ignoreThreshold && hErr->GetBinContent( hErr->GetMaximumBin() ) < ignoreThreshold )
        {
            delete hErr;
            continue;
        }

        hists.push_back( hErr );
    }

    // Plot each of the fractional contributions from lateral errors
    vector<string> latNames = h->GetLatErrorBandNames();
    for ( unsigned int i = 0; i != latNames.size(); ++i )
    {
        const MnvLatErrorBand *errBand = h->GetLatErrorBand( latNames[i] );
        TH1D *hErr = (TH1D*)errBand->GetErrorBand( doFractional, covAreaNormalize ).Clone(Form("tmp_lat_error%d_%d", i, __LINE__) );
        hErr->SetTitle( latNames[i].c_str() );

        //is this histogram part of a group?
        bool inGroup = false;
        for (ErrorSummaryGroupMap::const_iterator itGroup = error_summary_group_map.begin(); itGroup != error_summary_group_map.end(); ++itGroup )
        {
            const string& errName = itGroup->first;
            const vector<string>& histNames = itGroup->second;

            //if this histogram is in the group add it to the group
            if ( find( histNames.begin(), histNames.end(), latNames[i]) != histNames.end() )
            {
                map<string,TH1D*>::iterator itGroupHist = errGroupHists.find(errName);

                //if this group has no histogram yet, add it
                //otherwise, add in quadrature
                if ( errGroupHists.end() == itGroupHist )
                    errGroupHists[ errName ] = hErr;
                else
                {
                    MnvHist::AddInQuadrature( itGroupHist->second, hErr );
                    delete hErr;
                }

                inGroup = true;
                break;
            }
        }
        //if we added to the group, move on
        if ( inGroup )
            continue;

        //delete and move on if too small
        if ( 0 < ignoreThreshold && hErr->GetBinContent( hErr->GetMaximumBin() ) < ignoreThreshold )
        {
            delete hErr;
            continue;
        }

        hists.push_back( hErr );
    }


    // Plot each of the fractional contributions from uncorrelated sys errors
    vector<string> uncorrNames = h->GetUncorrErrorNames();
    for ( unsigned int i = 0; i != uncorrNames.size(); ++i )
    {
        TH1D *hErr = (TH1D*)h->GetUncorrErrorAsHist( uncorrNames[i], doFractional ).Clone(Form("tmp_uncorr_error%d_%d", i, __LINE__) );
        hErr->SetTitle( uncorrNames[i].c_str() );

        //is this histogram part of a group?
        bool inGroup = false;
        for (ErrorSummaryGroupMap::const_iterator itGroup = error_summary_group_map.begin(); itGroup != error_summary_group_map.end(); ++itGroup )
        {
            const string& errName = itGroup->first;
            const vector<string>& histNames = itGroup->second;

            //if this histogram is in the group add it to the group
            if ( find( histNames.begin(), histNames.end(), uncorrNames[i]) != histNames.end() )
            {
                map<string,TH1D*>::iterator itGroupHist = errGroupHists.find(errName);

                //if this group has no histogram yet, add it
                //otherwise, add in quadrature
                if ( errGroupHists.end() == itGroupHist )
                    errGroupHists[ errName ] = hErr;
                else
                {
                    MnvHist::AddInQuadrature( itGroupHist->second, hErr );
                    delete hErr;
                }

                inGroup = true;
                break;
            }
        }
        //if we added to the group, move on
        if ( inGroup )
            continue;

        //delete and move on if too small
        if ( 0 < ignoreThreshold && hErr->GetBinContent( hErr->GetMaximumBin() ) < ignoreThreshold )
        {
            delete hErr;
            continue;
        }

        hists.push_back( hErr );
    }//end of uncorr errors

    //add error groups
    for ( map<string,TH1D*>::iterator itGroup = errGroupHists.begin(); itGroup != errGroupHists.end(); ++itGroup )
    {
        TH1* hist = itGroup->second;

        if ( 0 < ignoreThreshold && hist->GetBinContent( hist->GetMaximumBin() ) < ignoreThreshold )
        {
            //throw away an error group if it is too small, or it will leak
            delete hist;
            continue;
        }
        hist->SetTitle( itGroup->first.c_str() );

        hists.push_back( hist );
    }

    //YOU OWN ALL THESE OBJECTS!
    return hists;
}

//=======================================
// Draw a summary of MC only errors
//=======================================
bool MnvPlotter::DrawErrorSummary(
        MnvH1D* h,
        const std::string& legPos  /* = "TR"    */,
        const bool   includeStat     /* = true    */,
        const bool   solidLinesOnly  /* = true    */,
        const double ignoreThreshold /* = 0.00001 */,
        const bool covAreaNormalize/* = false*/,
        const std::string& errorGroupName /* = "" */,
        const bool  asfrac  /* = false */,
        const std::string &Ytitle,
        bool ignoreUngrouped,
        const std::string& histDrawOption
        )
{
    if ( ! h )
    {
        Error("DrawErrorSummary", "You passed me a NULL MnvH1D.  Nothing to do.");
        return false;
    }

    //set max digits to the default, since we almost never want scientific notation for errors.
    //restore the setting before returning.
    const int oldMaxDigits = TGaxis::GetMaxDigits();
    TGaxis::SetMaxDigits( axis_max_digits );




    // Store the pieces for a legend
    vector<TH1*>   hists;
    vector<string> names;
    vector<string> opts;

    bool useDifferentLineStyles = !solidLinesOnly;

    //! Get the total error and apply styles
    TH1D *hTotalErr = (TH1D*)h->GetTotalError( includeStat, asfrac, covAreaNormalize ).Clone( Form("h_total_err_errSum_%d", __LINE__) );
    AddToTmp( hTotalErr );
    ApplyAxisStyle(hTotalErr);

    //Error("DrawDataMCErrorSummary", Form("Total Err Max: %.2f",hTotalErr->GetMaximum() ) );

    ApplyNextLineStyle( hTotalErr, true, useDifferentLineStyles );

    //respect max/min setting the user may have used
    if ( MnvHist::IsAutoAxisLimit( axis_minimum ) )
        hTotalErr->SetMinimum( 0. );
    else
        hTotalErr->SetMinimum( axis_minimum );

    if ( MnvHist::IsAutoAxisLimit( axis_maximum ) )
        hTotalErr->SetMaximum( headroom * hTotalErr->GetMaximum() );
    else
        hTotalErr->SetMaximum( axis_maximum );

    hTotalErr->GetYaxis()->SetTitle( "Fractional Uncertainty" );
    if (!asfrac )hTotalErr->GetYaxis()->SetTitle( Ytitle.c_str() );
    if (errorGroupName == "") {
        if (!asfrac) hTotalErr->Scale( hTotalErr->GetBinWidth(1), "width" );
        hTotalErr->Draw( histDrawOption.c_str() );
        const string totalName = ( includeStat ? "Total Uncertainty" : "Total Sys. Uncertainty" );
        hists.push_back( hTotalErr );
        names.push_back( totalName );
        opts.push_back( "l" );
    }

    if ( includeStat && errorGroupName == "")
    {
        TH1D *statErr = (TH1D*)h->GetStatError(asfrac).Clone( Form("this_stat_err_%d", __LINE__) );
        AddToTmp( statErr );

        statErr->SetLineColor( 12 );//dark gray
        statErr->SetLineStyle( 2 ); //dashed
        statErr->SetLineWidth( 3 );
        statErr->Draw((histDrawOption + " SAME").c_str());
        hists.push_back( statErr );
        names.push_back( stat_error_name );
        opts.push_back( "l" );
    }

    TH1D *hTmpErr = (TH1D*)hTotalErr->Clone( Form("h_tmp_err_errSum_%d", __LINE__) );
    hTmpErr->Reset();
    map<string,TH1D*> errGroupHists;

    // plot each of the fractional contributions from all the errors.
    // first, we make a list of error bands to plot...
    bool drawn_already = (errorGroupName == "") ? true : false;
    vector<string> errNames = h->GetVertErrorBandNames();
    vector<string> otherNames = h->GetLatErrorBandNames();
    errNames.insert( errNames.end(), otherNames.begin(), otherNames.end() );
    otherNames = h->GetUncorrErrorNames();
    errNames.insert( errNames.end(), otherNames.begin(), otherNames.end() );
    for ( vector<string>::const_iterator it_name = errNames.begin();
            it_name != errNames.end();
            ++it_name)
    {
      TH1D * hErr = NULL;

      if (h->HasVertErrorBand(*it_name))
        hErr = dynamic_cast<TH1D*>(h->GetVertErrorBand(*it_name)->GetErrorBand( asfrac, covAreaNormalize ).Clone( Form("tmp_vertError_%s", (*it_name).c_str()) ));
      else if (h->HasLatErrorBand(*it_name))
        hErr = dynamic_cast<TH1D*>(h->GetLatErrorBand(*it_name)->GetErrorBand( asfrac, covAreaNormalize ).Clone( Form("tmp_latError_%s", (*it_name).c_str()) ));
      else if (h->HasUncorrError(*it_name))
        hErr = dynamic_cast<TH1D*>(h->GetUncorrErrorAsHist( *it_name, asfrac ).Clone(Form("tmp_uncorr_error_%s", (*it_name).c_str()) ));
      else
        throw std::runtime_error( Form("MnvPlotter::DrawErrorSummary(): Couldn't determine error band type for error name '%s'", (*it_name).c_str()) );

      //is this histogram part of a group?
      bool inGroup = false;
      for (ErrorSummaryGroupMap::const_iterator itGroup = error_summary_group_map.begin();
              itGroup != error_summary_group_map.end(); ++itGroup ) {
        const string& errName           = itGroup->first;
        const vector<string>& histNames = itGroup->second;

        //if this histogram is not in the group we're considering, skip to the next one
        if ( find( histNames.begin(), histNames.end(), *it_name) == histNames.end() )
          continue;
        //std::cout << " MnvPlotter found " << errName << " " << *it_name << std::endl;
        // if plotting the errors from only one group,
        // we don't want to "sub-group" them any further.
        // therefore we don't do any adding of histograms.
        if (errorGroupName==errName) {
          inGroup=true;
          break;
        }
        // otherwise, if no group was specifically chosen,
        // then (since this error band is already known to be in this group)
        // we need to add it into the histogram for the group
        // (or create that histogram if it doesn't exist yet).
        else if (errorGroupName == "") {
          map<string,TH1D*>::iterator itGroupHist = errGroupHists.find(errName);

          if ( errGroupHists.end() == itGroupHist ) {
            errGroupHists[ errName ] = hErr;
          }
          else {
            MnvHist::AddInQuadrature( itGroupHist->second, hErr );
            delete hErr;
          }

          inGroup = true;
          break;
        }
      }

      // if we haven't selected a group whose constituents we want to see...
      if ( errorGroupName=="" ) {
        // we never want to show the individual plots
        // when a histogram was included in a group:
        // the sums will be drawn later.
        if (inGroup)
          continue;

        // ... when we are ignoring ungrouped errors,
        // then NOTHING gets drawn here:
        // the grouped errors were added to the group histogram
        // above (and so we don't want them here),
        // and the ungrouped ones we're ignoring altogether.
        if (ignoreUngrouped)
          continue;
      }
      // if we DID select a group to draw the constituents of,
      // and this histogram isn't part of it,
      // then we also don't want to see it.
      else if ( errorGroupName != "" && !inGroup)
        continue;

      AddToTmp( hErr );

      // don't draw any errors that have no bins above threshold
      if ( 0 < ignoreThreshold && hErr->GetBinContent( hErr->GetMaximumBin() ) < ignoreThreshold )
        continue;

      ApplyNextLineStyle( hErr, false, useDifferentLineStyles);

      hErr->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());

      map<string,int>::const_iterator itErrCol = error_color_map.find( *it_name );
      if ( error_color_map.end() != itErrCol )
        hErr->SetLineColor( itErrCol->second );

      hErr->SetLineWidth( mc_line_width );
      if (drawn_already) {
        hErr->Draw( (histDrawOption + " SAME").c_str() );
      }
      else {
        drawn_already = true;
        hTmpErr->GetYaxis()->SetTitle( "Fractional Uncertainty" );
        if (!asfrac) hTmpErr->GetYaxis()->SetTitle( Ytitle.c_str() );
        if (!asfrac) hTmpErr->Scale(hTmpErr->GetBinWidth(1),"width");
        hTmpErr->Draw(histDrawOption.c_str());
	
        ApplyAxisStyle(hErr);
        ////respect max/min setting the user may have used
        if ( MnvHist::IsAutoAxisLimit( axis_minimum ) )
	  hErr->SetMinimum( 0. );
        else
	  hErr->SetMinimum( axis_minimum );
	
	//if (inGroup) {
	//hErr->SetMaximum( headroom * axis_maximum_group );
	//}
        if ( MnvHist::IsAutoAxisLimit( axis_maximum ) )
	  hErr->SetMaximum( headroom * hTotalErr->GetMaximum() );
        else
	  hErr->SetMaximum( axis_maximum);

        hErr->GetYaxis()->SetTitle( "Fractional Uncertainty" );
        if (!asfrac) hErr->GetYaxis()->SetTitle( Ytitle.c_str() );
        //if (!asfrac) hErr->Scale(hErr->GetBinWidth(1),"width");
        hErr->Draw((histDrawOption + " SAME").c_str());
      }
      hists.push_back( hErr );
      names.push_back( *it_name );
      opts.push_back( "l" );
    }

    //add error groups
    for ( map<string,TH1D*>::iterator itGroup = errGroupHists.begin(); itGroup != errGroupHists.end(); ++itGroup )
    {
        //   std::cout << "  (plot " << h->GetName() << ") drawing error group '" << itGroup->first << "'" << std::endl;
        const string& name = itGroup->first;
        TH1* hist = itGroup->second;

        if ( 0 < ignoreThreshold && hist->GetBinContent( hist->GetMaximumBin() ) < ignoreThreshold )
        {
            //     std::cout << "    (... ignored because its maximum (" << hist->GetBinContent( hist->GetMaximumBin() ) << ") was below threshold)" << std::endl;
            continue;
        }

        ApplyNextLineStyle( hist, false, useDifferentLineStyles);

        map<string,int>::const_iterator itErrCol = error_color_map.find( name );
        if ( error_color_map.end() != itErrCol )
            hist->SetLineColor( itErrCol->second );

        hist->SetLineWidth( mc_line_width );
        if (!asfrac)hist->Scale(hist->GetBinWidth(1),"width");
        hist->Draw( (histDrawOption + " SAME").c_str() );
        hists.push_back( hist );
        names.push_back( name );
        opts.push_back( "l" );
    }

    if ( legPos != "N" )
    {
        size_t legendWidth = GetLegendWidthInLetters( names );
        double x1,y1,x2,y2;
        DecodeLegendPosition( x1, y1, x2, y2, legPos, hists.size(), legendWidth );
        AddPlotLegend( hists, names, opts, x1, y1, x2-x1, y2-y1, legend_text_size );
    }

    gPad->RedrawAxis();
    gPad->Update();

    TGaxis::SetMaxDigits( oldMaxDigits );

    return true;
}

//==============================================================
// DrawStackedMC - vector of hists to turn into a stack
//==============================================================
void MnvPlotter::DrawStackedMC(
        const TObjArray* mcHists,
        const Double_t mcScale,
        const std::string& legPos,
        const Int_t mcBaseColor,               // Color of first stacked histo.
        const Int_t mcColorOffset,             // 2nd histo in stack has color mcBaseColor+mcColorOffset, etc.
        const Int_t mcFillStyle,               // Fill style for histograms only.
        const char* xaxislabel,
        const char* yaxislabel
        )
{
    if (!gPad)
      throw std::runtime_error("MnvPlotter requires a TCanvas. Please make one first.");

    //----------------------
    // start creating legend
    //----------------------

    unsigned int nHists = mcHists->GetEntries();
    if (nHists==0) return;

    //find the longest title in each column
    vector<string> titles;
    for ( unsigned int i = 0; i != nHists; ++i )
    {
        MnvH1D *mnv = (MnvH1D*)mcHists->At(i);
        titles.push_back( mnv->GetTitle() );
    }

    size_t legendWidth = GetLegendWidthInLetters( titles );
    double x1, x2, y1, y2;
    DecodeLegendPosition( x1, y1, x2, y2, legPos, nHists+1, legendWidth, legend_text_size );

    TLegend *leg  = new TLegend(x1, y1, x2, y2);
    leg->SetNColumns( legend_n_columns );
    AddToTmp(leg);

    leg->SetNColumns( legend_n_columns );
    leg->SetBorderSize( legend_border_size );
    leg->SetFillColor( legend_fill_color );
    if ( legend_fill_color < 0 )
        leg->SetFillStyle(0);
    leg->SetTextSize( legend_text_size );
    leg->SetTextFont( legend_text_font );

    THStack *hs  = new THStack("hs", "Stacked 1D histograms");
    AddToTmp(hs);

    TH1     *hst(NULL);
    MnvH1D  *mnvhst(NULL);

    int fillCol = mcBaseColor;
    std::string first_xaxis_label, first_yaxis_label;
    for ( unsigned int i = 0; i < nHists; i++ )
    {
        fillCol += i*mcColorOffset;
        mnvhst = (MnvH1D*)mcHists->At(i);

        if ( draw_normalized_to_bin_width )
            hst    = (TH1*)mnvhst->GetBinNormalizedCopy().GetCVHistoWithError().Clone( Form( "tmpMC_%04d_%d", i, __LINE__) );
        else
            hst    = (TH1*)mnvhst->GetCVHistoWithError().Clone( Form( "tmpMC_%04d_%d", i, __LINE__) );
        AddToTmp( hst );

        if ( mcBaseColor > 0 )
            hst->SetFillColor( fillCol );
        if ( mcFillStyle > 0 )
            hst->SetFillStyle( mcFillStyle );
        hst->SetLineWidth(mc_line_width);
        hst->Scale(mcScale);
        hs->Add( hst );

        if (i==0) {
            first_xaxis_label = hst->GetXaxis()->GetTitle();
            first_yaxis_label = hst->GetYaxis()->GetTitle();
        }
    }

    for ( unsigned int i = 0; i != nHists; ++i )
    {

        mnvhst = (MnvH1D*)mcHists->At(nHists-1-i);
        hst    = (TH1*)hs->GetHists()->At(nHists-1-i);
        leg->AddEntry( hst, mnvhst->GetTitle(), "f" );
    }

    //HACK -need to draw a smaller histogram first to set the x-axis range correctly
    //for stacks with variable-width bins
    MnvH1D* tmp_mnv = (MnvH1D*)hs->GetHists()->At(0);
    tmp_mnv->GetXaxis()->SetTitle( xaxislabel );
    tmp_mnv->GetYaxis()->SetTitle( yaxislabel );
    tmp_mnv->GetXaxis()->SetTitleFont(axis_title_font_x);
    tmp_mnv->GetYaxis()->SetTitleFont(axis_title_font_y);
    tmp_mnv->GetXaxis()->SetTitleSize(axis_title_size_x);
    tmp_mnv->GetYaxis()->SetTitleSize(axis_title_size_y);
    tmp_mnv->GetXaxis()->SetLabelFont(axis_label_font);
    tmp_mnv->GetYaxis()->SetLabelFont(axis_label_font);
    tmp_mnv->GetXaxis()->SetLabelSize(axis_label_size);
    tmp_mnv->GetYaxis()->SetLabelSize(axis_label_size);
    tmp_mnv->GetXaxis()->CenterTitle(kTRUE);
    if ( strlen(xaxislabel) == 0 )
        tmp_mnv->GetXaxis()->SetTitle( first_xaxis_label.c_str() );
    if ( strlen(yaxislabel) == 0 )
        tmp_mnv->GetYaxis()->SetTitle( first_yaxis_label.c_str() );

    //respect max/min setting the user may have used
    if ( MnvHist::IsAutoAxisLimit( axis_minimum ) )
    {
        if ( hist_min_zero && !gPad->GetLogy() )
            tmp_mnv->SetMinimum( 0. );
        else
            tmp_mnv->SetMinimum( footroom * hs->GetMinimum() );
    }
    else
        tmp_mnv->SetMinimum( axis_minimum );

    if ( MnvHist::IsAutoAxisLimit( axis_maximum ) )
        tmp_mnv->SetMaximum( headroom * hs->GetMaximum() );
    else
        tmp_mnv->SetMaximum( axis_maximum );


    tmp_mnv->Draw("HIST");
    hs->Draw("SAME HIST");
    hs->Draw("SAME AXIS");

    if ( legPos != "N" )
        leg->Draw();

    gPad->Update();
}




//==============================================================
// DrawDataStackedMC - vector of hists to turn into a stack & "data"
//==============================================================

void MnvPlotter::DrawDataStackedMC(
        const MnvH1D* dataHist,
        const TObjArray* mcHists,
        const Double_t mcScale,
        const std::string& legPos,
        const std::string& dataName,
        const Int_t mcBaseColor,               // Color of first stacked histo.
        const Int_t mcColorOffset,             // 2nd histo in stack has color mcBaseColor+mcColorOffset, etc.
        const Int_t mcFillStyle,               // Fill style for histograms only.
        const char* xaxislabel,
        const char* yaxislabel,
        bool cov_area_normalize
        )
{
    if (!gPad)
      throw std::runtime_error("MnvPlotter requires a TCanvas. Please make one first.");

    //----------------------
    // start creating legend
    //----------------------

    unsigned int nHists = mcHists->GetEntries();
    vector<string> titles;
    titles.push_back( dataName );
    for ( unsigned int i = 0; i != nHists; ++i )
    {
        MnvH1D *mnv = (MnvH1D*)mcHists->At(i);
        // Let's just make double sure you didn't pass an array of MnvH2D.
        // (Because if you did you'll get an indecipherable seg fault)
        if (std::string(typeid(*mnv).name()).find(std::string("MnvH1D"))
                == std::string::npos) {
          throw std::runtime_error("The TObjArray you passed to DrawDataStackedMC does not contain all MnvH1Ds!");
        }
        titles.push_back( mnv->GetTitle() );
    }

    size_t legendWidth = GetLegendWidthInLetters( titles );
    double x1, x2, y1, y2;
    DecodeLegendPosition( x1, y1, x2, y2, legPos, nHists+1, legendWidth, legend_text_size );

    TLegend *leg  = new TLegend(x1, y1, x2, y2);
    leg->SetNColumns( legend_n_columns );
    AddToTmp(leg);

    leg->SetBorderSize( legend_border_size );
    leg->SetFillColor( legend_fill_color );
    if ( legend_fill_color < 0 )
        leg->SetFillStyle(0);
    leg->SetTextSize( legend_text_size );
    leg->SetTextFont( legend_text_font );

    THStack *hs  = new THStack("hs", "Stacked 1D histograms");
    AddToTmp(hs);

    TH1     *hst(NULL);
    MnvH1D  *mnvhst(NULL);

    // Note, Sumw2 is enforced automatically for MnvH1D's.
    TH1* tmpData(0);
    if ( draw_normalized_to_bin_width ) {
      tmpData = (TH1*)dataHist->GetCVHistoWithError(true, cov_area_normalize).Clone( Form("tmp_data_%d", __LINE__) );
      if (dataHist->GetNormBinWidth() > 0 ) tmpData->Scale( dataHist->GetNormBinWidth(), "width" );
    }
    else
      tmpData = (TH1*)dataHist->GetCVHistoWithError().Clone( Form("tmp_data_%d", __LINE__) );

    AddToTmp( tmpData );

    //style data
    tmpData->SetMarkerStyle(data_marker);
    tmpData->SetMarkerSize(data_marker_size);
    tmpData->SetLineWidth(data_line_width);
    tmpData->SetLineStyle(data_line_style);
    tmpData->SetLineColor(data_color);

    leg->AddEntry(tmpData,dataName.c_str(),"ple");
    int fillCol = mcBaseColor;
    for ( unsigned int i = 0; i < nHists; i++ )
    {
      fillCol += i*mcColorOffset;
      mnvhst = (MnvH1D*)mcHists->At(i);

      if ( draw_normalized_to_bin_width )
        hst    = (TH1*)mnvhst->GetBinNormalizedCopy().GetCVHistoWithError().Clone( Form( "tmpMC_%04d_%d", i, __LINE__) );
      else
        hst    = (TH1*)mnvhst->GetCVHistoWithError().Clone( Form( "tmpMC_%04d_%d", i, __LINE__) );
      AddToTmp( hst );

      if ( mcBaseColor > 0 )
        hst->SetFillColor( fillCol );
      if ( mcFillStyle > 0 )
        hst->SetFillStyle( mcFillStyle );
      hst->SetLineWidth(mc_line_width);
      hst->Scale(mcScale);
      hs->Add( hst );
    }

    //add the legend in reverse order so vertical alignment is same as stack
    for ( unsigned int i = 0; i != nHists; ++i )
    {
        mnvhst = (MnvH1D*)mcHists->At(nHists-1-i);
        hst    = (TH1*)hs->GetHists()->At(nHists-1-i);
        leg->AddEntry( hst, mnvhst->GetTitle(), "f" );
    }

    //HACK -need to draw a smaller histogram first to set the x-axis range correctly
    //for stacks with variable-width bins
    MnvH1D* tmp_mnv = (MnvH1D*)hs->GetHists()->At(0);
    tmp_mnv->GetXaxis()->SetRange(  tmpData->GetXaxis()->GetFirst(), tmpData->GetXaxis()->GetLast() );
    tmp_mnv->GetXaxis()->SetTitle( xaxislabel );
    tmp_mnv->GetYaxis()->SetTitle( yaxislabel );
    tmp_mnv->GetXaxis()->SetTitleFont(axis_title_font_x);
    tmp_mnv->GetYaxis()->SetTitleFont(axis_title_font_y);
    tmp_mnv->GetXaxis()->SetTitleSize(axis_title_size_x);
    tmp_mnv->GetYaxis()->SetTitleSize(axis_title_size_y);
    tmp_mnv->GetXaxis()->SetTitleOffset(axis_title_offset_x);
    tmp_mnv->GetYaxis()->SetTitleOffset(axis_title_offset_y);
    tmp_mnv->GetXaxis()->SetLabelFont(axis_label_font);
    tmp_mnv->GetYaxis()->SetLabelFont(axis_label_font);
    tmp_mnv->GetXaxis()->SetLabelSize(axis_label_size);
    tmp_mnv->GetYaxis()->SetLabelSize(axis_label_size);
    tmp_mnv->GetXaxis()->CenterTitle(kTRUE);
    if ( strlen(xaxislabel) == 0 )
        tmp_mnv->GetXaxis()->SetTitle( dataHist->GetXaxis()->GetTitle() );
    if ( strlen(yaxislabel) == 0 )
        tmp_mnv->GetYaxis()->SetTitle( dataHist->GetYaxis()->GetTitle() );


    //respect max/min setting the user may have used
    if ( MnvHist::IsAutoAxisLimit( axis_minimum ) )
    {
      if (!gPad)
        throw std::runtime_error("MnvPlotter requires a TCanvas. Please make one first.");

      if ( hist_min_zero && !gPad->GetLogy() )
        tmp_mnv->SetMinimum( 0. );
      else
        tmp_mnv->SetMinimum( footroom * std::min( tmpData->GetMinimum(), hs->GetMinimum() ) );
    }
    else
      tmp_mnv->SetMinimum( axis_minimum );

    if ( MnvHist::IsAutoAxisLimit( axis_maximum ) )
        tmp_mnv->SetMaximum( headroom * std::max( tmpData->GetMaximum(), hs->GetMaximum() ) );
    else
        tmp_mnv->SetMaximum( axis_maximum );



    tmp_mnv->Draw("HIST");
    hs->Draw("HIST SAME");


    tmpData->SetMarkerStyle(data_marker);
    tmpData->SetMarkerSize(data_marker_size);
    tmpData->SetLineWidth(data_line_width);
    tmpData->SetLineStyle(data_line_style);
    tmpData->SetLineColor(data_color);
    tmpData->Draw("SAME E1 X0");
    tmp_mnv->Draw("AXIS SAME");

    if ( legPos != "N" )
        leg->Draw();



    gPad->Update();
}

//==============================================================
// DrawDataStackedMC - vector of hists to turn into a stack & "data"
//==============================================================

void MnvPlotter::DrawDataStackedMC(
        const MnvH1D* dataHist,
        const TObjArray* mcHists,
        const Int_t* mcColors,
        const Double_t mcScale,
        const std::string& legPos,
        const std::string& dataName,
        const Int_t mcFillStyle,               // Fill style for histograms only.
        const char* xaxislabel,
        const char* yaxislabel,
        bool cov_area_normalize
        )
{
    if (!gPad)
      throw std::runtime_error("MnvPlotter requires a TCanvas. Please make one first.");

    //----------------------
    // start creating legend
    //----------------------

    unsigned int nHists = mcHists->GetEntries();
    vector<string> titles;
    titles.push_back( dataName );
    for ( unsigned int i = 0; i != nHists; ++i )
    {
        MnvH1D *mnv = (MnvH1D*)mcHists->At(i);
        // Let's just make double sure you didn't pass an array of MnvH2D.
        // (Because if you did you'll get an indecipherable seg fault)
        if (std::string(typeid(*mnv).name()).find(std::string("MnvH1D"))
                == std::string::npos) {
          throw std::runtime_error("The TObjArray you passed to DrawDataStackedMC does not contain all MnvH1Ds!");
        }
        titles.push_back( mnv->GetTitle() );
    }

    size_t legendWidth = GetLegendWidthInLetters( titles );
    double x1, x2, y1, y2;
    DecodeLegendPosition( x1, y1, x2, y2, legPos, nHists+1, legendWidth, legend_text_size );

    TLegend *leg  = new TLegend(x1, y1, x2, y2);
    leg->SetNColumns( legend_n_columns );
    AddToTmp(leg);

    leg->SetBorderSize( legend_border_size );
    leg->SetFillColor( legend_fill_color );
    if ( legend_fill_color < 0 )
        leg->SetFillStyle(0);
    leg->SetTextSize( legend_text_size );
    leg->SetTextFont( legend_text_font );

    THStack *hs  = new THStack("hs", "Stacked 1D histograms");
    AddToTmp(hs);

    TH1     *hst(NULL);
    MnvH1D  *mnvhst(NULL);

    // Note, Sumw2 is enforced automatically for MnvH1D's.
    TH1* tmpData(0);
    if ( draw_normalized_to_bin_width ) {
      tmpData = (TH1*)dataHist->GetCVHistoWithError(true, cov_area_normalize).Clone( Form("tmp_data_%d", __LINE__) );
      if (dataHist->GetNormBinWidth() > 0 ) tmpData->Scale( dataHist->GetNormBinWidth(), "width" );
    }
    else
      tmpData = (TH1*)dataHist->GetCVHistoWithError().Clone( Form("tmp_data_%d", __LINE__) );

    AddToTmp( tmpData );

    //style data
    tmpData->SetMarkerStyle(data_marker);
    tmpData->SetMarkerSize(data_marker_size);
    tmpData->SetLineWidth(data_line_width);
    tmpData->SetLineStyle(data_line_style);
    tmpData->SetLineColor(data_color);

    leg->AddEntry(tmpData,dataName.c_str(),"ple");
    for ( unsigned int i = 0; i < nHists; i++ )
    {
      mnvhst = (MnvH1D*)mcHists->At(i);

      if ( draw_normalized_to_bin_width )
        hst    = (TH1*)mnvhst->GetBinNormalizedCopy().GetCVHistoWithError().Clone( Form( "tmpMC_%04d_%d", i, __LINE__) );
      else
        hst    = (TH1*)mnvhst->GetCVHistoWithError().Clone( Form( "tmpMC_%04d_%d", i, __LINE__) );
      AddToTmp( hst );

      if ( mcColors )
        hst->SetFillColor( mcColors[i] );
      if ( mcFillStyle > 0 )
        hst->SetFillStyle( mcFillStyle );
      hst->SetLineWidth(mc_line_width);
      hst->Scale(mcScale);
      hs->Add( hst );
    }

    //add the legend in reverse order so vertical alignment is same as stack
    for ( unsigned int i = 0; i != nHists; ++i )
    {
        mnvhst = (MnvH1D*)mcHists->At(nHists-1-i);
        hst    = (TH1*)hs->GetHists()->At(nHists-1-i);
        leg->AddEntry( hst, mnvhst->GetTitle(), "f" );
    }

    //HACK -need to draw a smaller histogram first to set the x-axis range correctly
    //for stacks with variable-width bins
    MnvH1D* tmp_mnv = (MnvH1D*)hs->GetHists()->At(0);
    tmp_mnv->GetXaxis()->SetRange(  tmpData->GetXaxis()->GetFirst(), tmpData->GetXaxis()->GetLast() );
    tmp_mnv->GetXaxis()->SetTitle( xaxislabel );
    tmp_mnv->GetYaxis()->SetTitle( yaxislabel );
    tmp_mnv->GetXaxis()->SetTitleFont(axis_title_font_x);
    tmp_mnv->GetYaxis()->SetTitleFont(axis_title_font_y);
    tmp_mnv->GetXaxis()->SetTitleSize(axis_title_size_x);
    tmp_mnv->GetYaxis()->SetTitleSize(axis_title_size_y);
    tmp_mnv->GetXaxis()->SetTitleOffset(axis_title_offset_x);
    tmp_mnv->GetYaxis()->SetTitleOffset(axis_title_offset_x);
    tmp_mnv->GetXaxis()->SetLabelFont(axis_label_font);
    tmp_mnv->GetYaxis()->SetLabelFont(axis_label_font);
    tmp_mnv->GetXaxis()->SetLabelSize(axis_label_size);
    tmp_mnv->GetYaxis()->SetLabelSize(axis_label_size);
    tmp_mnv->GetXaxis()->CenterTitle(kTRUE);
    if ( strlen(xaxislabel) == 0 )
        tmp_mnv->GetXaxis()->SetTitle( dataHist->GetXaxis()->GetTitle() );
    if ( strlen(yaxislabel) == 0 )
        tmp_mnv->GetYaxis()->SetTitle( dataHist->GetYaxis()->GetTitle() );


    //respect max/min setting the user may have used
    if ( MnvHist::IsAutoAxisLimit( axis_minimum ) )
    {
      if (!gPad)
        throw std::runtime_error("MnvPlotter requires a TCanvas. Please make one first.");

      if ( hist_min_zero && !gPad->GetLogy() )
        tmp_mnv->SetMinimum( 0. );
      else
        tmp_mnv->SetMinimum( footroom * std::min( tmpData->GetMinimum(), hs->GetMinimum() ) );
    }
    else
      tmp_mnv->SetMinimum( axis_minimum );

    if ( MnvHist::IsAutoAxisLimit( axis_maximum ) )
        tmp_mnv->SetMaximum( headroom * std::max( tmpData->GetMaximum(), hs->GetMaximum() ) );
    else
        tmp_mnv->SetMaximum( axis_maximum );



    tmp_mnv->Draw("HIST");
    hs->Draw("HIST SAME");


    tmpData->SetMarkerStyle(data_marker);
    tmpData->SetMarkerSize(data_marker_size);
    tmpData->SetLineWidth(data_line_width);
    tmpData->SetLineStyle(data_line_style);
    tmpData->SetLineColor(data_color);
    tmpData->Draw("SAME E1 X0");
    tmp_mnv->Draw("AXIS SAME");

    if ( legPos != "N" )
        leg->Draw();



    gPad->Update();
}


//==============================================================
// DrawDataStackedMCWithErrorBand - vector of hists to turn into a stack & "data"
//==============================================================
void MnvPlotter::DrawDataStackedMCWithErrorBand( )
{
    std::cout << "Implement me by editing DrawDataStackedMC()!" << std::endl;
}

void MnvPlotter::DrawNormalizedMigrationHistogram(
        const TH2D* h_migration,
        const bool drawAsMatrix,
        const bool coarseContours, /* = false */
        const bool includeFlows, /* = true */
        const bool noText /* = false */
        )
{

    int first_bin = includeFlows ? 0 : 1;
    int last_bin = includeFlows ? h_migration->GetNbinsX()+1 : h_migration->GetNbinsX();
    Int_t nbins = includeFlows ? h_migration->GetNbinsX()+2 : h_migration->GetNbinsX();

    TMatrixD m_migration(nbins, nbins);
    TH2D tmp(*h_migration);
    tmp.Reset();
    for (int y = first_bin; y <= last_bin; ++y) {
        Double_t norm = 0.;
        for (int x = first_bin; x <= last_bin; ++x)
            norm += h_migration->GetBinContent(x,y);

        if ( fabs(norm) > 1E-8) {
            for (int x = first_bin; x <= last_bin; ++x) {
                double percentage =  100 * h_migration->GetBinContent(x,y) / norm;
                if (includeFlows) {
                    m_migration[y][x] = percentage; //yeah that's right  y/x
                }else{
                    m_migration[y-1][x-1] = percentage; //yeah that's right  y/x
                }
                tmp.SetBinContent( x, y, percentage);
            }
        }
    }

    if ( drawAsMatrix ) {
      tmp = TH2D( m_migration );
      tmp.GetXaxis()->SetTitle( "Reco Bins" );
      tmp.GetYaxis()->SetTitle( "True Bins" );
    }

    tmp.GetXaxis()->SetTitleOffset( axis_title_offset_x );
    tmp.GetYaxis()->SetTitleOffset( axis_title_offset_y );
    tmp.GetZaxis()->SetTitleOffset( axis_title_offset_z );
    tmp.GetZaxis()->SetTitle( "Fraction of Row in Cell" );

    if ( coarseContours ) {
        //set a low,ok,high,too high color scale
        double contours[20] = {
            0, 0.00001, 0.0002, 0.0003, 0.0004,
            25, 25.00001, 25.0002, 25.0003, 25.0004,
            50, 50.00001, 50.0002, 50.0003, 50.0004,
            75, 75.00001, 75.0002, 75.0003, 75.0004};
        tmp.SetContour( 20, contours );
    }

    gStyle->SetPaintTextFormat("2.0f");
    tmp.SetMarkerSize(2);
    if (noText) {
        tmp.DrawCopy("colz");
    } else {
        tmp.DrawCopy("colz text");
    }
}

void MnvPlotter::DrawAllUniverses(
        const MnvH1D *h,
        const bool covAreaNormalize, /* = false */
        const bool binWidthNormalize /* = true */
        )
{
    if (!gPad)
      throw std::runtime_error("MnvPlotter requires a TCanvas. Please make one first.");

    //! Figure out how to divide the canvas to fit all plots
    const int nplots   = h->GetErrorBandNames().size();
    const int nrows    = (int)( sqrt( nplots ) );
    int ncolumns = (int)( nplots/nrows );
    if ( nplots%nrows!=0 )
        ++ncolumns;

    //have to cast gPad to TVirtualPad to divide it (not sure why)
    TVirtualPad *thePad = gPad;
    thePad->Divide( ncolumns, nrows );

    //! Plot all vertical error bands
    int nplot = 0;
    vector<string> names;
    names = h->GetVertErrorBandNames();
    for ( vector<string>::iterator name=names.begin(); name!=names.end(); ++name ) {
        const MnvVertErrorBand *error = h->GetVertErrorBand( *name );
        ++nplot;
        thePad->cd( nplot );
        //! copy the MnvH1D axes attributes to the error band, so they are propagates to the draw
        TAxis mx = *error->GetXaxis();
        TAxis my = *error->GetYaxis();
        h->GetXaxis()->Copy( mx );
        h->GetYaxis()->Copy( my );
        error->DrawAll( "HIST", true/*draw CV*/, covAreaNormalize, binWidthNormalize ? h->GetNormBinWidth() : 0.0 );//0.0 means do not bin width normalize
        AddHistoTitle( (*name).c_str() );
    }

    //! Plot all vertical error bands
    names = h->GetLatErrorBandNames();
    for ( vector<string>::iterator name=names.begin(); name!=names.end(); ++name ) {
        const MnvLatErrorBand *error = h->GetLatErrorBand( *name );
        ++nplot;
        thePad->cd(nplot);
        //! copy the MnvH1D axes attributes to the error band, so they are propagates to the draw
        TAxis mx = *error->GetXaxis();
        TAxis my = *error->GetYaxis();

        h->GetXaxis()->Copy( mx);
        h->GetYaxis()->Copy( my);
        error->DrawAll( "HIST", true/*draw CV*/, covAreaNormalize, binWidthNormalize ? h->GetNormBinWidth() : 0.0 );//0.0 means do not bin width normalize
        AddHistoTitle( (*name).c_str() );
    }

}

void MnvPlotter::DrawErrorMatrices(
        const MnvH1D *h,
        const bool area_norm,
        const bool asCorr,
        const bool asFrac
        )
{

    if (!gPad)
      throw std::runtime_error("MnvPlotter requires a TCanvas. Please make one first.");

    // Figure out how to divide the canvas to fit all plots
    // one plot for each systematic, one for stat, one for total
    const int nplots   = h->GetErrorBandNames().size() + 2;
    const int nrows    = (int)( sqrt( nplots ) );
    int ncolumns = (int)( nplots/nrows );
    if ( nplots%nrows!=0 )
        ++ncolumns;

    //have to cast gPad to TVirtualPad to divide it (not sure why)
    TVirtualPad *thePad = gPad;
    thePad->Divide( ncolumns, nrows );

    // draw systematic error matrices
    int nplot = 0;
    vector<string> names;
    names = h->GetSysErrorMatricesNames();
    for ( vector<string>::iterator name=names.begin(); name!=names.end(); ++name )
    {

        ++nplot;
        thePad->cd(nplot);

        TMatrixD sysMatrix = asCorr ? h->GetSysCorrelationMatrix( *name, area_norm ) : h->GetSysErrorMatrix( *name, asFrac, area_norm );
        DrawErrorMatrix( sysMatrix, h->GetXaxis() );

        if ( asCorr )
        {
            this->AddHistoTitle( Form( "%s Cor. Matrix", (*name).c_str() ) );
        }
        else{
            this->AddHistoTitle( Form( "%s Cov. Matrix", (*name).c_str() ) );
        }

    }

    string name = "";

    // draw the statistical error matrix
    // no correlation for stat. matrix
    if ( ! asCorr )
    {
        name = "Stat.";
        ++nplot;
        thePad->cd(nplot);
        TMatrixD statMatrix =  h->GetStatErrorMatrix( asFrac );
        DrawErrorMatrix( statMatrix , h->GetXaxis() );
        this->AddHistoTitle( Form( "%s Cov. Matrix", name.c_str() ) );
    }

    // draw the total error matrix
    name = "Total";
    ++nplot;
    thePad->cd(nplot);
    TMatrixD totalMatrix = asCorr ? h->GetTotalCorrelationMatrix( area_norm ) : h->GetTotalErrorMatrix( true, asFrac, area_norm );
    DrawErrorMatrix( totalMatrix, h->GetXaxis() );
    if ( asCorr )
        this->AddHistoTitle( Form( "%s Cor. Matrix", name.c_str() ) );
    else
        this->AddHistoTitle( Form( "%s Cov. Matrix", name.c_str() ) );
}

void MnvPlotter::DrawErrorMatrix(
        const TMatrixD &matrix,
        const TAxis* axis,
        const double maximum /*= -1.0*/ )
{

    // create a 2D histogram with the matrix elements
    TH2D *h2D = new TH2D( "h_matrix",
            Form( "matrix;%s;%s", axis->GetTitle(), axis->GetTitle()),
            axis->GetNbins(), axis->GetXbins()->GetArray(),
            axis->GetNbins(), axis->GetXbins()->GetArray() );

    for ( int i = 0; i < matrix.GetNrows(); i++ )
    {
        for ( int j = 0; j < matrix.GetNcols(); j++ )
        {
            //      int sign = matrix[i][j]>0 ? 1 : -1;
            h2D->SetBinContent( h2D->GetBin(i,j), matrix[i][j] );
        }
    }

    h2D->GetXaxis()->SetRange( axis->GetFirst(), axis->GetLast() );
    h2D->GetYaxis()->SetRange( axis->GetFirst(), axis->GetLast() );
    h2D->SetTitleSize( axis->GetTitleSize(), "XY");
    h2D->GetXaxis()->SetNdivisions( axis->GetNdivisions() );
    h2D->GetYaxis()->SetNdivisions( axis->GetNdivisions() );

    // set minimum and maximum by type
    double min = h2D->GetBinContent( h2D->GetMinimumBin() );
    double max = h2D->GetBinContent( h2D->GetMaximumBin() );

    // if maximum is supplied, use it, and assume symmetric
    if ( maximum > 0.0 ) {
        max = maximum;
        min = -maximum;
    }

    // for correlation matrix fix maximum to 1.
    // this test usually means this is a correlation matrix.
    if ( draw_corr_max1 && 0 <= min && min < 1. && 0 < max && max < 1. )
    {
        h2D->SetMinimum( floor( min*10. )/10. );
        h2D->SetMaximum( 1.0 );
    }
    else
    {

        // separate minimum and maximum if they are very close
        if ( (max-min)<0.1 )
        {
            int middle = floor( 0.5+(min+max)/2.0 );
            h2D->SetMinimum( middle-0.5 );
            h2D->SetMaximum( middle+0.5 );
        }

        // if minimum is negative and maximum is positive
        // center z axis in zero and make limits symetric
        else if ( min*max<0 )
        {
            const double absmax = std::max( fabs(min), fabs(max) );
            h2D->SetMinimum( -absmax );
            h2D->SetMaximum(  absmax );
        }

        //
        // I don't understand why one would want to always separate by at least 1,
        // but I don't want to force the change.
        // I think the symmetric min/max is what we normally would want
        if ( draw_corr_red_blue )
        {
            const double absmax = std::max( fabs(min), fabs(max) );
            h2D->SetMinimum( -absmax );
            h2D->SetMaximum(  absmax );

        }

    }

    // draw the 2d histogram
    if ( draw_corr_red_blue )
        SetCorrelationPalette();

    // make the number of sigfigs on the z axis be the same, why this is not the default I do not know
    h2D->GetZaxis()->SetDecimals(true);

    h2D->DrawCopy( "colz" );

    //if ( draw_corr_red_blue )
    //  gStyle->SetPalette( palette_style );

    // clean up
    delete h2D;

}

void MnvPlotter::DrawDoubleGausFit(
        const  TH1D *h,
        double lowFitBound,
        double highFitBound,
        const char* legPos,
        double* parameters,
        double* errors,
        double  chisquared,
        int     ndf )
{

    TH1D *htmp = (TH1D*)h->Clone( "htmp" );

    // if the boundaries are not set, use the histogram to get the initial guess parameters and bin edges
    double lowbin  = lowFitBound;
    double highbin = highFitBound;

    if ( lowFitBound == highFitBound )
    {
        lowbin  = htmp->GetBinLowEdge( htmp->GetMaximumBin() - 4 );
        highbin = htmp->GetBinLowEdge( htmp->GetMaximumBin() + 4 );
    }

    // get an initial guess at the fit parameters
    TF1* tmp_fit = new TF1("tmp_fit","gaus",lowbin,highbin);
    htmp->Fit(tmp_fit,"QR");

    double tmp_pars[3];
    tmp_fit->GetParameters(tmp_pars);

    double minbin    = lowFitBound == highFitBound ? htmp->GetXaxis()->GetXmin() : lowbin;
    double maxbin    = lowFitBound == highFitBound ? htmp->GetXaxis()->GetXmax() : highbin;
    double maxvalue  = htmp->GetMaximum();
    double rmsvalue  = tmp_pars[2];
    double meanvalue = tmp_pars[1];
    delete tmp_fit;

    // fit function
    TF1* fit = new TF1("fit","gaus(0)+gaus(3)",minbin,maxbin);
    fit->SetLineColor(data_color);

    // set initial guess parameters
    fit->SetParameters(maxvalue,meanvalue,rmsvalue,maxvalue/4.0,meanvalue,2*rmsvalue);

    // fit the histogram
    htmp->Fit(fit,"QR");

    // get the parameters
    double pars[6];
    fit->GetParameters(pars);

    if ( parameters != (double*)NULL )
    {
        for (int i = 0; i < 6; ++i)
            parameters[i] = pars[i];
    }

    const double* errs = fit->GetParErrors();

    if ( errors != (double*)NULL )
    {
        for (int i = 0; i < 6; ++i)
            errors[i] = fabs(errs[i]);
    }

    chisquared = fit->GetChisquare();
    ndf = fit->GetNDF();

    // set attributes
    if ( MnvHist::IsAutoAxisLimit( axis_maximum ) )
        htmp->SetMaximum( headroom * htmp->GetMaximum() );
    else
        htmp->SetMaximum( axis_maximum );

    htmp->GetXaxis()->SetTitleFont(axis_title_font_x);
    htmp->GetYaxis()->SetTitleFont(axis_title_font_y);
    htmp->GetXaxis()->SetTitleSize(axis_title_size_x);
    htmp->GetYaxis()->SetTitleSize(axis_title_size_y);
    htmp->GetXaxis()->SetLabelFont(axis_label_font);
    htmp->GetYaxis()->SetLabelFont(axis_label_font);
    htmp->GetXaxis()->SetLabelSize(axis_label_size);
    htmp->GetYaxis()->SetLabelSize(axis_label_size);
    htmp->GetXaxis()->CenterTitle(kTRUE);

    htmp->SetLineColor(data_color);
    htmp->SetLineWidth(data_line_width);
    htmp->SetLineStyle(data_line_style);
    htmp->SetMarkerStyle(data_marker);
    htmp->SetMarkerSize(data_marker_size);
    htmp->SetMarkerColor(data_color);

    int gaus1_color = kRed;
    int gaus2_color = kBlue;

    int line_style  = 1;

    // draw histogram
    htmp->DrawCopy("p e1 x0");

    // draw curves
    double min2 = htmp->GetXaxis()->GetXmin();
    double max2 = htmp->GetXaxis()->GetXmax();

    double min1 = min2/2.0;
    double max1 = max2/2.0;

    TF1 *f1 = new TF1("f1","gaus",min1,max1);
    TF1 *f2 = new TF1("f2","gaus",min2,max2);

    f1->SetLineColor(gaus1_color);
    f1->SetLineStyle(line_style);
    f1->SetLineWidth(mc_line_width);

    f2->SetLineColor(gaus2_color);
    f2->SetLineStyle(line_style);
    f2->SetLineWidth(mc_line_width);

    f1->SetParameters(pars);
    f2->SetParameters(&pars[3]);

    f1->Draw("same");
    f2->Draw("same");

    // print parameters information
    int text_align = 11;
    int index[4]   = { 1, 2, 4, 5 };

    double x1, y1, x2, y2;
    DecodeLegendPosition( x1, y1, x2, y2, legPos, 1,10 );

    for (int i = 0; i < 4; ++i)
    {
        int    text_color = i < 2 ? gaus1_color : gaus2_color;
        string text_line  = Form("%s = %.2f#pm%.2f",i==0?"#mu":"#sigma",pars[index[i]],fabs(errs[index[i]]));
        AddPlotLabel(text_line.c_str(),x1,y1,legend_text_size,text_color,legend_text_font,text_align);
        y1 -= i == 1 ? 0.055 : 0.037;
    }

    y1 -= 0.05;

    string text_line = Form("#chi^{2} / ndf = %.2f / %d",fit->GetChisquare(),fit->GetNDF());
    AddPlotLabel(text_line.c_str(),x1,y1,legend_text_size,kBlack,legend_text_font,text_align);

    gPad->Update();

    // clean up
    delete htmp;
    delete fit;
}

//Draw Model comparisons and return chi2 for all objects in the array (no 2D chi2 yet.)
vector<double> MnvPlotter::draw2DXSecModelComparisons(TObjArray allHists,vector<string> histnames, int nModels, double mcScale, string legPos, bool areaNormalize, string savePath, bool drawChi2, string XVar, string YVar, string target,bool drawSmooth,bool genie_ratio) {
  cout << "Drawing model comparisons with mc scaled by " << mcScale << " and area normalized = " << areaNormalize << ". These plots will be written to " << savePath << endl;
  vector<double> chi2;
  TCanvas *c1 = new TCanvas();
  string drawstyle = "SAMEHIST][";
  if (drawSmooth) drawstyle = "SAMEHISTL][";
  //First two elements of the allHists are the 2D data and mc. The next N entries are the model 2D histograms.
  //Grab the data 2D histogram to get labels, n-bins.
  TObject *dataObj = (TObject*)allHists.At(0);
  MnvH2D *tmpdata = (MnvH2D*)dataObj->Clone("tmpdata2d_plotting");
  string xAxisBaseLabel = tmpdata->GetXaxis()->GetTitle();
  string yAxisBaseLabel = tmpdata->GetYaxis()->GetTitle();
  string zAxisBaseLabel = tmpdata->GetZaxis()->GetTitle();
  const int nXBins = tmpdata->GetNbinsX();
  const int nYBins = tmpdata->GetNbinsY();
  //legend options
  vector<string> options;
  options.push_back("l");//data
  options.push_back("l");//genie
  for (int i=0;i<nModels;i++) {
    options.push_back("l");//models
  }
  Double_t edgesvecX[nXBins+1];
  Double_t edgesvecY[nYBins+1];
  for (int i=1;i<nXBins+2;i++) edgesvecX[i-1]=tmpdata->GetXaxis()->GetBinLowEdge(i);
  for (int i=1;i<nYBins+2;i++) edgesvecY[i-1]=tmpdata->GetYaxis()->GetBinLowEdge(i);

  double binXErrorsZero[nXBins];
  double binYErrorsZero[nYBins];

  for (int i=0;i<nXBins;i++) binXErrorsZero[i]=0;
  for (int i=0;i<nYBins;i++) binYErrorsZero[i]=0;

  int start = 0;
  int end = start+2+nModels;
  int ndf = 0;
  //Doing Chi2 for 2D!!
  for (int i=start;i<end;i++) {
    cout << ((MnvH2D*)allHists.At(i))->GetName() << endl;
    string histname = ((MnvH2D*)allHists.At(i))->GetName();
    chi2.push_back(MnvPlotter::Chi2DataMC((MnvH2D*)allHists.At(start),(MnvH2D*)allHists.At(i),ndf,mcScale,true,areaNormalize));
    if (i!=start) {
      TH2D *tmp_data = new TH2D(((MnvH2D*)allHists.At(0)->Clone("data_hist_ratio"))->GetCVHistoWithStatError());
      TH2D *tmp_model = new TH2D(((MnvH2D*)allHists.At(i)->Clone("data_hist_ratio"))->GetCVHistoWithStatError());
      tmp_data->Divide(tmp_model);
      tmp_data->GetXaxis()->SetRangeUser(1.5,20);
      tmp_data->GetYaxis()->SetRangeUser(0,2.5);
      tmp_data->GetZaxis()->SetTitle("");
      tmp_data->GetZaxis()->SetRangeUser(0.5,1.5);
      tmp_data->Draw("COLZ");
      string name = savePath+Form("_data_%s_ratio",histname.c_str());
      MnvPlotter::MultiPrint(c1,name);
      delete tmp_data;
      delete tmp_model;
    }
    else{
      TH2D mat = ((MnvH2D*)allHists.At(0))->GetTotalError(true,false,areaNormalize);
      TH2D *data = new TH2D(((MnvH2D*)allHists.At(0))->GetCVHistoWithError(true,areaNormalize));
      mat.GetZaxis()->SetTitle("");
      mat.Divide(data);
      mat.Draw("COLZ");
      string name = savePath+"_data_total_error";
      MnvPlotter::MultiPrint(c1,name);
      delete data;
    }
  }

  //Now let's draw full phase space projections
  //X projection
  start = 2+nModels;//Skip 2D plots
  end = start+2+nModels;//data,mc,nmodels
  vector<TH1*> hists;
  vector<string> mod_histnames;

  for (int i=start;i<end;i++) {
    cout << "Doing projection X << \t" << ((MnvH1D*)allHists.At(i))->GetName() << endl;
    if (i == start) {//This is the data histogram
      MnvPlotter::ApplyNextLineStyle((MnvH1D*)allHists.At(i) ,true,false);
      string label = "#frac{d#sigma}{"+XVar+"} (cm^{2}/GeV/"+target+")";
      string xlabel = XVar+ "(GeV)";
      xlabel.erase(xlabel.begin());
      ((MnvH1D*)allHists.At(i))->GetYaxis()->SetTitle(label.c_str());
      ((MnvH1D*)allHists.At(i))->GetXaxis()->SetTitle(xlabel.c_str());
      ((MnvH1D*)allHists.At(i))->SetTitle("");
      hists.push_back((TH1*)allHists.At(i));
      if (genie_ratio) {
	TH1D *tmp = new TH1D(((MnvH1D*)allHists.At(start+1))->GetCVHistoWithStatError());
	for (int b = 0;b<tmp->GetNbinsX();b++) tmp->SetBinError(b,0);//Zero out MC stat errors
	((MnvH1D*)allHists.At(i))->GetYaxis()->SetTitle("Ratio to GENIE");
	((MnvH1D*)allHists.At(i))->DivideSingle((MnvH1D*)allHists.At(i),tmp);//divide data by genie before bin width because genie hasn't been bin width normalized yet...
	((MnvH1D*)allHists.At(i))->SetMaximum(2);
	((MnvH1D*)allHists.At(i))->SetMinimum(0.4);
      }
      if (!genie_ratio) {
	((MnvH1D*)allHists.At(i))->Scale(1.0,"width");
	((MnvH1D*)allHists.At(i))->SetMaximum(1.5*((MnvH1D*)allHists.At(i))->GetBinContent(((MnvH1D*)allHists.At(i))->GetMaximumBin()));
      }
      ((MnvH1D*)allHists.At(i))->GetCVHistoWithError(true,areaNormalize).DrawClone("PE");
      if (!genie_ratio)chi2.push_back(MnvPlotter::Chi2DataMC((MnvH1D*)allHists.At(start),(MnvH1D*)allHists.At(i),ndf,mcScale,true,areaNormalize));

    }
    else{
      MnvPlotter::ApplyNextLineStyle((MnvH1D*)allHists.At(i) ,false,false);
      string label = "#frac{d#sigma}{"+XVar+"} (cm^{2}/GeV/"+target+")";
      ((MnvH1D*)allHists.At(i))->GetYaxis()->SetTitle(label.c_str());
      hists.push_back((TH1*)allHists.At(i));
      MnvH1D *tmp = (MnvH1D*)allHists.At(i)->Clone(Form("htmp_%d",i));
      MnvH1D *tmp2= (MnvH1D*)tmp->Rebin(nXBins,Form("htmp_%d_rebinned",i),edgesvecX);
      tmp2->SetError(binXErrorsZero);


      if (genie_ratio) {
	tmp2->Divide(tmp2,(MnvH1D*)allHists.At(start+1));//divide genie by genie after width and models by genie after bin width correction
	tmp2->DrawClone(drawstyle.c_str());
      }
      if (!genie_ratio) {
	tmp2->Scale(1.0,"width");
	((MnvH1D*)allHists.At(i))->Scale(1.0,"width");
	allHists.At(i)->DrawClone(drawstyle.c_str());
	chi2.push_back(MnvPlotter::Chi2DataMC((MnvH1D*)allHists.At(start),tmp2,ndf,mcScale,true,areaNormalize));
      }
    }
  }
  ((MnvH1D*)allHists.At(start))->GetCVHistoWithError(true,areaNormalize).DrawClone("SAMEPE1");
  ((MnvH1D*)allHists.At(start))->GetCVHistoWithStatError().DrawClone("SAMEPE1");
  MnvPlotter::MultiPrint(c1,savePath+"_FullProjX");

  //Y projection
  start = end;
  end = start+2+nModels;
  hists.clear();
  for (int i=start;i<end;i++) {
    cout << "Doing projection Y << \t" << ((MnvH1D*)allHists.At(i))->GetName() << endl;
    if (i == start) {//This is the data histogram
      MnvPlotter::ApplyNextLineStyle((MnvH1D*)allHists.At(i) ,true,false);
      string label = "#frac{d#sigma}{"+YVar+"} (cm^{2}/GeV/"+target+")";
      string xlabel = YVar+ "(GeV)";
      xlabel.erase(xlabel.begin());
      ((MnvH1D*)allHists.At(i))->GetYaxis()->SetTitle(label.c_str());
      ((MnvH1D*)allHists.At(i))->GetXaxis()->SetTitle(xlabel.c_str());
      ((MnvH1D*)allHists.At(i))->SetTitle("");
      hists.push_back((TH1*)allHists.At(i));
      if (genie_ratio) {
	TH1D *tmp = new TH1D(((MnvH1D*)allHists.At(start+1))->GetCVHistoWithStatError());
	for (int b = 0;b<tmp->GetNbinsX();b++) tmp->SetBinError(b,0);//Zero out MC stat errors
	((MnvH1D*)allHists.At(i))->GetYaxis()->SetTitle("Ratio to GENIE");
	((MnvH1D*)allHists.At(i))->DivideSingle((MnvH1D*)allHists.At(i),tmp);//divide data by genie before bin width because genie hasn't been bin width normalized yet...
	((MnvH1D*)allHists.At(i))->SetMaximum(2);
	((MnvH1D*)allHists.At(i))->SetMinimum(0.4);
      }
      if (!genie_ratio) {
	((MnvH1D*)allHists.At(i))->Scale(1.0,"width");
	((MnvH1D*)allHists.At(i))->SetMaximum(1.5*((MnvH1D*)allHists.At(i))->GetBinContent(((MnvH1D*)allHists.At(i))->GetMaximumBin()));
	chi2.push_back(MnvPlotter::Chi2DataMC((MnvH1D*)allHists.At(start),(MnvH1D*)allHists.At(i),ndf,mcScale,true,areaNormalize));
      }
      ((MnvH1D*)allHists.At(i))->GetCVHistoWithError().DrawClone("PE");

    }
    else{
      MnvPlotter::ApplyNextLineStyle((MnvH1D*)allHists.At(i) ,false,false);
      string label = "#frac{d#sigma}{"+YVar+"} (cm^{2}/GeV/"+target+")";
      ((MnvH1D*)allHists.At(i))->GetYaxis()->SetTitle(label.c_str());
      hists.push_back((TH1*)allHists.At(i));
      MnvH1D *tmp = (MnvH1D*)allHists.At(i)->Clone(Form("htmp_%d",i));
      MnvH1D *tmp2= (MnvH1D*)tmp->Rebin(nYBins,Form("htmp_%d_rebinned",i),edgesvecY);
      tmp2->SetError(binYErrorsZero);


      if (genie_ratio) {
	tmp2->Divide(tmp2,(MnvH1D*)allHists.At(start+1));//divide genie by genie after width and models by genie after bin width correction
	tmp2->DrawClone(drawstyle.c_str());
      }
      if (!genie_ratio) {
	//tmp2->Scale(1.0,"width");
	((MnvH1D*)allHists.At(i))->Scale(1.0,"width");
	allHists.At(i)->DrawClone(drawstyle.c_str());
	chi2.push_back(MnvPlotter::Chi2DataMC((MnvH1D*)allHists.At(start),tmp2,ndf,mcScale,true,areaNormalize));
      }
    }
  }
  ((MnvH1D*)allHists.At(start))->GetCVHistoWithError(true,areaNormalize).DrawClone("SAMEPE1");
  ((MnvH1D*)allHists.At(start))->GetCVHistoWithStatError().DrawClone("SAMEPE1");
  MnvPlotter::MultiPrint(c1,savePath+"_FullProjY");


  //Now lets do bin by bin
  //X projections (uses y binning)

  for (int j=0;j<nYBins;j++) {
    start = end;
    end = start+2+nModels;
    hists.clear();
    string output = savePath+Form("_ProjX_bin_%d",j);
    for (int i=start;i<end;i++) {
      cout << "Doing projdection X << \t" << ((MnvH1D*)allHists.At(i))->GetName() << "\t slice " << j<< endl;
      if (i == start) {//This is the data histogram
	MnvPlotter::ApplyNextLineStyle((MnvH1D*)allHists.At(i) ,true,false);
	string label = "#frac{d#sigma}{"+XVar+YVar+"} (cm^{2}/GeV^{2}/"+target+")";
	string xlabel = XVar+ "(GeV)";
	string ylabel = YVar;
	xlabel.erase(xlabel.begin());
	ylabel.erase(ylabel.begin());
	((MnvH1D*)allHists.At(i))->GetYaxis()->SetTitle(label.c_str());
	((MnvH1D*)allHists.At(i))->GetXaxis()->SetTitle(xlabel.c_str());
	((MnvH1D*)allHists.At(i))->SetTitle(Form("%13.2f GeV < %s #leq %2.2f GeV",edgesvecY[j], ylabel.c_str(),edgesvecY[j+1]));
	hists.push_back((TH1*)allHists.At(i));
	if (genie_ratio) {
	  TH1D *tmp = new TH1D(((MnvH1D*)allHists.At(start+1))->GetCVHistoWithStatError());//Get GENIE cv+stat only (other errors lurking around in the error band
	  for (int b = 0;b<tmp->GetNbinsX();b++) tmp->SetBinError(b,0);//Zero out MC stat errors
	  ((MnvH1D*)allHists.At(i))->GetYaxis()->SetTitle("Ratio to GENIE");
	  ((MnvH1D*)allHists.At(i))->DivideSingle((MnvH1D*)allHists.At(i),tmp);//divide data by genie before bin width because genie hasn't been bin width normalized yet...
	  // ((MnvH1D*)allHists.At(i))->GetYaxis()->SetTitle("Ratio to GENIE");
	  // ((MnvH1D*)allHists.At(i))->Divide((MnvH1D*)allHists.At(i),(MnvH1D*)allHists.At(start+1));//divide data by genie before bin width because genie hasn't been bin width normalized yet... Also DO NOT bin width normalize
	  ((MnvH1D*)allHists.At(i))->SetMaximum(2);
	  ((MnvH1D*)allHists.At(i))->SetMinimum(0.4);
	}
	if (!genie_ratio) {
	  ((MnvH1D*)allHists.At(i))->Scale(1.0,"width");
	  ((MnvH1D*)allHists.At(i))->Scale(1/(edgesvecY[j+1]- edgesvecY[j]));//Scale by perpendicular bin width
	  ((MnvH1D*)allHists.At(i))->SetMaximum(1.5*((MnvH1D*)allHists.At(i))->GetBinContent(((MnvH1D*)allHists.At(i))->GetMaximumBin()));
	}
	((MnvH1D*)allHists.At(i))->GetCVHistoWithError(true,areaNormalize).DrawClone("PE");
	if (!genie_ratio)chi2.push_back(MnvPlotter::Chi2DataMC((MnvH1D*)allHists.At(start),(MnvH1D*)allHists.At(i),ndf,mcScale,true,areaNormalize));
      }
      else{
	MnvPlotter::ApplyNextLineStyle((MnvH1D*)allHists.At(i) ,false,false);
	string label = "#frac{d#sigma}{"+XVar+YVar+"} (cm^{2}/GeV^{2}/"+target+")";
	((MnvH1D*)allHists.At(i))->GetYaxis()->SetTitle(label.c_str());
	hists.push_back((TH1*)allHists.At(i));
	MnvH1D *tmp = (MnvH1D*)allHists.At(i)->Clone(Form("htmp_%d_%d",i,j));
	MnvH1D *tmp2= (MnvH1D*)tmp->Rebin(nXBins,Form("htmp_%d_%d_rebinned",i,j),edgesvecX);
	tmp2->SetError(binYErrorsZero);


	if (genie_ratio) {
	  tmp2->Divide(tmp2,(MnvH1D*)allHists.At(start+1));//divide genie by genie after width and models by genie after bin width correction
	  tmp2->DrawClone(drawstyle.c_str());
	}
	if (!genie_ratio) {
	  tmp2->Scale(1.0,"width");
	  tmp2->Scale(1/(edgesvecY[j+1]- edgesvecY[j]));//Scale by perpendicular bin width
	  ((MnvH1D*)allHists.At(i))->Scale(1.0,"width");
	  ((MnvH1D*)allHists.At(i))->Scale(1/(edgesvecY[j+1]- edgesvecY[j]));//Scale by perpendicular bin width
	  ((MnvH1D*)allHists.At(i))->DrawClone(drawstyle.c_str());
	  chi2.push_back(MnvPlotter::Chi2DataMC((MnvH1D*)allHists.At(start),tmp2,ndf,mcScale,true,areaNormalize));
	}
      }
    }
    ((MnvH1D*)allHists.At(start))->GetCVHistoWithError(true,areaNormalize).DrawClone("SAMEPE1");
    ((MnvH1D*)allHists.At(start))->GetCVHistoWithStatError().DrawClone("SAMEPE1");
    MnvPlotter::MultiPrint(c1,output);
  }

  //Y projections (uses x binning)

  for (int j=0;j<nXBins;j++) {
    hists.clear();
    start = end;
    end = start+2+nModels;
    string output = savePath+Form("_ProjY_bin_%d",j);
    for (int i=start;i<end;i++) {
      cout << "Doing projdection Y << \t" << ((MnvH1D*)allHists.At(i))->GetName() << "\t slice " << j<< endl;
      if (i == start) {//This is the data histogram
	MnvPlotter::ApplyNextLineStyle((MnvH1D*)allHists.At(i) ,true,false);
	string label = "#frac{d#sigma}{"+XVar+YVar+"} (cm^{2}/GeV^{2}/"+target+")";
	string xlabel = YVar+ "(GeV)";
	string ylabel = XVar;
	xlabel.erase(xlabel.begin());
	ylabel.erase(ylabel.begin());
	((MnvH1D*)allHists.At(i))->GetYaxis()->SetTitle(label.c_str());
	((MnvH1D*)allHists.At(i))->GetXaxis()->SetTitle(xlabel.c_str());
	((MnvH1D*)allHists.At(i))->SetTitle(Form("%13.2f GeV < %s #leq %2.2f GeV",edgesvecX[j], ylabel.c_str() ,edgesvecX[j+1]));
	hists.push_back((TH1*)allHists.At(i));
	if (genie_ratio) {
	  TH1D *tmp = new TH1D(((MnvH1D*)allHists.At(start+1))->GetCVHistoWithError());//Get GENIE cv+stat only (other errors lurking around in the error band
	  for (int b = 0;b<tmp->GetNbinsX();b++) tmp->SetBinError(b,0);// Zero out MC stats?
	  ((MnvH1D*)allHists.At(i))->GetYaxis()->SetTitle("Ratio to GENIE");
	  ((MnvH1D*)allHists.At(i))->DivideSingle((MnvH1D*)allHists.At(i),tmp);//divide data by genie before bin width because genie hasn't been bin width normalized yet...
	  // ((MnvH1D*)allHists.At(i))->Divide(((MnvH1D*)allHists.At(i)),(MnvH1D*)allHists.At(start+1));//divide data by genie before bin width because genie hasn't been bin width normalized yet...
	  // ((MnvH1D*)allHists.At(i))->GetYaxis()->SetTitle("Ratio to GENIE");
	  ((MnvH1D*)allHists.At(i))->SetMaximum(2);
	  ((MnvH1D*)allHists.At(i))->SetMinimum(0.4);
	}

	if (!genie_ratio) {
	  ((MnvH1D*)allHists.At(i))->Scale(1.0,"width");
	  ((MnvH1D*)allHists.At(i))->Scale(1/(edgesvecX[j+1]- edgesvecX[j]));//Scale by perpendicular bin width
	  ((MnvH1D*)allHists.At(i))->SetMaximum(1.5*((MnvH1D*)allHists.At(i))->GetBinContent(((MnvH1D*)allHists.At(i))->GetMaximumBin()));
	  chi2.push_back(MnvPlotter::Chi2DataMC((MnvH1D*)allHists.At(start),(MnvH1D*)allHists.At(i),ndf,mcScale,true,areaNormalize));
	}
	((MnvH1D*)allHists.At(i))->GetCVHistoWithError(true,areaNormalize).DrawClone("PE");

      }
      else{
	MnvPlotter::ApplyNextLineStyle((MnvH1D*)allHists.At(i) ,false,false);
	string label = "#frac{d#sigma}{"+XVar+YVar+"} (cm^{2}/GeV^{2}/"+target+")";
	((MnvH1D*)allHists.At(i))->GetYaxis()->SetTitle(label.c_str());
	hists.push_back((TH1*)allHists.At(i));
	MnvH1D *tmp = (MnvH1D*)allHists.At(i)->Clone(Form("htmp_%d_%d",i,j));
	MnvH1D *tmp2= (MnvH1D*)tmp->Rebin(nYBins,Form("htmp_%d_%d_rebinned",i,j),edgesvecY);
	tmp2->SetError(binYErrorsZero);
	if (genie_ratio) {
	  tmp2->Divide(tmp2,(MnvH1D*)allHists.At(start+1));//divide genie by genie after width and models by genie after bin width correction
	  tmp2->DrawClone(drawstyle.c_str());
	}
	if (!genie_ratio) {
	  tmp2->Scale(1.0,"width");
	  tmp2->Scale( 1/(edgesvecX[j+1] - edgesvecX[j] ));//Scale by perpendicular bin width
	  ((MnvH1D*)allHists.At(i))->Scale(1.0,"width");
	  ((MnvH1D*)allHists.At(i))->Scale( 1/(edgesvecX[j+1] - edgesvecX[j] ));//Scale by perpendicular bin width
	  ((MnvH1D*)allHists.At(i))->DrawClone(drawstyle.c_str());
	  chi2.push_back(MnvPlotter::Chi2DataMC((MnvH1D*)allHists.At(start),tmp2,ndf,mcScale,true,areaNormalize));
	}
      }
    }
    ((MnvH1D*)allHists.At(start))->GetCVHistoWithError(true,areaNormalize).DrawClone("SAMEPE1");
    ((MnvH1D*)allHists.At(start))->GetCVHistoWithStatError().DrawClone("SAMEPE1");
    MnvPlotter::MultiPrint(c1,output);
  }
  TLegend leg(0,0,1,1);
  for (unsigned int i=0;i<hists.size();i++) {
    leg.AddEntry(hists[i],histnames[i].c_str(),options[i].c_str());
  }
  c1->Clear();
  leg.Draw();
  string output = savePath+"Legend";
  MnvPlotter::MultiPrint(c1,output);



  //Use MnvPlotter::Chi2DataMC to get chi. Add to plots!!


  //Use MnvPlotter::MultiPrint(c,name) include full path in name so /path/to/write/area/name the .X is controlled by the function. Might need to add C types.

  return chi2;

}



//#####################################################################################
//
// END
//
//#####################################################################################

#endif
