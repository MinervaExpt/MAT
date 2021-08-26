#ifndef MNV_MnvPlotter_h
#define MNV_MnvPlotter_h 1

#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"

#include "TDirectory.h"
#include "TClass.h"
#include "TObjArray.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TDecompSVD.h"

#include <iostream>
#include <vector>
#include <map>
#include <string>

namespace PlotUtils
{
  typedef TH1* TH1Ptr;
  //! Quickset style options
  enum t_PlotStyle
  {
    kDefaultStyle                = 0, //!< The Default style
    kCompactStyle                = 1, //!< Same as the default style but with less whitespace
    kNukeCCStyle                 = 2, //!< Similar to compact but with NukeCC-specifics
    kNukeCCPrintStyle            = 3,
    kCCNuPionIncStyle            = 4, //!< compact with CCNuPionInc color scheme
    kCCCohStyle                  = 5,  //!< Coherent rocks!!!
    kCCQENuStyle                 = 6,  //!< CCQENu
    kCCQENuInclusiveStyle        = 7,  //!< CCQENuInclusive
    kCCQEAntiNuStyle             = 8,  //!< CCQEAntiNu 
    kCCInclusiveHeliumStyle      = 9,  //!< cc-helium	
    kIMDStyle                    =10,  //!< CCQENuIMD
    kCCQENuTransverseStyle       =11,  //!< CCQENuTransverse
    kCCPi0AnaStyle               =12,  //!< CCPi0Ana style

 };

  //! enums for the available colour palettes, copied from ROOT 6's TColor. Usable via MnvPlotter::SetROOT6Palette()
#ifndef MNVROOT6
  enum EColorPalette {kDeepSea=51,          kGreyScale=52,    kDarkBodyRadiator=53,
                      kBlueYellow= 54, /* PAR removed kRainbow */   kInvertedDarkBodyRadiator=56,
                      kBird=57,             kCubehelix=58,    kGreenRedViolet=59,
                      kBlueRedYellow=60,    kOcean=61,        kColorPrintableOnGrey=62,
                      kAlpine=63,           kAquamarine=64,   kArmy=65,
                      kAtlantic=66,         kAurora=67,       kAvocado=68,
                      kBeach=69,            kBlackBody=70,    kBlueGreenYellow=71,
                      kBrownCyan=72,        kCMYK=73,         kCandy=74,
                      kCherry=75,           kCoffee=76,       /* PAR removed kDarkRainbow */
                      kDarkTerrain=78,      kFall=79,         kFruitPunch=80,
                      kFuchsia=81,          kGreyYellow=82,   kGreenBrownTerrain=83,
                      kGreenPink=84,        kIsland=85,       kLake=86,
                      kLightTemperature=87, kLightTerrain=88, kMint=89,
                      kNeon=90,             kPastel=91,       kPearl=92,
                      kPigeon=93,           kPlum=94,         kRedBlue=95,
                      kRose=96,             kRust=97,         kSandyTerrain=98,
                      kSienna=99,           kSolar=100,       kSouthWest=101,
                      kStarryNight=102,     kSunset=103,      kTemperatureMap=104,
                      kThermometer=105,     kValentine=106,   kVisibleSpectrum=107,
                      kWaterMelon=108,      kCool=109,        kCopper=110,
                      kGistEarth=111,       kViridis=112};
#endif

  class MnvPlotter
  {
    public:
      //! Default Constructor
      MnvPlotter();

      //! Constructor with a style
      explicit MnvPlotter( t_PlotStyle style );

      //! Default Destructor
      virtual ~MnvPlotter();


    //=======================================================
    // Garbage collection stuff
    //=======================================================
    private:
      //! Vector of temp objects to be garbage collected
      TObjArray fTmpObjects;

      //! Add a new object to tmpObjects
      void AddToTmp( TObject* obj );


    public:
      //! Delete the objects held in tmpObjects
      void CleanTmp();


    //=========================================
    // STYLE Variables
    //=========================================
    /** @name Public Style Variables
      These variables are used to control the style of plots.
      They are publicy accessible, so you can change them with 'mnvPlotter.var_name = X;'

    @{*/
    public:
      //-- chi2 calculation
      bool chi2_use_overflow_err; ///< Do you want to consider correlations to under/overflow bins in chi2 calculation?

      //! do you want to draw histograms normalized to bin width
      bool draw_normalized_to_bin_width;

      //-- marker settings
      int data_marker;
      int ratio_marker;
      double data_marker_size;
      double ratio_marker_size;

      //-- line settings
      int data_line_width;
      int data_line_style;
      int mc_line_width;
      int mc_line_style;
      int ratio_line_width;

      //-- cut arrow settings
      int arrow_line_width;
      int arrow_line_style;
      int arrow_line_color;
      float arrow_size;
      std::string arrow_type;

      //-- color settings
      int palette_style;

      int data_color;
      int mc_color;
      int mc_error_color;
      int mc_error_style;
      int ratio_color;

      int mc_bkgd_color;
      int mc_bkgd_width;
      int mc_bkgd_line_color;
      int mc_bkgd_style;

      int data_bkgd_color;
      int data_bkgd_style;
      double data_bkgd_size;

      //-- correlation/covariance
      bool draw_corr_max1;
      bool draw_corr_red_blue;
      unsigned int n_color_contours;

      //-- title settings
      int title_font;
      double title_size;

      //-- axis options
      bool   hist_min_zero;        ///< argument to SetHistMinimumZero()
      bool   axis_draw_grid_x;     ///< Set to true to get grid drawn for X axis
      bool   axis_draw_grid_y;     ///< Set to true to get grid drawn for Y axis
      int    axis_max_digits;
      int    axis_title_font_x;
      int    axis_title_font_y;
      int    axis_title_font_z;
      double axis_title_offset_x;
      double axis_title_offset_y;
      double axis_title_offset_z;
      double axis_title_size_x;
      double axis_title_size_y;
      double axis_title_size_z;
      double axis_minimum;        ///< argument to SetMinimum() use -1111 for automatic scale
      double axis_maximum;        ///< argument to SetMinimum() use -1111 for automatic scale
      double axis_maximum_group;  ///< argument to SetMinimum() use -1111 for automatic scale

      //-- axis label options
      int    axis_label_font;
      double axis_label_size;

      //-- margins
      double extra_top_margin;
      double extra_bottom_margin;
      double extra_left_margin;
      double extra_right_margin;

      //-- layout
      double headroom;
      double footroom; //because legroom made be think of room for a TLegend...

      //-- legend
      double height_nspaces_per_hist;
      double width_xspace_per_letter;
      int    legend_text_font;
      int    legend_fill_color;
      int    legend_border_size;
      double legend_text_size;
      double legend_offset_x;
      double legend_offset_y;
      unsigned int legend_n_columns;

      //--stat error name
      std::string stat_error_name;

      //-- colors
      typedef std::map<std::string, std::vector<std::string> > ErrorSummaryGroupMap;
      //! Vector of good ROOT colors
      std::vector<int> good_colors;
      //! Map from sys error name to ROOT color
      std::map<std::string,int> error_color_map;
      //! Name of a group of systematics that will be added in quadrature for the summary
      ErrorSummaryGroupMap error_summary_group_map;

      //! Vector of print formats
      std::vector<std::string> print_formats;
      //! Top directory for printed plots
      std::string print_topdir;

      //@}


    //=======================================================
    // Public functions
    //=======================================================
    public:
      //=========================================
      // BEGIN - Style functions
      //=========================================

        /** @name Style Functions
          These functions change the class's or ROOT's style setting.

          @{*/

        //! utility to set up basic root environment
        void SetRootEnv();

        //! A colour palette that goes blue->white->red, useful for correlation matrices
        void SetCorrelationPalette() const;
        //! A colour palette that goes white->red->redder, can be easier to read than raibow
        void SetRedHeatPalette() const;
        //! A rainbow colour palette that goes to white at 0, makes more sense than purple=0
        void SetWhiteRainbowPalette() const;
        //! A blackbody color palette
        void SetBlackbodyPalette() const;

        //! Set the global color palette to any ROOT 6 option. See: https://root.cern.ch/doc/master/classTColor.html#C06
        //! You can use any of the enums from EColorPalette, eg MnvPlotter::SetROOT6Palette(kRust)
        static void SetROOT6Palette(Int_t palette);

        //! utility to apply one of the PlotStyles
        void ApplyStyle(
            PlotUtils::t_PlotStyle style = PlotUtils::kDefaultStyle
            );

        //! Set the number of columns to be used for legends
        void SetLegendNColumns(int user_legend_n_columns);

        //! Set the style of all axes to this histogram (using the setting of this class)
        bool ApplyAxisStyle(
            TH1 *h,
            bool centerXTitle = true,
            bool centerYTitle = false,
            bool centerZTitle = false
            ) const;

        //! Set the style of all axes in this TDirectory (using the setting of this class)
        bool ApplyAxisStyle(
            TDirectory *td,
            bool centerXTitle = true,
            bool centerYTitle = false,
            bool centerZTitle = false
            ) const;

        //! Set the style of all axes in gDirectory
        bool ApplyAxisStyle(
            bool centerXTitle = true,
            bool centerYTitle = false,
            bool centerZTitle = false
            ) const;

        //! Set the axis limits to auto-detect mode
        void UseAutoAxisLimits();

        //! Reset the vector of good colors to the default
        //void ResetGoodColors();

        //!Add a cut arrow to the current pad
        void AddCutArrow(
            const double cut_location,
            const double y1,
            const double y2,
            const double arrow_length,
            const std::string& arrow_direction
            ) const;

        //! Get a distinct line style
        void ApplyNextLineStyle(
            TH1* h,
            bool startOver = false,
            bool changeStyle = true
            ) const;

        //@}

      //=========================================
      // END - Style functions
      //=========================================


      //==========================================================
      // Write Specific Information to Plots
      //==========================================================

        //! easily add Latex labels on plots
        void AddPlotLabel(
            const char* label,
            const double x,
            const double y,
            const double size = 0.05,
            const int color = 1,
            const int font = 62,
            const int align = 22,
            const double angle = 0
            );

        //! Obtain a better estimate for the length of a string in a legend by ignoring some special latex characters
        size_t GetLegendEntrySize(
            const std::string& title
            ) const;

        //! How many letters of width do we need for a legend with these titles
        //! @param[in] titles Vector of histogram titles
        //! return Sum of the size of the longest title in each column of the legend
        size_t GetLegendWidthInLetters(
            const std::vector<std::string>& titles
            ) const;

        //! decodes a string to determine location, alignment of plot label
        void DecodePosition(
            const std::string& opts,
            double size,
            int &align,
            double &xLabel,
            double &yLabel
            ) const;

        /*! decodes a string to determine the 4 corners of a TLegend
          @param[out] x1 position of left side of the legend
          @param[out] y1 position of bottom side of the legend
          @param[out] x2 position of right side of the legend
          @param[out] y2 position of top side of the legend
          @param[in] opts String description of position (see DecodePosition)
          @param[in] nHists Number of histograms in the legend
          @param[in] longestTitleSize Number of characters in the longest title entered into the legend
          @param[in] size Desired size of letters in the legend (not exact)
          */
        void DecodeLegendPosition(
            double &x1,
            double &y1,
            double &x2,
            double &y2,
            const std::string& opts,
            const size_t nHists,
            const size_t longestTitleSize = 14,
            double size = .04
            ) const;

        //! adds a title above the plot
        void AddHistoTitle (
            const char* title,
            double titleSize,
            int titleFont
            );

        void AddHistoTitle (
            const char* title,
            double titleSize
            )
        {
          AddHistoTitle( title, titleSize, title_font );
        }

        void AddHistoTitle( const char* title )
        {
          AddHistoTitle( title, title_size, title_font );
        }



        //! writes the chi2/ndf between two histograms on the plot
        void AddChi2Label(
            const TH1* dataHist,
            const TH1* mcHist,
            const std::string& opts
            );

        //! writes the chi2/ndf between two histograms on the plot - MnvH1D overload
        void AddChi2Label(
            const MnvH1D* dataHist,
            const MnvH1D* mcHist,
            const std::string& opts,
            const bool useDataErrorMatrix = false,
            const bool useOnlyShapeErrors = false
            );

        //! writes the chi2/ndf between two histograms on the plot
        void AddChi2Label(
            const TH1* dataHist,
            const TH1* mcHist,
            Double_t mcScale = 1.0,
            const std::string& opts = "BC",
            double size = 0.04,
            double yOffset = 0.0
            );

        //! writes the chi2/ndf between two histograms on the plot - MnvH1D overload
        void AddChi2Label(
            const MnvH1D* dataHist,
            const MnvH1D* mcHist,
            Double_t mcScale = 1.0,
            const std::string& opts = "BC",
            double size = 0.04,
            double yOffset = 0.0,
            const bool useDataErrorMatrix = false,
            const bool useOnlyShapeErrors = false,
            const std::string& pre_tag = ""
            );

        //! add 'preliminary' to the plot
        void WritePreliminary(
            const std::string& opts = "TC",
            double size = 0.035,
            double yOffset = 0,
            double xOffset = 0,
            bool isWorkInProgress = false
            );

        //! add 'preliminary' to the plot
        void WritePreliminary(
            double x,
            double y,
            double size = 0.035,
            bool isWorkInProgress = false
            );

        //! add 'normalization' label to the plot
        void WriteNorm(
            const char *norm,
            const std::string& opts = "TR",
            double size = 0.045,
            double yOffset = 0,
            double xOffset = 0,
            double pot = 0
            );

        //! add 'normalization' label to the plot
        void WriteNorm(
            const char *norm,
            double x,
            double y,
            double size = 0.045,
            double pot = 0
            );

        //! calculate the chi2 between two histograms
        //! Bin range comes from mcHist
        Double_t Chi2DataMC(
            const TH1* dataHist,
            const TH1 *mcHist,
            int& ndf,
            const Double_t mcScale = 1.0
            );

        //! calculate the chi2 between two histograms
        //! using the total Error Matrix
        //! Bin range comes from mcHist
        Double_t Chi2DataMC(
            const MnvH2D* dataHist,
            const MnvH2D* mcHist,
            Int_t& ndf,
            const Double_t mcScale = 1.0,
            const bool useDataErrorMatrix = false,
            const bool useOnlyShapeErrors = false,
            const bool useModelStat = true,
            TMatrixD *Chi2ByBin = NULL
            );

        //! calculate the chi2 between two histograms
        //! using the total Error Matrix
        //! Bin range comes from mcHist
        Double_t Chi2DataMC(
            const MnvH1D* dataHist,
            const MnvH1D* mcHist,
            int& ndf,
            const Double_t mcScale = 1.0,
            const bool useDataErrorMatrix = false,
            const bool useOnlyShapeErrors = false,
            const bool useModelStat = true,
            TMatrixD *Chi2ByBin = NULL
            );

        //! calculate the chi2 between two histograms
        //! Bin range comes from mcHist
        Double_t Chi2DataMC(
            const TH1* dataHist,
            const TH1 *mcHist,
            const Double_t mcScale = 1.0
            );

        //! calculate the chi2 between two histograms
        //! using the total Error Matrix
        //! Bin range comes from mcHist
        Double_t Chi2DataMC(
            const MnvH2D* dataHist,
            const MnvH2D* mcHist,
            const Double_t mcScale = 1.0,
            const bool useDataErrorMatrix = false,
            const bool useOnlyShapeErrors = false,
            const bool useModelStat = true,
            TMatrixD *Chi2ByBin = NULL
            );

        //! calculate the chi2 between two histograms
        //! using the total Error Matrix
        //! Bin range comes from mcHist
        Double_t Chi2DataMC(
            const MnvH1D* dataHist,
            const MnvH1D* mcHist,
            const Double_t mcScale = 1.0,
            const bool useDataErrorMatrix = false,
            const bool useOnlyShapeErrors = false,
            const bool useModelStat = true,
            TMatrixD *Chi2ByBin = NULL
            );

        //! calculate the chi2 between two histograms
        //! each with a full covariance matrix
        //! using the total Error Matrix
        //! taking into account hist-hist correlations
        //! Bin range comes from histA
        Double_t Chi2MCMC(
            const MnvH1D* histA,
            const MnvH1D* histB,
            int& ndf,
            const Double_t bScale = 1.0,
            const bool useOnlyShapeErrors = false
            );

        void AddPOTNormBox(
            const double dataPOT,
            const double mcPOT,
            const double xLeft,
            const double yTop,
            const double size = 0.03
            );

        void AddAreaNormBox(
            const double dataScale,
            const double mcScale,
            const double xLeft,
            const double yTop,
            const double size = 0.03
            );

      //==========================================================
      // draw specific histogram comparisons
      //==========================================================

        //! draw MC data points on current Pad
        void DrawDataMC(
            const TH1* dataHist,
            const TH1* mcHist,
            const Double_t mcScale = 1.0, /*= 1.0*/
            const std::string& legPos = "L",
            const bool useHistTitles = false  //default is to use Data/MC/Background
            );

        //! draw MC with error band + data points on current Pad
        void DrawDataMCWithErrorBand(
            const TH1* dataHist,
            const TH1* mcHist,
            const Double_t mcScale = 1.0,
            const std::string& legPos = "L",
            const bool useHistTitles = false,  //default is to use Data/MC/Background
            const TH1* bkgdHist = NULL,
            const TH1* dataBkgdHist = NULL,
            const bool isSmooth = false
            );


        //! draw MC with error band + data points on current Pad - MnvH1D overload.
        void DrawDataMCWithErrorBand(
            const MnvH1D* dataHist,
            const MnvH1D* mcHist,
            const Double_t mcScale = 1.0,
            const std::string& legPos = "L",
            const bool useHistTitles = false,  //default is to use Data/MC/Background
            const MnvH1D* bkgdHist = NULL,
            const MnvH1D* dataBkgdHist = NULL,
            const bool covAreaNormalize = false,
            const bool statPlusSys = false,
            const bool isSmooth = false
            );

        //! draw data/MC ratio in current pad
        void DrawMCWithErrorBand(
            const TH1 *mcHist,
            const Double_t mcScale = 1.0,
            const TH1* bkgdHist = NULL,
            const bool isSmooth = false
            );

        bool AddSysError(
            TH1D *h,
            const TH1D *hErr,
            bool sumSquares = true
            ) const;

        void DrawErrorBand(
            const TH1 *h,
            const int color
            );

        //! draw MC with error band + data points on current Pad
        void DrawDataMCRatio(
            const TH1* dataHist,
            const TH1 *mcHist,
            const Double_t mcScale = 1.0,
            const bool drawOneLine = true,
            const double plotMin = -1.0,    // Default to setting plot min and max automatically.
            const double plotMax = -1.0,
            const char* yaxisLabel = "Data / MC"
            );  // Default to setting plot min and max automatically.

        //! Draw data/MC ratio in current pad with a sysematic error envelope
        void DrawDataMCRatio(
            const MnvH1D *dataHist,
            const MnvH1D *mcHist,
            const Double_t mcScale = 1.0,
            const bool drawSysLines = true,
            const bool drawOneLine  = true,
            const double plotMin = -1.0,    // Default to setting plot min and max automatically.
            const double plotMax = -1.0,
            const char* yaxisLabel = "Data / MC",
            const bool covAreaNormalize = false
            );  // Default to setting plot min and max automatically.

        /*! Draw data and a bunch of MC, or generically a histogram with points and other histograms as lines.
          Axis labels come from dataHist.
          @param[in] dataHist the data histogram
          @param[in] mcHists list of MC histograms
          @param[in] mcScale each MC histogram is scaled by this amount
          @param[in] legPos legend placement code
          @param[in] dataAsPoints draw data as points.  Alternative is a solid black line.
          @param[in] allSolidLines Do not apply different line styles.  Make them all solid.
          @param[in] leaveStyleAlone Do not apply any line/marker styles or colors
          */
        void DrawDataMCVariations(
            const MnvH1D* dataHist,
            const TObjArray* mcHists,
            const Double_t mcScale = 1.0,
            const std::string& legPos = "TR",
            const bool dataAsPoints = true,
            const bool allSolidLines = false,
            const bool leaveStyleAlone = false,
            const bool covAreaNormalize = false
            );

        //! adds a legend to the plot
        void AddPlotLegend(
            const std::vector< TH1* >       & hists,
            const std::vector< std::string >& names,
            const std::vector< std::string >& opts ,
            const std::string& legPos ="R"
            );

        //! adds a legend to the plot
        void AddPlotLegend(
            const std::vector< TH1* >       & hists,
            const std::vector< std::string >& names,
            const std::vector< std::string >& opts ,
            const double x,
            const double y,
            const double x_width = 0.25,
            const double y_width = 0.15,
            const double textSize = 0.02
            );


        //! print a canvas in multiple formats
        void MultiPrint(
            TCanvas *c,
            const std::string& name = ""
            ) const;

        //! supply a comma-separated list of formats you want to print
        void MultiPrint(
            TCanvas *c,
            const std::string& name,
            const std::string& typeStr
            ) const;


        //! get the mean of a histogram in some restricted range
        double GetHistoMean(
            const TH1 * hist,
            const int minBin = -1,
            const int maxBin = -1
            ) const;

        double GetHistoMean(
            const TH1 * hist,
            double & err,
            const int minBin = -1,
            const int maxBin = -1
            ) const;

        //! Draw a hexagon
        void DrawHex(
            const double apothem,
            const double xCenter = 0.0,
            const double yCenter = 0.0
            ) const;


      //===========================================================================
      // Histogram axis gymnastics
      //===========================================================================

        //! Reverse the axis
        void ReverseXAxis( TH1 *h );

        //! Reverse X  axis of a 2D hist
        void ReverseXAxis( TH2 *h );

        //! Create bins suitable for plotting in log10 scale
        static bool SetLogBins(
            const int nbins,
            const double min,
            const double max,
            double *bins
            );

        /*! What are the highest and lowest bins with BinContent greater than some value?
          @param[in] h The histogram
          @param[out] lowBin The lowest bin with content above minContentLow
          @param[out] highBin The highest bin with content above minContentHigh
          @param[in] minContentLow BinContent must be greater than this to be considered for the lowest bin with content
          @param[in] minContentHigh BinContent must be greater than this to be considered for the highest bin with content
          @return true if everything worked.  false if there was a problem.
          */
        bool GetBinRangeWithMinimumBinContent(
            const TH1 *h,
            int &lowbin,
            int &highbin,
            double minContentLow = 0.,
            double minContentHigh = 0.
            ) const;

        //! Cal GetBinRangeWithMinimumBinContent and require bin content to be greater than 0.
        bool GetNonZeroBinRange(
            const TH1 *h,
            int &lowbin,
            int &highbin
            ) const;

        /*! The min and max values among the histograms
          @param[in] hists Array of histograms
          @param[out] minmin Stores the minimum value found in the histograms
          @param[out] maxmax Stores the maximum value found in the histograms
          @return true if everything worked
          */
        bool GetHistsRange(
            const TObjArray& hists,
            double& minmin, double& maxmax
            ) const;


      //=======================================
      // Get histograms of errors using user defined error groups
      //=======================================
        /*! Get histograms of errors using user defined error groups
          @return Vector of histogram clones that the caller owns (delete them!)
          @param[in] h The MnvH1D whose errors you want
          @param[in] doFractional Do you want fraction errors
          @param[in] covAreaNormalize Do you want to get shape only errors?
          @param[in] ignoreThreshold Histograms with max total/fractional error below this value will be omitted
          */
        std::vector<TH1*> GetSysErrorGroupHists(
            MnvH1D* h,
            const bool doFractional = false,
            const bool covAreaNormalize = false,
            const double ignoreThreshold = 0.00001
            ) const;


      //====================================================================================
      // Draw stuff using Mnv classes
      //====================================================================================
        /*! Draw all error sources and their sum in quadrature on a single plot
          @param[in] h The MnvH1D whose errors you want to draw
          @param[in] legPos Specifies the placement of the legend
          @param[in] includeStat Do you want to include statistical errors when calculating the total error?
          @param[in] solidLinesOnly Draw all lines as solid instead of varying the line style
          @param[in] ignoreThreshold If an error source never contributes more than this fractional error, don't draw it
          @param[in] covAreaNormalize Do you want to compute the covariance on the shape only instead of the full covariance?
          @param[in] errorGroupName Draw an error summary using only errors in this group.  Default is to draw all errors
          @param[in] asFrac Do you want the errors to be shown as fractional (relative to the CV) instead of absolute?
          @param[in] Ytitle Y-axis title (default: "Fractional Uncertainty")
          @param[in] ignoreUngrouped  Do you want to show only errors that are in error groups?
          */
        bool DrawErrorSummary(
            MnvH1D* h,
            const std::string& legPos   = "TR",
            const bool         includeStat     = true,
            const bool         solidLinesOnly  = true,
            const double       ignoreThreshold = 0.00001,
            const bool         covAreaNormalize = false,
            const std::string& errorGroupName = "",
            const bool         asfrac = true,
            const std::string  &Ytitle = "",
            bool               ignoreUngrouped = false
            );

        /*! Draw all error sources and their sum in quadrature on a single plot
          @param[in] hMC The MonteCarlo MnvH1D whose systematic (and maybe statistical) errors you want to draw
          @param[in] hData The Data MnvH1D whose statistical errors you want to draw
          @param[in] legPos Specifies the placement of the legend
          @param[in] includeMCStat Do you want to include MC statistical errors in the total errors and statistical errors?
          @param[in] solidLinesOnly Draw all lines as solid instead of varying the line style
          @param[in] ignoreThreshold If an error source never contributes more than this fractional error, don't draw it
          @param[in] yaxisMin If this value is equal to yaxisMax, set bounds automatically. Otherwise use this value.
          @param[in] yaxisMax If this value is equal to yaxisMin, set bounds automatically. Otherwise use this value.
          @param[in] errorGroupName Draw an error summary using only errors in this group.  Default is to draw all errors
          */
        /* Deprecated see cxx for reason. DGR
           bool DrawDataMCErrorSummary(
           MnvH1D* hMC,
           MnvH1D* hData,
           const std::string& legPos   = "TR",
           const bool    includeMCStat   = true,
           const bool    solidLinesOnly  = true,
           const double  ignoreThreshold = 0.00001,
           const bool covAreaNormalize = false,
           const std::string& errorGroupName = ""
           );
           */
        //!Draw a mc stack
        void DrawStackedMC(
            const TObjArray* mcHists,
            const Double_t mcScale = 1.0,
            const std::string& legPos = "L",
            const Int_t mcBaseColor = 2,                  // Color of first stacked histo. // enter -1 to use the existing fill colors
            const Int_t mcColorOffset = 1,                // 2nd histo in stack has color mcBaseColor+mcColorOffset, etc.
            const Int_t mcFillStyle = 3001,               // Fill style for histograms only // enter -1 to use the existing fill style
            const char* xaxislabel = "",
            const char* yaxislabel = ""
            );

        //! Draw an overlay plot with a stack and a "data" histogram.
        void DrawDataStackedMC(
            const MnvH1D* dataHist,
            const TObjArray* mcHists,
            const Double_t mcScale = 1.0,
            const std::string& legPos = "L",
            const std::string& dataName = "Data",
            const Int_t mcBaseColor = 2,                  // Color of first stacked histo. // enter -1 to use the existing fill colors
            const Int_t mcColorOffset = 1,                // 2nd histo in stack has color mcBaseColor+mcColorOffset, etc.
            const Int_t mcFillStyle = 3001,               // Fill style for histograms only // enter -1 to use the existing fill style
            const char* xaxislabel = "",
            const char* yaxislabel = "",
            bool cov_area_normalize = false
            );

        //! Draw an overlay plot with a stack and a "data" histogram using user color setting
        void DrawDataStackedMC(
            const MnvH1D* dataHist,
            const TObjArray* mcHists,
            const Int_t* mcColors,                  // Array of Color for mc histograms
            const Double_t mcScale = 1.0,
            const std::string& legPos = "L",
            const std::string& dataName = "Data",
            const Int_t mcFillStyle = 3001,               // Fill style for histograms only // enter -1 to use the existing fill style
            const char* xaxislabel = "",
            const char* yaxislabel = "",
            bool cov_area_normalize = false
            );

        //! Draw an overlay plot with a stack and a "data" histogram.
        void DrawDataStackedMCWithErrorBand( );

        void DrawNormalizedMigrationHistogram(
            const TH2D* h_migration,
            const bool drawAsMatrix = false,
            const bool coarseContours = false,
            const bool includeFlows = true,
            const bool no_text = false
            );

        void DrawAllUniverses(
            const MnvH1D* h,
            const bool covAreaNormalize = false,
            const bool binWidthNormalize = true
            );

        void DrawErrorMatrices(
            const MnvH1D *h,
            const bool area_norm,
            const bool asCorr,
            const bool asFrac
            );

        void DrawErrorMatrix(
            const TMatrixD &matrix,
            const TAxis *axis,
            const double maximum = -1.0
            );

        /*! Draw a histogram fitted to a double Gaussian and returns both Gaussian fit parameters
          @param[in] h The histogram to fit with a double Gaussian
          @param[in] lowFitBound The bin low edge for fitting the histogram
          @param[in] highFitBound The bin high edge for fitting the histogram
          @param[in] legPos The location of the printed parameters and errors
          @param[out] parameters The double Gaussian parameters
          @param[out] errors The double Gaussian errors
          @param[out] chisquared the chisquared from the fit
          @param[out] ndf the degrees of freedom from the fit
          */
        void DrawDoubleGausFit(
            const   TH1D *h,
            double  lowFitBound  = -9999,
            double  highFitBound = -9999,
            const char* legPos   = "TR",
            double* parameters   = (double*)NULL,
            double* errors       = (double*)NULL,
            double  chisquared   = 1.0,
            int     ndf          = 1
            );
        std::vector< double > draw2DXSecModelComparisons(
            TObjArray allHists,
            std::vector <std::string> histnames,
            int nModels          = 0,
            double mcScale       = 1.0,
            std::string legPos   = "TR",
            bool areaNormalize   = false,
            std::string savePath = "",
            bool drawChi2        = false,
            std::string XVar = "",
            std::string YVar = "",
            std::string target = "",
            bool drawSmooth      = false,
            bool genie_ratio     = false
            );

  }; //end of MnvPlotter class

}//close MnvPlot namespace

#endif
